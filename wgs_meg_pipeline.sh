#!/usr/bin/env bash

#################
#################
##
##
##
##
##
##
##
##
##
#################

if [ ! "$BASH_VERSION" ] ; then
    echo "Please do not use sh to run this script ($0), just execute it directly." 1>&2
    exit 1
fi

shopt -s extglob


########
# Help #
########

#Display help menu
display_help () {
    echo "
    Usage: wgs_meg_pipeline.sh [options] -1 forward_reads.fastq -2 reverse_reads.fastq
    
        -h | --help                 help info
        
    Pipeline Options:
        -1 | --i1   FILE            Forward input fastq file
        -2 | --i2   FILE            Reverse input fastq file
        -o | --output DIR           Directory for output of important files (main output dir)
        -s | --sample_name  STR     Main file name for this sample throughout the pipeline
        -spp| --species STR         Species identifier string (see documentation for details)
        -t | --threads  INT         Threads to use where applicable
        -td| --temp_dir DIR         Temporary directory for intermediate dirs/files

"
}


###############
## Variables ##
###############
## Paths to program executables or jar files where noted
bbmap="/s/angus/index/common/tools/BBMap_35.85/bbmap/bbmerge.sh"
metamos="/s/angus/index/common/tools/metAMOS-1.5rc3"
nucmer="/s/angus/index/common/tools/MUMmer3.23/nucmer"
cisa="/s/angus/index/common/tools/CISA1.3"
blastdb="/usr/bin/makeblastdb"
blastn="/usr/bin/blastn"
prokka="/s/angus/index/common/tools/prokka/prokka-1.11/bin/prokka"

## Paths to output directories
temp_dir=""
output_dir=""

## Flags and variables used in the pipeline
sample_name=""
spp_pipeline="Lmonocytogenes"
best_assembly=""  # used in the event of less than 3 assemblies from metamos
insert=0  # insert size as determined by bbmerge
kmer=0  # k-mer size determined by specialk in metamos
genome_size=0  # genome size for this particular organism
threads=1  # threads to use where applicable, recommend ~ 20


###############
## Functions ##
###############
## Check that the program and directory paths exist
validate_paths() {
    echo -e "\n@.@\nMEG WGS Pipeline started on `date`\n" >> WGS_LabNotebook.txt

    local missing=""

    echo "
    Sample name this run:     ${sample_name}
    Species for this run:     ${spp_pipeline}
    
    Paths currently set:
    BBmap bbmerge:            ${bbmap}
    iMetAmos:                 ${metamos}
    nucmer:                   ${nucmer}
    CISA:                     ${cisa}
    makeblastdb:              ${blastdb}
    blastn:                   ${blastn}
    prokka:                   ${prokka}

    
    Directories selected:
    Temporary directory:      ${temp_dir}
    Output directory:         ${output_dir}
   	
   	" >> WGS_LabNotebook.txt
   	
   	if [ ! -e "${bbmap}" ]; then
        local missing="$missing:BBMap;"
    fi
    
    if [ ! -e "${metamos}" ]; then
    	local missing="$missing:iMetAmos;"
    fi
    
    if [ ! -e "${nucmer}" ]; then
    	local missing="$missing:nucmer;"
    fi
    
    if [ ! -e "${cisa}" ]; then
    	local missing="$missing:CISA;"
    fi
    
    if [ ! -e "${blastdb}" ]; then
    	local missing="$missing:makeblastdb;"
    fi
    
    if [ ! -e "${blastn}" ]; then
    	local missing="$missing:blastn;"
    fi
    
    if [ ! -e "${prokka}" ]; then
    	local missing="$missing:prokka;"
    fi
    
    #--------------------------------
    
    if [ "${sample_name}" == "" ]; then
        local missing="$missing:SampleName;"
    fi
    
    if [ "${spp_pipeline}" == "" ]; then
        local missing="$missing:SpeciesPipeline;"
    fi
    
    #--------------------------------
    
    if [ ! -d "${temp_dir}" ]; then
        local missing="$missing:TemporaryDirectory;"
    fi
    
    if [ ! -d "${output_dir}" ]; then
        local missing="$missing:OutputDirectory;"
    fi
    
    ## Check if components are missing
    if [ ! -e $missing ]; then
        echo -e "Directories or Paths are missing, check the following:\n"
        echo -e $missing | sed 's/;/\n/g'
        exit 1
    else
        echo -e "\nPaths and Directories are valid. Proceeding..."
    fi
}


validate_inputs() {
	filename1=$( basename $1 )
	filename2=$( basename $2 )
	if [ ! -e "$1" ] || [ ! -e "$2" ]; then
		echo -e "$1 or $2 does not exist"
		exit 1
	elif [ ! "${spp_pipeline}" == "Lmonocytogenes" ]; then
	    echo -e "${spp_pipeline} is not in the list of accepted species identifiers. See documentation"
	    echo -e "${spp_pipeline} is not in the list of accepted species identifiers. See documentation" >> WGS_LabNotebook.txt
	    exit 1
	else
		echo -e "Begin WGS Pipeline for $filename1 and $filename2..."
		echo -e "\nBegin WGS Pipeline for $filename1 and $filename2..." >> WGS_LabNotebook.txt
	fi
}


get_versions() {
	echo -e "Software versions:" >> WGS_LabNotebook.txt
	$bbmap | grep "^BBMerge" >> WGS_LabNotebook.txt
	head -n 1 ${metamos}/README.md | grep -Po "MetAMOS.*\" " >> WGS_LabNotebook.txt
	$nucmer -version 3>&1 1>&2 2>&3 3>&- | grep "NUCmer" >> WGS_LabNotebook.txt
	echo $cisa | grep -Po "CISA[0-9].[0-9]" >> WGS_LabNotebook.txt
	$blastdb -version | head -n 1 >> WGS_LabNotebook.txt
	$blastn -version | head -n 1 >> WGS_LabNotebook.txt
	$prokka --version 3>&1 1>&2 2>&3 3>&- >> WGS_LabNotebook.txt
}



bbmap_insert_size() {
	## We need to capture stderr from bbmerge to regex capture the insert range.
	## Below is the way to do that by juggling output streams: 3>&1 1>&2 2>&3 3>&-
	insert=$( $bbmap in1="$1" in2="$2" out=/dev/null 3>&1 1>&2 2>&3 3>&- | grep "^Insert range" | grep -Po "[0-9]{1,4} - [0-9]{1,4}" | sed 's/ - /:/' )
	
	## Check to make sure we succeeded
	if [ "$insert" == "0" ]; then
	    echo -e "Error: insert size detection with BBMerge failed, insert size is still zero."
	    echo -e "Error: insert size detection with BBMerge failed, insert size is still zero." >> WGS_LabNotebook.txt
	    exit 1
	fi
}


imetamos_run() {
    ## Run the iMetAMOS pipeline for assembly.  This produces four files that we need, one for each assembler
    ## We then record the location of these and general other information for input to CISA later.
    ## We will create BOTH CISA input files within this function, then check to see they are of valid
    ## format to run with CISA.
    ## Note that this function utilizes global variables.
    
    #if [ ! -f "${output_dir}"/velvet_assembly.contig] || [ ! -f "${output_dir}"/abyss_assembly.contig ] || [ ! -f "${output_dir}"/spades_assembly.contig ] || [ ! -f "${output_dir}"/idba_assembly.contig ] || [ ! -f "${output_dir}"/edena_assembly.contig ]; then
        echo -e "Beginning MetAMOS assembly using idba-ud, velvet, spades, and abyss..."
        echo -e "Beginning MetAMOS assembly using idba-ud, velvet, spades, and abyss..." >> WGS_LabNotebook.txt
        #"${metamos}"/initPipeline -q -1 $forward -2 $reverse -d "${temp_dir}/${sample_name}" -i $insert -W iMetAMOS
        #kmer=$( "${metamos}"/runPipeline -p $threads -t eautils -q -a velvet,spades,idba-ud,abyss,edena -b -z genus -d "${temp_dir}/${sample_name}" | grep "Selected kmer size" | grep -Po "[0-9]{1,4}" )
        
        #echo "Kmer test: $kmer" >> kmer_test.txt  # debugging
        kmer="71"
        
        valid_assemblies=0
        which_assemblies=()
        
        if [ -f "${temp_dir}/${sample_name}/Assemble/out/idba-ud.${kmer}.seq100.contig" ]; then
            ((valid_assemblies++))
            which_assemblies=(${which_assemblies[@]} ${output_dir}/idba_assembly.contig)
            cp "${temp_dir}/${sample_name}/Assemble/out/idba-ud.${kmer}.seq100.contig" "${output_dir}"/idba_assembly.contig
        fi
        
        if [ -f "${temp_dir}/${sample_name}/Assemble/out/spades.${kmer}.seq100.contig" ]; then
            ((valid_assemblies++))
            which_assemblies=(${which_assemblies[@]} ${output_dir}/spades_assembly.contig)
            cp "${temp_dir}/${sample_name}/Assemble/out/spades.${kmer}.seq100.contig" "${output_dir}"/spades_assembly.contig
        fi
        
        if [ -f "${temp_dir}/${sample_name}/Assemble/out/velvet.${kmer}.seq100.contig" ]; then
            ((valid_assemblies++))
            which_assemblies=(${which_assemblies[@]} ${output_dir}/velvet_assembly.contig)
            cp "${temp_dir}/${sample_name}/Assemble/out/velvet.${kmer}.seq100.contig" "${output_dir}"/velvet_assembly.contig
        fi
        
        if [ -f "${temp_dir}/${sample_name}/Assemble/out/abyss.${kmer}.seq100.contig" ]; then
            ((valid_assemblies++))
            which_assemblies=(${which_assemblies[@]} ${output_dir}/abyss_assembly.contig)
            cp "${temp_dir}/${sample_name}/Assemble/out/abyss.${kmer}.seq100.contig" "${output_dir}"/abyss_assembly.contig
        fi
        
        if [ -f "${temp_dir}/${sample_name}/Assemble/out/edena.${kmer}.seq100.contig" ]; then
            ((valid_assemblies++))
            which_assemblies=(${which_assemblies[@]} ${output_dir}/edena_assembly.contig)
            cp "${temp_dir}/${sample_name}/Assemble/out/edena.${kmer}.seq100.contig" "${output_dir}"/edena_assembly.contig
        fi
    #else
    #    echo -e "MetAMOS assembly files already detected for sample, proceeding..."
    #    echo -e "MetAMOS assembly files already detected for sample, proceeding..." >> WGS_LabNotebook.txt
    #fi
    
    ## This note is for when we run into bugs with bad assemblies (and we will at some point):
    ## Make sure this section is assigning variables correctly for best assembly
    if [ "$valid_assemblies" == "0" ]; then
        echo -e "Error: Zero valid assemblies detected.  Check error logs."
        echo -e "Error: Zero valid assemblies detected.  Check error logs." >> WGS_LabNotebook.txt
    elif [ "$valid_assemblies" -lt "3" ]; then
        echo -e "Less than three valid assemblies detected, skipping CISA and continuing with annotation..."
        echo $which_assemblies | sed 's/;/\n/g'
        echo -e "Less than three valid assemblies detected, skipping CISA and continuing with annotation..." >> WGS_LabNotebook.txt
        echo $which_assemblies | sed 's/;/\n/g' >> WGS_LabNotebook.txt
        best_assembly=$( head -n 1 "${temp_dir}/${sample_name}/Postprocess/out/best.asm" | sed 's/\n//' | sed 's/\..*//' | sed 's/-ud$//' )
        best_assembly="${output_dir}/${best_assembly}_assembly.contig"
    else
        ## Proceed with CISA config file creation
        
        echo -e "Creating CISA config files for ${spp_pipeline} organism..."
        echo -e "Creating CISA config files for ${spp_pipeline} organism..." >> WGS_LabNotebook.txt
        
        if [ "${spp_pipeline}" == "Lmonocytogenes" ]; then
            genome_size="2944528"
        fi
        
        echo -e "count=${valid_assemblies}" >> "${temp_dir}/Merge.config"
        for i in "${which_assemblies[@]}"; do
            fullpath=$( readlink -f  "$i" )
            assembler=$( echo $i | grep -Po "(abyss|idba|velvet|spades|edena)" )
            echo -e "data=${fullpath},title=${assembler}" >> "${temp_dir}/Merge.config"
        done
        echo -e "Master_file=${output_dir}/${sample_name}_merged_contigs.fa" >> "${temp_dir}/Merge.config"
        
        
        echo -e "genome=${genome_size}\ninfile=${output_dir}/${sample_name}_merged_contigs.fa\noutfile=${output_dir}/${sample_name}_cisa_integrated.fa" >> "${temp_dir}/CISA.config"
        echo -e "nucmer=${nucmer}\nR2_gap=0.95\nCISA=${cisa}\nmakeblastdb=${blastdb}\nblastn=${blastn}" >> "${temp_dir}/CISA.config"
    fi
    
    rm flowchart.svg
}


cisa_run() {
    ## CISA takes the assemblies from iMetAMOS and integrates them into a single file
    ## It then takes that file and computes a combined set of contigs if possible
    $cisa/Merge.py "${temp_dir}/Merge.config"
    echo "y" | $cisa/CISA.py "${temp_dir}/CISA.config"
    best_assembly="${output_dir}/${sample_name}_cisa_integrated.fa"
    
    if [ ! -f "$best_assembly" ]; then
        echo -e "CISA failed to produce the expected merged file.  Please check the logs."
        echo -e "CISA failed to produce the expected merged file.  Please check the logs." >> WGS_LabNotebook.txt
    fi
    
    for i in ${output_dir}/*; do
        if [[ $i =~ (.p.fa) ]]; then
            rm $i
        fi
    done
    
    rm -rf CISA1 CISA2 CISA3 CISA4
    rm CISA2_thr.out
    rm info info1 info2 Merge_info
    rm R1_Contigs.fa
    rm R2_Contigs.fa
}


prokka_annotate() {
    ## Prokka uses several methods to annotate the best assembly based on all known sources of information
    ## We need to pass some genus/species specific information to Prokka for this to work correctly
    ## so we need to detect the organism pipeline here
    if [ "${spp_pipeline}" == "Lmonocytogenes" ]; then
        genus="Listeria"
        species="monocytogenes"
    fi
    
    echo -e "Beginning annotation with prokka..."
    
    $prokka --genus $genus --species $species --usegenus --addgenes --cpus $threads --prefix "${sample_name}_prokka" ${best_assembly}
    
}


##########
## Main ##
##########
while [[ "${1+defined}"  ]]; do
    case "$1" in
		-1 | --i1)
			forward=$2
			shift 2
			;;
	    -2 | --i2)
	    	reverse=$2
	    	shift 2
	    	;;
        -h | --help)
            display_help
            exit 0
            ;;
	    -o | --output)
	        output_dir=$2
	        shift 2
	        ;;
	    -s | --sample_name)
	        sample_name=$2
	        shift 2
	        ;;
	    -spp | --species)
	        spp_pipeline=$2
	        shift 2
	        ;;
	    -t | --threads)
	    	threads=$2
	    	shift 2
	    	;;
	   	-td | --temp_dir)
	   		temp_dir=$2
	   		shift 2
	   		;;
        --) #End of options
            shift 1
            break
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
        *) #No more options
            break
            ;;
    esac
done

## Make temp_dir, since it will be deleted later (careful with modifying deletion code in shell scripts)
## If you aren't careful with explicit passing of variables, you can end up deleting unintended things
if [ ! -d "$temp_dir" ]; then
	mkdir "$temp_dir"
fi

if [ ! -d "$output_dir" ]; then
    mkdir "$output_dir"
fi

## Make sure everything is set up properly
validate_paths
get_versions
validate_inputs $forward $reverse


##############
## Assembly ##
##############
## Find average insert size with BBMerge from BBMap
bbmap_insert_size $forward $reverse

## Assemble using MetAmos: velvet, spades, idba, abyss
imetamos_run

## Merge contigs with CISA if applicable
if [ "${best_assembly}" == "" ]; then
    cisa_run
fi

## Annotate the best assembly, either CISA or the single best from MetAMOS
prokka_annotate



















exit 0

