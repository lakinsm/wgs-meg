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

## Paths to output directories
temp_dir=""
output_dir=""

## Flags and variables used in the pipeline
sample_name=""
insert=0  # insert size as determined by bbmerge
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
    
    Paths currently set:
    BBmap bbmerge:            ${bbmap}
    iMetAmos:                 ${metamos}

    
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
    
    #--------------------------------
    
    if [ "${sample_name}" == "" ]; then
        local missing="$missing:SampleName;"
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
	else
		echo -e "Begin WGS Pipeline for $filename1 and $filename2..."
		echo -e "\nBegin WGS Pipeline for $filename1 and $filename2..." >> WGS_LabNotebook.txt
	fi
}


get_versions() {
	echo -e "Software versions:" >> WGS_LabNotebook.txt
	$bbmap | grep "^BBMerge" >> WGS_LabNotebook.txt
	head -n 1 ${metamos}/README.md | grep -Po "MetAMOS.*\" " >> WGS_LabNotebook.txt
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
    
    "${metamos}"/initPipeline -q -1 $forward -2 $reverse -d "${temp_dir}/${sample_name}" -i $insert -W iMetAMOS
    kmer=$( "${metamos}"/runPipeline -p $threads -t eautils -q -a velvet,spades,idba-ud,abyss -b -z genus -d "${temp_dir}/${sample_name}" | grep "Selected kmer size" | grep -Po "[0-9]{1,4}" )
    echo "Kmer size is: $kmer"
    
    #if [ ! -f "${metamos}/${sample_name}"/Assemble/out/idba-ud
    
    
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

## 



















exit 0

