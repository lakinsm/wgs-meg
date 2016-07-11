#usr/bin/env bash

if [ ! "$BASH_VERSION" ] ; then
    echo "Please do not use sh to run this script ($0), just execute it directly." 1>&2
    exit 1
fi

RELPATH="${BASH_SOURCE%/*}"

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
        -a | --assembly             Flag: run assembly pipeline on this sample (slower)
        -l | --last                 Flag: is this the last sample? If so, run kSNP
        -o | --output DIR           Directory for output of important files (main output dir)
        -s | --sample_name  STR     Main file name for this sample throughout the pipeline
        -spp| --species STR         Species identifier string (see documentation for details)
        -t | --threads  INT         Threads to use where applicable
        -td| --temp_dir DIR         Temporary directory for intermediate dirs/files
        -T | --threshold INT        Quality score threshold for N-masking in kSNP pipeline

"
}


###############
## Variables ##
###############
## Paths to program executables or jar files where noted
bbmap="/usr/local/bin/bbmerge.sh"
metamos="/home/chris_dean/software/metAMOS-1.5rc3"
nucmer="/usr/bin/nucmer"
cisa="/home/chris_dean/software/CISA1.3"
blastdb="/usr/local/bin/makeblastdb"
blastn="/usr/local/bin/blastn"
prokka="/usr/local/bin/prokka"
trimmomatic="/usr/local/Trimmomatic-0.36"
bwa="/usr/local/bin/bwa"
amrdb="/home/chris_dean/sequence_data/amrdb/mmarc.fasta"
vfdb="/home/chris_dean/sequence_data/vfdb/VFDB.fa"
plasmiddb="/home/chris_dean/sequence_data/plasmid/plasmid_seqs.fasta"
csa="/usr/local/bin/csa"
samtools="/usr/local/bin/samtools"
freebayes="/usr/local/bin/freebayes"
bcftools="/usr/local/bin/bcftools"
vcfutils="/usr/local/bin/vcfutils.pl"
PERL5LIB="/home/chris_dean/software/vcftools-vcftools-4491144/src/perl"
tabix="/home/chris_dean/software/tabix-0.2.6"
ksnp="/usr/local/kSNP3/kSNP3"
ksnp_repo="/home/chris_dean/sequence_data/ksnp_repository"
integration_limbo="/home/chris_dean/sequence_data/ksnp_repository/awaiting_integration"
Lmono="/home/chris_dean/sequence_data/list/Listeria.fa"

## Paths to output directories
temp_dir="/home/chris_dean/wgs_pipeline/wgs-meg/temp_dir"
output_dir="/home/chris_dean/wgs_pipeline/wgs-meg/output_dir"

## Flags and variables used in the pipeline
sample_name="Listeria"
spp_pipeline="Lmonocytogenes"
run_assembly=0
best_assembly=""  # used in the event of less than 3 assemblies from metamos
consensus_file=""  # used to verify that a proper consensus was created
insert=0  # insert size as determined by bbmerge
kmer=0  # k-mer size determined by specialk in metamos
genome_size=0  # genome size for this particular organism
threshold=25  # quality score threshold for the N_masking step in kSNP
chosen_k=0  # k-value chosen by kchooser
last_sample=0  # is this the last sample in the pipline?
threads=1  # threads to use where applicable, recommend ~ 20

export PERL5LIB="/home/chris_dean/software/vcftools-vcftools-4491144/src/perl"
#export PATH="${PATH}:/s/angus/index/common/tools/tabix"

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
    Trimmomatic:              ${trimmomatic}
    BWA:                      ${bwa}
    SAMTools:                 ${samtools}
    BCFTools:                 ${bcftools}
    VCFUtils:                 ${vcfutils}
    PERL5LIB:                 ${PERL5LIB}
    Tabix:                    ${tabix}
    kSNP3:                    ${ksnp}
    
    AMR Database:             ${amrdb}
    Virulence Database:       ${vfdb}
    Plasmid Database:         ${plasmiddb}

    
    Directories selected:
    Temporary directory:      ${temp_dir}
    Output directory:         ${output_dir}
    kSNP repository:          ${ksnp_repo}
    kSNP integration dir:     ${integration_limbo}
   	
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
    
    if [ ! -e "${trimmomatic}" ]; then
    	local missing="$missing:Trimmomatic;"
    fi
    
    if [ ! -e "${bwa}" ]; then
    	local missing="$missing:BWA;"
    fi

    if [ ! -e "${samtools}" ]; then
    	local missing="$missing:SAMTools;"
    fi
    
    if [ ! -e "${bcftools}" ]; then
    	local missing="$missing:BCFTools;"
    fi
    
#    if [ ! -e "${vcfutils}" ]; then
#    	local missing="$missing:VCFUtils;"
#    fi
    
    if [ ! -d "${PERL5LIB}" ]; then
    	local missing="$missing:PERL5LIB;"
    fi
    
#    if [ ! -d "${tabix}" ]; then
#    	local missing="$missing:TabixPath;"
#    fi
    
    if [ ! -e "${ksnp}" ]; then
    	local missing="$missing:kSNP3;"
    fi
    
    if [ ! -e "${amrdb}" ]; then
    	local missing="$missing:AMRDatabase;"
    fi
    
    if [ ! -e "${vfdb}" ]; then
    	local missing="$missing:VirulenceDatabase;"
    fi
    
    if [ ! -e "${plasmiddb}" ]; then
    	local missing="$missing:PlasmidDatabase;"
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

#    if [ ! -d "${ksnp_repo}" ]; then
#	local missing="$missing:kSNPRepository"
#    fi
#
#    if [ ! -d "${integration_limbo}" ]; then
#        local missing="$missing:kSNPLimboDirectory"
#    fi
    
    ## Check if components are missing
    if [ ! -e $missing ]; then
        echo -e "Directories or Paths are missing, check the following:\n"
        echo -e $missing | sed 's/;/\n/g'
        exit 1
    else
        bgzip="${tabix}/bgzip"
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
	$samtools 3>&1 1>&2 2>&3 3>&- | grep "Version" | sed 's/^/SamTools /' >> WGS_LabNotebook.txt
}


bbmap_insert_size() {
	## We need to capture stderr from bbmerge to regex capture the insert range.
	## Below is the way to do that by juggling output streams: 3>&1 1>&2 2>&3 3>&-
	insert=$( $bbmap in1="$1" in2="$2" out=/dev/null 3>&1 1>&2 2>&3 3>&- | grep "^Insert range" | grep -Po "[0-9]{1,4} - [0-9]{1,4}" | sed 's/ - /:/' )
	
	## Check to make sure we succeeded
	if [ "$insert" == "0" ]; then
	    echo -e "Error: insert size detection with BBMerge failed, insert size is still zero."
	    echo -e "\tError: insert size detection with BBMerge failed, insert size is still zero." >> WGS_LabNotebook.txt
	    exit 1
	fi
}


imetamos_run() {
    ## Run the iMetAMOS pipeline for assembly.  This produces four files that we need, one for each assembler
    ## We then record the location of these and general other information for input to CISA later.
    ## We will create BOTH CISA input files within this function, then check to see they are of valid
    ## format to run with CISA.
    ## Note that this function utilizes global variables.

    echo -e "Beginning MetAMOS assembly using idba-ud, velvet, spades, abyss, and edena..."
    echo -e "\tBeginning MetAMOS assembly using idba-ud, velvet, spades, abyss, and edena..." >> WGS_LabNotebook.txt
    "${metamos}"/initPipeline -q -1 $forward -2 $reverse -d "${temp_dir}/${sample_name}" -i $insert -W iMetAMOS
    kmer=$( "${metamos}"/runPipeline -p $threads -t eautils -q -a velvet,spades,idba-ud,abyss,edena -b -z genus -d "${temp_dir}/${sample_name}" | grep "Selected kmer size" | grep -Po "[0-9]{1,4}" )
    
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
    
    ## This note is for when we run into bugs with bad assemblies (and we will at some point):
    ## Make sure this section is assigning variables correctly for best assembly
    if [ "$valid_assemblies" == "0" ]; then
        echo -e "Error: Zero valid assemblies detected.  Check error logs."
        echo -e "\tError: Zero valid assemblies detected.  Check error logs." >> WGS_LabNotebook.txt
    elif [ "$valid_assemblies" -lt "3" ]; then
        echo -e "Less than three valid assemblies detected, skipping CISA and continuing with annotation..."
        echo $which_assemblies | sed 's/;/\n/g'
        echo -e "\tLess than three valid assemblies detected, skipping CISA and continuing with annotation..." >> WGS_LabNotebook.txt
        echo $which_assemblies | sed 's/;/\n/g' >> WGS_LabNotebook.txt
        best_assembly=$( head -n 1 "${temp_dir}/${sample_name}/Postprocess/out/best.asm" | sed 's/\n//' | sed 's/\..*//' | sed 's/-ud$//' )
        best_assembly="${output_dir}/${best_assembly}_assembly.contig"
    else
        ## Proceed with CISA config file creation
        
        echo -e "Creating CISA config files for ${spp_pipeline} organism..."
        echo -e "\tCreating CISA config files for ${spp_pipeline} organism..." >> WGS_LabNotebook.txt
        
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
    
    assembler=""
    valid_assemblies=""
    which_assemblies=()
    rm flowchart.svg
}


cisa_run() {
    ## CISA takes the assemblies from iMetAMOS and integrates them into a single file
    ## It then takes that file and computes a combined set of contigs if possible
    echo -e "Running CISA for contig integration..."
    echo -e "\tRunning CISA for contig integration..." >> WGS_LabNotebook.txt
    
    $cisa/Merge.py "${temp_dir}/Merge.config"
    echo "y" | $cisa/CISA.py "${temp_dir}/CISA.config"
    best_assembly="${output_dir}/${sample_name}_cisa_integrated.fa"
    
    if [ ! -f "$best_assembly" ]; then
        echo -e "CISA failed to produce the expected merged file.  Please check the logs."
        echo -e "\tCISA failed to produce the expected merged file.  Please check the logs." >> WGS_LabNotebook.txt
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
    echo -e "\tBeginning annotation with prokka..." >> WGS_LabNotebook.txt
    
    $prokka --genus $genus --species $species --usegenus --addgenes --cpus $threads --prefix "${sample_name}_prokka_${1}" $2
    
}


trim_reads() {
    if [ ! -f "${temp_dir}/1p.fastq" ] || [ ! -f "${temp_dir}/2p.fastq" ]; then
        echo -e "Beginning read trimming..."
        echo -e "\tBeginning read trimming..." >> WGS_LabNotebook.txt
        java -jar "${trimmomatic}/trimmomatic-0.36.jar" PE -threads $threads -phred33 $forward $reverse "${temp_dir}/1p.fastq" "${temp_dir}/1u.fastq" "${temp_dir}/2p.fastq" "${temp_dir}/2u.fastq" ILLUMINACLIP:"${trimmomatic}/adapters/TruSeq3-PE.fa:2:30:10:3:TRUE" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    else
        echo -e "Trimmed reads already detected, proceeding with alignments..."
        echo -e "\tTrimmed reads already detected, proceeding with alignments..." >> WGS_LabNotebook.txt
    fi
}


amr_align() {
    if [ ! -f "${output_dir}/${sample_name}_amr_parsed.tabular" ]; then
        echo -e "Aligning to AMR database..."
        echo -e "\tAligning to AMR database..." >> WGS_LabNotebook.txt
        $bwa mem -t $threads ${amrdb} ${temp_dir}/1p.fastq ${temp_dir}/2p.fastq > ${temp_dir}/amr.sam
        $csa -ref_fp ${amrdb} -sam_fp ${temp_dir}/amr.sam -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp "${output_dir}/${sample_name}_amr_parsed.tabular"
        
        rm ${temp_dir}/amr.sam
    else
        echo -e "AMR parsed file detected, proceeding..."
        echo -e "\tAMR parsed file detected, proceeding..." >> WGS_LabNotebook.txt
    fi
}


vfdb_align() {
    if [ ! -f "${output_dir}/${sample_name}_vfdb_parsed.tabular" ]; then
        echo -e "Aligning to Virulence database..."
        echo -e "\tAligning to Virulence database..." >> WGS_LabNotebook.txt
        $bwa mem -t $threads ${vfdb} ${temp_dir}/1p.fastq ${temp_dir}/2p.fastq > ${temp_dir}/vfdb.sam
        $csa -ref_fp ${vfdb} -sam_fp ${temp_dir}/vfdb.sam -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp "${output_dir}/${sample_name}_vfdb_parsed.tabular"
        
        rm ${temp_dir}/vfdb.sam
    else
        echo -e "VFDB parsed file detected, proceeding..."
        echo -e "\tVFDB parsed file detected, proceeding..." >> WGS_LabNotebook.txt
    fi
}


plasmid_align() {
    if [ ! -f "${output_dir}/${sample_name}_plasmid_parsed.tabular" ]; then
        echo -e "Aligning to Plasmid database..."
        echo -e "\tAligning to Plasmid database..." >> WGS_LabNotebook.txt
        $bwa mem -t $threads ${temp_dir}/1p.fastq ${temp_dir}/2p.fastq > ${temp_dir}/plasmid.sam
        $csa -ref_fp ${plasmiddb} -sam_fp ${temp_dir}/plasmid.sam -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp "${output_dir}/${sample_name}_plasmid_parsed.tabular"
        
        rm ${temp_dir}/plasmid.sam
    else
        echo -e "Plasmid parsed file detected, proceeding..."
        echo -e "\tPlasmid parsed file detected, proceeding..." >> WGS_LabNotebook.txt
    fi
}


choose_reference() {
    ## Pick the appropriate reference based on the input species
    if [ "${spp_pipeline}" == "Lmonocytogenes" ]; then
        refgenome="${Lmono}"
    fi
}


ref_align() {
    ## Align to the reference genome chosen in the previous steps
    ## We will use BWA mem here, since whole genome alignment is what it was designed to do
    $bwa mem -t $threads ${refgenome} ${temp_dir}/1p.fastq ${temp_dir}/2p.fastq > ${temp_dir}/ref.sam
    $samtools view -hbS ${temp_dir}/ref.sam > ${temp_dir}/ref.bam
    $samtools sort ${temp_dir}/ref.bam -o ${output_dir}/${sample_name}_ref_sorted.bam
    $samtools index ${output_dir}/${sample_name}_ref_sorted.bam
    
    ## We also want to keep reads that are unmapped to reference
    ## To do this, we can use the filter_sam.py script with a filter threshold of 70bp read length
    ## since we don't want to keep short read fragments that are unlikely to be useful in BLAST
    cat ${temp_dir}/ref.sam | python3 ${RELPATH}/filter_sam.py - -l 70 > ${output_dir}/${sample_name}_reads_unmapped_to_ref.fasta
    
    ## Now we mpileup to create a consensus
    #$samtools mpileup -uD -f ${refgenome} ${output_dir}/${sample_name}_ref_sorted.bam | $bcftools view - > ${temp_dir}/ref_raw.bcf
    $freebayes -p 1 -f ${refgenome} ${output_dir}/${sample_name}_ref_sorted.bam | $bcftools view - > ${temp_dir}/ref_raw.bcf    
    $bcftools view ${temp_dir}/ref_raw.bcf | $bgzip -c > ${output_dir}/${sample_name}_ref_snps.vcf.gz
    ${tabix}/tabix -p vcf ${output_dir}/${sample_name}_ref_snps.vcf.gz
    cat ${refgenome} | ${PERL5LIB}/vcf-consensus ${output_dir}/${sample_name}_ref_snps.vcf.gz > ${output_dir}/${sample_name}_consensus.fa
    consensus_file="${output_dir}/${sample_name}_consensus.fa"
}


n_mask() {
    ## For KChooser, we need to first mask low-quality bases with N and then concatenate the resulting files
    ## into a single FASTA
    python3 ${RELPATH}/dna_threshold.py ${temp_dir}/1p.fastq $threshold > ${temp_dir}/masked_reads.fasta
    python3 ${RELPATH}/dna_threshold.py ${temp_dir}/2p.fastq $threshold >> ${temp_dir}/masked_reads.fasta
}


ksnp_build() {
	## kSNP3 generates a distance matrix and Phylip tree output
	## We will store those outputs as well as turn the Phylip into strict Phylip for further use
	if [ ! -d ${sample_name}_ksnp ]; then
	    mkdir ${sample_name}_ksnp
	fi
	
	## Need to make KSNP input file here, but we need to agree on our method for this first
	if [ ! -d "$ksnp_repo" ]; then
	    mkdir "$ksnp_repo"
        fi
	if [ ! -d "$integration_limbo" ]; then
	    mkdir "$integration_limbo"
	fi
        ## The above ksnp_repo is hardcoded into this script.  This is intentional, as that directory
	## should be a read-only repository of the fasta files we have collected and n_masked to date.
	## Likewise, the "integration_repo" is a location where, after the pipeline has been run, new
	## genomes are concatenated to the old list and await approval before final integration into
	## the master genome list.  This pipeline will, for every run, copy THE MASTER FILE ONLY,
	## and it will append all of the current genomes being run and output that file to limbo.
	## Integration into the master file must be done using a separate script or manually.
	## DO THIS BETWEEN PIPELINE RUNS; the above precautions are to prevent destruction of all
	## previous work.
	mv "${temp_dir}/masked_reads.fasta" "${ksnp_repo}/${sample_name}.fasta"
	fullpath=$( readlink -f "${ksnp_repo}/${sample_name}.fasta" )
	if [ -f "${integration_limbo}/current_ksnp_run.infile" ]; then
	    echo -e "${fullpath}\t${sample_name}" >> "${integration_limbo}/current_ksnp_run.infile"
	elif [ -f "${ksnp_repo}/master_paths.infile" ]; then
	    cat "${ksnp_rpo}/master_paths.infile" > "${integration_limbo}/current_ksnp_run.infile"
	    echo -e "${fullpath}\t${sample_name}" >> "${integration_limbo}/current_ksnp_run.infile"
	else
	    echo -e "${fullpath}\t${sample_name}" >> "${integration_limbo}/current_ksnp_run.infile"
	fi
	
	## If this is the last sample, then we run ksnp
	if [ "$last_sample" == "1" ]; then
	    chosen_k=19  # this is temporary while we debug kchooser
	    echo -e "${refgenome}\tListeria" >> "${integration_limbo}/current_ksnp_run.infile"
	    $ksnp -in ${integration_limbo}/current_ksnp_run.infile -outdir ${sample_name}_ksnp -k ${chosen_k} -CPU $threads -core -ML
	fi
}


cleanup() {
    rm ${temp_dir}/1p.fastq ${temp_dir}/2p.fastq ${temp_dir}/1u.fastq ${temp_dir}/2u.fastq ${temp_dir}/amr.sam
    echo "Files cleaned up.  Temporary directory ${temp_dir} can optionally be deleted by the user."
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
	    -a | --assembly)
	        run_assembly=1
	        shift 1
	        ;;
        -h | --help)
            display_help
            exit 0
            ;;
	-l | --last)
	    last_sample=1
	    shift 1
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
	    -T | --threshold)
	        threshold=$2
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
if [ "${run_assembly}" == "1" ]; then
    echo -e "Beginning Assembly pipeline..."
    echo -e "Beginning Assembly pipeline..." >> WGS_LabNotebook.txt
    
    ## Find average insert size with BBMerge from BBMap
    bbmap_insert_size $forward $reverse

    ## Assemble using MetAmos: velvet, spades, idba, abyss
    imetamos_run

    ## Merge contigs with CISA if applicable
    if [ "${best_assembly}" == "" ]; then
        cisa_run
    fi
    
    ## Truncate headers of best assembly for prokka
    mv "${best_assembly}" "${best_assembly}_temp"
    python3 ${RELPATH}/truncate_headers.py "${best_assembly}_temp" > "${best_assembly}"
    rm "${best_assembly}_temp"

    ## Annotate the best assembly, either CISA or the single best from MetAMOS
    prokka_annotate "assembly" "${best_assembly}"
    
    best_assembly=""

    echo -e "Assembly pipeline complete."
    echo -e "Assembly pipeline complete." >> WGS_LabNotebook.txt
fi


###############
## Alignment ##
###############
echo -e "Beginning Alignment pipeline..."
echo -e "Beginning Alignment pipeline..." >> WGS_LabNotebook.txt
## Preprocess the reads
trim_reads

## Align to the AMR database, VFDB, and Plasmid database
amr_align
vfdb_align
plasmid_align

## Pick the correct reference genome
choose_reference

## Align to reference genome
ref_align

# Annotate the consensus file created in the reference alignment step
if [ "${consensus_file}" == "" ]; then
    echo -e "Consensus file creation failed; check the logs."
    echo -e "\tConsensus file creation failed; check the logs." >> WGS_LabNotebook.txt
else
    mv "${consensus_file}" "${consensus_file}_temp"
    python3 ${RELPATH}/truncate_headers.py "${consensus_file}_temp" > "${consensus_file}"
    rm "${consensus_file}_temp"
    prokka_annotate "align_consensus" "${consensus_file}"
fi

echo -e "Alignment pipeline complete."
echo -e "Alignment pipeline complete." >> WGS_LabNotebook.txt

###################
## kSNP Distance ##
###################
echo -e "Beginning kSNP pipeline..."
echo -e "Beginning kSNP pipeline..." >> WGS_LabNotebook.txt
## Perform N-masking on the trimmed reads from the Trimmomatic step
## Pass these into the kSNP pipeline
n_mask

## Run KChooser (we're still figuring out the debugging on this, to be implemented later)

## Run pre-reqs to kSNP3
ksnp_build

echo -e "kSNP prep complete."
echo -e "kSNP prep complete." >> WGS_LabNotebook.txt

cleanup

exit 0
