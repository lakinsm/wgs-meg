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
        -i | --input                A DIRECTORY with ONLY the files to analyze, MUST have _R1_ _R2_ identifiers
        -a | --assembly             Flag: run assembly pipeline on this sample (slower)
        -o | --output DIR           Directory for output of important files (main output dir)
        -spp| --species STR         Species identifier string (see documentation for details)
        -t | --threads  INT         Threads to use where applicable
        -td| --temp_dir DIR         Temporary directory for intermediate dirs/files
        -T | --threshold INT        Quality score threshold for N-masking in kSNP pipeline

"
}


###############
## Variables ##
###############
## Paths to output directories
temp_dir=""
output_dir=""

## Flags and variables used in the pipeline
spp_pipeline="Lmonocytogenes"
run_assembly=0
treshold=25  # quality score threshold for the N_masking step in kSNP
threads=1  # threads to use where applicable, recommend ~ 20


##########
## Main ##
##########
while [[ "${1+defined}"  ]]; do
    case "$1" in
		-i | --input)
			indir=$2
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
	    -o | --output)
	        output_dir=$2
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

declare -a files
files=($indir/*)
pos=$(( ${#files[*]} - 1 ))

for i in `seq 0 $pos`; do
	echo "${files[$i]}"
done





















exit 0

