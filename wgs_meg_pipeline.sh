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
    echo "Please do not use sh to run this script ($0), just execute it directly:\n./wgs_meg_pipeline.sh [args]" 1>&2
    exit 1
fi

shopt -s extglob


###############
## Variables ##
###############
## Paths to program executables or jar files where noted
bbmap="/s/angus/index/common/tools/BBMap_35.85/bbmap/bbmerge.sh"

## Flags and variables used in the pipeline
insert=0  # insert size as determined by bbmerge


###############
## Functions ##
###############
## Check that the program and directory paths exist
validate_paths() {
    echo -e "\n@.@\nMEG WGS Pipeline started on `date`\n" >> WGS_LabNotebook.txt

    local missing=""

    echo "
    Paths currently set:
    BBmap bbmerge:            ${bbmap}

    
    Directories selected:
    Temporary directory:      ${temp_dir}
   	
   	" >> WGS_LabNotebook.txt
   	
   	if [ ! -e "${bbmap}" ]; then
        local missing="$missing:BBMap;"
    fi
    
    #--------------------------------
    
    if [ ! -d "${temp_dir}" ]; then
        local missing="$missing:TemporaryDirectory;"
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
		echo -e "Begin WGS Pipeline for $filename1 and $filename2..." >> WGS_LabNotebook.txt
	fi
}


get_versions() {
	echo -e "Software versions:"
	$bbmap | grep "^BBMerge" >> WGS_LabNotebook.txt
}



bbmap_insert_size() {
	## We need to capture stderr from bbmerge to regex capture the avg insert size.
	## Below is the way to do that by juggling output streams: 3>&1 1>&2 2>&3 3>&-
	insert=$( $bbmap in1="$1" in2="$2" out=/dev/null 3>&1 1>&2 2>&3 3>&- | grep "^Avg Insert" | grep -Po "[0-9]{1,4}\.[0-9]{1}" )
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

## Make sure everything is set up properly
validate_paths
get_versions
validate_inputs $forward $reverse


# Find optimal k-mer length with BBMerge from BBMap
bbmap_insert_size $forward $reverse



















exit 0

