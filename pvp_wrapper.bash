#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
condapath=$DIR/.mc/bin
export PATH=$condapath:$PATH
force=false
help=false
threads=1


TEMP=`getopt -o q:fr:c:ho:t:: --long fastq-list,force-overwrite,reference:,host-reference:,help,output-directory:,threads:: -n 'bact_vcp' -- "$@"`
eval set -- "$TEMP"

while true ; do
    case "$1" in
        -q | --fastq-list ) fastqs="$2"; shift 2 ;;
        -b | --bloomfilter ) bf="$2"; shift 2 ;;
        -f | --force-overwrite ) force=true; shift ;;
        -r | --reference ) ref="$2"; shift 2 ;;
        -c | --host-reference ) host_ref="$2"; shift 2 ;;
        -h | --help ) help=true; shift ;;
        -o | --output-directory ) out="$2"; shift 2 ;;
        -t | --threads ) 
            case "$2" in
                "") threads=1; shift 2 ;;
                *) threads="$2"; shift 2 ;;
            esac ;;
        -- ) shift; break;;
        * ) exit 1 ;;
    esac
done

if [ fastqs == "" ]; then
    echo "Required option: -q | --fastq-list not given."
	exit 1
fi

if [ ref == "" ]; then
    echo "Required option: -r | --reference not given."
	exit 1
fi

if [ host_ref == "" ]; then
    echo "Required option: -c | --cattle host reference not given."
	exit 1
fi

if [ out == "" ]; then
    echo "Required option: -o | --output-directory not given."
	exit 1
fi


if [ $help == true ]; then
    echo "Another bacterial variant calling pipeline:   
-q | --fastq-list    A comma separated list containing fastq file paths.

-b | --bloomfilter    Path to bloom filter file (.bf) of the host genome.

-f | --force-overwrite    Overwrite files in the output directories.

-r | --reference    Path to a reference fasta file

-c | --host-reference    Path to cattle host reference fasta file

-h | --help    Prints this help output.

-o | --output-directory    Path to the output directory. A directory will be created if one does not exist.

-t | --threads Number of threads to use for multiprocessing. Must be written without space between flag and value. i.e. for 4 threads argument is -t4 NOT -t 4. Defaults to 1.

"
exit 0
fi

$DIR/pvp.py $fastqs $ref $out $force $threads $host_ref