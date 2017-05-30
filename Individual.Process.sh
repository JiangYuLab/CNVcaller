#!/bin/sh
# Samtools executables must be on your path.
which samtools > /dev/null || exit 1
# CNVcaller is the installation directory for CNVcaller - it must be an exported environment variable
export CNVcaller=`pwd`
echo "CNVcaller install directory $CNVcaller"
if [ ! -f "$CNVcaller/bin/1.1.CNVprocess.pl" ]; then echo "You should set CNVcaller installation directory"; exit 1; fi
##############
usage() {
  echo "Usage: $0 -b <bam> -h <header> -d <dup> -s <sex_chromosome>"
  echo "required options:"
    echo "-b|--bam           alignment file in BAM format"
    echo "-h|--header        header of bam file, the prefix of output file [same with SM tag of input BAM file]"
    echo "-d|--dup           duplicated window record file used for absolute copy number correction"
    echo "-s|--sex           the name of sex chromosome"
  1>&2; exit 1;
}

##############
OPTS=`getopt -o b:h:d:s: --long bam:,header:,dup:,sex: -- "$@"`
if [ $? != 0 ]; then usage;  exit 1; fi
eval set -- "$OPTS"
while true ; do
  case "$1" in
    -b|--bam ) 
        bam=$2 ; shift 2 ;;
    -h|--header ) 
        header=$2 ; shift 2 ;;
    -d|--dup )
        dup=$2 ; shift 2 ;;
    -s|--sex ) 
        sex=$2 ; shift 2 ;;
    --) shift ; break ;;
    *) echo "Option error!" ; exit 1 ;;
  esac
done
if [ -z "$header" ] || [ -z "$dup" ] || [ -z "$bam" ] || [ -z "$sex" ] ; then usage; exit 1; fi
refdb=$(ls `pwd`/referenceDB*)
windowsize=$(echo $refdb | grep -oP "referenceDB.\d+" | grep -oP "\d+")
echo "refdb $refdb";
echo "header $header";
echo "duplicate file $dup";
echo "bam $bam";
echo "sex $sex";
echo "windowsize $windowsize";
mkdir -p RD_raw
cd RD_raw
perl $CNVcaller/bin/1.1.CNVprocess.pl $refdb $bam -k1 $windowsize
echo "raw reads count finished!"
cd ..
mkdir -p RD_absolute
perl $CNVcaller/bin/1.2.CNVcorrect.pl RD_raw/$header\_raw $dup RD_absolute/$header
echo "absolute correct finished!"
mkdir -p RD_normalized
cd RD_normalized
perl $CNVcaller/bin/1.3.CNVnormalize.pl ../RD_absolute/$header\_link -k1 $windowsize -k2 $sex
cd ..
echo "normalization finished!"
