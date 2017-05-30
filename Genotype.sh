#!/bin/sh
# CNVcaller installation directory
export CNVcaller=`pwd`
echo "CNVcaller install directory $CNVcaller"
if [ ! -f "$CNVcaller/bin/2.1.CNVDiscoveryMerge.pl" ]; then echo "You should set CNVcaller installation directory"; exit 1; fi
##############
usage() {
  echo "Usage: $0 -i <input> -n <sampleNumber> -o <output>"
  echo "required options:"
    
    echo "-i|--input          merged CNVR file"

    echo "-n|--sampleNumber   number of samples"
    
    echo "-o|--output         vcf output"
    1>&2; exit 1;
}

##############
OPTS=`getopt -o i:n:o: --long input:,sampleNumber:,output -- "$@"`
if [ $? != 0 ]; then usage;  exit 1; fi
eval set -- "$OPTS"
while true ; do
  case "$1" in
    -i|--input )
        input=$2 ; shift 2 ;;
    -n|--sampleNumber )
        sampleNumber=$2 ; shift 2 ;;
    -o|--output )
        output=$2 ; shift 2 ;;
    --) shift ; break ;;
    *) echo "$1 Internal error!" ; exit 1 ;;
  esac
done
if [ -z "$input" ] || [ -z "$sampleNumber" ] || [ -z "$output" ]; then usage; exit 1; fi
echo "merge CNVR $input";
echo "sampleNumber $sampleNumber";
echo "output $output";

perl $CNVcaller/bin/CNVGenotyping.pl $input $sampleNumber >Genotype_Constrained
Rscript $CNVcaller/bin/mclust.R $input
perl $CNVcaller/bin/Genotype_to_VCF.pl -gc Genotype_Constrained -gu Genotype_Unsupervised -n $sampleNumber -o $output
echo "CNVR genotyping finished!";
