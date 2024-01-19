# CNVcaller
`CNVcaller` is designed to detect copy number variation using sequencing data from populations.

## Overview
`CNVcaller` is a program for detecting the integrated copy number veriation regions (CNVRs) using population sequencing data. The high-confidence CNVRs are discovered and refined by both individual and population criteria. The result is a VCF format genotype file which can be used in GWAS/QLT research. Based on our validation, CNVcaller can report CNVRs from large populations with more than 1000 individuals within one week on one computational node. It can be applied to complicated genomes such as wheat and pan-genome.

## Installation and requirements

### Requirements
The following software must be installed on your machine:
* Perl5+ : tested with version 5.10.1
* samtools : tested with version 1.3 (using htslib 1.3)
* Python : tested with version 3.6
* blasr : tested with version 5.2+ (Optinal, in order to generate your own duplicated window record file)

### Installation
Clone the `CNVcaller` git:</br>
    `git clone https://github.com/JiangYuLab/CNVcaller.git`</br>
Go to `CNVcaller` directory</br>
    `cd CNVcaller`

### Install blasr
Since blasr has stopped updating, it is highly recommended to create a separate conda virtual environment for installation:
````
conda create -n blasr blasr=5.3.2
conda activate blasr
````

## Test with example:
To grab sample data and test `CNVcaller`, please download it from [here](http://animal.nwsuaf.edu.cn/code/source/download/CNVcaller/demo.tar.gz).

## Running the program
CNVcaller contains four steps consisting one perl script, two bash scripts and one python script. You need to set `CNVcaller` variables in the three bash scripts based on your environment.

### 1. Indexing Reference Genome
The reference genome is segmented into overlapping sliding windows. The windows are indexed to form a reference database used in all samples. This commend will create the file `referenceDB.windowsize` in current directory by default.
````
$ perl CNVReferenceDB.pl <ref>
Required arguments
 <ref> Reference sequence
Optional arguments
 -w the window size (bp) for all samples [default=800]
 -l the lower limit of GC content [default=0.2]
 -u the upper limit of GC content [default=0.7]
 -g the upper limit of gap content [default=0.5]
````
* Argument details</br>
 `-w` We recommend 400-1000bp window size for >10X coverage sequencing data, 1000-2000 window size for <10 X coverage sequencing data. Increasing the window size will reduce the noise at the cost of sensitivity.
 
 ### 2. Individual RD processing 
Count the reads of each window across genome from BAM file and generate a comparable read depth (RD) file of each individual. `referenceDB.windowsize` must be placed in current directory.</br>

Three default directories `RD_raw`  `RD_absolute`  `RD_normalized` will be created in current directory in order, containing the raw read depth, read depth after absolute copy number correction and the final GC corrected normalized read depth of each sample. The name of the normalized RD file indicates the average RD (mean), STDEV of the RD and the gender (1=XX/ZZ, 2=XY/ZW) of this sample. The final read depths are normalized to one.</br>

This step consumes about 500 MB for each individual, multiple tasks can be run in parallel. Shell script `Individual.Process.sh` is provided to complete these procedures.</br>
````
$ bash Individual.Process.sh -b <bam> -h <header> -d <dup> -s <sex_chromosome>
Required arguments 
 -b|--bam      alignment file in BAM format
 -h|--header   header of bam file, the prefix of output file
 -d|--dup      duplicated window record file used for absolute copy number correction
 -s|--sex      the name of sex chromosome
````
* Argument details</br>
 `-dup` The duplicated window record files. We provide duplicated window record files for different species, such as [human](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/human.tar.gz), [goat](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/goat.tar.gz), [sheep](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/sheep.tar.gz), [pig](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/pig.tar.gz), [cattle](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/cattle.tar.gz), [chicken](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/chicken.tar.gz), [maize](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/maize.tar.gz), [wheat](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/wheat.tar.gz), and [soybean](http://animal.omics.pro/code/source/download/CNVcaller/database/duplicated_window_record_file/soybean.tar.gz). If you work with other organisms, you will want to create duplicated window record file in order to use absolute copy number correction function of CNVcaller. Follow the [instruction](https://github.com/JiangYuLab/CNVcaller/tree/master#generate-your-own-duplicated-window-record-file).</br>
 `-s` The gender of this individual will be determines by the ratio of RD of the given sex chromosome and the RD of the other autosomes. The name of X or Z chromosome should be given for the XY or ZW genomes.</br>

* Example, to convert ERR340328.bam to normalized copy number using 1000bp window size.</br>
  `bash Individual.Process.sh -b ERR340328.bam -h ERR340328 -d dupfile -s X `

### 3. CNVR detection
The normalized RD files of all samples are piled up into a two-dimensional population RD file. The integrated CNVR are detected by scanning the population RD file with aberrantly RD, CNV allele frequency and significantly correlation with adjacent windows. The adjacent candidate windows showing high correlation will be further merged.
````
$ bash CNV.Discovery.sh -l <RDFileList> -e <excludedFileList> -f <frequency> -h <homozygous> -r <pearsonCorrelation> -p <primaryCNVR> -m <mergedCNVR>
Required arguments: 
-l|--RDFileList          individual normalized read depth file list
-e|--excludedFileList    list of samples exclude from CNVR detection
-f|--frequency           minimum frequency of gain/loss individuals for candidate CNV window definition
                         [recommend 0.1]
-h|--homozygous          number of homozygous gain/loss individuals for candidate CNV window definition 
                         [recommend 3]
-r|--pearsonCorrelation  minimum of Pearson’s correlation coefficient between the two adjacent non-overlapping windows
                         0.5 for sample size (0, 30]
                         0.4 for sample size (30, 50]
                         0.3 for sample size (50, 100]
                         0.2 for sample size (100, 200]
                         0.15 for sample size (200, 500]
                         0.1 for sample size (500,+∞)
-p|--primaryCNVR         primary CNVR result 
-m|--mergedCNVR          merged CNVR result
````
* Argument details</br>
 `-e`  The samples in this list will be exclude from CNVR detection, and their copy numbers are deduced based on the CNVR boundaries defined by other samples. This option is applicable to the outgroup or the poor quality precious samples. An empty file means all individuals are included in the CNVR detection.</br>
 `-f/-h` Windows satisfied any of this two conditions will be selected as candidate CNV windows.</br>
 `-r` The adjacent windows with significant correlation will be merged in to one call. The recommend value is significant at p=0.01 level. Raise this index will increase the detection accuracy with a decrease of sensitivity.</br>

* Example, run `CNV.Discovery.sh` on all your individual normalized RD files for discovering CNV.</br>
 An example of normalized read depth file list `-l`:
````
  RD_normalized/ERR340328_mean_70.81_SD_10.84_sex_1
  RD_normalized/ERR340329_mean_62.00_SD_10.52_sex_1
  RD_normalized/ERR340330_mean_135.66_SD_13.96_sex_1
  RD_normalized/ERR340331_mean_128.76_SD_15.27_sex_1
  RD_normalized/ERR340333_mean_69.30_SD_10.19_sex_1
  RD_normalized/ERR340334_mean_132.30_SD_14.59_sex_1
  RD_normalized/ERR340335_mean_73.50_SD_10.16_sex_1
  RD_normalized/ERR340336_mean_72.52_SD_10.03_sex_1
  RD_normalized/ERR340338_mean_124.12_SD_13.24_sex_1
  RD_normalized/ERR340340_mean_131.00_SD_14.74_sex_1
````
  `bash CNV.Discovery.sh -l list -e exclude_list -f 0.1 -h 3 -r 0.5 -p primaryCNVR -m mergeCNVR`

### 4. Genotyping
Clustering the input samples into genotypes uses Gaussian mixture modes. The output contain a genotype VCF - a VCF format file containing the input site descriptions, additional site-specific information and a called genotype for each input sample.
````
$ python Genotype.py --cnvfile <input> --outprefix <outfile prefix>
Required arguments:
--cnvfile          merged CNVR file
--outprefix        prefix of out files
Optional arguments:
--nproc            number of process will be used, default is one.
````
* Example</br>
 `python Genotype.py --cnvfile mergeCNVR --outprefix Genotype`

## Generate your own duplicated window record file
It is highly recommended to generate the corresponding file for each chromosome separately and then concatenate the output files in the correct order.

Step 1: Split genome into short kmer sequences.
````
$ python 0.1.Kmer_Generate.py [OPTIONS] FAFILE WINSIZE OUTFILE

 Required arguments
 <FAFILE> Reference sequence in FASTA format
 <WINSIZE> The size of the window to use for CNV calling
 <OUTFILE> Output kmer file in FASTA format
````
 `Example: python 0.1.Kmer_Generate.py reference.fa 800 kmer.fa`

Step 2: Align the kmer FASTA (from step 1) to reference genome using blasr sawriter and in the conda `blasr` environment.

Create the reference.fa.sa file: 
`Example: sawriter reference.fasta.sa reference.fasta`

Align the kmer FASTA to reference genome: 
`Example: blasr kmer.fa reference.fa --sa reference.fa.sa --out kmer.aln -m 5 --noSplitSubreads --minMatch 15 --maxMatch 20 --advanceHalf --advanceExactMatches 10 --fastMaxInterval --fastSDP --aggressiveIntervalCut --bestn 10`

Step 3: Generate duplicated window record file.
````
$ python 0.2.Kmer_Link.py [OPTIONS] BLASR WINSIZE OUTFILE

 Required arguments
 <BLASR> blasr results (-m 5 format)
 <WINSIZE> The size of the window to use for CNV calling
 <OUTFILE> Output genome duplicated window record file
````
 `Example: python 0.2.Kmer_Link.py kmer.aln 800 window.link`
 
## Contact 
Any questions, bug reports and suggestions can be posted to Email:</br>
yu.jiang@nwafu.edu.cn

## CNVcaller has a paper
https://academic.oup.com/gigascience/advance-article/doi/10.1093/gigascience/gix115/4689116

## Citation
Xihong Wang, Zhuqing Zheng, Yudong Cai, Ting Chen, Chao Li, Weiwei Fu, Yu Jiang; CNVcaller: Highly Efficient and Widely Applicable Software for Detecting Copy Number Variations in Large Populations, GigaScience, , gix115, https://doi.org/10.1093/gigascience/gix115
