#!/bin/bash
#$ -N fastq2sam
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 1-10

## load modules
module load java/1.8.0.111      # language of gatk & picard-tools
module load gatk/4.1.2.0        # includes picard tools
module load bwa/0.7.8           # aligner
module load samtools/1.9        # filter reads by quality score, convert sam2bam

## to run a script on hpc use `qsub path/to/filename.sh`
    ## `qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/fastq2vcf-gatk-v0.1.0.sh`
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`

######################
## create directories
######################

## before running this script
## make and move to this directory 
## /dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## make a directory for generated genomic files
mkdir ../genomic/


####################
## create variables
####################

## create absolute path2datadir
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'

## create path2fastq holding directory variable
path2fastq=${path2datadir}'raw/nguyen_nc_2018/ind1/'

# create path to deposit data
path2genomic=${path2datadir}'gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/genomic/'

## create varibles for lists of files
r1list=${path2datadir}'gen/nguyen_nc_2018/20190809-r1r2lists-i1-rename/r1list.txt'
r2list=${path2datadir}'gen/nguyen_nc_2018/20190809-r1r2lists-i1-rename/r2list.txt'

## create a name variable
## note this name has '-R#.fastq.gz' extension included
## e.g. i1-lib001-L001-A01-TAAGGCGA-GCGTAAGA-R1.fastq.gz
name1wext=`head -n $SGE_TASK_ID $r1list | tail -n 1 | cut -f1`
name2wext=`head -n $SGE_TASK_ID $r2list | tail -n 1 | cut -f1`

## remove R# & ext from name
## e.g. i1-lib001-L001-A01-TAAGGCGA-GCGTAAGA
name=${name1wext%'-R'[1-2].fastq.gz}

## assign variables to every piece of info in name
ind=${name:0:2}
lib=${name:3:6}
lane=${name:10:4}
well=${name:15:3}
barcodes=${name:19}
bar1=${name:19:8}
bar2=${name:28}

## create output path
path2output=${path2datadir}'gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/genomic/'

## create filename extenstion variables
samext='.sam'
bamext='.bam'
txtext='.txt'
statext='.stat.txt'
table='table.csv'
ptable='posttable.csv'
vcfext='.vcf'

## file prefixes
aligned='aligned-'
reord='reordered-'
sort='sorted-'
realign='realigned-'
dupm='dupmark-'
filter='filtered-'
recal='recal-'
proper='proper-'
unpair='unpaired-'

## define readgroupinfo variable
## bwa uses RGinfo to label things as normal or tumor
## A readgroup is a set of reads that were generated 
    ## from a single run of a sequencing instrument.
    ## '\t' inserts tabs between variables 
## ID is readgroup ID (one for each illumina run)
## ID in this case is a library
## PU is {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
## PU takes precidence over ID for base recalibration, if present
## PU is not required by GATK
## SM is sample name
## PL is platform/technology used to sequence
readgroupinfo='@RG\tID:'${lib}'\tPU:'${lib}.${lane}.${bar1}${bar2}'\tSM:'${name}'\tPL:Illumina\tLB:'${lib}


##########################
## create sub-directories
##########################

## make directory for each set of files
mkdir $path2genomic$name

## make path variable for file deposition
path2output=$path2genomic$name'/'


######################
## gatk pre-processing
######################

## come back to this later
## Make uBAM with Picard FastqToSam
## https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam
# java -jar picard.jar FastqToSam \
#     F1=$path2fastq$name1wext \
#     O=


## repeat make uBAM for read2 files
# java -jar picard.jar FastqToSam \
#     F1=$path2fastq$name1wext \
#     O=


## align reads to human reference
## `bwa mem [flags] ref.fa.gz r1.fq.gz r2.fq.gz > aligned-name.sam`
## note: ref.fasta.fai and ref.dict can't be zipped!
## http://bio-bwa.sourceforge.net/bwa.shtml
## -M Mark shorter split hits as secondary (for Picard compatibility)
## -t Number of threads
## -R readgroup
bwa mem -M -t 32 \
    -R $readgroupinfo \
    $path2datadir'ref/broad/bundles/b37/human_g1k_v37.fasta.gz' \
    $path2fastq$name1wext $path2fastq$name2wext \
    > $path2output$aligned$name$samext


## create bam and filter reads by quality & location
## `samtools view [options] input.sam [region]`
## `-b` output files in BAM format
## `-q` skip alignments with MAPQ score smaller than INT
    ## - default [0]
    ## use 20 according to Xu et al., eLife, 2019 
## `-o FILE` sends output to FILE 
## `chrM` output all alignments mapped to the reference sequence named `chrM` (i.e. @SQ SN:chrM)
samtools view -b -q 20 \
    -o $path$filter$aligned$name$bamext \
    $path2output$aligned$name$samext chrM 


## sort bam


## create bam index



## mark duplicates







