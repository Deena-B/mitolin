#!/bin/bash
#$ -N fastq2abam
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 1-10

## load modules
module load java/1.8.0.111      # language of gatk & picard-tools
module load gatk/4.1.2.0        # includes picard tools
module load bwa/0.7.8           # aligner


######################
## create directories
######################

## before running this script
## make and move to this directory 
## /dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/20190718-fastq2vcf/ind1/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## make a directory for generated genomic files
mkdir ../genomic/


####################
## create variables
####################

## create absolute path2data
path2data='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'

## create path2fastq holding directory
path2fastq=${path2data}'raw/nguyen_nc_2018/ind1/'

## create varibles for lists of files
r1list=${path2data}'gen/nguyen_nc_2018/20190715-pairr1r2-renamed/ind1/r1_list.txt'
r2list=${path2data}'gen/nguyen_nc_2018/20190715-pairr1r2-renamed/ind1/r2_list.txt'

## create a name variable
## note this name has '-R#.fastq.gz' extension included
## e.g. L001-A02-CGTACTAG-GCGTAAGA-R1.fastq.gz
namewfqgzext=`head -n $SGE_TASK_ID $r1list | tail -n 1 | cut -f1`

## remove R# & ext from name
## e.g. L001-A02-CGTACTAG-GCGTAAGA
name=${namewfqgzext%'-R'[0-2].fastq.gz}
platewell=${name:8}
plate=${name:0:4}
well=${name:5:3}
barcodes=${name:9}
bar1=${name:9:8}
bar2=${name:18:8}

## create output path
path2output=${path2data}'gen/nguyen_nc_2018/20190718-fastq2vcf/ind1/genomic/'

## create filename extenstion variables
samext='.sam'
bamext='.bam'
txtext='.txt'
statext='.stat.txt'
table='table.csv'
ptable='posttable.csv'
vcfext='.vcf'

## file prefixes
reord='reordered-'
sort='sorted-'
realign='realigned-'
dupm='dupmark-'
filter='filtered-'
recal='recal-'
proper='proper-'
unpair='unpaired-'

## assign a line of text to the readgroupinfo variable
## bwa uses RGinfo to label things as normal or tumor
## A readgroup is a set of reads that were generated 
## from a single run of a sequencing instrument.
## Read as: ID = EGA, SM = name, PL = Illumina
## with '\t' = tabs between them 
## ID = readgroup ID (one for each illumina run)
## Used platename, e.g.: "L001" (probably not correct)
## SM = sample name
## PU = {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
## PU takes precidence over ID for base recalibration, if present
## PU is not required by GATK
## PL = platform/technology used to sequence
readgroupinfo='@RG\tID:'${plate}'\tSM:'${name}'\tPL:Illumina'


##########################
## create sub-directories
##########################

## dirname = name - 
dir=${}

