#!/bin/bash
#$ -N markdups
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 1-118


#############
## Docstring
#############

## before running this script
## make and move to this directory 
## /dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/ \
    ## 20190909-markdups-DRB/ind1/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## to run a script on hpc use `qsub path/to/filename.sh`
    ## `qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/gatkpp-3-markdup2bqsr-v0.1.0.sh`
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`


################################################
## load modules, create variables & directories
################################################

## add packages to PATH
module load java/1.8.0.111          # language of gatk & picard-tools
module load gatk/4.1.3.0            # includes picard tools
module load samtools/1.9            # view (filter/convert) sort index

# create path variables to access data 
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'
# the path below is to samples that were only run in a single lane
    # those listed in 'unpairedlist.txt'
path2filuamgbam=${path2datadir}'gen/nguyen_nc_2018/20190809-gatkpp1-fastq2uamgfil-DRB/ind1/genomic/4-filuamg'
# the path below is to lane-merged bam & bai files
path2luamgfilbam=${path2datadir}'gen/nguyen_nc_2018/20190821-gatkpp2-lanesmerged-DRB/ind1/genomic/luamgfil/'
# the path below has lists: unpairedlist.txt & lmglist.txt
path2lists=${path2datadir}'gen/nguyen_nc_2018/20190906-celllist-DRB/'

# create paths to deposit data
path2genomic=${path2datadir}'gen/nguyen_nc_2018/20190909-markdups-DRB/ind1/genomic/'
path2dupsmarked=${path2genomic}'dupsmarked/'
path2dupmetrics=${path2genomic}'dupmetrics/'

## make directories for each 'deposit data' path above
mkdir $path2genomic
mkdir $path2dupsmarked
mkdir $path2dupmetrics

## create variables for lists of bam files
unpairedlist=${path2lists}'unpairedlist.txt'
lmglist=${path2lists}'lmglist.txt'

## create a name variable
## note this name has '.bam' extension included
## e.g. i1-lib001-A01.bam
namewext=`head -n $SGE_TASK_ID $unpairedlist | tail -n 1 | cut -f1`

## remove ext from name
## e.g. name=i1-lib001-A01
name=${namewext%'.bam'}

## create filename extenstion variables
samext='.sam'
bamext='.bam'
txtext='.txt'
statext='.stat.txt'
table='table.csv'
ptable='posttable.csv'
vcfext='.vcf'


#############################
## gatk pre-processing con't
#############################


## mark duplicates using gatk MarkDuplicates
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/picard_sam_markduplicates_MarkDuplicates.php
    ## All this does is mark duplicates with flags and generate metrics
    ## You could choose to remove duplicates
    ## I think marking dups is enough for them to be exculded from vcf
    ## Inputs must be coordinate sorted

# below currently only processes luamgfilbams
gatk MarkDuplicates \
    -I $path2filuamgbam$name$bamext \
    -O $path2dupsmarked$name$bamext \
    -M $path2dupmetrics$name'.txt' \
    --TAG_DUPLICATE_SET_MEMBERS true \
    --CREATE_INDEX 


## What goes next?





