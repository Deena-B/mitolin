#!/bin/bash
#$ -N fastq2abam
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 1-10

## load modules
module load gatk/4...
module load java...

## make and move to this directory before running this script
## ~/deepcelllineage/data/gen/nguyen_nc_2018/20190712-fastq2abam/ind2/erroroutput
## run this script from this ^ directory
## so e & o files get deposited there
## on UCI hpc: '~/'' == '/data/users/dalawson/''

## make a directory for generated genomic files
mkdir ../genomic

####################
## create variables
####################

## create path2fastq data variable
path2fastq='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/raw/nguyen_nc_2018/ind2'

## create variables for lists of files
read1list='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/20190527-pairr1r2/ind2/r1_list_pairs.txt'
read2list='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/20190527-pairr1r2/ind2/r2_list_pairs.txt'

## create a name with fastq extension variable
namewfqext=`head -n $SGE_TASK_ID $read1list | tail -n 1 | cut -f1`

