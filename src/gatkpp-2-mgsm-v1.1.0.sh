#!/bin/bash
#$ -N mgsm
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 2-4


#############
## Docstring
#############

## before running this script
## make and move to this directory 
## /dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/20191029-mgsm/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## to run a script on hpc use `qsub path/to/filename.sh`
    ## `qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/gatkpp-2-mgsm-v1.1.0.sh`
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`


################################################
## load modules, create variables & directories
################################################

## add packages to PATH
module load java/1.8.0.111          # language of gatk & picard-tools
module load gatk/4.1.3.0            # includes picard tools

# create path variables to access and deposit data 
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'
path2gen_nguyen18=${path2datadir}'gen/nguyen_nc_2018/'
path2bam=${path2gen_nguyen18}'20191025-fastq2ummg/output/3-ummg/'
path2output=${path2gen_nguyen18}'20191029-mgsm/output'
path2lummg=${path2output}'lummg/'

## make directories for each path above
mkdir $path2output
mkdir $path2lummg

## make .keep files so folders are tracked by git
touch $path2output/.keep
touch $path2lummg/.keep 

## create path variables to access lists of bam files
path2lists=${path2gen_nguyen18}'20190906-pairlanelists-DRB/'

## create variables for lists of bam files
## generated using ipynb *pairl1l2*
l1list=${path2lists}'lane1list-paired.txt'
l2list=${path2lists}'lane2list-paired.txt'
singlerunlist=${path2lists}'unpairedlist.txt'

## create a name variable
## note this name has '-L00#.bam' extension included
## e.g. i1-lib001-A01-L001.bam
name1wext=`head -n $SGE_TASK_ID $l1list | tail -n 1 | cut -f1`
name2wext=`head -n $SGE_TASK_ID $l2list | tail -n 1 | cut -f1`

## remove Lane # & ext from name
## e.g. i1-lib001-A01
name=${name1wext%'-L00'[1-2].bam}


#############################
## gatk pre-processing con't
#############################


## merge aligned bams from the same cell 
    ## a cell is identified by its individual, library, and well, e.g.: 'i1-lib001-A01'
    ## equivalent cells should also have the same @RG sample name 'SM' 
    ## lists of bam files to be merged were generated using the nb '20190820-pairl1l2-DRB.ipynb'
        ## MergeSamFiles - use to combine SAM and/or BAM files from different runs (Lanes) (same as samtools merge)
            ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/picard_sam_MergeSamFiles.php
            ## https://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles

gatk MergeSamFiles \
    -I ${path2bam}${name1wext} \
    -I ${path2bam}${name2wext} \
    -O ${path2lummg}${name}'.bam' \
    --CREATE_INDEX true

## next gatkpp3-markdup








