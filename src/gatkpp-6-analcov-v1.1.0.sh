#!/bin/bash
#$ -N analcov
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 2-3


#############
## Docstring
#############

## before running this script
## make and move to this directory 
## /dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/ \
    ## 20191031-bqsr/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## change -t ^ to the samples that you want to run the script on
    ## for testing, use `-t 2-3`
    ## for all, count the number of lines: 
    ## `wc -l dupsmarkedlist.txt'

## to run a script on hpc use `qsub path/to/filename.sh`
    ## `qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/gatkpp-6-analcov-v1.1.0.sh`
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`


################################################
## load modules, create variables & directories
################################################

## add packages to PATH
module load java/1.8.0.111          # language of gatk & picard-tools
module load gatk/4.1.3.0            # includes picard tools
module load R/3.6.0                 # allows analcov to generate pdf 

## create parent path variables
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'
path2gen_nguyen2018=${path2datadir}'gen/nguyen_nc_2018/'

## create path variable to access list: dupsmarkedlist.txt
path2markdupout=${path2gen_nguyen2018}'20191030-markdup/output/'
path2dupsmarked=${path2markdupout}'dupsmarked/'
path2list=${path2dupsmarked}

## these directories were made by gatkpp-4-baserecal
## these paths are needed to access & deposit data
path2output=${path2gen_nguyen2018}'20191031-bqsr/output/'
path2tables=${path2bqsrout}'baserecaltables/'

## create new paths to deposit data
path2analcovs=${path2output}'analcovs/'

## make directories for each new 'deposit data' path above
mkdir $path2analcovs

## make .keep files so new directories are tracked by git
touch $path2analcovs/.keep

## create variable for list of table files
dupsmarkedlist=${path2list}'dupsmarkedlist.txt'

## create a name variable
## note this name has '.bam' extension included
## e.g. i1-lib001-A01.bam
namewext=`head -n $SGE_TASK_ID $dupsmarkedlist | tail -n 1 | cut -f1`

## remove ext from name
## e.g. name=i1-lib001-A01
name=${namewext%'.bam'}


#############################
## gatk pre-processing con't
#############################


## apply base quality score recalibration to bams
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_AnalyzeCovariates.php
    ## generates plots to assess the quality of a recalibration

gatk AnalyzeCovariates \
    -bqsr $path2tables$name'.table' \
    -plots $path2analcovs$name'.pdf' \
    -csv $path2analcovs$name'.csv'


## Next run gatksomsnv-1-mutect2-v0.1.0.sh


