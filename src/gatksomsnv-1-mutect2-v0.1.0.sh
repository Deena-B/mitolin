#!/bin/bash
#$ -N mutect2
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
    ## 20190917-mutect2-DRB/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## to run a script on hpc use `qsub path/to/filename.sh`
    ## `qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/gatksomsnv-1-mutect2-v0.1.0.sh`
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`


################################################
## load modules, create variables & directories
################################################

## add packages to PATH
module load java/1.8.0.111          # language of gatk & picard-tools
module load gatk/4.1.3.0            # includes picard tools

# create path variables to access data 
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'
# path to bam files generated by gatkpp-5-applybqsr*.sh
path2bqsrbams=${path2datadir}'gen/nguyen_nc_2018/20190911-bqsr-DRB/genomic/bqsrbams/'
# the path below has list: bqsrbamlist.txt
path2list=${path2bqsrbams}

# create path variables (for data deposit)
path2genomic=${path2datadir}'gen/nguyen_nc_2018/20190911-bqsr-DRB/genomic/'
path2vcfs=${path2genomic}'vcfs/'

## make directories for each path to a directory that doesn't yet exist (in list above above)
mkdir $path2vcfs

## create variable for list of table files
bqsrbamslist=${path2list}'bqsrbamlist.txt'

## create a name variable
## note this name has an extension included (e.g. '.bam')
## e.g. i1-lib001-A01.bam
namewext=`head -n $SGE_TASK_ID $bqsrbamslist | tail -n 1 | cut -f1`

## create name variable with ext removed
## e.g. name=i1-lib001-A01
name=${namewext%'.bam'}


############################################
## gatk somatic snv & indel variant calling
############################################

## make vcf with Mutect2 in mito mode
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php
    ## Call somatic short diffs (SNVs & indels)
    ## Blog with details of Mutect2 mitochondrial mode
        ## https://software.broadinstitute.org/gatk/blog?id=23598
    ## Mitochondria mode sets (lod=='limit of detection')
        ## --initial-tumor-lod to 0
        ## --tumor-lod-to-emit to 0
        ## --af-of-alleles-not-in-resource to 4e-3
        ## --pruning-lod-threshold to -4
    ## Mito mode only accepts a single sample

gatk Mutect2 \
    -R $path2datadir'ref/broad/bundles/b37/human_g1k_v37.fasta' \
    -L MT \
    --mitochondria-mode true \
    -I $path2bqsrbams$name.bam \
    -O $path2mutect2$name.vcf.gz