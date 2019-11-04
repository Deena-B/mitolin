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
    ## 20191101-somsnv/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## change -t ^ to the samples that you want to run the script on
    ## for testing, use `-t 2-3`
    ## for all, count the number of lines in list
    ## `cd DATE-markdup/output/dupsmarked/`
    ## `wc -l dupsmarkedlist.txt'

## to run a script on hpc use `qsub path/to/filename.sh`
    ## `qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/gatksomsnv-1-mutect2-v1.2.0.sh`
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`


################################################
## load modules, create variables & directories
################################################

## add packages to PATH
module load java/1.8.0.111          # language of gatk & picard-tools
module load gatk/4.1.3.0            # includes picard tools

# create parent path variables
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'
path2ref=${path2datadir}'ref/broad/bundles/b38/v0/'
path2gen_nguyen2018=${path2datadir}'gen/nguyen_nc_2018/'

## create path variable to access list: dupsmarkedlist.txt
path2markdupout=${path2gen_nguyen2018}'20191030-markdup/output/'
path2dupsmarked=${path2markdupout}'dupsmarked/'
path2list=${path2dupsmarked}

## create variable for list of file names
dupsmarkedlist=${path2list}'dupsmarkedlist.txt'

## these directories were made by gatkpp-4/5-baserecal/applybqsr
## these paths are needed to access & deposit data
path2bqsroutput=${path2gen_nguyen2018}'20191031-bqsr/output/'
path2bqsrbams=${path2bqsroutput}'bqsrbams/'

## create new paths to deposit data
path2somsnv=${path2gen_nguyen2018}'20191101-somsnv/'
path2mutectoutput=${path2somsnv}'output/'
path2mutect2=${path2mutectoutput}'mutect2/'
# path2pileup=${path2genomic}'pileup/'

## make directories for each new 'deposit data' path above
mkdir $path2mutectoutput
mkdir $path2mutect2
# mkdir $path2pileup

## make .keep files so new directories are tracked by git
touch $path2somsnv/.keep
touch $path2mutectoutput/.keep
touch $path2mutect2/.keep
# touch $path2pileup/.keep

## create a name variable
## note this name has an extension included (e.g. '.bam')
## e.g. i1-lib001-A01.bam
namewext=`head -n $SGE_TASK_ID $dupsmarkedlist | tail -n 1 | cut -f1`

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
    -R $path2ref'Homo_sapiens_assembly38.fasta' \
    -L chrM \
    --mitochondria-mode true \
    -I $path2bqsrbams$name'.bam' \
    -O $path2mutect2$name'.vcf.gz'


## get pileup summaries
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php
    ## Tabulates pileup metrics for inferring contamination

# gatk GetPileupSummaries \
#     -I $path2bqsrbams$name'.bam' \
#     -V $path2mutect2$name'.vcf.gz' \
#     -L $path2mutect2$name'.vcf.gz' \
#     -O $path2pileup$name'.table'