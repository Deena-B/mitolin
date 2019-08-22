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
## /dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## to run a script on hpc use `qsub path/to/filename.sh`
    ## `qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/uamgfil2bqsrecal-gatk-v0.1.0.sh`
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`


################################################
## load modules, create variables & directories
################################################

## add packages to PATH
module load java/1.8.0.111          # language of gatk & picard-tools
# module load picard-tools/1.96     # for generating uBAM & merging uBAM w aBAM
module load gatk/4.1.3.0            # includes picard tools
module load samtools/1.9            # view (filter/convert) sort index

# create path variables to access and deposit data 
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'
path2bam=${path2datadir}'gen/nguyen_nc_2018/20190809-fastq2uamgfil/ind1/genomic/uamgfil/'
path2genomic=${path2datadir}'gen/nguyen_nc_2018/20190821-uamgfil2bqsrecal-DRB/ind1/genomic/'
path2lanemerged=${path2genomic}'lanemerged/'
path2cells=${path2genomic}'cells/'
path2cell=${path2cells}${name}'/'

## make directories for each path above
mkdir $path2genomic
mkdir $path2lanemerged
mkdir $path2cells

## create variables for lists of bam files
l1list=${path2datadir}'gen/nguyen_nc_2018/20190820-pairlanelists-DRB/lane1list-paired.txt'
l2list=${path2datadir}'gen/nguyen_nc_2018/20190820-pairlanelists-DRB/lane2list-paired.txt'
singlerunlist=${path2datadir}'gen/nguyen_nc_2018/20190820-pairlanelists-DRB/unpairedlist.txt'

## create a name variable
## note this name has '-L00#.bam' extension included
## e.g. filtered-uamerged-aligned-i1-lib001-A01-L001.bam
name1wext=`head -n $SGE_TASK_ID $l1list | tail -n 1 | cut -f1`
name2wext=`head -n $SGE_TASK_ID $l2list | tail -n 1 | cut -f1`

## remove L# & ext from name
## e.g. name=filtered-uamerged-aligned-i1-lib001-A01
name=${name1wext%'-L00'[1-2].bam}

## make directory for individual cells
mkdir $path2cell

## create filename extenstion variables
samext='.sam'
bamext='.bam'
txtext='.txt'
statext='.stat.txt'
table='table.csv'
ptable='posttable.csv'
vcfext='.vcf'

## file prefixes
unaligned='unaligned-'
aligned='aligned-'
filter='filtered-'
uamerged='uamerged-'
lanesmerged='lanesmerged-'
sort='sorted-'
dupm='dupmark-'


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

# gatk MergeSamFiles \
#     -I ${path2bam}${name1wext} \
#     -I ${path2bam}${name2wext} \
#     -O ${path2lanemerged}${name}${bamext} \
#     --CREATE_INDEX true


## mark duplicates using gatk MarkDuplicatesSpark
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/org_broadinstitute_hellbender_tools_spark_transforms_markduplicates_MarkDuplicatesSpark.php
    ## 

# gatk MarkDuplicates \
#     -I $path2output$sort$filter$aligned$name$bamext \
#     -O $path2output$dupm$sort$filter$aligned$name$bamext \
#     -M $path2output'lib_complexity_metrics_from_markdups'$name'.txt'

## next add create index ^
## check remove duplicates = True as default ^



## output of MergeSamFiles is, by default, sorted by coordinate order, so I probably don't need 
## sort bam using gatk gatk SortSamSpark
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/org_broadinstitute_hellbender_tools_spark_pipelines_SortSamSpark.php
    ## sorts reads by coordinate order  

# gatk SortSamSpark \
#     -I $path2output$filter$aligned$name$bamext \
#     -O $path2output$sort$filter$aligned$name$bamext



