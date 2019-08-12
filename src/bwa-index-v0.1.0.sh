#!/bin/bash
#$ -N bwa-idx-bwtsw
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 1


#################################
## how to run this script on hpc
#################################

## before running this script
## move to the directory below:
    ## /dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/ref/broad/bundles/b37/
## run this script from this ^ directory

## to run this script on hpc use `qsub path/to/filename.sh`
    ## qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/bwa-index-v0.1.0.sh
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`


###############
## load modules
###############

module load bwa/0.7.8


###############
## run packages 
###############

## index the reference fasta
## -p prefix (allows you to rename the output file prefix)
    ## suffixes are .amb, .ann, .bwt, .pac, .sa
## `-a bwtsw` for long sequences 
bwa index -p human_g1k_v37_long -a bwtsw human_g1k_v37.fasta.gz

