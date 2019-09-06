#!/bin/bash
#$ -N fastq2markdups
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 2-4

#############
## Docstring
#############

## before running this script
## make and move to this directory 
## /dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/erroroutput/
## run this script from this ^ directory
## so e & o files get deposited there

## to run a script on hpc use `qsub path/to/filename.sh`
    ## `qsub /dfs3/som/dalawson/drb/deepcelllineage/mitolin/src/fastq2vcf-gatk-v0.1.0.sh`
## to check on the script's status: `qstat -u dalawson`
## to stop a submitted project by job-ID number: `qdel job-ID-number`

## Picard-tools, which are wrapped in gatk, are Java Archive (JAR) file tools 
## the following syntax should be followed
## java -Xmx2g -jar JARfilename toolname
    ## the -Xmx2g flag assigns Java Virtual Machine (JVM) 2GB of memory, this should work for most tools
    ## the -jar flag indicates that a JAR file will follow
    ## JARfilename for picard-tools is picard.jar
    ## toolnames vary by task, flags for specific tools are covered below


################################################
## load modules, create variables & directories
################################################

## on hpc add packages to PATH
module load java/1.8.0.111          # language of gatk & picard-tools
# module load picard-tools/1.96     # for generating uBAM & merging uBAM w aBAM
module load gatk/4.1.2.0            # includes picard tools
module load bwa/0.7.17-5            # aligner
module load samtools/1.9            # view (filter/convert) sort index

# create path variables to access and deposit data 
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'
path2fastq=${path2datadir}'raw/nguyen_nc_2018/ind1/'
path2genomic=${path2datadir}'gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/genomic/'
path2ubams=${path2genomic}'1-ubams/'
path2aligned=${path2genomic}'2-aligned/'
path2uamg=${path2genomic}'3-uamg/'
path2filuamg=${path2genomic}'4-filuamg/'

## make directories for each path above
mkdir $path2genomic
mkdir $path2ubams
mkdir $path2aligned
mkdir $path2uamg
mkdir $path2filuamg

## make .keep files so folders are tracked by git
touch $path2genomic/.keep
touch $path2ubams/.keep
touch $path2aligned/.keep
touch $path2uamg/.keep
touch $path2filuamg/.keep 

## create varibles for lists of fastq files
r1list=${path2datadir}'gen/nguyen_nc_2018/20190809-r1r2lists-i1-rename/r1list.txt'
r2list=${path2datadir}'gen/nguyen_nc_2018/20190809-r1r2lists-i1-rename/r2list.txt'

## create a name variable
## note this name has '-R#.fastq.gz' extension included
## e.g. i1-lib001-L001-A01-TAAGGCGA-GCGTAAGA-R1.fastq.gz
name1wext=`head -n $SGE_TASK_ID $r1list | tail -n 1 | cut -f1`
name2wext=`head -n $SGE_TASK_ID $r2list | tail -n 1 | cut -f1`

## remove R# & ext from name
## e.g. name=i1-lib001-L001-A01-TAAGGCGA-GCGTAAGA
name=${name1wext%'-R'[1-2].fastq.gz}

## assign variables to every piece of info in name
ind=${name:0:2}
lib=${name:3:6}
lane=${name:10:4}
well=${name:15:3}
barcodes=${name:19}
bar1=${name:19:8}
bar2=${name:28}
cell=${ind}'-'${lib}'-'${well}

## create filename extenstion variables
samext='.sam'
bamext='.bam'

## define readgroupinfo variables
    ## https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
    ## bwa uses RGinfo to label things as normal or tumor
    ## '\t' inserts tabs between variables 
    ## for ind1,2,3: lib.lane is used instead of flowcell name & barcode
        ## this is okay since this dataset doesn't have any runs that contain mixed libraries
    ## ID = readgroup ID (one number for each sequencing run)
        ## Illumina recommends ID = FLOWCELL_NAME.FLOWCELL_BC.LANEno
        ## e.g. lib001.L001
    ## PU = Platform Unit
        ## PU is {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
        ## PU takes precidence over ID for base recalibration, if present
        ## PU is not required by GATK
    ## SM is sample name
        ## e.g. i1-lib001-L001-A01-TAAGGCGA-GCGTAAGA
    ## PL is platform/technology used to sequence
    ## LB is DNA preparation library identifier 
readgroupinfo='@RG\tID:'${lib}.${lane}'\tPU:'${lib}.${lane}.${bar1}${bar2}'\tSM:'${cell}'\tPL:Illumina\tLB:'${lib}


######################
## gatk pre-processing
######################


## Make uBAM with picard-tools FastqToSam
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/picard_sam_FastqToSam.php
    ## https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam
    ## put all outputs into one folder, so they are easy to access later

gatk FastqToSam \
    -F1 $path2fastq$name1wext \
    -F2 $path2fastq$name2wext \
    -O $path2ubams$cell'-'$lane$bamext \
    -RG ${lib}'.'${lane} \
    -PU ${lib}.${lane}.${bar1}${bar2} \
    -SM $cell \
    -PL 'Illumina' \
    -LB $lib


## align reads to human reference
    ## `bwa mem [flags] ref.fa.gz r1.fq.gz r2.fq.gz > aligned-name.sam`
    ## note: ref.fasta.fai and ref.dict can't be zipped!
    ## http://bio-bwa.sourceforge.net/bwa.shtml
    ## -M Mark shorter split hits as secondary (for Picard compatibility)
    ## -t Number of threads
    ## -R readgroup
    ## Error file report starts with "[M::main_mem]"

bwa mem -M -t 32 \
    -R $readgroupinfo \
    $path2datadir'ref/broad/bundles/b37/human_g1k_v37.fasta.gz' \
    $path2fastq$name1wext $path2fastq$name2wext \
    > $path2aligned$cell'-'$lane$samext 


## merge bwa aligned, (optionally samtools filtered), sam or bam files with uBAM files
    ## MergeBamAlignment - merges aligned with unaligned to create unaligned bam
    ## produces a third SAM or BAM file of aligned and unaligned reads
        ## https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_MergeBamAlignment.php
        ## https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment

gatk MergeBamAlignment \
    -UNMAPPED $path2ubams$cell'-'$lane$bamext \
    -ALIGNED $path2aligned$cell'-'$lane$samext \
    -O $path2uamg$cell'-'$lane$bamext \
    -R $path2datadir'ref/broad/bundles/b37/human_g1k_v37.fasta.gz' \
    -SO coordinate \
    --CREATE_INDEX true \
    --ALIGNED_READS_ONLY true 


## filter by quality & location
    # input BAM file must have an index & be sorted in coordinate order 

samtools view -b -q 20 \
    -o $path2filuamg$cell'-'$lane$bamext \
    $path2uamg$cell'-'$lane$bamext MT 


## see ipynb *pairl1l2* for generation of lists of paired samples
## see script gatkpp-2* for next steps