#!/bin/bash
#$ -N fastq2markdups
#$ -ckpt restart
#$ -q som,pub64,free64,asom
#$ -pe make 64
#$ -t 1-4

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

## add packages to PATH
module load java/1.8.0.111      # language of gatk & picard-tools
module load gatk/4.1.3.0        # includes picard tools
module load bwa/0.7.8           # aligner
module load samtools/1.9        # filter reads by quality score, convert sam2bam

## make a directory for generated genomic files
mkdir ../genomic/
mkdir ../genomic/ubams/

## create absolute path2datadir
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'

## create path2fastq holding directory variable
path2fastq=${path2datadir}'raw/nguyen_nc_2018/ind1/'

# create path to deposit data
path2genomic=${path2datadir}'gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/genomic/'

## create varibles for lists of files
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

## create output path
path2output=${path2datadir}'gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/genomic/'
path2ubams=${path2output}'ubams/'

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
reord='reordered-'
sort='sorted-'
realign='realigned-'
dupm='dupmark-'
filter='filtered-'
recal='recal-'
proper='proper-'
unpair='unpaired-'

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


##########################
## create sub-directories
##########################

## make directory for each set of files
mkdir $path2genomic$name

## make path variable for file deposition
path2output=$path2genomic$name'/'


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
    -O $path2ubams$unaligned$cell'-'$lane$bamext \
    -SM $cell \
    -RG ${lib}'.'${lane}


## align reads to human reference
    ## `bwa mem [flags] ref.fa.gz r1.fq.gz r2.fq.gz > aligned-name.sam`
    ## note: ref.fasta.fai and ref.dict can't be zipped!
    ## http://bio-bwa.sourceforge.net/bwa.shtml
    ## -M Mark shorter split hits as secondary (for Picard compatibility)
    ## -t Number of threads
    ## -R readgroup

bwa mem -M -t 32 \
    -R $readgroupinfo \
    $path2datadir'ref/broad/bundles/b37/human_g1k_v37.fasta.gz' \
    $path2fastq$name1wext $path2fastq$name2wext \
    > $path2output$aligned$name$samext


## create bam and filter reads by quality & location using samtools view
    ## http://www.htslib.org/doc/samtools.html
    ## `samtools view [options] input.sam [region]`
    ## `-b` output files in BAM format
    ## `-q` skip alignments with MAPQ score smaller than INT
        ## default [0]
        ## use 20 according to Xu et al., eLife, 2019 
    ## `-o FILE` sends output to FILE 
    ## `chrM` output all alignments mapped to the reference sequence named `chrM` (i.e. @SQ SN:chrM)

samtools view -b -q 20 \
    -o $path2output$filter$aligned$name$bamext \
    $path2output$aligned$name$samext chrM 



## merge bwa aligned, samtools filtered, bam files from the same cell
    ## select all filenames with lane='L001' & store as L1_list, select all filenames with 'L002' & store as L2_list
    ## edit the filenames in the two lists to remove L1 or L2, assign to variables s1 and s2, respectively
    ## if s1=s2, merge bam files with filenames that include L001 and L002
    ## is there a gatk merge bam files tool?
        ## MergeBamAlignment - merges aligned with unaligned to create unaligned bam
            ## https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_MergeBamAlignment.php
        ## MergeSamFiles - use to combine SAM and/or BAM files from different runs (Lanes) (same as samtools merge)
            ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/picard_sam_MergeSamFiles.php
            ## https://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles





## sort bam using gatk gatk SortSamSpark
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/org_broadinstitute_hellbender_tools_spark_pipelines_SortSamSpark.php
    ## sorts reads by coordinate order  

gatk SortSamSpark \
    -I $path2output$filter$aligned$name$bamext \
    -O $path2output$sort$filter$aligned$name$bamext


## mark duplicates using gatk MarkDuplicatesSpark
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/org_broadinstitute_hellbender_tools_spark_transforms_markduplicates_MarkDuplicatesSpark.php
    ## 

gatk MarkDuplicates \
    -I $path2output$sort$filter$aligned$name$bamext \
    -O $path2output$dupm$sort$filter$aligned$name$bamext \
    -M $path2output'lib_complexity_metrics_from_markdups'$name'.txt'

# next add create index ^
# check remove duplicates = True as default ^





