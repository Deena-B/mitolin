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
# module load picard-tools/1.96   # for generating uBAM & merging uBAM w aBAM
module load gatk/4.1.2.0        # includes picard tools
module load bwa/0.7.8           # aligner
module load samtools/1.9        # filter reads by quality score, convert sam2bam

# create path variables to access and deposit data 
path2datadir='/dfs3/som/dalawson/drb/deepcelllineage/mitolin/data/'
path2fastq=${path2datadir}'raw/nguyen_nc_2018/ind1/'
path2genomic=${path2datadir}'gen/nguyen_nc_2018/20190809-fastq2vcf/ind1/genomic/'
path2ubams=${path2genomic}'ubams/'
path2aligned=${path2genomic}'aligned/'
path2filtered=${path2genomic}'filtered/'
path2uamerged=${path2genomic}'uamerged/'
path2uamgfil=${path2genomic}'uamgfil/'
path2lanemerged=${path2genomic}'lanemerged/'
path2cells=${path2genomic}'cells/'
path2cell=${path2cells}${cell}'/'

## make directories for each path above
mkdir $path2genomic
mkdir $path2ubams
mkdir $path2aligned
mkdir $path2filtered
mkdir $path2uamerged
mkdir $path2uamgfil
mkdir $path2lanemerged
mkdir $path2cells

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
readsmerged='readsmerged-'
sort='sorted-'
dupm='dupmark-'

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
    > $path2aligned$aligned$cell'-'$lane$samext


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
    -o $path2filtered$filter$aligned$cell'-'$lane$bamext \
    $path2aligned$aligned$cell'-'$lane$samext chrM 


## merge bwa aligned, samtools filtered, bam files with uBAM
    ## MergeBamAlignment - merges aligned with unaligned to create unaligned bam
    ## produces a third SAM or BAM file of aligned and unaligned reads
        ## https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_MergeBamAlignment.php
        ## https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment

gatk MergeBamAlignment \
      -ALIGNED $path2filtered$filter$aligned$cell'-'$lane$bamext \
      -UNMAPPED $path2ubams$unaligned$cell'-'$lane$bamext \
      -O $path2uamerged$uamerged$filter$aligned$cell'-'$lane$bamext \
      -R $path2datadir'ref/broad/bundles/b37/human_g1k_v37.fasta.gz' \
      -PAIRED_RUN TRUE
      -CREATE_INDEX TRUE


## create bam index using picard BuildBamIndex
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/picard_sam_BuildBamIndex.php
    ## generates a BAM index ".bai" file
    ## need this step because samtools view requires an index file
    ## input BAM file must be sorted in coordinate order

gatk BuildBamIndex \
    -I $path2uamerged$uamerged$filter$aligned$cell'-'$lane$bamext


## re-filter by quality & location
    # input BAM file must be sorted in coordinate order 

samtools view -b -q 20 \
    -o $path2uamgfil$filter$uamerged$aligned$cell'-'$lane$bamext \
    $path2uamerged$uamerged$filter$aligned$cell'-'$lane$bamext chrM 


## merge aligned bams that have the same 'SM' (i.e. they are from the same cell)
    ## select all filenames with lane='L001' & store as L1_list, select all filenames with 'L002' & store as L2_list
    ## edit the filenames in the two lists to remove L1 or L2, assign to variables s1 and s2, respectively
    ## if s1=s2, merge bam files with filenames that include L001 and L002
        ## MergeSamFiles - use to combine SAM and/or BAM files from different runs (Lanes) (same as samtools merge)
            ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/picard_sam_MergeSamFiles.php
            ## https://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles




## output of MergeBamAlignment is by default sorted by coordinate order, so I probably don't need 
## sort bam using gatk gatk SortSamSpark
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/org_broadinstitute_hellbender_tools_spark_pipelines_SortSamSpark.php
    ## sorts reads by coordinate order  

# gatk SortSamSpark \
#     -I $path2output$filter$aligned$name$bamext \
#     -O $path2output$sort$filter$aligned$name$bamext


## mark duplicates using gatk MarkDuplicatesSpark
    ## https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.3.0/org_broadinstitute_hellbender_tools_spark_transforms_markduplicates_MarkDuplicatesSpark.php
    ## 

# gatk MarkDuplicates \
#     -I $path2output$sort$filter$aligned$name$bamext \
#     -O $path2output$dupm$sort$filter$aligned$name$bamext \
#     -M $path2output'lib_complexity_metrics_from_markdups'$name'.txt'

# next add create index ^
# check remove duplicates = True as default ^





