





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





