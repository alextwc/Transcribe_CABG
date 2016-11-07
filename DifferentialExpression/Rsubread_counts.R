library(Rsubread)
library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
bam_file      <- args[1]
count_file    <- args[2]
gtf_file      <- args[3]
count_feature <- args[4]
multi_map     <- as.logical(args[5])

main <- function(){
  count_object <- create_count_object() 
  save(count_object, file = str_c(count_file, '.rda'))
  count_object %>%   
    use_series(counts) %>% 
    data.frame %>% 
    add_rownames('rowname') %>% 
    write_tsv(count_file, col_names = F)
}

create_count_object <- function(){
  count_object <- featureCounts(
    files                  = bam_file,
    annot.ext              = gtf_file,
    isGTFAnnotationFile    = TRUE,
    GTF.attrType           = count_feature,
    GTF.featureType        = "exon",
    isPairedEnd            = TRUE,
    countMultiMappingReads = multi_map,
    fraction               = multi_map)
    
    
}



main()
