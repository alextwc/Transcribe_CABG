library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(magrittr)

  

# -----------------------------------------------------------------------------


args <- commandArgs(trailingOnly = TRUE)
input_file        <- args[1]
output_file       <- args[2]
count_method      <- args[3]

main <- function(){
  autosome_feature_list <- create_autosome_feature_list(count_method)
  count_table           <- read_tsv(input_file, col_names = F)
  autosome_count_table  <- subset_count_table(count_method, count_table, autosome_feature_list)
  write_tsv(autosome_count_table, output_file, col_names = F)
}


create_autosome_feature_list <- function(count_method){
  if (count_method == 'rsubread_gene') file <- '/udd/reala/reference_files/DE/gene_autosome_only_features.tsv'
  if (count_method == 'rsubread_exon') file <- '/udd/reala/reference_files/DE/exon_autosome_only_features.tsv'
  if (count_method == 'htseq_gene')   file <- '/udd/reala/reference_files/DE/gene_autosome_only_features.tsv'
  if (count_method == 'htseq_exon')   file <- '/udd/reala/reference_files/DE/gene_autosome_only_features.tsv'
  autosome_feature_list <- file %>% 
    read_tsv %>% 
    use_series('feature')
}

subset_count_table <- function(count_method, count_table, autosome_feature_list){
  if (count_method == 'rsubread_gene') return(subset_count_table_standard(count_table, autosome_feature_list)) 
  if (count_method == 'rsubread_exon') return(subset_count_table_standard(count_table, autosome_feature_list))
  if (count_method == 'htseq_gene')   return(subset_count_table_standard(count_table, autosome_feature_list))
  if (count_method == 'htseq_exon')   return(subset_count_table_htseq_exon(count_table, autosome_feature_list))
}

subset_count_table_standard <- function(count_table, autosome_feature_list){
  filter(count_table, (X1 %in% autosome_feature_list))
} 

subset_count_table_htseq_exon <- function(count_table, autosome_feature_list){
  count_table %>% 
    rowwise %>% 
    mutate(gene = str_split(X1, ':')[[1]][[1]]) %>% 
    filter(gene %in% autosome_feature_list) %>% 
    select(-gene)
} 

main()
