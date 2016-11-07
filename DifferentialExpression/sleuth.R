library(sleuth)
library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(magrittr)



#------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input_folder        <- args[1]
output_folder       <- args[2]
pheno_file          <- args[3]
design_formula      <- as.formula(args[4])


main <- function(){
  sample_table  <- create_sample_table()
  sleuth_object <- create_sleuth_object(sample_table)
  write_sleuth_table(sleuth_object)
}

create_sample_table <- function(){
  sample_table  <- pheno_file %>% 
    read.table(sep = "\t", header = T, stringsAsFactors = F) %>% 
    select(-fileName) %>% 
    rowwise %>% 
    mutate(path = str_c(input_folder, sampleName)) %>% 
    rename(sample = sampleName)
}

create_sleuth_object <- function(sample_table){
  gene_table    <- read_tsv('/udd/reala/reference_files/DE/transcript_to_gene_table.tsv')
  sleuth_object <- sample_table %>% 
    sleuth_prep(design_formula, target_mapping = gene_table, max_bootstrap = 25) %>% 
    sleuth_fit 
}

write_sleuth_table <- function(sleuth_object){
  covariate <- get_covariate(sleuth_object)
  sleuth_object %>% 
    sleuth_wt(covariate) %>% 
    sleuth_results(covariate) %>%
    select(target_id, ens_gene, ext_gene, everything()) %>% 
    write_tsv(str_c(output_folder, 'results.tsv'))
}

get_covariate <- function(sleuth_object){
  index <- length(colnames(sleuth_object$design_matrix))
  return(colnames(sleuth_object$design_matrix)[index])
}



main()



