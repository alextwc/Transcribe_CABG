library(sva)
library(edgeR)
library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(magrittr)




# -----------------------------------------------------------------------------




args <- commandArgs(trailingOnly = TRUE)
input_folder        <- args[1]
output_folder       <- args[2]
pheno_file          <- args[3]
design_string       <- args[4]


main <- function(){
  pheno_table            <- create_pheno_table(pheno_file)
  counts_matrix          <- create_counts_matrix(pheno_table)
  filtered_counts_matrix <- filter_count_matrix(counts_matrix)
  design_string          <- create_design_string(design_string)
  full_model_matrix      <- model.matrix(as.formula(design_string), data = pheno_table)
  reduced_model_matrix   <- model.matrix(create_reduced_formula(design_string), data = pheno_table)
  svseq_obj              <- svaseq(filtered_counts_matrix, full_model_matrix, reduced_model_matrix)
  saveRDS(svseq_obj, str_c(output_folder, 'sva_object.rds'))
  write_tsv(pheno_table, (str_c(output_folder,'pheno_file.tsv')))
}


create_pheno_table <- function(pheno_file){
  pheno_table <- read.table(pheno_file, sep = "\t", header = T, stringsAsFactors = F)
  if ('case.control' %in% colnames(pheno_table))
    pheno_table <- mutate(pheno_table, case.control = ifelse(case.control == 1, 'control', 'case'))
  if ('race.V2' %in% colnames(pheno_table))
    pheno_table <- mutate(pheno_table, race.V2 = ifelse(race.V2 == 1, 'NHW', 'AA'))
  if ('gender.V2' %in% colnames(pheno_table))
    pheno_table <- mutate(pheno_table,gender.V2 = ifelse(gender.V2 == 1, 'male', 'female'))
  if ('SmokCigNow.V2' %in% colnames(pheno_table))
    pheno_table <- mutate(pheno_table,SmokCigNow.V2 = ifelse(SmokCigNow.V2 == 1, 'current', 'former'))
  return(pheno_table)
}

create_counts_matrix <- function(pheno_table){
  counts_matrix <- pheno_table %>%
    use_series(fileName) %>%
    str_c(input_folder, .) %>%
    map(read_tsv, col_names = F) %>%
    map2(pheno_table$sampleName,  function(x,y) set_colnames(x, c('feature', y))) %>% 
    reduce(left_join) %>% 
    filter(!(str_detect(feature, '^_'))) %>% 
    set_rownames(.$feature) %>% 
    select(-feature) %>% 
    as.matrix
}

create_design_string <- function(design_string){
  design_string %>% 
    str_match(pattern = '^~([:print:]+)') %>% 
    extract2(2) %>% 
    str_c('~1+', .)
}


create_reduced_formula <- function(design_string){
  design_string %>%
    str_match(pattern = '^(~[:print:]+)\\+[:print:]+$') %>% 
    extract2(2) %>% 
    as.formula
}

filter_count_matrix <- function(counts_matrix){
  filtered_matrix <- counts_matrix[rowSums(cpm(counts_matrix) > 1) >= 2, ]
}

main()




