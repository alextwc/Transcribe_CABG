library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(magrittr)


setup_limma_app <- function(input_folder, output_folder, batch, mode = 'gene', 
                            script_dir = '/udd/reala/repos/prod/Differential_expression/'){
  data_folder <- str_c(output_folder, '/data/')
  batch_folder <- str_c(data_folder, batch)
  dir.create(output_folder)
  dir.create(data_folder)
  dir.create(batch_folder)
  copy_shiny_files(output_folder, mode, script_dir)
  link_data_folders(input_folder, batch_folder)
  if(mode == 'gene') file.symlink('/udd/reala/reference_files/DE/ensemble_to_hgnc_gene_table.tsv', data_folder)
  if(mode == 'exon') file.symlink('/udd/reala/reference_files/DE/exon_to_gene_table.tsv', data_folder)
}

copy_shiny_files <- function(output_folder, mode, script_dir){
  if(mode == 'gene') file.copy(str_c(script_dir, 'limma_gene_app.R'), str_c(output_folder, 'app.R'), overwrite = T)
  if(mode == 'exon') file.copy(str_c(script_dir, 'limma_exon_app.R'), str_c(output_folder, 'app.R'), overwrite = T)
  file.copy(str_c(script_dir, 'limma_app_functions.R'), str_c(output_folder, 'functions.R'), overwrite = T)
}

link_data_folders <- function(input_folder, batch_folder){
  limma_folders <- input_folder %>% 
    list.dirs(recursive = F) %>% 
    discard(str_detect(., 'rout')) %>% 
    map(function(x) file.symlink(x, batch_folder))
}

