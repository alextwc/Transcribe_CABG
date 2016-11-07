library(DEXSeq)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
input_folder        <- args[1]
output_folder       <- args[2]
pheno_file          <- args[3]
design_string       <- args[4]
add_hgnc            <- as.logical(args[5])
gff_file            <- args[6]


# -----------------------------------------------------------------------------

main <- function(){
  dexseq_object <- create_dexseq_object()
  saveRDS(dexseq_object, file =  str_c(output_folder, '/DEXSeq_object.rds'))
  results <- DEXSeqResults(dexseq_object)
  write_results_table(results, output_folder)
  write_count_table(results, output_folder)
}

create_dexseq_object <- function(){
  sample_table     <- read.table(pheno_file, row.names = NULL, header = T, stringsAsFactors = F) 
  dexseq_object    <- DEXSeqDataSetFromHTSeq(countfiles    = str_c(input_folder, sample_table$fileName),
                                             sampleData    = sample_table,
                                             flattenedfile = gff_file )
  dexseq_object <- dexseq_object %>% 
    estimateSizeFactors() %>% 
    estimateDispersions() %>% 
    testForDEU() %>% 
    estimateExonFoldChanges(fitExpToVar = "condition")
}

write_results_table <- function(results, output_folder){
  results_table <- results %>%
    data.frame %>%
    add_rownames("ensemble")
  if(add_hgnc) hgnc_table <- read_tsv('/udd/reala/reference_files/DE/ensemble_to_hgnc_exon_table.tsv')
  if(add_hgnc) results_table <- left_join(results_table, hgnc_table)
  write_tsv(results_table, str_c(output_folder, '/DEXSeq_results.tsv'))
}

write_count_table <- function(results, output_folder){
  count_table <- results %>%
    counts %>%
    data.frame %>%
    set_colnames(results@sampleData@listData$sampleName) %>%
    data.frame %>%
    add_rownames("ensemble")
  if(add_hgnc) hgnc_table <- read_tsv('/udd/reala/reference_files/DE/ensemble_to_hgnc_exon_table.tsv')
  if(add_hgnc) count_table <- left_join(count_table, hgnc_table)
  write_tsv(count_table, str_c(output_folder, '/DEXSeq_counts.tsv'))
}

main()