library(limma)
library(edgeR)
library(sva)
library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(magrittr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_folder        <- args[1]
output_folder       <- args[2]
pheno_file          <- args[3]
design_string       <- args[4]
add_hgnc            <- as.logical(args[5])
min_samples_filter  <- args[6]
do_sva              <- as.logical(args[7])


main <- function(){
  pheno_table   <- as_data_frame(fread(pheno_file))
  counts_matrix <- create_counts_matrix(pheno_table)
  sva_object    <- create_sva_object(do_sva, counts_matrix, pheno_table)
  model_matrix  <- create_model_matrix(pheno_table, sva_object)
  DGE_object    <- create_DGE_exon_object(counts_matrix, pheno_table$sampleName)
  feature_table <- create_feature_table(rownames(DGE_object$counts))
  voom_object   <- voom(DGE_object, model_matrix, save.plot = T)
  limma_object  <- create_limma_exon_object(voom_object, model_matrix, feature_table$ensembl_exon_id, feature_table$external_gene_name)
  limma_table   <- topSplice(limma_object, coef = ncol(model_matrix), test = "t", number = Inf)
  # write output
  create_MV_plot(voom_object)
  create_MDS_object(voom_object)
  saveRDS(limma_object, file = (str_c(output_folder, 'limma_object.rds')))
  write_count_tables(limma_table, DGE_object, counts_matrix, feature_table)
  if(do_sva) saveRDS(sva_object, str_c(output_folder, 'sva_object.rds'))
}


# counts matrix ---------------------------------------------------------------
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

# sva -------------------------------------------------------------------------
create_sva_object <- function(do_sva, counts_matrix, pheno_table){
  if(!do_sva) return(F)
  filtered_counts_matrix <- filter_count_matrix(counts_matrix)
  design_string          <- create_design_string(design_string)
  full_model_matrix      <- model.matrix(as.formula(design_string), data = pheno_table)
  reduced_model_matrix   <- model.matrix(create_reduced_formula(design_string), data = pheno_table)
  svseq_obj              <- svaseq(filtered_counts_matrix, full_model_matrix, reduced_model_matrix)
}

filter_count_matrix <- function(counts_matrix){
  filtered_matrix <- counts_matrix[rowSums(cpm(counts_matrix) > 1) >= min_samples_filter, ]
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

create_model_matrix <- function(pheno_table, sva_object){
  design_formula <- as.formula(design_string)
  if (do_sva) 
    model_matrix <- create_sva_model_matrix(pheno_table, sva_object, design_formula)
  if (!do_sva) 
    model_matrix <- model.matrix(design_formula, data = pheno_table)
  return(model_matrix)
}

create_sva_model_matrix <- function(pheno_table, sva_object, design_formula){
  model_matrix <- cbind(sva_object$sv, 
                        model.matrix(design_formula, data = pheno_table))
}


create_feature_table <- function(exon_names){
  exon_to_gene_table <- read_tsv('/udd/reala/reference_files/DE/exon_to_gene_table.tsv')
  feature_table <- exon_names %>% 
    data_frame(ensembl_exon_id = .) %>% 
    left_join(exon_to_gene_table)
}

create_DGE_exon_object <- function(counts_matrix, sample_names){
  DGE_object <- counts_matrix %>%
    DGEList(group = sample_names) %>% 
    filter_by_cpm %>% 
    calcNormFactors()
}

filter_by_cpm <- function(DGE_object){
  keep <- rowSums(cpm(DGE_object) > 1) >= min_samples_filter
  DGE_object <- DGE_object[keep, , keep.lib.sizes=FALSE]
  return(DGE_object)
}

create_limma_exon_object <- function(voom_object, model_matrix, exon_names, gene_names){
  limma_object <- voom_object %>%  
    lmFit(model_matrix) %>%
    eBayes %>% 
    diffSplice(geneid = gene_names, exonid = exon_names)
}

create_MV_plot <- function(voom_object){
  png(str_c(output_folder,"MV_plot.png"), 
      width = 4, height = 4, units = "in", res = 300)
  plot(voom_object$voom.xy, pch = 16, cex = 0.25)
  title("voom: Mean-variance trend")
  lines(voom_object$voom.line, col = "red")
  dev.off()
}

create_MDS_object <- function(voom_object){
  voom_object %>%
    plotMDS %>%
    saveRDS(str_c(output_folder, 'MDS_object.rds'))
}

# count tables ----------------------------------------------------------------

write_count_tables <- function(limma_table, DGE_object, counts_matrix, feature_table){
  pval_order_table <- limma_table %>% 
    arrange(P.Value) %>% 
    select(ExonID) %>% 
    set_colnames('ensembl_exon_id')
  write_normalized_count_table(DGE_object, pval_order_table, feature_table)
  write_count_table(counts_matrix, pval_order_table, 'raw_count', feature_table)
}

write_normalized_count_table <- function(DGE_object, pval_order_table, feature_table){
  lib_size         <- with(DGE_object$samples, lib.size * norm.factors)
  norm_count_table <- DGE_object %>% 
    use_series(counts) %>% 
    add(0.5) %>% 
    divide_by(lib_size + 1) %>% 
    multiply_by(1e+06) %>% 
    write_count_table(pval_order_table, 'norm_count', feature_table)
}

write_count_table <- function(counts_matrix, pval_order_table, prefix, feature_table){
  count_table <- counts_matrix %>% 
    data.frame %>% 
    add_rownames('ensembl_exon_id') %>% 
    left_join(pval_order_table, .) %>% 
    write_table(prefix, feature_table)
}



write_table <- function(table, prefix, feature_table){
  if(add_hgnc) table <- left_join(table, feature_table)
  if(add_hgnc) table <- select(table, ensembl_exon_id, ensembl_gene_id, external_gene_name, everything())
  write_tsv(table, str_c(output_folder, prefix, '_table.tsv'))
}


main()
