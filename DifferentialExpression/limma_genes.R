library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(magrittr)
library(data.table)
library(tibble)
library(limma)
library(edgeR)
library(sva)
library(goseq)
library(gage)
library(org.Hs.eg.db)
library(GO.db)
library(KEGG.db)




# -----------------------------------------------------------------------------




args <- commandArgs(trailingOnly = TRUE)
input_folder        <- args[1]
output_folder       <- args[2]
pheno_file          <- args[3]
design_string       <- args[4]
add_hgnc            <- as.logical(args[5])
min_samples_filter  <- args[6]
do_sva              <- as.logical(args[7])
do_roast            <- as.logical(args[8])
do_go               <- as.logical(args[9])
do_gage             <- as.logical(args[10])


main <- function(){
  # create objects
  pheno_table   <- as_data_frame(fread(pheno_file))
  counts_table  <- create_counts_table(pheno_table)
  counts_matrix <- table_to_matrix(counts_table)
  sva_object    <- create_sva_object(do_sva, counts_matrix, pheno_table)
  model_matrix  <- create_model_matrix(pheno_table, sva_object)
  DGE_object    <- create_DGE_object(counts_matrix)
  voom_object   <- voom(DGE_object, model_matrix, save.plot = T)
  limma_object  <- create_limma_object(voom_object, model_matrix)
  limma_table   <- create_limma_table(limma_object, ncol(model_matrix))
  go_table      <- create_go_table(do_go, limma_table)
  if(do_gage) gage_tables  <- create_gage_tables(limma_table)
  if(do_roast) roast_table <- create_roast_table(counts_table, model_matrix)
  
  # write output
  create_MV_plot(voom_object)
  create_MDS_object(voom_object)
  saveRDS(limma_object, file = (str_c(output_folder, 'limma_object.rds')))
  write_count_tables(limma_table, DGE_object, counts_matrix)
  write_tsv(limma_table, str_c(output_folder, 'limma_table.tsv'))
  if(do_go) write_tsv(go_table, str_c(output_folder, 'go_table.tsv'))
  if(do_sva) saveRDS(sva_object, str_c(output_folder, 'sva_object.rds'))
  if(do_gage) map2(gage_tables, 
                  c('kegg_up_table.tsv', 'kegg_down_table.tsv', 'go_up_table.tsv', 'go_down_table.tsv'),
                  function(x,y) write_tsv(x, str_c(output_folder, y)))
  if(do_roast) write_tsv(roast_table, str_c(output_folder, 'roast_table.tsv'))
  
  
}


# counts table ---------------------------------------------------------------
create_counts_table <- function(pheno_table){
  counts_matrix <- pheno_table %>%
    use_series(fileName) %>%
    str_c(input_folder, .) %>%
    map(read_tsv, col_names = F) %>%
    map2(pheno_table$sampleName,  function(x,y) set_colnames(x, c('feature', y))) %>% 
    purrr::reduce(left_join) %>% 
    filter(!(str_detect(feature, '^_')))
}


# counts matrix ---------------------------------------------------------------
table_to_matrix <- function(table){
  table %>% 
    as.data.frame %>% 
    set_rownames(.$feature) %>% 
    dplyr::select(-feature) %>% 
    as.matrix
}


# sva object ------------------------------------------------------------------
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

# model matrix ----------------------------------------------------------------
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

# dge object ------------------------------------------------------------------
create_DGE_object <- function(counts_matrix){
  DGE_object <- counts_matrix %>% 
    DGEList %>% 
    filter_by_cpm %>% 
    calcNormFactors
}

filter_by_cpm <- function(DGE_object){
  keep <- rowSums(cpm(DGE_object) > 1) >= min_samples_filter
  DGE_object <- DGE_object[keep, , keep.lib.sizes=FALSE]
  return(DGE_object)
}

# voom object -----------------------------------------------------------------
create_voom_object <- function(counts_matrix, model_matrix){
  voom_object <- counts_matrix %>%
    create_DGE_object %>%
    voom(model_matrix, save.plot = T)
}

# limma object ----------------------------------------------------------------
create_limma_object <- function(voom_object, model_matrix){
  limma_object <- voom_object %>%  
    lmFit(model_matrix) %>%
    eBayes
}

# limma table -----------------------------------------------------------------

create_limma_table <- function(limma_object, coef){
  limma_table <- limma_object %>% 
    topTable(coef = coef, number = Inf, sort = "p") %>% 
    add_rownames('ensemble')
  if(add_hgnc) hgnc_file <-'/udd/reala/reference_files/DE/ensemble_to_hgnc_gene_table.tsv'
  if(add_hgnc) limma_table <- left_join(limma_table, read_tsv(hgnc_file))
}

# go object -------------------------------------------------------------------

create_go_table <- function(do_go, limma_table){
  if(!do_go) return(F)
  glt <- fread('/udd/reala/reference_files/DE/go_level_table.tsv') %>% 
    as_data_frame
  weight_table <- create_weight_table(limma_table)
  go_results   <- goseq(weight_table, "hg38","geneSymbol", test.cats = c("GO:CC", "GO:BP", "GO:MF")) %>% 
    inset('DB', value = 'GO') %>% 
    left_join(glt)
  kegg_results <- goseq(weight_table, "hg38","geneSymbol", test.cats = c("KEGG")) %>% 
    left_join(create_kegg_table()) %>% 
    inset('DB', value = 'KEGG')
 bind_rows(go_results, kegg_results)
}

create_weight_table <- function(limma_table){
  hgnc_table <- filter(limma_table, !is.na(hgnc_symbol))
  sig_genes <- hgnc_table  %>% 
    filter(adj.P.Val < 0.05) %>% 
    dplyr::select(hgnc_symbol) %>% 
    distinct %>% 
    use_series(hgnc_symbol)
  all_genes <- hgnc_table  %>% 
    dplyr::select(hgnc_symbol) %>% 
    distinct %>% 
    use_series(hgnc_symbol)
  weight_table <- all_genes %in% sig_genes %>% 
    set_names(all_genes) %>% 
    nullp("hg38","geneSymbol", plot.fit = F) 
}

create_kegg_table <- function(){
  lines <- readLines(
    "http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext")
  pathways <- do.call(
    rbind,
    str_split( grep( "^[ABCD]\\s+\\d{5}\\s+.*?$", lines, value=TRUE ), "\\s{2,}" )
  )
  pathways %>% 
    as.data.frame %>% 
    as_data_frame %>% 
    dplyr::select(V2, V3) %>% 
    set_colnames(c('category', 'term')) %>% 
    inset('ontology', value = '')
}

# gage table ------------------------------------------------------------------

create_gage_tables <- function(limma_table){
  data(kegg.gs)
  data(go.gs)
  entrez_table <- add_entrez_names(limma_table)
  fold_changes <- limma_table$logFC %>% 
    set_names(entrez_table$entrezgene)
  kegg_object <- gage(fold_changes, gsets = kegg.gs, ref = NULL, samp = NULL)
  kegg_up_table <- create_table(kegg_object$greater, 'pathway')
  kegg_down_table <- create_table(kegg_object$less, 'pathway') 
  go_object <- gage(fold_changes, gsets = go.gs, ref = NULL, samp = NULL)
  go_up_table <- create_table(go_object$greater, 'pathway')
  go_down_table <- create_table(go_object$less, 'pathway') 
  return(list(kegg_up_table, kegg_down_table, go_up_table, go_down_table))
}

create_table <- function(matrix, rowname){
  matrix %>% 
    as.data.frame %>% 
    rownames_to_column(rowname) %>% 
    as_data_frame
}

add_entrez_names <- function(limma_table){
  entrez_table <- limma_table %>% 
    left_join(create_entrez_table(), by = c('ensemble' = 'ensembl_gene_id')) %>% 
    filter(!is.na(entrezgene))
}

create_entrez_table <- function(){
  '/udd/reala/reference_files/DE/gene_to_entrez_table.tsv' %>% 
    fread %>% 
    as_data_frame
}

# roast table -----------------------------------------------------------------

create_roast_table <- function(counts_table, model_matrix){
  entrez_counts_matrix <- create_entrez_matrix(counts_table)
  entrez_voom_object <- entrez_counts_matrix %>%
    create_DGE_entrez_object(as.character(rownames(.))) %>% 
    voom(model_matrix, save.plot = T)
  mroast_kegg_table <- create_mroast_table(entrez_voom_object, entrez_counts_matrix, model_matrix, 'KEGG')
  mroast_go_table <- create_mroast_table(entrez_voom_object, entrez_counts_matrix, model_matrix, 'GO')
  mroast_table <- bind_rows(mroast_go_table, mroast_kegg_table)  
}

create_entrez_matrix <- function(counts_table){
  entrez_counts_matrix <- counts_table %>% 
    dplyr::rename(ensembl_gene_id = feature) %>% 
    inner_join(create_entrez_table()) %>% 
    dplyr::select(-ensembl_gene_id) %>% 
    group_by(entrezgene) %>% 
    summarise_all(sum) %>% 
    dplyr::rename(feature = entrezgene) %>% 
    table_to_matrix 
}

create_DGE_entrez_object <- function(counts_matrix, names){
  DGE_object <- counts_matrix %>% 
    DGEList(genes = names) %>% 
    filter_by_cpm %>% 
    calcNormFactors
}

create_mroast_table  <- function(entrez_voom_object, entrez_counts_matrix, model_matrix, db){
  mroast_index <- get_mroast_index(rownames(entrez_counts_matrix), entrez_voom_object$genes$genes, db)
  mroast_table <- mroast(entrez_voom_object, mroast_index, model_matrix, nrot = 99999) %>% 
    as_data_frame %>% 
    rownames_to_column('gene_set_id') %>% 
    inset('db', value = db)
  description_table <- create_description_table(mroast_table, db)
  left_join(mroast_table, description_table)
}

get_mroast_index <- function(unfiltered_gene_names, filtered_gene_names, db){
  db_list     <- create_db_list(db)
  db_vector   <- unlist(db_list)
  db_lengths  <- sapply(db_list, length)
  db_index <- match(db_vector, unfiltered_gene_names) %>% 
    split(rep(names(db_list), db_lengths)) %>% 
    lapply(function(x) x[!is.na(x)]) %>% 
    .[sapply(., length)  > 10] %>% 
    map(function(x) filtered_gene_names %in% x) %>% 
    .[map_int(., sum) > 0]
}

create_db_list <- function(db){
  if(db == 'GO') db_list <- as.list(org.Hs.egGO2EG)
  if(db == 'KEGG') {
    db_list <- as.list(org.Hs.egPATH2EG)
    names <- names(db_list) %>% 
      map(as.integer) %>% 
      map(as.character)
    names(db_list) <- names
  }
  return(db_list)
}

create_description_table <- function(mroast_table, db){
  if(db == 'GO')  description_table <- select(GO.db, keys=mroast_table$gene_set_id,
                                              columns = c("GOID","TERM"), 
                                              keytype = "GOID")
  if(db == 'KEGG') description_table <- as.data.frame(KEGGPATHID2NAME) %>% 
      rowwise %>% 
      mutate(path_id = as.character(as.integer(path_id)))
  description_table <- set_colnames(description_table, c('gene_set_id', 'description'))
  return(description_table)
}


# output ----------------------------------------------------------------------
create_MV_plot <- function(voom_object){
  png(str_c(output_folder,"MV_plot.png"), 
      width = 4, height = 4, units = "in", res = 300)
  plot(voom_object$voom.xy, pch = 16, cex = 0.25)
  title("voom: Mean-variance trend")
  lines(voom_object$voom.line, col = "red")
  dev.off()
}

create_MDS_object <- function(voom_object){
  dev.new()
  voom_object %>%
    plotMDS %>%
    saveRDS(str_c(output_folder, 'MDS_object.rds'))
  dev.off()
}

write_count_tables <- function(limma_table, DGE_object, counts_matrix){
  pval_order_table <- limma_table %>% 
    arrange(P.Value) %>% 
    dplyr::select(ensemble)
  write_normalized_count_table(DGE_object, pval_order_table)
  write_count_table(counts_matrix, pval_order_table, 'raw_count')
}

write_normalized_count_table <- function(DGE_object, pval_order_table){
  lib_size         <- with(DGE_object$samples, lib.size * norm.factors)
  norm_count_table <- DGE_object %>% 
    use_series(counts) %>% 
    add(0.5) %>% 
    divide_by(lib_size + 1) %>% 
    multiply_by(1e+06) %>% 
    write_count_table(pval_order_table, 'norm_count')
}

write_count_table <- function(counts_matrix, pval_order_table, prefix){
  count_table <- counts_matrix %>% 
    data.frame %>% 
    add_rownames('ensemble') %>% 
    left_join(pval_order_table, .) %>% 
    write_table(prefix)
}

write_table <- function(table, prefix){
  if(add_hgnc) hgnc_file <-'/udd/reala/reference_files/DE/ensemble_to_hgnc_gene_table.tsv'
  if(add_hgnc) table <- left_join(table, read_tsv(hgnc_file))
  write_tsv(table, str_c(output_folder, prefix, '_table.tsv'))
}

write_model_table <- function(model_matrix){
  model_matrix %>% 
    data.frame %>% 
    write_tsv(str_c(output_folder, 'model_table.tsv'))
}


main()
