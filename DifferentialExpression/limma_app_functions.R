
make_analyses_list <- function(){
    analyses <- list.dirs(path = "data", full.names = T, recursive = F) %>% 
        map(function(dir) list.dirs(dir, full.names = T, recursive = F)) %>% 
        flatten_chr %>% 
        data_frame %>% 
        set_colnames('folder_name') %>% 
        rowwise %>% 
        mutate(ui_name = basename(folder_name))
}

read_rds_object <- function(analysis_folder, object_string){
    object <- analysis_folder %>% 
        str_c(object_string) %>% 
        readRDS
}

# -----------------------------------------------------------------------------

get_analysis_folder <- function(analyses, analysis_name){
    analysis_folder <- analyses %>% 
        filter(ui_name == analysis_name) %>% 
        use_series(folder_name)
}

get_limma_covariates <- function(limma_object){
    limma_covariates <- limma_object %>% 
        use_series(coefficients) %>% 
        colnames %>% 
        rev
}

make_limma_table <- function(limma_object, covariate, hgnc_table){
    limma_object %>% 
        topTable(coef = covariate, number = Inf) %>%
        add_rownames("ensemble") %>% 
        left_join(hgnc_table)
}

make_limma_exon_table <- function(limma_object, covariate){
    limma_object %>% 
        topSplice(coef = covariate, test = "t", number = Inf) 
}

get_pheno_covariates <- function(analysis_folder){
    analysis_folder %>% 
        str_c('/pheno_file.tsv') %>% 
        read.table(sep = "\t", header = T, stringsAsFactors = F) %>% 
        colnames %>% 
        .[-(1:2)]
}



make_exon_count_table <- function(analysis_folder, table_name){
    count_table <- analysis_folder %>% 
        str_c(table_name) %>% 
        read_tsv %>% 
        arrange(external_gene_name) %>% 
        set_colnames(map(colnames(.), rename_sample))
}



make_exon_count_vector <- function(count_table, index){
    count_table[as.numeric(index), -(1:3)] %>% 
        unlist %>% 
        unname
}




filter_limma_table <- function(limma_table, max_pvalue){
    limma_table %>% 
        filter(adj.P.Val <= max_pvalue) %>% 
        select(ensemble, hgnc_symbol, everything()) %>% 
        arrange(adj.P.Val) %>% 
        round_limma_table
}

filter_limma_exon_table <- function(limma_table, max_FDR){
    limma_table %>% 
        filter(FDR <= max_FDR) %>% 
        arrange(FDR)
}

round_limma_table <- function(table){
    table %>% 
        mutate(logFC     = round(logFC,     3)) %>% 
        mutate(AveExpr   = round(AveExpr,   3)) %>% 
        mutate(t         = round(t,         3)) %>% 
        mutate(P.Value   = round(P.Value,   8)) %>% 
        mutate(adj.P.Val = round(adj.P.Val, 8)) %>% 
        mutate(B         = round(B,         3)) 
}


make_pheno_vector <- function(analysis_folder, covariate){
    analysis_folder %>% 
        str_c('/pheno_file.tsv') %>% 
        read.table(sep = "\t", header = T, stringsAsFactors = T) %>% 
        select_(as.character(covariate)) %>% 
        extract2(1)
}


make_exon_covariate_row <- function(model_table, covariate){
    covariate_row <- model_table %>% 
        select_(as.character(covariate))
}

make_covariate_vector <- function(limma_object, covariate){
    covariate_row <- limma_object %>% 
        use_series('design') %>% 
        data.frame %>% 
        select_(as.character(covariate)) %>% 
        extract2(1)
}

# -----------------------------------------------------------------------------

make_raw_count_table <- function(analysis_folder, table_name, norm_count_table, pval_order_table){
    raw_count_table <- make_count_table(analysis_folder, table_name, pval_order_table) %>% 
        .[.$ensemble %in% norm_count_table$ensemble, ]
}

make_count_table <- function(analysis_folder, table_name, pval_order_table){
    count_table <- analysis_folder %>% 
        str_c(table_name) %>% 
        read_tsv %>% 
        select(ensemble, hgnc_symbol, everything()) %>% 
        arrange(ensemble) %>% 
        left_join(pval_order_table) %>% 
        arrange(P.Value) %>% 
        select(-P.Value) %>% 
        set_colnames(map(colnames(.), rename_sample))
}

make_exon_raw_count_table <- function(analysis_folder, table_name, norm_count_table, pval_order_table){
    raw_count_table <- make_exon_count_table(analysis_folder, table_name, pval_order_table) %>% 
        .[.$ensembl_exon_id %in% norm_count_table$ensembl_exon_id, ]
}

make_exon_count_table <- function(analysis_folder, table_name, pval_order_table){
    count_table <- analysis_folder %>% 
        str_c(table_name) %>% 
        read_tsv %>% 
        left_join(pval_order_table, by = c('ensembl_exon_id' = 'ExonID')) %>% 
        arrange(FDR) %>% 
        select(-FDR, -ensembl_gene_id) %>% 
        set_colnames(map(colnames(.), rename_sample))
}

make_count_vector <- function(count_table, index){
    count_table[as.numeric(index), -(1:2)] %>% 
        unlist %>% 
        unname
}

rename_sample <- function(name){
    if(!str_detect(name, 'S.[:digit:]+_[:print:]+')) return(name)
    name %>% 
        str_split('_') %>% 
        extract2(1) %>% 
        extract(1)
}


# -----------------------------------------------------------------------------

create_mds_plot <- function(mds_obj, labels){
    plot_df <- data.frame('coordinate1' = unname(mds_obj$x), 
                          'coordinate2' = unname(mds_obj$y), 
                          'groups'      = labels)
    ggplot(plot_df, aes(coordinate1, coordinate2, label = labels)) +
        geom_point(aes(color = labels)) 
}

create_exon_plot <- function(limma_object, covariate, gene_index, filtered_table){
    gene_name <- filtered_table %>% 
        use_series('GeneID') %>% 
        extract(as.numeric(gene_index))
    plotSplice(limma_object, coef = covariate,  geneid = gene_name)
}

create_qq_plot <- function(p_vals){
    plot_df <- data.frame('p_vals' = p_vals)
    ggplot(plot_df, aes(sample = p_vals)) +
        stat_qq()
}

create_histogram <- function(p_vals){
    plot_df <- data.frame('p_vals' = p_vals)
    ggplot(plot_df, aes(p_vals)) +
        geom_histogram()
}


create_count_plot <- function(covariate, counts, xlab){
    if(is.factor(covariate)) 
        return(create_count_boxplot(covariate, counts, xlab))
    else
        return(create_count_scatterplot(covariate, counts, xlab))
}

create_count_boxplot <- function(factor, counts, xlab){
    plot_df <- data.frame('factor' = factor, 'counts' = counts)
    ggplot(plot_df, aes(factor, counts)) +
        geom_boxplot(aes(fill = factor(factor)), outlier.colour = NA) + 
        geom_jitter(size = 4, width = 0.2) +
        labs(x = xlab)
}

create_count_scatterplot <- function(value, counts, xlab){
    plot_df <- data.frame('value' = value, 'counts' = counts)
    ggplot(plot_df, aes(value, counts)) +
        geom_point() +
        labs(x = xlab)
}


create_MA_plot <- function(limma_table, cutoff){
    limma_table <- mutate(limma_table, sig = adj.P.Val < cutoff)
    ggplot(limma_table) + 
        geom_point(aes(AveExpr, logFC, color = sig), alpha = 1/10) +
        geom_point(data = subset(limma_table, sig == T),
                   aes(AveExpr, logFC, color = sig)) +
        xlab("log(Average Expression)") +
        ylab("log(Fold Change)")
}

create_sva_plot <- function(pheno_table, sva_object, grouping_variable, sv){
    covs <- pheno_table %>% 
        select_(as.character(grouping_variable)) %>% 
        extract2(1)
    values <- sva_object$sv[,as.numeric(sv)]
    plot_df <- data_frame(covs, values)
    if(is.character(covs)) 
        return(create_sva_boxplot(plot_df, grouping_variable))
    else
        return(create_sva_scatterplot(plot_df, grouping_variable))
}

create_sva_boxplot <- function(plot_df, xlab){
    ggplot(plot_df, aes(covs, values)) +
        geom_boxplot(aes(fill = factor(covs)), outlier.colour = NA) + 
        geom_jitter(size = 4, width = 0.2) +
        labs(x = xlab)
}

create_sva_scatterplot <- function(plot_df, xlab){
    ggplot(plot_df, aes(covs, values)) +
        geom_point() +
        labs(x = xlab)
}


create_go_table <- function(analysis, db, ont, min_genes, max_genes,  filter, go_min, go_max){
    table <- analysis %>% 
        str_c(., '/go_table.tsv') %>% 
        read_tsv
    table <- filter(table, numDEInCat >= min_genes)
    table <- filter(table, numDEInCat <= max_genes)
    if(db != 'all') table <- filter(table, DB == db)
    if(ont != 'all') table <- filter(table, ontology == ont)
    if(filter) table <- filter(table, level >= go_min)
    if(filter) table <- filter(table, level <= go_max)
    table <- select(table, -level) %>% 
        distinct
    return(table)
}

create_gage_table <- function(analysis, input_table, min_genes, max_genes){
    input_table <- analysis %>% 
        str_c(., '/', input_table, '_table.tsv') %>% 
        read_tsv %>% 
        filter(set.size >= min_genes) %>% 
        filter(set.size <= max_genes)
}

create_roast_table <- function(analysis, roast_db, min_genes, max_genes){
    table <- analysis %>% 
        str_c(., '/roast_table.tsv') %>% 
        read_tsv %>% 
        filter(NGenes >= min_genes) %>% 
        filter(NGenes <= max_genes)
    if(roast_db != 'all') table <- filter(table, db == roast_db)
    return(table)
}
