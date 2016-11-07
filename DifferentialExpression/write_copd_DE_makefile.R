library(stringr)
library(dplyr)
library(purrr)
library(readr)
library(magrittr)
library(testthat)



write_copd_makefile <- function(
  # parameters
  batch_name, 
  conversion_file,
  kallisto_bootstrap       = 100,
  kallisto_index_file      = '/udd/reala/reference_files/kallisto/Homo_sapiens.GRCh38.rel79.cdna.all.idx',
  script_dir               = '/udd/reala/repos/prod/Differential_expression/',
  gtf_file                 = '/udd/reala/reference_files/DE/Homo_sapiens.GRCh38.81.gtf',
  gff_file                 = '/udd/reala/reference_files/DE/Homo_sapiens.GRCh38.81.gff'){
  
  # folders
  counts_folder      <- '/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/QCed_data/Counts/'
  bam_folder         <- '/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/QCed_data/BAM_files/BAM_files/'
  sorted_bam_folder  <- '/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/QCed_data/BAM_files/sorted_BAM_files/'
  abundances_folder  <- '/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/QCed_data/Counts/Kallisto_transcripts/Grch38/'
  input_bam_folder   <- str_c('/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/archivalroot/', batch_name, '/Alignment/bam_files/')
  input_fastq_folder <- str_c('/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/archivalroot/', batch_name, '/Alignment/fastq_files/')
  makefile_folder    <- str_c('/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/QCed_data/Makefiles/', batch_name, '/')
  
  conversion_table       <- make_conversion_table(conversion_file, batch_name)
  
  # makefile line
  convert_bam_file_lines <- make_convert_bam_file_lines(conversion_table, input_bam_folder, bam_folder)
  sort_bam_file_lines    <- make_lines_by_row(conversion_table, bam_folder, sorted_bam_folder, 'samout/', make_sort_bam_file_lines_by_sample)
  Rsubread_count_lines   <- make_rsubread_count_lines(conversion_table, bam_folder, counts_folder, gtf_file, script_dir)
  HTSeq_count_lines      <- make_HTSeq_count_lines(conversion_table, bam_folder, sorted_bam_folder, counts_folder, gtf_file, gff_file, script_dir)
  autosomal_count_lines  <- make_autosomal_count_lines(conversion_table, counts_folder, gtf_file, gff_file, script_dir)
  kallisto_lines         <- make_kallisto_lines_by_row(conversion_table, input_fastq_folder, abundances_folder, 'kallisto_out/', kallisto_bootstrap, kallisto_index_file)
  header_lines           <- make_header_lines(counts_folder, conversion_table)
  
  create_directories(conversion_table, makefile_folder, abundances_folder, batch_name)
  
  writeLines(c(header_lines, autosomal_count_lines, Rsubread_count_lines, HTSeq_count_lines, sort_bam_file_lines, kallisto_lines, convert_bam_file_lines), 
             'makefile')



}

create_directories <- function(conversion_table, makefile_folder, abundances_folder, batch_name){
  conversion_table$COPDgeneID %>% 
    str_c(abundances_folder, ., '_', batch_name) %>% 
    map(dir.create)
  dir.create(makefile_folder) 
  setwd(makefile_folder)
  dir.create('rout')
  dir.create('samout')
  dir.create('htseq_out')
  dir.create('kallisto_out')
}

# conversion table ------------------------------------------------------------

make_conversion_table <- function(conversion_file, batch_name){
  conversion_file %>% 
    read.table(sep = "\t", header = T, stringsAsFactors = F) %>% 
    filter(action != 'drop') %>% 
    rowwise %>% 
    mutate(BAM_name = str_c(COPDgeneID, '_', batch_name, '.bam')) %>% 
    mutate(Count_name = str_c(COPDgeneID, '_', batch_name, '_counts.tsv')) %>% 
    mutate(fastq_files = get_fastq_names(BAM_files)) %>% 
    mutate(abundance_name = str_c(COPDgeneID, '_', batch_name, '/abundance.h5')) %>% 
    mutate(kallisto_folder = str_c(COPDgeneID, '_', batch_name))
}

get_fastq_names <- function(bam_names){
  bam_names %>% 
    str_split(';') %>% 
    extract2(1) %>% 
    map(get_fastq_names_per_bam) %>% 
    str_c(collapse = ';')
  
}

get_fastq_names_per_bam <- function(bam_name){
  bam_name %>% 
    str_match('([:print:]+)_accepted_hits.sorted.bam') %>% 
    extract2(2) %>% 
    str_c(c('_R1.fastq', '_R2.fastq'), collapse = ';')
}

# sam tools -------------------------------------------------------------------

make_convert_bam_file_lines <- function(conversion_table, input_folder, output_folder){
  link_bam_file_lines    <- make_lines_by_row(filter(conversion_table, action == 'link'), input_folder, output_folder, 'samout/', make_link_bam_file_lines_by_sample)
  merge_bam_file_lines   <- make_lines_by_row(filter(conversion_table, action == 'merge'), input_folder, output_folder, 'samout/', make_merge_bam_file_lines_by_sample)
  convert_bam_file_lines <- c(link_bam_file_lines, merge_bam_file_lines, '')
}

make_lines_by_row <-  function(table, input_folder, output_folder, error_folder, fun){
  table %>% 
    split(.$COPDgeneID) %>% 
    map(fun, input_folder, output_folder, error_folder) %>% 
    flatten_chr
}

make_link_bam_file_lines_by_sample <- function(table, input_folder, output_folder, error_folder){
  input  <- str_c(input_folder, table$BAM_files)
  output <- str_c(output_folder, table$BAM_name)
  command <- str_c('ln -s', input, output, sep = ' ')
  link_bam_file_lines <- create_make_rule_lines(input, output, command)
}

make_merge_bam_file_lines_by_sample <- function(table, input_folder, output_folder, error_folder){
  input  <- table$BAM_files %>% 
    str_split(';') %>% 
    unlist %>% 
    str_c(input_folder, .) %>% 
    str_c(collapse = ' ')
  output <- str_c(output_folder, table$BAM_name)
  error <- str_c(error_folder, 'merge_', basename(output), '_error.txt')
  command <- str_c('samtools merge', output, input, sep = " ")
  merge_bam_file_lines <- create_make_rule_lines(input, output, command)
}

make_sort_bam_file_lines_by_sample <- function(table, input_folder, output_folder, error_folder){
  input  <- str_c(input_folder, table$BAM_name)
  output <- str_c(output_folder, table$BAM_name)
  output_prefix <- output %>% 
    str_split('.bam') %>% 
    extract2(1) %>% 
    extract(1)
  command <- str_c('samtools sort -n', input, output_prefix, sep = " ")
  sort_bam_file_lines <- create_make_rule_lines(input, output, command)
}

test_samtools <- function(){
  expect_equal(make_sort_bam_file_lines_by_sample(data_frame('BAM_name' = 'one.bam'),
                                                'input_folder/', 'output_folder/', 'error_folder/'),
               c('output_folder/one.bam : input_folder/one.bam', 
                 '\tsamtools sort -n input_folder/one.bam output_folder/one', ''))
  
  expect_equal(make_merge_bam_file_lines_by_sample(data_frame('BAM_files' = 'input1.bam;input2.bam', 'BAM_name' = 'output.bam'),
                                                  'input_folder/', 'output_folder/', 'error_folder/'),
               c('output_folder/output.bam : input_folder/input1.bam input_folder/input2.bam', 
                 '\tsamtools merge output_folder/output.bam input_folder/input1.bam input_folder/input2.bam',
                 ''))
}

# Rsubread_count lines --------------------------------------------------------

make_rsubread_count_lines <- function(conversion_table, bam_folder, counts_folder, gtf_file, script_dir){
  Rsubread_gene_count_lines          <- make_R_lines_by_row(conversion_table, bam_folder, str_c(counts_folder, 'Rsubread_gene/Grch38/'), make_rsubread_count_lines_by_sample, c(gtf_file, 'gene_id', 'F'), script_dir, 'rout/', 'gene_')
  Rsubread_exon_count_lines          <- make_R_lines_by_row(conversion_table, bam_folder, str_c(counts_folder, 'Rsubread_exon/Grch38/'), make_rsubread_count_lines_by_sample, c(gtf_file, 'exon_id', 'F'), script_dir, 'rout/', 'exon_')
  Rsubread_multimap_gene_count_lines <- make_R_lines_by_row(conversion_table, bam_folder, str_c(counts_folder, 'Rsubread_gene/Grch38_multi/'), make_rsubread_count_lines_by_sample, c(gtf_file, 'gene_id', 'T'), script_dir, 'rout/', 'gene_multi_')
  Rsubread_multimap_exon_count_lines <- make_R_lines_by_row(conversion_table, bam_folder, str_c(counts_folder, 'Rsubread_exon/Grch38_multi/'), make_rsubread_count_lines_by_sample, c(gtf_file, 'exon_id', 'T'), script_dir, 'rout/', 'exon_multi_')
  Rsubread_count_lines <- c(Rsubread_gene_count_lines, Rsubread_exon_count_lines, Rsubread_multimap_gene_count_lines, Rsubread_multimap_exon_count_lines)
}

make_R_lines_by_row <-  function(table, input_folder, output_folder, fun, args = '', script_dir, rout = 'rout/', rout_mod = ''){
  table %>% 
    split(.$COPDgeneID) %>% 
    map(fun, input_folder, output_folder, args, script_dir, rout, rout_mod) %>% 
    flatten_chr
}

make_rsubread_count_lines_by_sample <- function(table, input_folder, output_folder, args, script_dir, rout, rout_mod = ''){
  input   <- str_c(input_folder, table$BAM_name)
  output  <- str_c(output_folder, table$Count_name)
  command <- create_R_command(input, output, args, script_dir, 'Rsubread_counts.R', rout, rout_mod)
  rsubread_count_lines <- create_make_rule_lines(input, output, command, request = '32')
}

# HTSeq_count lines -----------------------------------------------------------

make_HTSeq_count_lines <- function(conversion_table, bam_folder, sorted_bam_folder, counts_folder, gtf_file, gff_file, script_dir){
  HTSeq_gene_count_lines <- make_HTSeq_lines_by_row(conversion_table, sorted_bam_folder, str_c(counts_folder, 'HTSeq_gene/Grch38/'), 'htseq_out/', make_HTSeq_count_lines_by_sample, script_dir, '/htseq_gene.sh', gtf_file)
  HTSeq_exon_count_lines <- make_HTSeq_lines_by_row(conversion_table, bam_folder, str_c(counts_folder, 'HTSeq_exon/Grch38/'), 'htseq_out/', make_HTSeq_count_lines_by_sample, script_dir, '/htseq_exon.sh', gff_file)
  HTSeq_count_lines <- c(HTSeq_gene_count_lines, HTSeq_exon_count_lines)
}

make_HTSeq_lines_by_row <- function(table, input_folder, output_folder, error_folder, fun, script_dir, script, genome_file){
  table %>% 
    split(.$COPDgeneID) %>% 
    map(fun, input_folder, output_folder, error_folder, script_dir, script, genome_file) %>% 
    flatten_chr
}

make_HTSeq_count_lines_by_sample <- function(table, input_folder, output_folder, error_folder, script_dir, script, genome_file){
  input   <- str_c(input_folder, table$BAM_name)
  output  <- str_c(output_folder, table$Count_name)
  error   <- str_c(error_folder, script, "_", basename(output), '_error.txt')
  command <- str_c(str_c(script_dir, script), input, genome_file, output, '&>', error, sep = " ")
  HTSeq_count_lines <- create_make_rule_lines(input, output, command)
}

# Autosome count lines --------------------------------------------------------

make_autosomal_count_lines <- function(conversion_table, counts_folder, gtf_file, gff_file, script_dir){
  autosomal_Rsubread_gene_count_lines          <- make_R_lines_by_row(conversion_table, str_c(counts_folder, 'Rsubread_gene/Grch38/'),       str_c(counts_folder, 'Rsubread_gene/Grch38_autosome/'), make_autosomal_count_lines_by_sample, 'rsubread_gene', script_dir, 'rout/', 'rsubread_gene_')
  autosomal_Rsubread_exon_count_lines          <- make_R_lines_by_row(conversion_table, str_c(counts_folder, 'Rsubread_exon/Grch38/'),       str_c(counts_folder, 'Rsubread_exon/Grch38_autosome/'), make_autosomal_count_lines_by_sample, 'rsubread_exon', script_dir, 'rout/', 'rsubread_exon_')
  autosomal_Rsubread_multimap_gene_count_lines <- make_R_lines_by_row(conversion_table, str_c(counts_folder, 'Rsubread_gene/Grch38_multi/'), str_c(counts_folder, 'Rsubread_gene/Grch38_multi_autosome/'), make_autosomal_count_lines_by_sample, 'rsubread_gene', script_dir, 'rout/', 'rsubread_gene_multi_')
  autosomal_Rsubread_multimap_exon_count_lines <- make_R_lines_by_row(conversion_table, str_c(counts_folder, 'Rsubread_exon/Grch38_multi/'), str_c(counts_folder, 'Rsubread_exon/Grch38_multi_autosome/'), make_autosomal_count_lines_by_sample, 'rsubread_exon', script_dir, 'rout/', 'rsubread_exon_multi_')
  autosomal_HTSeq_gene_count_lines             <- make_R_lines_by_row(conversion_table, str_c(counts_folder, 'HTSeq_gene/Grch38/'),          str_c(counts_folder, 'HTSeq_gene/Grch38_autosome/'), make_autosomal_count_lines_by_sample, 'htseq_gene', script_dir, 'rout/', 'htseq_gene_')
  autosomal_HTSeq_exon_count_lines             <- make_R_lines_by_row(conversion_table, str_c(counts_folder, 'HTSeq_exon/Grch38/'),          str_c(counts_folder, 'HTSeq_exon/Grch38_autosome/'), make_autosomal_count_lines_by_sample, 'htseq_exon', script_dir, 'rout/', 'htseq_exon_')
  autosomal_count_lines <- c(autosomal_Rsubread_gene_count_lines, autosomal_Rsubread_exon_count_lines, autosomal_Rsubread_multimap_gene_count_lines, autosomal_Rsubread_multimap_exon_count_lines, autosomal_HTSeq_gene_count_lines, autosomal_HTSeq_exon_count_lines)
}

make_autosomal_count_lines_by_sample <- function(table, input_folder, output_folder, args, script_dir, rout, rout_mod = ''){
  input   <- str_c(input_folder, table$Count_name)
  output  <- str_c(output_folder, table$Count_name)
  command <- create_R_command(input, output, args, script_dir, 'remove_non_autosomes.R', rout, rout_mod)
  autosome_count_lines <- create_make_rule_lines(input, output, command)
}

# Kallisto lines --------------------------------------------------------------

make_kallisto_lines_by_row <-  function(table, input_folder, output_folder, error_folder, bootstrap, index_file){
  table %>% 
    split(.$COPDgeneID) %>% 
    map(make_kallisto_lines_by_sample, input_folder, output_folder, error_folder, bootstrap, index_file) %>% 
    flatten_chr
}

make_kallisto_lines_by_sample <- function(table, input_folder, output_folder, error_folder, bootstrap, index_file){
  input  <- table$fastq_files %>% 
    str_split(';') %>% 
    unlist %>% 
    str_c(input_folder, .) %>% 
    str_c(collapse = ' ') 
  output1 <- str_c(output_folder, table$abundance_name)
  output2 <- str_c(output_folder, table$kallisto_folder)
  error <- str_c(error_folder, 'kallisto_', basename(output2), '_error.txt')
  command <- str_c('kallisto quant -i', index_file, '-o', output2, '-b', bootstrap, input, '&>', error, sep = ' ')
  kallisto_lines <- create_make_rule_lines(input, output1, command)
}

# Header lines ----------------------------------------------------------------

make_header_lines <- function(counts_folder, conversion_table){
  header <- str_c('all :', 
                  str_c(counts_folder, 'Rsubread_gene/Grch38_autosome/', conversion_table$Count_name, collapse = ' '),
                  str_c(counts_folder, 'Rsubread_exon/Grch38_autosome/', conversion_table$Count_name, collapse = ' '),
                  str_c(counts_folder, 'Rsubread_gene/Grch38_multi_autosome/', conversion_table$Count_name, collapse = ' '),
                  str_c(counts_folder, 'Rsubread_exon/Grch38_multi_autosome/', conversion_table$Count_name, collapse = ' '),
                  str_c(counts_folder, 'HTSeq_gene/Grch38_autosome/', conversion_table$Count_name, collapse = ' '),
                  str_c(counts_folder, 'HTSeq_exon/Grch38_autosome/', conversion_table$Count_name, collapse = ' '),
                  str_c(counts_folder, 'Kallisto_transcripts/Grch38/', conversion_table$abundance_name, collapse = ' '),
                  sep = ' ')
  header_lines <- c(header, '')
}

# misc ------------------------------------------------------------------------

create_R_command <- function(input, output, args, script_dir, script, rout_dir, rout_mod = ''){
  command <- str_c("R CMD BATCH '--args", input, output, 
                   str_c(args, collapse = ' '),
                   "'", 
                   str_c(script_dir, script),
                   str_c(rout_dir, script, '_', rout_mod, basename(input), '.Rout'), sep = ' ')
}

test_create_R_command <- function(){
  expect_equal(create_R_command('input.tsv', 'output.tsv', 'arg1', 'script_dir/', 'script.R', 'rout_dir/'),
               "R CMD BATCH '--args input.tsv output.tsv arg1 ' script_dir/script.R rout_dir/script.R_input.tsv.Rout")
  expect_equal(create_R_command('input.tsv', 'output.tsv', 'arg1', 'script_dir/', 'script.R', 'rout_dir/', 'option'),
               "R CMD BATCH '--args input.tsv output.tsv arg1 ' script_dir/script.R rout_dir/script.R_optioninput.tsv.Rout")
  
}

create_make_rule_lines <- function(input, output, command, request = ''){
  target_line <- str_c(output, ' : ', input)
  if (request == '') recipe_line <- str_c('\t', command)
  else recipe_line <- str_c('\tSGE_RREQ="-l lx,virtual_free=', request, '" ', command)
  make_rule_lines <- c(target_line, recipe_line, '')
}

test_create_make_rule_lines <- function(){
  expect_equal(create_make_rule_lines('input.tsv', 'output.tsv', 'script.sh input.tsv > output,tsv'),
               c('output.tsv : input.tsv', '\tscript.sh input.tsv > output,tsv', ''))
  expect_equal(create_make_rule_lines('input.tsv', 'output.tsv', 'script.sh input.tsv > output,tsv', '32G'),
               c('output.tsv : input.tsv', '\tSGE_RREQ="-l lx,virtual_free=32G" script.sh input.tsv > output,tsv', ''))
  
}

# run tests -------------------------------------------------------------------

test_samtools()
test_create_R_command()
test_create_make_rule_lines()
