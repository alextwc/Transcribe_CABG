source("/udd/reala/repos/prod/kallisto/write_kallisto_makefile.R")
# write_copd_kallisto_makefile('/proj/regeps/regep00/studies/COPDGene/analyses/rerpc/petepilotPrimary/pilot/', '417590015_4')
# write_copd_kallisto_makefile('/proj/regeps/regep00/studies/COPDGene/analyses/rerpc/september32015Primary/september3/', '1300881131_8')
# write_copd_kallisto_makefile('/proj/regeps/regep00/studies/COPDGene/analyses/rerpc/taffetaRNAseq20160223/february2016/', '3924469546_4', merge_files = F )



# input_folder <- '/proj/regeps/regep00/studies/COPDGene/analyses/rerpc/taffetaRNAseq20160223/february2016/'
# batch_number = '417590015_4'
# merge_files = F
# input_mode  = 'nested'
# index_file  = '/udd/reala/reference_files/kallisto/Homo_sapiens.GRCh38.rel79.cdna.all.idx'
# bootstrap   = 100
# counts_folder  = '/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/archivalroot/Counts/Kallisto_transcripts/Grch38/'

write_copd_kallisto_makefile <- function(
  input_folder, batch_number,
  index_file  = '/udd/reala/reference_files/kallisto/Homo_sapiens.GRCh38.rel79.cdna.all.idx',
  counts_folder  = '/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/archivalroot/Counts/Kallisto_transcripts/Grch38/',
  bootstrap   = 100, input_mode  = 'nested', sample_file = 'none',
  merge_files = T){
  
  
  sample_table   <- create_sample_table(input_folder, merge_files, input_mode, sample_file)
  kallisto_lines <- write_copd_kallisto_lines(sample_table, counts_folder, index_file, bootstrap)
  header_lines   <- write_copd_header_lines(sample_table, counts_folder)
  create_makefile(batch_number, header_lines, kallisto_lines)
}

create_makefile <- function(batch_number, header_lines, kallisto_lines){
  kallisto_folder <- str_c('/proj/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/archivalroot/',
                           batch_number, '/Alignment/kallisto') 
  dir.create(kallisto_folder)
  setwd(kallisto_folder)
  writeLines(c(header_lines, kallisto_lines), 'makefile')
}

write_copd_kallisto_lines <- function(sample_table, counts_folder, index_file, bootstrap){
  kallisto_lines <- sample_table %>% 
    split(.$sample_name) %>% 
    map(write_copd_kallisto_lines_by_sample, counts_folder, index_file, bootstrap) %>% 
    flatten_chr()
}

write_copd_kallisto_lines_by_sample <- function(table, counts_folder, index_file, bootstrap){
  sample_file   <- str_c(counts_folder, table$sample_name[[1]])
  output        <- str_c(sample_file, '/abundance.h5')
  input         <- map(table$folder_path, get_fastq_files) %>% 
    flatten_chr %>% 
    str_c(collapse = ' ')
  target_line   <- str_c(output, ' : ', input)
  rule_line     <- str_c('\tkallisto quant -i', index_file, '-o', sample_file, '-b', bootstrap, input, sep = ' ')
  lines         <- c(target_line, rule_line, '')
}

write_copd_header_lines <- function(sample_table, counts_folder){
  header_lines <- sample_table %>% 
    use_series(sample_name) %>% 
    unique %>% 
    str_c(counts_folder, ., '/abundance.h5') %>% 
    str_c(collapse = ' ') %>% 
    str_c('all : ', .) %>% 
    c('')
}




