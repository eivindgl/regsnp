
# Performs counting. and writes csv file. 
# Data analysis done elsewhere.
#
pacman::p_load(
  tidyverse,
  stringr,
  rtracklayer,
  GenomicRanges
)

  
# raggr webservice computes all proxy snps from a list of tag snps
# read csv file
# http://raggr.usc.edu/
read_raggr_csv <- function(path) {
  read_csv(path,
           col_types = cols(`SNP1 Chr` = col_character(),
                            `SNP2 Chr` = col_character()))
}  

raggr_df_to_genomic_ranges <- function(snps) {
  snps %>% 
    dplyr::select(chrom = `SNP2 Chr`,
                  start = `SNP2 Pos`,
                  tag_snp = `SNP1 Name`) %>% 
    mutate(end = start + 1,
           tag_snp = str_extract(tag_snp, '^[^:]+')) %>% 
    makeGRangesFromDataFrame(ignore.strand = TRUE, keep.extra.columns = TRUE)
}

total_MB_coverage <- function(gr) {
  gr %>% 
    coverage() %>% 
    sum() %>% 
    as.numeric() %>% 
    sum() %>% 
    `/`(1e6)
}

extract_name <- function(path) {
  name <- basename(path)
  str_replace(name, '.narrowPeak.gz', '')
}

read_raggr <- purrr::compose(raggr_df_to_genomic_ranges, read_raggr_csv)

snps <- read_raggr('input_data/external_static/CeD_proxy-SNPs_08_CEU.csv')
seqlevelsStyle(snps) <- 'UCSC'


paths <- list.files('input_data/external/encode_dnase', full.names = TRUE)
sample_names <- paths %>% map_chr(extract_name)
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
dfs <- paths %>% 
  map(partial(import.bed, extraCols = extraCols_narrowPeak)) %>% 
  set_names(nm = sample_names)

#
# Weighted count
# 
proxies_per_tag <- tibble(tag_snp = snps$tag_snp) %>% 
  group_by(tag_snp) %>% 
  summarize(n_proxies = length(tag_snp))

tagSNPs_per_experiment <- function(gr, experiment_name, snps, proxies_per_tag) {
  hits <- gr %>% 
    findOverlaps(snps)
  tibble(experiment_name,
         tag_snp=snps$tag_snp[to(hits)]) %>% 
    group_by(experiment_name, tag_snp) %>% 
    summarize(n_overlapping = length(tag_snp)) %>% 
    ungroup() %>% 
    inner_join(proxies_per_tag, by = 'tag_snp')
}

snps_weighted_per_exp <-  local({
  exp_names <- names(dfs)
  f <- partial(tagSNPs_per_experiment, 
               snps = snps, proxies_per_tag = proxies_per_tag)
  map2_df(dfs, exp_names, f) 
})

dir.create('out/process/encode', recursive = TRUE, showWarnings = FALSE)
snps_weighted_per_exp %>% 
  write_csv('out/process/encode/experiment_CeD-SNP_overlap.csv')

dnase_exp_cov <- map_dbl(dfs, total_MB_coverage)
tibble(experiment_name = names(dnase_exp_cov), coverage_MB = dnase_exp_cov) %>% 
  write_csv('out/process/encode/cell_type_coverage.csv')
