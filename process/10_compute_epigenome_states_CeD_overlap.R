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

coverage_by_state <- function(gr) {
  xs <- gr %>% 
    split(gr$name) %>% 
    as.list() %>% 
    map_dbl(total_MB_coverage)
   tibble(state = names(xs),
          coverage_MB = xs) %>% 
     arrange(coverage_MB)
}


read_raggr <- purrr::compose(raggr_df_to_genomic_ranges, read_raggr_csv)

snps <- read_raggr('input_data/external_static/CeD_proxy-SNPs_08_CEU.csv')
seqlevelsStyle(snps) <- 'UCSC'

states_desc <- read_csv('input_data/external_static/metadata/epigenome_roadmap/chromatin_states_description.csv')

extract_name <- function(path) {
  name = basename(path)
  str_extract(name, '[^_]+')
}

paths <- list.files('input_data/external/epigenome_roadmap/states', full.names = TRUE)
sample_names <- paths %>% map_chr(extract_name)
dfs <- paths %>% 
  map(import.bed) %>% 
  set_names(nm = sample_names)

SNPs_per_state <- function(x, snps) {
  x %>% 
    split(x$name) %>% 
    countOverlaps(snps)
}

# Weighted count
# 
proxies_per_tag <- tibble(tag_snp = snps$tag_snp) %>% 
  group_by(tag_snp) %>% 
  summarize(n_proxies = length(tag_snp))

tagSNPs_per_state <- function(gr, experiment_name, snps, proxies_per_tag) {
  by_state <- gr %>% 
    split(gr$name)
  hits <- by_state %>% 
    findOverlaps(snps)
  tibble(state=names(by_state)[from(hits)],
         tag_snp=snps$tag_snp[to(hits)]) %>% 
    group_by(state, tag_snp) %>% 
    summarize(n_overlapping = length(tag_snp)) %>% 
    ungroup() %>% 
    inner_join(proxies_per_tag, by = 'tag_snp') %>% 
    mutate(eid = experiment_name) %>% 
    dplyr::select(eid, everything())
}

snps_weighted_per_state <-  local({
  exp_names <- names(dfs)
  f <- partial(tagSNPs_per_state, 
               snps = snps, proxies_per_tag = proxies_per_tag)
  map2_df(dfs, exp_names, f) 
})
sample_meta <- read_csv('input_data/external_static/metadata/epigenome_roadmap/chromatin_state_samples_meta.csv')
snps_weighted_per_state
sdf <- sample_meta %>% 
  inner_join(snps_weighted_per_state) %>% 
  mutate(source = 'epigenome') %>% 
  dplyr::select(sample=epigenome_Mnemonic, group, n_overlapping, n_proxies, tag_snp, state, source, eid)
dir.create('out/process/epigenome', recursive = TRUE, showWarnings = FALSE)

filter_enh <- function(x) {
  x %>% 
    filter(state == '7_Enh')
  
}

sdf <- sdf %>% 
  mutate(cell_category = ifelse(group == 'Blood & T-cell', 'T-cell', 'unspecified'))

sdf %>% 
  write_csv('out/process/epigenome/state_CeD-SNP_overlap.csv')
x <- sdf %>% 
  filter_enh() %>% 
  group_by(sample, tag_snp) %>% 
  summarize(n_overlapping = sum(n_overlapping))
sdf %>% 
  dplyr::select(-state, -n_overlapping) %>% 
  distinct() %>% 
  inner_join(x) %>% 
  write_csv('out/process/epigenome/Enhancer_CeD-SNP_overlap.csv')

covdf <- map2_df(dfs, names(dfs), 
        ~ .x %>% 
          coverage_by_state %>% 
          mutate(eid = .y) %>% 
          dplyr::select(eid, everything()))
covdf <- covdf %>% 
  inner_join(sample_meta) %>% 
  mutate(source = 'epigenome') %>% 
  dplyr::select(sample=epigenome_Mnemonic, coverage_MB, group, state, source, eid) %>% 
  mutate(cell_category = ifelse(group == 'Blood & T-cell', 'T-cell', 'unspecified'))

covdf %>%  
  write_csv('out/process/epigenome/cell_type_state_coverage.csv')

enh_covdf <- covdf %>% 
  filter_enh() %>% 
  group_by(sample, group, source, eid, cell_category) %>% 
  summarise(coverage_MB = sum(coverage_MB))
enh_covdf %>% 
  write_csv('out/process/epigenome/cell_type_Enhancer_coverage.csv')
