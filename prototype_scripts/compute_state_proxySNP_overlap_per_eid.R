# Performs counting. and writes csv file. 
# Data analysis done elsewhere.
#
pacman::p_load(
  tidyverse,
  stringr,
  rtracklayer,
  GenomicRanges
)
e001 <- import.bed('input_data/external/epigenome_roadmap/states/E001_15_coreMarks_dense.bed.gz')
e129 <- import.bed('input_data/external/epigenome_roadmap/states/E129_15_coreMarks_dense.bed.gz')

  
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

str_detect_any_of <- function(needles, haystack) {
  haystack <- tolower(haystack)
  any(needles %>% map_lgl(~ str_detect(haystack, .)))
}
states_desc <- read_csv('input_data/external_static/metadata/epigenome_roadmap/chromatin_states_description.csv')
states_desc

is_active_chrom <- partial(str_detect_any_of, c('tss', 'enh', 'transcr'))
open_states <- states_desc %>% 
  filter(map_lgl(description, is_active_chrom)) %>% 
  filter(!str_detect(description, 'Weak transcription')) %>% 
  mutate(id_str = str_c(state_number, name, sep = '_')) %>% 
  `[[`('id_str')
e001 %>%
  subset(e001$name %in% open_states) %>% 
  total_MB_coverage()
e129 %>%
  subset(e129$name %in% open_states) %>% 
  total_MB_coverage()

e001 %>% 
  coverage_by_state()

e129 %>% 
  coverage_by_state()

  
x <- e001 #%>% 
  head(1000)
x %>% 
  split(x$name) %>% 
  countOverlaps(snps) %>% 
  sum()

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

# compute snps overlap with different state types.
# does not take SNP LD into account.
# so tag SNPs with lots of proxy SNPs will
# contribute more to total score than other sites.
# perhaps divide SNP score by number of proxy SNPs?
snps_per_state <-  local({
  f <- partial(SNPs_per_state, snps = snps)
  x <- dfs %>% 
    map(f) 
  x %>% 
    as_tibble() %>% 
    mutate(state = names(x[[1]]),
           state = str_replace(state, '^\\d+_', '')) %>% 
    gather(sample, nSNPs, -state) %>% 
    spread(state, nSNPs)
})

#
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
    dplyr::select(experiment_name, everything())
}

snps_weighted_per_state <-  local({
  exp_names <- names(dfs)
  f <- partial(tagSNPs_per_state, 
               snps = snps, proxies_per_tag = proxies_per_tag)
  map2_df(dfs, exp_names, f) 
})

dir.create('out/prototype_scripts', recursive = TRUE, showWarnings = FALSE)
snps_weighted_per_state %>% 
  write_csv('out/prototype_scripts/snp_state_count_by_exp.csv')

