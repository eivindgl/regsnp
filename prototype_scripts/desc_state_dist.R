# How constant is state coverage between samples?
# It looks to be "quite" stable. E.g. coverage of Enh sites is pretty stable (Enh [30MB to 120MB])
#
pacman::p_load(
  tidyverse,
  stringr
)

read_bedfile <- function(path) {
}

extract_name <- function(path) {
  name = basename(path)
  str_extract(name, '[^_]+')
}

compute_bed_coverage <- function(path) {
  read_tsv(path, col_names = c('chrom', 'start', 'end', 'state')) %>% 
  mutate(size = end - start) %>% 
  group_by(state) %>%
  summarize(coverage = sum(as.numeric(size)) / 1e6) %>% 
  mutate(sample = extract_name(path))
}
  
paths <- list.files('input_data/external/epigenome_roadmap/states', full.names = TRUE)
sample_names <- paths %>% map_chr(extract_name)

df <- paths %>% map_df(compute_bed_coverage)
df <- map2_df(
  map(paths, extract_name),
  dfs,
  compute_bed_coverage
)

df %>% 
  ggplot(aes(state, log10(coverage + 1))) +
  geom_boxplot()
  
df %>% 
  filter(str_detect(state, 'Enh')) %>% 
  group_by(state) %>% 
  summarise(
    min = min(coverage),
    quantile_1 = quantile(coverage, 0.25),
    median = median(coverage), 
    quantile_3 = quantile(coverage, 0.75),
    max = max(coverage)
  )
    
  
