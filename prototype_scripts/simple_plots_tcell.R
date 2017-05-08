pacman::p_load(
  tidyverse,
  stringr,
  forcats
)
snps_weighted_per_state <- read_csv('out/prototype_scripts/snp_state_count_by_exp.csv')
sample_meta <- read_csv('input_data/external_static/metadata/epigenome_roadmap/chromatin_state_samples_meta.csv')

eid_counts <- snps_weighted_per_state %>% 
  filter(str_detect(state, 'EnhG?$')) %>% 
  distinct(eid, tag_snp) %>% 
  group_by(eid) %>% 
  count() %>% 
  inner_join(sample_meta) %>% 
  arrange(desc(n)) %>% 
  mutate(group= fct_reorder(group, n)) %>% 
  filter(group == 'Blood & T-cell')
eid_counts %>% 
  dplyr::select(-name,-edacc9_name) %>% 
  arrange(n)


edf <- local({
  x <- snps_weighted_per_state %>% 
    filter(str_detect(state, 'EnhG?$')) %>%
    filter(eid %in% eid_counts$eid) %>% 
    group_by(eid, tag_snp, n_proxies) %>% 
    summarise(n_overlapping = sum(n_overlapping))
  x <- sample_meta %>% 
    dplyr::rename(idname=epigenome_Mnemonic) %>% 
    distinct(eid, idname) %>% 
    inner_join(x, by = 'eid')
  # select max number of SNPs overlapped by any track per tag_snp
  # best_overlap <- x %>% 
  #   group_by(tag_snp, idname) %>% 
  #   summarise(n_overlapping = median(n_overlapping)) %>% 
  #   summarize(max_overlap = max(n_overlapping))
  x %>% 
    # inner_join(best_overlap, by = 'tag_snp') %>% 
    mutate(score = n_overlapping / n_proxies)  
})

edf %>% 
  ggplot() +
  geom_jitter(aes(eid, score, color = idname))

edf %>% 
  ggplot() +
  geom_tile(aes(tag_snp, idname, fill = score, color = 'white')) +
  scale_fill_gradient(low = 'white', high = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
