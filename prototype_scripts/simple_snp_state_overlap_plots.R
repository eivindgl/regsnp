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
  mutate(group= fct_reorder(group, n))
group_order <- levels(eid_counts$group)
eid_counts %>% 
  ggplot() +
  geom_boxplot(aes(group, n)) +
  coord_flip() +
  labs(title = 'Number of risk loci SNPs overlapped by at least 1 enhancer',
       subtitle = 'No track coverage normalization or scoring of amount of overlap per locus')
dir.create('out/prototype_scripts/plots', recursive = TRUE, showWarnings = FALSE)
ggsave('out/prototype_scripts/plots/number_of_risk_loci_covered_by_enhancers.png')


edf <- snps_weighted_per_state %>% 
  filter(str_detect(state, 'EnhG?$')) %>%
  group_by(eid, tag_snp, n_proxies) %>% 
  summarise(n_overlapping = sum(n_overlapping))

# select max number of SNPs overlapped by any track per tag_snp
best_overlap <- edf %>% 
  group_by(tag_snp) %>% 
  summarize(max_overlap = max(n_overlapping))

edf <- sample_meta %>% 
  distinct(eid, group) %>% 
  inner_join(edf, by = 'eid') %>% 
  inner_join(best_overlap, by = 'tag_snp')
edf %>% 
  mutate(group = factor(group, levels = group_order)) %>% 
  group_by(group, tag_snp, n_proxies, max_overlap) %>% 
  summarise(n_overlapping = median(n_overlapping)) %>% 
  mutate(overlap_score = n_overlapping / max_overlap) %>% 
  ggplot() +
  geom_tile(aes(tag_snp, group, fill = overlap_score, color = 'white')) +
  scale_fill_gradient(low = 'white', high = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('out/prototype_scripts/plots/tag_snps_scored_by_category.png')
