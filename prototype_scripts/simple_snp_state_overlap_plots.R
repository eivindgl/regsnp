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


edf <- local({
  x <- snps_weighted_per_state %>% 
    filter(str_detect(state, 'EnhG?$')) %>%
    group_by(eid, tag_snp, n_proxies) %>% 
    summarise(n_overlapping = sum(n_overlapping))
  x <- sample_meta %>% 
    distinct(eid, group) %>% 
    inner_join(x, by = 'eid')
  # select max number of SNPs overlapped by any track per tag_snp
  best_overlap <- x %>% 
    group_by(tag_snp, group) %>% 
    summarise(n_overlapping = median(n_overlapping)) %>% 
    summarize(max_overlap = max(n_overlapping))
  x %>% 
    inner_join(best_overlap, by = 'tag_snp') %>% 
    mutate(score = n_overlapping / max_overlap,
           group = factor(group, levels = group_order)) %>% 
    dplyr::select(-max_overlap)
})


edf %>% 
  group_by(group, tag_snp) %>% 
  summarise(score = median(score)) %>% 
  ggplot() +
  geom_tile(aes(tag_snp, group, fill = score, color = 'white')) +
  scale_fill_gradient(low = 'white', high = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('out/prototype_scripts/plots/tag_snps_scored_by_category.png')
ggsave('out/prototype_scripts/plots/tag_snps_scored_by_category.pdf') # pdf so I can copy SNP ids


# One group has score 1.0 per tag SNP.
# Drop Blood & T-cell group and select SNPs with max score below e.g. 0.5 (because then T-cell had 1.0)
tag_snps <- read_tsv('input_data/external_static/CeD_tag_SNPs.bed', col_names = c('chrom', 'start', 'end', 'tag_snp'))
#except T-cells
top_group_per_tag_snp <- edf %>% 
  group_by(group, tag_snp) %>% 
  summarise(score = median(score)) %>% 
  filter(!str_detect(group, 'T-cell')) %>% 
#  filter(!str_detect(group, 'Thymus')) %>% # did not improve overlap with ASE effects
  ungroup() %>% 
  group_by(tag_snp) %>% 
  summarise(n_over_05 = sum(score > 0.5)) %>% 
  filter(n_over_05 < 2) %>%
  dplyr::select(tag_snp) 
  # slice(which.max(score)) %>% 
  # ungroup()

# It is potentially very interesting that T-cell specific tag_snps coincide
# with the strongest ASE effects.
# However, it is perhaps not too surprising that there is some overlap.
#   - We have in total ~60 tag snps
#   - We pick up ASE effects in 24 cis genes of 13 different tag_snps.
#   - We find 12 tag_snps to be highly t-cell specific.
#   - 5 tag snps overlap covering 9 genes.
#     - 5 which are clearly heterozygous. 
#     - 3 out of 4 remaining are obviously interesting.


tcell_snps <- top_group_per_tag_snp %>% 
  inner_join(tag_snps, by = 'tag_snp') %>% 
  dplyr::select(chrom, start, end, everything())
tcell_snps
tcell_snps %>% 
  write_tsv('out/prototype_scripts/T-cell_specific_tag_snps.bed', col_names = FALSE)
tcell_ase <- read_csv('input_data/external_static/potential_ASE_by_tag_snp.csv')
tcell_ase %>% 
  inner_join(tcell_snps)
tcell_snps
tcell_ase %>% 
  mutate(ratio  = heterozygous_tag_SNP / (heterozygous_tag_SNP + homozygous_tag_SNP)) %>% 
  print(n = Inf)

