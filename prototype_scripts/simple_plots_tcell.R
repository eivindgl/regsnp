pacman::p_load(
  tidyverse,
  stringr,
  forcats
)
snps_weighted <- read_csv('out/process/20_combined_counts.csv')
cell_coverage <- read_csv('out/process/20_combined_coverage.csv')

proxy_counts <- snps_weighted %>% 
  group_by(sample) %>% 
  summarise(nSNPs = sum(n_overlapping)) %>% 
  arrange(desc(nSNPs))

tag_counts <- snps_weighted %>% 
  distinct(sample, tag_snp) %>% 
  count(sample, sort = TRUE) %>% 
  dplyr::rename(nTagSNPs = n)

counts_and_cov <- cell_coverage %>% 
  inner_join(tag_counts) %>% 
  mutate(tag_per_mb = nTagSNPs / coverage_MB)

proxy_counts_and_cov <- cell_coverage %>% 
  inner_join(proxy_counts)

counts_and_cov %>% 
  ggplot(aes(coverage_MB, nTagSNPs)) +
  geom_jitter(aes(color = source), size = 3)

#
# Plot number of proxy snippe overlapped by track by track size.
# Regress on T-cells and compute R^2 for noverlaps given T-cell genome coverage.
#
pm <- proxy_counts_and_cov %>% 
  filter(cell_category == 'T-cell') %>% 
  lm(nSNPs ~ coverage_MB, data = .)

p <- proxy_counts_and_cov %>% 
  #filter(source != 'epigenome') %>% 
  # filter(!str_detect(group, '[TB]-cell')) %>%
  # mutate(is_tcell = str_detect(group, 'T-cell')) %>% 
  ggplot(aes(coverage_MB, nSNPs)) +
  geom_point(aes(shape = source, color = cell_category), size = 2) +
  geom_smooth(method = 'lm', se = FALSE, data = filter(proxy_counts_and_cov, cell_category == 'T-cell'))
p + 
  geom_text(aes(50, 210), check_overlap = TRUE, label = sprintf('R^2 is %.3f', summary(pm)$r.square)) +
  labs(title = 'CeD proxy SNP overlap by track coverage',
       x = 'Coverage (MB)',
       y = '#CeD risk SNPs')
  ggsave('out/prototype_scripts/proxy_SNPs_overlap_by_coverage.png')

#
# gsTCC display a good ratio of SNPs vs genome coverage in the plot above,
# but it looks to be explained by a simple regression nSNPs by coverage on t-cells.
#

# tdf <- snps_weighted %>% 
#   filter(cell_category == 'T-cell') %>% 
#   filter(source == 'encode')
# 
# tdf <- tdf %>% group_by(tag_snp) %>% 
#   summarise(max_overlapping = max(n_overlapping)) %>% 
#   inner_join(tdf) %>% 
#   mutate(score = if_else(max_overlapping == 0, 0, n_overlapping / max_overlapping))
# 
# tdf %>% 
#   ggplot() +
#   geom_tile(aes(tag_snp, sample, fill = score, color = 'white')) +
#   scale_fill_gradient(low = 'white', high = 'steelblue') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
