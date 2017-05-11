pacman::p_load(
  assertthat,
  tidyverse,
  stringr
)
paths <- list(
  gsTCC = list(counts = 'out/process/gsTCC/experiment_CeD-SNP_overlap.csv',
               coverage = 'out/process/gsTCC/cell_type_coverage.csv'),
  encode = list(counts = 'out/process/encode/experiment_CeD-SNP_overlap.csv',
                coverage = 'out/process/encode/cell_type_coverage.csv'),
  epigenome = list(counts = 'out/process/epigenome/Enhancer_CeD-SNP_overlap.csv',
                   coverage = 'out/process/epigenome/cell_type_Enhancer_coverage.csv',
                   meta = 'input_data/external_static/metadata/epigenome_roadmap/chromatin_state_samples_meta.csv'),
  greenleaf_2015 = list(counts = 'out/process/greenleaf_2015/experiment_CeD-SNP_overlap.csv',
                        coverage = 'out/process/greenleaf_2015/cell_type_coverage.csv')
)
#
# Combine tag SNP count files
#
paths %>% 
  unlist() %>% #use.names = FALSE) %>% 
  map_lgl(file.exists) %>% 
  all() %>% 
  assert_that()

df <- paths%>% 
  map(~ read_csv(.x$counts)) %>% 
  bind_rows() %>% 
  dplyr::select(-eid)

df_cov <- paths%>% 
  map(~ read_csv(.x$coverage)) %>% 
  bind_rows()  %>% 
  dplyr::select(-eid)

df %>% 
  write_csv('out/process/20_combined_counts.csv')
df_cov %>% 
  write_csv('out/process/20_combined_coverage.csv')
