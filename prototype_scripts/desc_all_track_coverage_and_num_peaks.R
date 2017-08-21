pacman::p_load(
  tidyverse,
  stringr,
  rtracklayer,
  GenomicRanges
)

peaks_per_sample <- function(grs) {
  grs %>% 
    map_int(length)
}

coverage_per_sample <- function(grs) {
  grs %>% 
    map_dbl(compose(sum, as.numeric, width))
}

median_coverage_per_sample <- function(grs) {
  grs %>% 
    map_dbl(compose(median, width))
}

common_df <- function(sample_names, dfs, project) {
  tibble(
    project = project,
    sample = sample_names,
    num_peaks = peaks_per_sample(dfs),
    covMB = coverage_per_sample(dfs) / 1e6,
    covMedian = median_coverage_per_sample(dfs)
  )
}

gsTCC_meta <- function(paths){
  extract_name <- function(path) {
    name <- basename(path)
    str_replace(name, '.bed', '')
  }
  sample_names<- paths %>% map_chr(extract_name)
  
  dfs <- paths %>% 
    map(import.bed) %>% 
    set_names(nm = sample_names)
  common_df(sample_names, dfs, 'gsTCC')
}

greenleaf_meta <- function(paths) {
  extract_name <- function(path) {
    name <- basename(path)
    str_replace(name, '.bed.gz', '')
  }

  sample_names <- paths %>% 
    map_chr(extract_name)
  # I don't know what the extra columns encode, 95 conf int is just a guess
  extraCols_greenleaf <- c(low96 = 'numeric', high95 = 'numeric')
  dfs <- paths %>% 
    map(partial(import.bed, extraCols = extraCols_greenleaf)) %>% 
    set_names(nm = sample_names)
  common_df(sample_names, dfs, 'greenleaf_primary_ATACseq')
}

encode_meta <- function(paths) {
  extract_name <- function(path) {
    name <- basename(path)
    str_replace(name, '.narrowPeak.gz', '')
  }
  sample_names <- paths %>% map_chr(extract_name)
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
   sample_names <- paths %>% 
     map_chr(extract_name)
  dfs <- paths %>% 
    map(partial(import.bed, extraCols = extraCols_narrowPeak)) %>% 
    set_names(nm = sample_names)
  common_df(sample_names, dfs, 'encode')
}

epigenome_meta <- function(paths) {
  extract_name <- function(path) {
    name = basename(path)
    str_extract(name, '[^_]+')
  }
  sample_names <- paths %>% 
    map_chr(extract_name)
  dfs <- paths %>% 
    map(import.bed) %>% 
    set_names(nm = sample_names) %>% 
    map(function(x) x[x$name == '7_Enh', ])
  common_df(sample_names, dfs, 'epigenome_pred_enh')
}

list.files('input_data/external/epigenome_roadmap/states', full.names = TRUE) %>% 
  epigenome_meta() ->
  epigenome_df
  
list.files('input_data/external/encode_dnase', full.names = TRUE) %>% 
  encode_meta() ->
  encode_df

list.files('input_data/external/greenleaf_2015', full.names = TRUE) %>% 
  greenleaf_meta() ->
  greenleaf_df
  
list.files('input_data/external/gsTCC_DNase', full.names = TRUE) %>% 
  gsTCC_meta() ->
  gsTCC_df

df <- bind_rows(
  gsTCC_df,
  greenleaf_df,
  encode_df,
  epigenome_df
) %>% 
  mutate(project = fct_recode(project,
    "Encode" = "encode",
    "T-cell biopsy Atac-seq (Qu et al 2015)" = "greenleaf_primary_ATACseq",
    "Roadmap Epigenomics predicted enhancers" = "epigenome_pred_enh"
  ))


df %>% 
  write_csv('out/process/30_peak_desc_data.csv')

encode_df %>% 
  mutate(peaks_per_MB = num_peaks / covMB) %>% 
  summarise(m)
df %>% 
  group_by(project) %>% 
  #sample_n(23) %>% 
  ggplot() +
  geom_point(aes(covMB, num_peaks, color = project)) +
  labs(
    x = 'Total coverage in MB',
    y = 'Number of peaks',
  #  caption = "Encode peaks are with a few exceptions 150bp",
    title = "Data source greatly influence genomic properties"
    ) +
  scale_y_continuous(labels = scales::comma)
  # scale_y_log10() +
  # scale_x_log10()

# roadmap epigenomics uses chromHMM, which has a 200bp resolution
# encode has mostly just 150bp peaks
# greenleaf has 25bp intervals, but mean size converges
df %>% 
  ggplot() +
  geom_jitter(aes(project, covMedian, color = project)) +
  labs(
    caption = "Encode peaks are with a few exceptions 150bp",
    title = "Genome coverage and "
    ) +
  scale_y_continuous(labels = scales::comma)
  # scale_y_log10() +
  # scale_x_log10()
df %>% 
  group_by(project) %>% 
  sample_n(5)

df %>% 
  filter(project == 'epigenome_pred_enh')
