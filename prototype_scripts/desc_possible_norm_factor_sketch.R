#
# What properties of tracks makes good normalization factors?
# - Promoters have poor overlap with enhancers, but perhaps an indication regardless?
#
pacman::p_load(
  tidyverse,
  stringr,
  rtracklayer,
  EnsDb.Hsapiens.v75
)

total_MB_coverage <- function(gr) {
  gr %>% 
    coverage() %>% 
    sum() %>% 
    as.numeric()
    sum() %>% 
    `/`(1e6)
}

promoter_regions <- promoters(EnsDb.Hsapiens.v75, downstream = 500, upstream = 3000) %>% 
  reduce() # remove overlaps

# promoter genome coverage
promoter_regions %>% 
  total_MB_coverage()

x <- import.bed('input_data/external/epigenome_roadmap/states/E001_15_coreMarks_dense.bed.gz')
seqlevelsStyle(x) <- 'NCBI'

enhancers <- x[str_detect(mcols(x)$name, 'Enh'), ]
enhancers %>% 
  total_MB_coverage()

enhancers %>% 
  intersect(promoter_regions, ignore.strand = TRUE) %>% 
  total_MB_coverage()
x
promoter_regions
