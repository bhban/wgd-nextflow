#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(GENESPACE)
})

option_list <- list(
  make_option("--genomes-tsv", type="character", default=NULL),
  make_option("--raw-genomerepo", type="character", default=NULL),
  make_option("--genespace-wd", type="character", default=NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

req <- c("genomes-tsv", "raw-genomerepo", "genespace-wd")
missing <- req[vapply(req, function(k) is.null(opt[[k]]) || opt[[k]] == "", logical(1))]

if (length(missing) > 0) {
  stop("Missing required option(s): --", paste(missing, collapse=" --"))
}

df <- readr::read_tsv(
  opt$`genomes-tsv`,
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    genome_source = tolower(genome_source),
    ploidy = as.integer(ploidy)
  )

stopifnot(all(df$genome_source %in% c("ncbi", "ensembl", "phytozome")))
stopifnot(!any(is.na(df$ploidy)))

run_one <- function(sub, ..., label) {
  if (nrow(sub) == 0) {
    message("No genomes for source: ", label, " (skipping)")
    return(invisible(NULL))
  }

  genomeDirs <- sub$genome_id
  genomeIDs  <- sub$genome_id

  message("Running parse_annotations for ", label, " with ", length(genomeDirs), " genomes")
  parse_annotations(
    rawGenomeRepo = opt$`raw-genomerepo`,
    genomeDirs = genomeDirs,
    genomeIDs = genomeIDs,
    gffString = "final.gff3",
    faString = "final.pep",
    genespaceWd = opt$`genespace-wd`,
    presets = "none",
    gffIdColumn = "ID",
    ...,
    convertSpecialCharacters = "_",
    troubleShoot = FALSE,
    overwrite = TRUE
  )
}

# ----------------------------
# NCBI
# ----------------------------
df_ncbi <- df %>% filter(genome_source == "ncbi")
run_one(
  df_ncbi,
  headerEntryIndex = 1,
  headerSep = " ",
  gffStripText = "cds-",
  headerStripText = "",
  label = "ncbi"
)

# ----------------------------
# Ensembl
# ----------------------------
df_ens <- df %>% filter(genome_source == "ensembl")
run_one(
  df_ens,
  headerEntryIndex = 1,
  headerSep = " ",
  gffStripText = "gene:",
  headerStripText = "",
  label = "ensembl"
)

# ----------------------------
# Phytozome
# ----------------------------
df_phy <- df %>% filter(genome_source == "phytozome")
run_one(
  df_phy,
  headerEntryIndex = 1,
  headerSep = " ",
  gffStripText = "",
  headerStripText = "",
  label = "phytozome"
)

message("All parse_annotations calls complete")
