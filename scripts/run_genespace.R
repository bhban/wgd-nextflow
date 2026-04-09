#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(GENESPACE)
})

option_list <- list(
  make_option("--config", type = "character"),
  make_option("--genespace-wd", type = "character"),
  make_option("--orthofinder-dir", type = "character", default = NULL),
  make_option("--genomes-tsv", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

need <- c("config", "genespace-wd", "genomes-tsv")
missing <- need[vapply(need, function(x) is.null(opt[[x]]) || opt[[x]] == "", logical(1))]
if (length(missing) > 0) {
  stop("Missing required option(s): ", paste0("--", missing, collapse = ", "))
}

cfg <- yaml::read_yaml(opt$config)
gs  <- cfg$genespace

# ----------------------------
# OrthoFinder directory
# ----------------------------
rawOrthofinderDir <- gs$rawOrthofinderDir
cli_orthofinder_dir <- opt$`orthofinder-dir`

if (!is.null(cli_orthofinder_dir) && cli_orthofinder_dir != "") {
  results_txt <- file.path(cli_orthofinder_dir, "results_dir.txt")

  if (!file.exists(results_txt)) {
    stop(paste0(
      "--orthofinder-dir was provided, but results_dir.txt was not found: ",
      results_txt,
      "\nThis should be created by orthofinder_or_skip.py"
    ))
  }

  rawOrthofinderDir <- readLines(results_txt, warn = FALSE)[1]
  rawOrthofinderDir <- trimws(rawOrthofinderDir)

  if (is.na(rawOrthofinderDir) || rawOrthofinderDir == "") {
    stop(paste0("results_dir.txt was empty: ", results_txt))
  }

  if (!dir.exists(rawOrthofinderDir)) {
    stop(paste0("Results dir from results_dir.txt does not exist: ", rawOrthofinderDir))
  }

  message("Using OrthoFinder results derived from --orthofinder-dir: ", cli_orthofinder_dir)
} else {
  if (is.null(rawOrthofinderDir) || rawOrthofinderDir == "") {
    stop(paste0(
      "No OrthoFinder directory available.\n",
      "Either provide --orthofinder-dir or set genespace.rawOrthofinderDir in config.yaml."
    ))
  }

  if (!dir.exists(rawOrthofinderDir)) {
    stop(paste0("genespace.rawOrthofinderDir does not exist: ", rawOrthofinderDir))
  }

  message("Using OrthoFinder results from genespace.rawOrthofinderDir in config.yaml")
}

# ----------------------------
# Genomes TSV (genomeIDs + ploidy)
# ----------------------------
genomes_tsv <- opt$`genomes-tsv`
if (is.null(genomes_tsv) || genomes_tsv == "") {
  stop("Missing required option: --genomes-tsv")
}
if (!file.exists(genomes_tsv)) {
  stop(paste0("genomes.tsv not found: ", genomes_tsv))
}

genomes_df <- read.delim(
  genomes_tsv,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t",
  check.names = FALSE
)

need_cols <- c("genome_id", "genome_source", "ploidy")
missing_cols <- setdiff(need_cols, colnames(genomes_df))
if (length(missing_cols) > 0) {
  stop("genomes.tsv missing required column(s): ", paste(missing_cols, collapse = ", "))
}

genomes_df$genome_source <- tolower(trimws(genomes_df$genome_source))
genomes_df$ploidy <- suppressWarnings(as.integer(genomes_df$ploidy))

if (any(is.na(genomes_df$ploidy))) {
  bad <- genomes_df$genome_id[is.na(genomes_df$ploidy)]
  stop("Ploidy could not be parsed as integers for genome_id: ", paste(bad, collapse = ", "))
}

genomeIDs <- genomes_df$genome_id
ploidy <- genomes_df$ploidy
names(ploidy) <- genomeIDs

# ----------------------------
# Parameters
# ----------------------------
blkSize <- as.integer(gs$blkSize)
if (is.na(blkSize)) {
  stop("genespace.blkSize must be an integer")
}

path2mcscanx <- gs$path2mcscanx
path2orthofinder <- gs$path2orthofinder
path2diamond <- gs$path2diamond

if (is.null(path2mcscanx) || path2mcscanx == "") {
  stop("genespace.path2mcscanx not set in config.yaml")
}
if (is.null(path2orthofinder) || path2orthofinder == "") {
  stop("genespace.path2orthofinder not set in config.yaml")
}
if (is.null(path2diamond) || path2diamond == "") {
  stop("genespace.path2diamond not set in config.yaml")
}

message("Using wd = ", opt$`genespace-wd`)
message("Using rawOrthofinderDir = ", rawOrthofinderDir)
message("Using blkSize = ", blkSize)

gpar <- init_genespace(
  wd = opt$`genespace-wd`,
  rawOrthofinderDir = rawOrthofinderDir,
  orthofinderInBlk = isTRUE(gs$orthofinder_in_blk),
  genomeIDs = genomeIDs,
  ploidy = ploidy,
  blkSize = blkSize,
  path2mcscanx = path2mcscanx,
  path2orthofinder = path2orthofinder,
  path2diamond = path2diamond
)

out <- run_genespace(gpar, overwrite = TRUE)
cat("GENESPACE complete\n")
