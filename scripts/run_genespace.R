#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GENESPACE)
})

option_list <- list(
  make_option("--genespace-wd", type = "character"),
  make_option("--orthofinder-dir", type = "character", default = NULL),
  make_option("--genomes-tsv", type = "character"),
  make_option("--raw-orthofinder-dir", type = "character", default = ""),
  make_option("--blk-size", type = "integer"),
  make_option("--orthofinder-in-blk", type = "character", default = "true")
)

opt <- parse_args(OptionParser(option_list = option_list))

need <- c(
  "genespace-wd",
  "genomes-tsv",
  "blk-size"
)

missing <- need[vapply(need, function(x) is.null(opt[[x]]) || opt[[x]] == "", logical(1))]
if (length(missing) > 0) {
  stop("Missing required option(s): ", paste0("--", missing, collapse = ", "))
}

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "yes", "y")
}

# ----------------------------
# Tool paths from container PATH
# ----------------------------
path2mcscanx <- dirname(Sys.which("MCScanX"))
path2orthofinder <- Sys.which("orthofinder")
path2diamond <- Sys.which("diamond")

if (path2mcscanx == "") {
  stop("MCScanX not found on PATH inside container")
}
if (path2orthofinder == "") {
  stop("orthofinder not found on PATH inside container")
}
if (path2diamond == "") {
  stop("diamond not found on PATH inside container")
}

# ----------------------------
# OrthoFinder directory
# ----------------------------
rawOrthofinderDir <- opt$`raw-orthofinder-dir`
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

} else if (!is.null(rawOrthofinderDir) && rawOrthofinderDir != "") {

  if (!dir.exists(rawOrthofinderDir)) {
    stop(paste0("--raw-orthofinder-dir does not exist: ", rawOrthofinderDir))
  }

  message("Using OrthoFinder results from --raw-orthofinder-dir")

} else {
  rawOrthofinderDir <- NULL
  message("No OrthoFinder directory provided; proceeding without rawOrthofinderDir")
}

# ----------------------------
# Genomes TSV
# ----------------------------
genomes_tsv <- opt$`genomes-tsv`
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
blkSize <- as.integer(opt$`blk-size`)
if (is.na(blkSize)) {
  stop("--blk-size must be an integer")
}

orthofinderInBlk <- as_bool(opt$`orthofinder-in-blk`)

message("Using wd = ", opt$`genespace-wd`)
message("Using blkSize = ", blkSize)
message("Using MCScanX = ", path2mcscanx)
message("Using orthofinder = ", path2orthofinder)
message("Using diamond = ", path2diamond)

gs_args <- list(
  wd = opt$`genespace-wd`,
  orthofinderInBlk = orthofinderInBlk,
  genomeIDs = genomeIDs,
  ploidy = ploidy,
  blkSize = blkSize,
  path2mcscanx = path2mcscanx,
  path2orthofinder = path2orthofinder,
  path2diamond = path2diamond
)

if (!is.null(rawOrthofinderDir) && rawOrthofinderDir != "") {
  gs_args$rawOrthofinderDir <- rawOrthofinderDir
}

gpar <- do.call(init_genespace, gs_args)

out <- run_genespace(gpar, overwrite = TRUE)
cat("GENESPACE complete\n")
                       
