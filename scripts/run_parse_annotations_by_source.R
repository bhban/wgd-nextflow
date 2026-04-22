#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
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

stopifnot(all(df$genome_source %in% c("ncbi", "ensembl", "phytozome", "helixer")))
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

get_attr <- function(x, key) {
  out <- str_match(x, paste0("(?:^|;)", key, "=([^;]+)"))[, 2]
  out
}

fix_ncbi_bed_one_genome <- function(genome_id, raw_genomerepo, genespace_wd, bed_start_zero_based = TRUE) {
  gff_file <- file.path(raw_genomerepo, genome_id, paste0(genome_id, ".final.gff3"))
  bed_file <- file.path(genespace_wd, "bed", paste0(genome_id, ".bed"))

  if (!file.exists(gff_file)) {
    stop("NCBI fix failed: missing GFF file: ", gff_file)
  }
  if (!file.exists(bed_file)) {
    stop("NCBI fix failed: missing BED file: ", bed_file)
  }

  gff <- readr::read_tsv(
    gff_file,
    comment = "#",
    col_names = FALSE,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )

  if (ncol(gff) != 9) {
    stop("Unexpected GFF format in ", gff_file, ": expected 9 columns")
  }

  colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

  gff <- gff %>%
    mutate(
      start = as.integer(start),
      end = as.integer(end)
    )

  genes <- gff %>%
    filter(type == "gene") %>%
    transmute(
      gene_id = get_attr(attributes, "ID"),
      chr = seqid,
      gene_start = start,
      gene_end = end
    )

  mrna <- gff %>%
    filter(type %in% c("mRNA", "transcript")) %>%
    transmute(
      mrna_id = get_attr(attributes, "ID"),
      gene_id = get_attr(attributes, "Parent")
    )

  cds <- gff %>%
    filter(type == "CDS") %>%
    transmute(
      protein_id = get_attr(attributes, "protein_id"),
      mrna_id = get_attr(attributes, "Parent")
    ) %>%
    filter(!is.na(protein_id), !is.na(mrna_id))

  protein_to_gene <- cds %>%
    distinct(protein_id, mrna_id) %>%
    left_join(mrna, by = "mrna_id") %>%
    left_join(genes, by = "gene_id") %>%
    filter(!is.na(chr), !is.na(gene_start), !is.na(gene_end)) %>%
    distinct(protein_id, chr, gene_start, gene_end)

  dup_proteins <- protein_to_gene %>%
    count(protein_id) %>%
    filter(n > 1)

  if (nrow(dup_proteins) > 0) {
    stop(
      "NCBI fix failed for ", genome_id,
      ": some protein IDs map to multiple gene spans, e.g. ",
      paste(head(dup_proteins$protein_id, 5), collapse = ", ")
    )
  }

  bed <- readr::read_tsv(
    bed_file,
    col_names = FALSE,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )

  if (ncol(bed) < 4) {
    stop("Unexpected BED format in ", bed_file, ": expected at least 4 columns")
  }

  colnames(bed)[1:4] <- c("chr", "start", "end", "protein_id")

  bed_fixed <- bed %>%
    select(-chr, -start, -end) %>%
    left_join(protein_to_gene, by = "protein_id") %>%
    mutate(
      start = if (bed_start_zero_based) gene_start - 1L else gene_start,
      end = gene_end
    )

  missing_map <- bed_fixed %>%
    filter(is.na(chr) | is.na(start) | is.na(end))

  if (nrow(missing_map) > 0) {
    stop(
      "NCBI fix failed for ", genome_id,
      ": could not map some BED protein IDs back to genes, e.g. ",
      paste(head(missing_map$protein_id, 5), collapse = ", ")
    )
  }

  bed_fixed <- bed_fixed %>%
    select(chr, start, end, protein_id)

  readr::write_tsv(bed_fixed, bed_file, col_names = FALSE)
  message("Rewrote NCBI BED with gene spans: ", bed_file)
}

fix_ncbi_beds <- function(sub, raw_genomerepo, genespace_wd, bed_start_zero_based = TRUE) {
  if (nrow(sub) == 0) {
    return(invisible(NULL))
  }

  for (genome_id in sub$genome_id) {
    fix_ncbi_bed_one_genome(
      genome_id = genome_id,
      raw_genomerepo = raw_genomerepo,
      genespace_wd = genespace_wd,
      bed_start_zero_based = bed_start_zero_based
    )
  }
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

fix_ncbi_beds(
  df_ncbi,
  raw_genomerepo = opt$`raw-genomerepo`,
  genespace_wd = opt$`genespace-wd`,
  bed_start_zero_based = TRUE
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

# ----------------------------
# Helixer
# ----------------------------
df_hel <- df %>% filter(genome_source == "helixer")
run_one(
  df_hel,
  headerEntryIndex = 1,
  headerSep = " ",
  gffStripText = "",
  headerStripText = "",
  label = "helixer"
)

message("All parse_annotations calls complete")
