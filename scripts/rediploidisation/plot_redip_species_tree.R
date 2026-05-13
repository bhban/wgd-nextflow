#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(readr)
  library(dplyr)
  library(stringr)
})

option_list <- list(
  make_option("--species-tree", type = "character", default = NULL,
              help = "Input species tree in Newick format."),
  make_option("--branch-definitions", type = "character", default = NULL,
              help = "Branch definitions TSV. Expected columns: branch_id, species. Header optional."),
  make_option("--annotation", type = "character", default = NULL,
              help = "Optional annotation TSV with tip_label and Scientific_name columns."),
  make_option("--colors", type = "character", default = NULL,
              help = "Optional colours TSV with branch_id and color columns. Header optional."),
  make_option("--output", type = "character", default = NULL,
              help = "Output plot path. Extension determines device where possible."),
  make_option("--width", type = "double", default = 8,
              help = "Output plot width in inches."),
  make_option("--height", type = "double", default = 6,
              help = "Output plot height in inches."),
  make_option("--dpi", type = "integer", default = 300,
              help = "Output DPI for raster formats."),
  make_option("--circle-size", type = "double", default = 7,
              help = "Branch marker circle size."),
  make_option("--circle-stroke", type = "double", default = 0.5,
              help = "Branch marker circle stroke width."),
  make_option("--branch-number-size", type = "double", default = 4,
              help = "Branch number label size."),
  make_option("--tip-label-size", type = "double", default = 5,
              help = "Tip label text size."),
  make_option("--tip-label-offset", type = "double", default = 0.4,
              help = "Tip label offset."),
  make_option("--branch-label-offset", type = "double", default = 0.15,
              help = "Offset for branch ID labels."),
  make_option("--prune-to-branch-species", type = "logical", default = TRUE,
              help = "Prune tree to species present in the branch definitions file.")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("species-tree", "branch-definitions", "output")
missing <- required[vapply(required, function(x) is.null(opt[[x]]) || opt[[x]] == "", logical(1))]

if (length(missing) > 0) {
  stop("Missing required option(s): --", paste(missing, collapse = " --"))
}

branch_sort_key <- function(x) {
  ifelse(grepl("^[0-9]+$", x), sprintf("%010d", as.integer(x)), x)
}

circos_color_to_hex <- function(x) {
  x <- str_trim(x)

  is_hex <- grepl("^#[0-9A-Fa-f]{6}$", x)
  if (is_hex) {
    return(toupper(x))
  }

  is_rgb <- grepl("^color=[0-9]+,[0-9]+,[0-9]+$", x)
  if (!is_rgb) {
    stop(
      "Colour value is not recognised: '", x, "'. ",
      "Use '#RRGGBB' or 'color=R,G,B'."
    )
  }

  rgb_values <- sub("^color=", "", x)
  rgb_values <- as.integer(strsplit(rgb_values, ",", fixed = TRUE)[[1]])

  if (length(rgb_values) != 3 || any(is.na(rgb_values)) || any(rgb_values < 0) || any(rgb_values > 255)) {
    stop("Invalid RGB colour value: '", x, "'.")
  }

  rgb(rgb_values[1], rgb_values[2], rgb_values[3], maxColorValue = 255)
}

load_colors <- function(path = NULL, branch_ids) {
  full_palette <- c(
    "color=68,1,84",
    "color=72,40,120",
    "color=62,73,137",
    "color=49,104,142",
    "color=38,130,142",
    "color=31,158,137",
    "color=53,183,121",
    "color=109,205,89",
    "color=180,222,44",
    "color=253,231,37"
  )

  branch_ids <- as.character(branch_ids)

  if (is.null(path) || path == "") {
    sorted_branch_ids <- sort(branch_ids, method = "radix")
    sorted_branch_ids <- sorted_branch_ids[order(branch_sort_key(sorted_branch_ids))]

    n_branches <- length(sorted_branch_ids)
    n_colors <- length(full_palette)

    if (n_branches > n_colors) {
      stop(
        "No --colors file was supplied, but there are ",
        n_branches,
        " branch IDs and only ",
        n_colors,
        " default colours."
      )
    }

    if (n_branches == 1) {
      palette_indices <- floor(n_colors / 2) + 1
    } else {
      palette_indices <- round(seq(1, n_colors, length.out = n_branches))
    }

    out <- full_palette[palette_indices]
    out <- vapply(out, circos_color_to_hex, character(1))
    names(out) <- sorted_branch_ids

    return(out)
  }

  color_tbl <- read_tsv(
    path,
    col_names = FALSE,
    comment = "#",
    show_col_types = FALSE,
    trim_ws = TRUE
  )

  if (ncol(color_tbl) < 2) {
    stop("Colors file must contain at least two columns: branch_id and color.")
  }

  color_tbl <- color_tbl[, 1:2]
  names(color_tbl) <- c("branch_id", "color")

  first_row <- as.character(color_tbl[1, ])
  if (
    tolower(first_row[1]) %in% c("branch_id", "branch", "id") &&
    tolower(first_row[2]) %in% c("color", "colour")
  ) {
    color_tbl <- color_tbl[-1, , drop = FALSE]
  }

  color_tbl <- color_tbl %>%
    mutate(
      branch_id = as.character(branch_id),
      color = as.character(color)
    ) %>%
    filter(!is.na(branch_id), branch_id != "", !is.na(color), color != "")

  color_dict <- setNames(
    vapply(color_tbl$color, circos_color_to_hex, character(1)),
    color_tbl$branch_id
  )

  missing <- setdiff(branch_ids, names(color_dict))

  if (length(missing) > 0) {
    missing <- missing[order(branch_sort_key(missing))]
    stop(
      "Colors file does not define colors for these branch IDs: ",
      paste(missing, collapse = ", ")
    )
  }

  color_dict[branch_ids]
}

read_branch_definitions <- function(path) {
  raw <- read_tsv(
    path,
    col_names = FALSE,
    comment = "#",
    show_col_types = FALSE,
    trim_ws = TRUE
  )

  if (ncol(raw) < 2) {
    stop("Branch definitions file must contain at least two columns: branch_id and species.")
  }

  raw <- raw[, 1:2]
  names(raw) <- c("branch_id", "species")

  first_row <- as.character(raw[1, ])
  if (
    tolower(first_row[1]) %in% c("branch_id", "branch", "id") &&
    tolower(first_row[2]) %in% c("species", "species_id", "genome_id", "taxon")
  ) {
    raw <- raw[-1, , drop = FALSE]
  }

  raw %>%
    mutate(
      branch_id = as.character(branch_id),
      species = str_trim(as.character(species))
    ) %>%
    filter(!is.na(branch_id), branch_id != "", !is.na(species), species != "")
}

read_tip_annotation <- function(path) {
  if (is.null(path) || path == "") {
    return(NULL)
  }

  annot <- read_tsv(
    path,
    show_col_types = FALSE,
    trim_ws = TRUE
  )

  if (!all(c("tip_label", "Scientific_name") %in% names(annot))) {
    stop("Annotation file must contain columns named 'tip_label' and 'Scientific_name'.")
  }

  annot %>%
    mutate(
      tip_label = str_trim(as.character(tip_label)),
      Scientific_name = str_trim(as.character(Scientific_name))
    ) %>%
    filter(
      !is.na(tip_label), tip_label != "",
      !is.na(Scientific_name), Scientific_name != ""
    )
}

get_branch_node <- function(tree, species_vec) {
  species_vec <- unique(species_vec)

  if (length(species_vec) == 1) {
    tip_index <- which(tree$tip.label == species_vec)

    if (length(tip_index) != 1) {
      stop("Could not uniquely match tip: ", species_vec)
    }

    return(tip_index)
  }

  node <- getMRCA(tree, species_vec)

  if (is.null(node)) {
    stop(
      "Could not find MRCA for species set: ",
      paste(species_vec, collapse = ", ")
    )
  }

  node
}

tr <- read.tree(opt$`species-tree`)

branch_defs <- read_branch_definitions(opt$`branch-definitions`)

allowed_species <- unique(branch_defs$species)

missing_species <- setdiff(allowed_species, tr$tip.label)

if (length(missing_species) > 0) {
  stop(
    "These species from the branch-definition file were not found in the tree: ",
    paste(missing_species, collapse = ", ")
  )
}

if (isTRUE(opt$`prune-to-branch-species`)) {
  tr_plot <- keep.tip(tr, allowed_species)
} else {
  tr_plot <- tr
}

annot <- read_tip_annotation(opt$annotation)

if (is.null(annot)) {
  tip_map <- setNames(tr_plot$tip.label, tr_plot$tip.label)
} else {
  tip_map <- setNames(annot$Scientific_name, annot$tip_label)
}

branch_ids <- unique(as.character(branch_defs$branch_id))
branch_colors <- load_colors(opt$colors, branch_ids)

branch_nodes <- branch_defs %>%
  group_by(branch_id) %>%
  summarise(species = list(unique(species)), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    node = get_branch_node(tr_plot, unlist(species)),
    color = unname(branch_colors[as.character(branch_id)])
  ) %>%
  ungroup() %>%
  arrange(branch_sort_key(branch_id))

p_base <- suppressWarnings(ggtree(tr_plot, size = 0.8))

tree_df <- p_base$data %>%
  mutate(label = if_else(is.na(label), "", label)) %>%
  mutate(
    display_label = if_else(
      isTip & label %in% names(tip_map),
      unname(tip_map[label]),
      label
    )
  )

circle_df <- tree_df %>%
  select(node, x, y, label, isTip) %>%
  inner_join(
    branch_nodes %>% select(branch_id, node, color),
    by = "node"
  ) %>%
  mutate(
    label_x = x + opt$`branch-label-offset`,
    label_y = y
  )

p_redip <- suppressWarnings(
  ggtree(tr_plot, size = 0.8) %<+% tree_df +
    geom_tiplab(
      aes(label = display_label),
      align = TRUE,
      linesize = 0.3,
      size = opt$`tip-label-size`,
      fontface = 3,
      offset = opt$`tip-label-offset`
    ) +
    geom_point(
      data = circle_df,
      aes(x = x, y = y, fill = color),
      shape = 21,
      size = opt$`circle-size`,
      color = "black",
      stroke = opt$`circle-stroke`
    ) +
    geom_text(
      data = circle_df,
      aes(x = label_x, y = label_y, label = branch_id),
      size = opt$`branch-number-size`,
      hjust = 0,
      vjust = 0.35
    ) +
    scale_fill_identity() +
    theme_tree2() +
    theme(
      legend.position = "none",
      plot.margin = margin(20, 80, 20, 20)
    )
)

ggsave(
  filename = opt$output,
  plot = p_redip,
  width = opt$width,
  height = opt$height,
  dpi = opt$dpi,
  limitsize = FALSE
)

message("Wrote redip species tree plot: ", opt$output)
