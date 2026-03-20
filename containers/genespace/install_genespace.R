cran_repo <- "https://cloud.r-project.org"
options(repos = c(CRAN = cran_repo))

github_commit <- "7561036be821e333c1dfb52461ccec7222e95582"
github_tarball <- sprintf(
  "https://github.com/jtlovell/GENESPACE/archive/%s.tar.gz",
  github_commit
)

check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed.", pkg), call. = FALSE)
  }
}

append_collate_entry <- function(description_file, new_file) {
  desc <- readLines(description_file, warn = FALSE)

  collate_idx <- grep("^Collate:", desc)
  if (length(collate_idx) == 0) {
    return(invisible(FALSE))
  }

  start <- collate_idx[1]
  end <- start

  while (end < length(desc) && grepl("^[[:space:]]", desc[end + 1])) {
    end <- end + 1
  }

  collate_block <- desc[start:end]

  if (any(grepl(sprintf("'%s'", new_file), collate_block, fixed = TRUE))) {
    return(invisible(TRUE))
  }

  desc[end] <- paste0(desc[end], sprintf("\n        '%s'", new_file))
  writeLines(desc, description_file, useBytes = TRUE)
  invisible(TRUE)
}

verify_patch <- function() {
  library(GENESPACE)

  fn <- get("ofInBlk_engine", envir = asNamespace("GENESPACE"))
  fn_text <- paste(deparse(fn), collapse = "\n")

  required_snippets <- c(
    "sb01md <- sb01[, list(n = (uid1[1] + uid2[1])/2), by = \"rid\"]",
    "propPass <- sum(sb01md$n[sb01md$n >= 40])/sum(sb01md$n)"
  )

  missing <- required_snippets[
    !vapply(required_snippets, grepl, logical(1), x = fn_text, fixed = TRUE)
  ]

  if (length(missing) > 0) {
    stop(
      paste(
        "Installed GENESPACE does not appear to contain the patched ofInBlk_engine.",
        "Missing expected snippet(s):",
        paste(missing, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

check_pkg("Biostrings")
check_pkg("rtracklayer")
check_pkg("remotes")
check_pkg("igraph")
check_pkg("yaml")
check_pkg("optparse")
check_pkg("dplyr")
check_pkg("readr")

tmp_dir <- tempfile("genespace_src_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

tarball_path <- file.path(tmp_dir, "GENESPACE.tar.gz")
download.file(github_tarball, destfile = tarball_path, mode = "wb", quiet = FALSE)
utils::untar(tarball_path, exdir = tmp_dir)

src_dirs <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
src_dir <- src_dirs[grepl("^GENESPACE-", basename(src_dirs))]

if (length(src_dir) != 1) {
  stop("Could not uniquely identify unpacked GENESPACE source directory.", call. = FALSE)
}

patch_file <- "/tmp/ofInBlk_engine_patched.R"
if (!file.exists(patch_file)) {
  stop("Patch file /tmp/ofInBlk_engine_patched.R was not found.", call. = FALSE)
}

target_patch_name <- "zz_ofInBlk_engine_patch.R"
target_patch_path <- file.path(src_dir, "R", target_patch_name)

file.copy(patch_file, target_patch_path, overwrite = TRUE)

description_file <- file.path(src_dir, "DESCRIPTION")
if (file.exists(description_file)) {
  append_collate_entry(description_file, target_patch_name)
}

remotes::install_local(
  path = src_dir,
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)

verify_patch()
cat("GENESPACE installed successfully with patched ofInBlk_engine\n")
