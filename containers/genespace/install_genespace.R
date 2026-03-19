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

strip_strings_and_comments <- function(x) {
  x <- gsub('"([^"\\\\]|\\\\.)*"', '""', x, perl = TRUE)
  x <- gsub("'([^'\\\\]|\\\\.)*'", "''", x, perl = TRUE)
  x <- sub("#.*$", "", x, perl = TRUE)
  x
}

find_function_bounds <- function(lines, fun_name) {
  start_idx <- grep(sprintf("^\\s*%s\\s*<-\\s*function\\b", fun_name), lines)
  if (length(start_idx) != 1) {
    stop(sprintf("Expected exactly one definition of %s, found %d.", fun_name, length(start_idx)))
  }

  brace_depth <- 0L
  saw_open <- FALSE
  end_idx <- NA_integer_

  for (i in seq.int(start_idx, length(lines))) {
    line_clean <- strip_strings_and_comments(lines[i])
    chars <- strsplit(line_clean, "", fixed = TRUE)[[1]]

    for (ch in chars) {
      if (ch == "{") {
        brace_depth <- brace_depth + 1L
        saw_open <- TRUE
      } else if (ch == "}") {
        brace_depth <- brace_depth - 1L
      }

      if (saw_open && brace_depth == 0L) {
        end_idx <- i
        break
      }
    }

    if (!is.na(end_idx)) {
      break
    }
  }

  if (is.na(end_idx)) {
    stop(sprintf("Could not determine end of %s function.", fun_name))
  }

  list(start = start_idx, end = end_idx)
}

patch_function_in_file <- function(filepath, fun_name, replacement_lines) {
  lines <- readLines(filepath, warn = FALSE)
  bounds <- find_function_bounds(lines, fun_name)

  new_lines <- c(
    lines[seq_len(bounds$start - 1L)],
    replacement_lines,
    lines[seq(bounds$end + 1L, length(lines))]
  )

  writeLines(new_lines, filepath, useBytes = TRUE)
}

verify_patch <- function() {
  library(GENESPACE)

  fn <- getFromNamespace("ofInBlk_engine", "GENESPACE")
  fn_text <- paste(deparse(fn), collapse = "\n")

  required_snippets <- c(
    "sb01md <- sb01[, list(n = (uid1[1] + uid2[1])/2), by = \"rid\"]",
    "propPass <- sum(sb01md$n[sb01md$n >= 40])/sum(sb01md$n)"
  )

  missing <- required_snippets[!vapply(required_snippets, grepl, logical(1), x = fn_text, fixed = TRUE)]
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

patch_lines <- readLines(patch_file, warn = FALSE)

r_files <- list.files(file.path(src_dir, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
matching_files <- r_files[
  vapply(
    r_files,
    function(f) any(grepl("^\\s*ofInBlk_engine\\s*<-\\s*function\\b", readLines(f, warn = FALSE))),
    logical(1)
  )
]

if (length(matching_files) != 1) {
  stop(
    sprintf(
      "Expected exactly one source file containing ofInBlk_engine, found %d.",
      length(matching_files)
    ),
    call. = FALSE
  )
}

patch_function_in_file(
  filepath = matching_files,
  fun_name = "ofInBlk_engine",
  replacement_lines = patch_lines
)

remotes::install_local(
  path = src_dir,
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)

verify_patch()

cat("GENESPACE installed successfully with patched ofInBlk_engine\n")
