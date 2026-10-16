# Run source("book/build.R") from the package root.
build_guide <- function() {
  if (!file.exists("DESCRIPTION") || !file.exists("book/index.Rmd")) {
    stop("Run source('book/build.R') from the evalpheno package root.")
  }
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  old_libraries <- .libPaths()
  on.exit(.libPaths(old_libraries), add = TRUE)
  if (dir.exists(".local-library")) {
    .libPaths(c(normalizePath(".local-library"), .libPaths()))
  }

  for (package in c("bookdown", "rmarkdown", "evalpheno")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install ", package, " before building the guide.")
    }
  }

  # Let fresh chapter sessions use the same package libraries as this session.
  old_r_libs <- Sys.getenv("R_LIBS", unset = NA_character_)
  on.exit(if (is.na(old_r_libs)) Sys.unsetenv("R_LIBS") else
    Sys.setenv(R_LIBS = old_r_libs), add = TRUE)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

  if (!rmarkdown::pandoc_available() && .Platform$OS.type == "windows") {
    bundled_pandoc <- file.path(Sys.getenv("ProgramFiles"), "RStudio",
                               "resources/app/bin/quarto/bin/tools")
    if (file.exists(file.path(bundled_pandoc, "pandoc.exe"))) {
      old_pandoc <- Sys.getenv("RSTUDIO_PANDOC", unset = NA_character_)
      on.exit(if (is.na(old_pandoc)) Sys.unsetenv("RSTUDIO_PANDOC") else
        Sys.setenv(RSTUDIO_PANDOC = old_pandoc), add = TRUE)
      Sys.setenv(RSTUDIO_PANDOC = bundled_pandoc)
    }
  }
  if (!rmarkdown::pandoc_available()) {
    stop("Pandoc is required. Build from RStudio or install Pandoc.")
  }

  stage <- file.path(root, "book", "_build")
  dir.create(stage, recursive = TRUE, showWarnings = FALSE)
  assets <- c("index.Rmd", "references.Rmd", "_bookdown.yml", "_output.yml",
              "style.css")
  if (!all(file.copy(file.path(root, "book", assets), stage,
                     overwrite = TRUE))) {
    stop("Could not copy the book configuration.")
  }
  bibliography <- c("vignette_evalpheno.bib", "elsevier-harvard.csl")
  if (!all(file.copy(file.path(root, "vignettes", bibliography), stage,
                     overwrite = TRUE))) {
    stop("Could not copy the bibliography and citation style.")
  }

  # Generate chapter inputs without changing or duplicating vignette sources.
  chapters <- c(
    "getting-started" = "Getting started",
    "model-structures-and-submodels" = "Model structures and submodels",
    "model-calibration_special-cases" = "Model calibration",
    "bud-population" = "Bud populations"
  )
  for (slug in names(chapters)) {
    source <- file.path(root, "vignettes", paste0(slug, ".Rmd"))
    lines <- readLines(source, warn = FALSE, encoding = "UTF-8")
    delimiters <- which(trimws(lines) == "---")
    if (length(delimiters) < 2L || delimiters[1L] != 1L) {
      stop("Expected YAML front matter in ", source)
    }
    body <- lines[seq.int(delimiters[2L] + 1L, length(lines))]
    body <- body[!grepl("^## References[[:space:]]*$", body)]
    heading <- paste0("# ", chapters[[slug]], " {#", slug, "}")
    writeLines(c(heading, "", body), file.path(stage, paste0(slug, ".Rmd")),
               useBytes = TRUE)
  }

  old_directory <- setwd(stage)
  on.exit(setwd(old_directory), add = TRUE)
  bookdown::render_book("index.Rmd", output_format = "bookdown::bs4_book",
                        encoding = "UTF-8")
  file.create(file.path("_book", ".nojekyll"))
  message("Guide built at ", file.path(stage, "_book", "index.html"))
  invisible(file.path(stage, "_book", "index.html"))
}

build_guide()
