# Metadata shared by the guide's generated chapters and downloadable citations.
guide_metadata <- function(source, build_date) {
  metadata <- rmarkdown::yaml_front_matter(source)
  authors <- unname(unlist(metadata$author))
  if (!length(authors) || any(!nzchar(authors))) {
    stop("Expected chapter authors in ", source)
  }
  date <- metadata$date
  rendered <- is.null(date) || identical(date, "`r Sys.Date()`")
  if (rendered) {
    date <- build_date
  } else {
    date <- as.Date(date)
    if (length(date) != 1L || is.na(date)) {
      stop("Use an ISO date or `r Sys.Date()` in ", source)
    }
  }
  list(authors = authors, date = date,
       date_label = if (rendered) "Rendered" else "Date")
}

guide_byline <- function(metadata, version) {
  c("::: {.chapter-meta}",
    paste(metadata$authors, collapse = ", "), "",
    paste0(metadata$date_label, ": ", metadata$date,
           " · evalpheno ", version),
    ":::", "")
}

guide_citation <- function(title, metadata, version, url, key, file,
                           chapter = FALSE) {
  note <- paste0(if (chapter) "Chapter in The evalpheno guide. " else "",
                 "evalpheno version ", version, ". ", metadata$date_label,
                 " ", metadata$date)
  entry <- utils::bibentry(
    bibtype = "Manual", key = key, title = title,
    author = utils::as.person(metadata$authors),
    year = format(metadata$date, "%Y"), note = note, url = url
  )
  bibtex <- format(entry, style = "Bibtex")
  writeLines(bibtex, file, useBytes = TRUE)
  c(paste0("## Cite this ", if (chapter) "chapter" else "guide", " {-}"),
    "", "::: {.guide-citation}",
    format(entry, style = "text"), "",
    paste0("[Download BibTeX](citations/", basename(file), ")"), "",
    "<details>", "<summary>BibTeX</summary>", "",
    "```bibtex", bibtex, "```", "", "</details>", ":::")
}
