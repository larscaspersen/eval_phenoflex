# Build and publish the evalpheno guide

The guide uses bookdown's `bs4_book` format, with a chapter sidebar, search,
and a contents panel for the current chapter.

Edit chapter content in `vignettes/*.Rmd`. `book/build.R` generates book inputs
under `book/_build/`, adds chapter headings, and uses separate R sessions for
each chapter. The ordinary package vignettes remain available.

## Preview locally

From the package root in RStudio, install the tools once:

```r
install.packages(c("bookdown", "bslib", "downlit", "jsonlite", "xml2", "remotes"))
remotes::install_deps(dependencies = TRUE)
```

Install the current package after changing its R or C++ code, then build:

```r
install.packages(".", repos = NULL, type = "source")
source("book/build.R")
browseURL("book/_build/_book/index.html")
```

The build script also uses `.local-library/` when that project library exists.
In that case, install the current package with
`install.packages(".", lib = ".local-library", repos = NULL, type = "source")`
instead, so the guide uses the updated copy. Pandoc is included
with RStudio; command-line builds can use a separate Pandoc installation.

`bs4_book` requires `bslib`, `downlit`, `jsonlite`, and `xml2` in addition to
bookdown. These are optional dependencies of bookdown, so installing bookdown
alone may leave them missing. `Config/Needs/website` in `DESCRIPTION` declares
the website dependencies for GitHub Actions and the local build's dependency
checks.

`book/index.Rmd` is the landing page. The chapter order is set in
`book/_bookdown.yml`; when adding a vignette, also add its filename and title
to the chapter list in `book/build.R`. Appearance is configured in
`book/_output.yml` and `book/style.css`.

`split_bib: true` in `book/_output.yml` adds the cited references to the end of
each chapter. The final References chapter keeps the complete bibliography
of works cited in the guide.

## Publish with GitHub Pages

1. In the GitHub repository, open **Settings → Pages** and choose
   **GitHub Actions** under **Build and deployment → Source**.
2. Commit and push the book configuration and
   `.github/workflows/guide.yaml` to `main`.
3. Follow the **Build and publish guide** workflow in the **Actions** tab.
   You can also start it using **Run workflow**.
4. After deployment succeeds, open
   <https://larscaspersen.github.io/eval_phenoflex/>.

The workflow installs the package from the checkout and renders its vignettes
on Linux. This keeps the documentation build independent of local Windows
Rtools execution restrictions. A successful documentation build does not
replace Windows package checks.

Rendered HTML and generated chapter inputs are ignored by Git. GitHub Actions
uploads the rendered book as a Pages artifact; no `gh-pages` branch is needed.
Pull requests build the guide without deploying it.
