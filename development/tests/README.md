# Development tests

Package checks run only `tests/testthat/test-smoke.R`: eight test blocks with
50 expectations covering native calculations, parameter conversion, the four
model structures, invalid inputs, fitting, collections, populations, and a small
fixed-parameter cross-validation/ensemble workflow. Optimizer smoke checks use
one iteration and a small evaluation budget, without testing convergence.

The original comprehensive suite and its fixtures are preserved here. It includes
optimizer convergence, repeated cross-validation, plotting, population sampling,
and the frozen weather regressions. This directory is excluded from source
packages by the existing `^development$` entry in `.Rbuildignore`.

Run the smoke checks and full development suite from the package root:

```r
source("development/tests/run-tests.R")
```

Or from a terminal:

```sh
Rscript development/tests/run-tests.R
```

The runner uses `pkgload` to load the current package source, including internal
functions needed by the development tests. The weather regression scripts and
their separate instructions are in `development/tests/regression/`.
