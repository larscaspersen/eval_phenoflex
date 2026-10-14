The numerical fixture is generated from the colleague-supplied
`PhenoFlex_modular.cpp`, using its `PhenoFlex_modular_par_cpp()` entry point.
It covers all five chill selections, both heat scalings, early stopping and
full trajectories, and modified chill-hours thresholds and Utah zone weights.
The fixture stores the input weather, parameters, returned chill/heat trajectories,
the source filename and its MD5 digest. Tests do not need the external source file.

Run from the package root to regenerate with the supplied source:

```r
source("tests/testthat/fixtures/chill-reference/generate_reference.R")
generate_modular_chill_reference("path/to/PhenoFlex_modular.cpp")
```

The supplied engine overwrites `bloomindex` on every threshold crossing when
`stopatzc = FALSE`. This package consistently retains the first event, so the
fixture records both the reference's returned index and the first crossing
derived independently from its heat trajectory. Tests compare the first event.
For Utah, its raw signed `y` is compared with the package's `y_raw`; the package's
last chill column clips that total at zero for coupling. Dynamic Model internal
pool diagnostics have different row conventions, so parity compares cumulative
chill and effective heat rather than those pool columns.
