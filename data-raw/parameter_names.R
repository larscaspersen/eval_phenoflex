# Run from the package root to regenerate the data() compatibility interface.
# Keep the exported vectors in R/parameter_names.R as the source of truth.
parameter_data <- new.env(parent = baseenv())
sys.source("R/parameter_names.R", envir = parameter_data)
for (name in c("phenoflex_parnames_kinetic", "phenoflex_parnames_characteristic",
               "phenoflex_parnames_old", "phenoflex_parnames_new")) {
  save(list = name, envir = parameter_data,
       file = file.path("data", paste0(name, ".rda")),
       compress = "xz", version = 2)
}
