# Run from the package root with the supplied PhenoFlex_modular.cpp as input.
generate_modular_chill_reference <- function(source_path,
    output = "tests/testthat/fixtures/chill-reference/reference.rds") {
  reference <- new.env(parent = globalenv())
  Rcpp::sourceCpp(source_path, env = reference, showOutput = FALSE)
  weights <- c(0, 0.5, 1, 0.5, 0, -0.5, -1)
  variants <- list(
    hours = list(name = "chill_hours", code = 1, p = c(CH_lower = 0, CH_upper = 7.2)),
    dynamic = list(name = "dynamic", code = 2,
      p = c(A0 = 6319.5, A1 = 5.939917e13, E0 = 3372.8, E1 = 9900.3, Tf = 4, slope = 1.6)),
    utah = list(name = "utah", code = 3, p = setNames(weights, paste0("utah_w", 1:7))),
    positive = list(name = "positive_utah", code = 4,
      p = setNames(weights, paste0("pos_utah_w", 1:7))),
    none = list(name = "none", code = 5, p = c(heat_start_index = 25)),
    hours_custom = list(name = "chill_hours", code = 1, p = c(CH_lower = 2, CH_upper = 8.5)),
    utah_custom = list(name = "utah", code = 3,
      p = setNames(c(0.2, 0.4, 1.2, 0.7, 0.1, -0.75, -1.2), paste0("utah_w", 1:7))),
    positive_custom = list(name = "positive_utah", code = 4,
      p = setNames(c(0.2, 0.4, 1.2, 0.7, 0.1, -0.75, -1.2), paste0("pos_utah_w", 1:7))))
  weather <- data.frame(Temp = rep(c(20, 5, 8, 12, 20, 2, 5, 8, 12, 5, 8, 12), each = 24),
    Year = 2008, JDay = rep(1:12, each = 24), Hour = rep(0:23, 12))
  records <- list()
  for (id in names(variants)) for (scaling in c("unscaled", "scaled")) for (stop in c(TRUE, FALSE)) {
    variant <- variants[[id]]
    structure_parameters <- if (variant$name == "none") c(zc = 2) else c(yc = 1, zc = 2, s1 = 0.5)
    if (scaling == "scaled") structure_parameters["zc"] <- structure_parameters["zc"] * 21
    p <- c(structure_parameters, variant$p, Tb = 4, Tu = 25, Tc = 36)
    result <- reference$PhenoFlex_modular_par_cpp(weather$Temp, seq_len(nrow(weather)), p,
      chill_model = variant$code, heat_model = if (scaling == "scaled") 2L else 1L,
      stopatzc = stop, basic_output = FALSE, day_group = rep(1:12, each = 24))
    crossings <- which(result$z >= p[["zc"]])
    key <- paste(id, scaling, stop, sep = "_")
    records[[key]] <- list(chill = variant$name, scaling = scaling, parameters = p,
      stopatzc = stop, y = result$y, z = result$z,
      first_bloomindex = if (length(crossings)) crossings[1] else 0L,
      reference_bloomindex = result$bloomindex)
  }
  fixture <- list(source = basename(source_path), source_md5 = unname(tools::md5sum(source_path)),
    weather = weather, records = records)
  saveRDS(fixture, output, version = 2)
  invisible(fixture)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: Rscript generate_reference.R path/to/PhenoFlex_modular.cpp")
  generate_modular_chill_reference(args[1])
}
