cv_example <- function() {
  model <- pheno_model("sequential", chill_dynamic("kinetic"), parameters = c(yc = 0.1, zc = 5))
  years <- 2001:2008
  seasons <- lapply(years, function(y) data.frame(Temp = rep(8, 480), Year = y,
    JDay = rep(50:69, each = 24), Hour = rep(0:23, 20)))
  names(seasons) <- as.character(years)
  observed <- predict_phenology(model, seasons)
  list(model = model, seasons = seasons, years = years, observed = observed)
}
