# Explicit one-time snapshot command. Never invoked by the regression test runner.
source("tests/regression/cases.R")
destination <- "tests/testthat/fixtures/weather-regression"
baseline <- file.path(destination, "baseline.csv")
if (file.exists(baseline) && !"--overwrite" %in% commandArgs(TRUE))
  stop("Baseline exists. Review changes before deliberately using --overwrite.")
fixtures <- readRDS(file.path(destination, "inputs.rds"))
results <- run_regression_cases(fixtures)
stopifnot(all(is.finite(results$value[results$case_id == "ordinary"])),
          all(is.na(results$value[results$case_id == "no_bloom"])),
          all(c("sequential", "parallel", "partial_overlap", "combined", "three_stages",
                "calendar_dec31", "calendar_jan01", "calendar_feb29") %in% results$case_id))
write.csv(results, baseline, row.names = FALSE, na = "NA")
for (station in unique(results$station))
  write.csv(results[results$station == station, ], file.path(destination, paste0("baseline_", station, ".csv")),
            row.names = FALSE, na = "NA")
files <- c("tests/regression/cases.R", file.path(destination, "inputs.rds"),
           list.files("R", full.names = TRUE), list.files("src", pattern = "[.](cpp|h)$", full.names = TRUE))
write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))),
          file.path(destination, "baseline_sources.csv"), row.names = FALSE)
writeLines(c(paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("evalpheno:", packageVersion("evalpheno")), paste("chillR:", packageVersion("chillR")),
  paste("Package library:", find.package("evalpheno")), capture.output(sessionInfo())),
  file.path(destination, "session.txt"))
cat("Saved", nrow(results), "baseline values, with separate CSVs for each station.\n")
