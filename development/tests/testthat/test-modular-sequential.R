test_that("Dynamic and GDH modules have independently checked values and alignment", {
  A1 <- exp(2 / 281)
  chill <- calculate_chill_dynamic(rep(8, 5), 0:4, A0 = 10 * A1 * exp(-1 / 281),
                                    A1 = A1, E0 = 1, E1 = 2, Tf = 8)
  first_pool <- 5 * (1 - exp(-1))
  expect_identical(colnames(chill), c("x", "xs", "y"))
  expect_equal(unname(chill[1, ]), c(0, 0, 0))
  expect_equal(unname(chill[2, ]), c(first_pool, 10, first_pool), tolerance = 1e-12)
  expect_gt(first_pool, 1) # Residual x can still convert unless explicitly frozen.
  expect_gt(chill[3, "y"], chill[2, "y"])
  expect_equal(calculate_heat_gdh(c(4, 25, 36, 0), 0:3), c(0, 21, 0), tolerance = 1e-12)
  expect_equal(calculate_heat_gdh(c(4, 25, 36, 999), 0:3), c(0, 21, 0), tolerance = 1e-12)
})

test_that("sequential coupling uses the last column and keeps all chill pools updating", {
  # The first column would reach yc immediately if selected by mistake.
  chill <- cbind(x = c(100, 101, 102, 103), measure = c(0, 0.4, 1, 1.2))
  original <- chill
  result <- apply_sequential_structure(chill, c(2, 3, 5), yc=1, zc=4,
                                       stopatzc=FALSE, basic_output=FALSE)
  expect_equal(result$z, c(0, 0, 3, 8))
  expect_equal(result$bloomindex, 4)
  expect_equal(result$chill[, 1], c(100, 101, 102, 103))
  expect_equal(result$chill[, 2], c(0, 0.4, 1, 1.2))
  expect_identical(chill, original)
  expect_equal(apply_sequential_structure(chill[, 2, drop=FALSE], c(2, 3, 5),
                                          yc=1, zc=4)$bloomindex, 4)
  no_bloom <- apply_sequential_structure(chill, c(2, 3, 5), yc=2, zc=4,
                                         basic_output=FALSE)
  expect_equal(no_bloom$z, rep(0, 4))
  expect_equal(no_bloom$bloomindex, 0)
  expect_equal(no_bloom$chill, chill)
  initial <- apply_sequential_structure(matrix(c(1, 2, 3), ncol=1), c(2, 3),
                                        yc=1, zc=2, basic_output=FALSE)
  expect_equal(initial$bloomindex, 2)
  expect_equal(as.vector(initial$chill), c(1, 2, 0))
})

test_that("seq_model composes the modules and continues chill accumulation", {
  A1 <- exp(2 / 281)
  args <- list(temp=rep(8, 5), times=0:4, A0=10*A1*exp(-1/281), A1=A1,
               E0=1, E1=2, Tf=8)
  chill <- do.call(calculate_chill_dynamic, args)
  heat <- calculate_heat_gdh(args$temp, args$times)
  modular <- apply_sequential_structure(chill, heat, yc=1, zc=1e6,
                                         stopatzc=FALSE, basic_output=FALSE)
  legacy_interface <- do.call(seq_model, c(args, list(yc=1, zc=1e6,
                                                      stopatzc=FALSE, basic_output=FALSE)))
  expect_gt(chill[2, "x"], 1)
  expect_equal(legacy_interface$x, unname(chill[, "x"]))
  expect_equal(legacy_interface$y, unname(chill[, "y"]))
  expect_equal(legacy_interface$x, unname(modular$chill[, "x"]))
  expect_equal(legacy_interface$y, unname(modular$chill[, "y"]))
  expect_equal(legacy_interface$z, modular$z)
  expect_equal(legacy_interface$xs, c(chill[-1, "xs"], 0))
  expect_equal(legacy_interface$bloomindex, modular$bloomindex)
})

test_that("module entry points reject incompatible or unsafe inputs", {
  expect_error(calculate_chill_dynamic(numeric(), numeric()), "at least two")
  expect_error(calculate_chill_dynamic(c(8, 8), c(0, 2)), "hourly")
  expect_error(calculate_chill_dynamic(c(8, 8), 0:1, A0=0), "A0")
  expect_error(calculate_heat_gdh(c(8, 8), c(0, 2)), "hourly")
  expect_error(calculate_heat_gdh(c(8, 8), 0:1, Tb=25), "Tb < Tu < Tc")
  expect_error(apply_sequential_structure(matrix(0, 1, 1), numeric()), "at least two")
  expect_error(apply_sequential_structure(matrix(0, 3, 0), c(1, 1)), "one column")
  expect_error(apply_sequential_structure(matrix(0, 3, 1), 1), "nrow")
  expect_error(apply_sequential_structure(matrix(c(0, NA), 2), 1), "finite")
  expect_error(apply_sequential_structure(matrix(0, 2, 1), -1), "non-negative")
  expect_error(apply_sequential_structure(matrix(0, 2, 1), NA_real_), "finite")
  expect_error(apply_sequential_structure(matrix(0, 2, 1), 1, yc=0), "yc")
})

test_that("sequential onset remains enabled while chill continues changing", {
  chill <- matrix(c(0,1,0.5,2),ncol=1)
  out <- apply_sequential_structure(chill,c(1,2,3),yc=1,zc=100,
                                    stopatzc=FALSE,basic_output=FALSE)
  expect_equal(out$z,c(0,1,3,6))
  expect_identical(out$chill,chill)
})
