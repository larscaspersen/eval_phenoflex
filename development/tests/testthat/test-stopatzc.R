test_that("structure kernels preserve first bloom while finishing trajectories", {
  cases <- list(
    sequential = list(fn = seq_model, parameters = list(zc = 3), threshold = 3),
    parallel = list(fn = parallel_model, parameters = list(zc = 3), threshold = 3),
    partial_overlap = list(fn = po_model,
                           parameters = list(b1 = 1, b2 = 2, b3 = 0), threshold = 3)
  )
  for (case in cases) {
    args <- c(list(temp = rep(8, 200), times = seq_len(200), yc = 0.1), case$parameters)
    run <- function(stop, basic = FALSE)
      do.call(case$fn, c(args, list(stopatzc = stop, basic_output = basic)))
    early <- run(TRUE)
    full <- run(FALSE)
    first <- which(full$z >= case$threshold)[1]
    expect_true(first > 1 && first < length(args$temp))
    expect_equal(early$bloomindex, first)
    expect_equal(full$bloomindex, first)
    expect_equal(run(FALSE, TRUE)$bloomindex, first)
    expect_equal(run(TRUE, TRUE)$bloomindex, first)
    expect_equal(full$z[seq_len(first)], early$z[seq_len(first)])
    expect_gt(tail(full$z, 1), full$z[first])
    expect_true(all(early$z[(first + 1):length(args$temp)] == 0))

    # With a constant threshold, bloom on the last available row is retained.
    if (identical(case$fn, po_model)) {
      args$b1 <- tail(full$z, 1) - args$b2
    } else {
      args$zc <- tail(full$z, 1)
    }
    expect_equal(run(TRUE)$bloomindex, length(args$temp))
    expect_equal(run(FALSE)$bloomindex, length(args$temp))

    if (identical(case$fn, po_model)) args$b1 <- 1e12 else args$zc <- 1e12
    expect_equal(run(TRUE)$bloomindex, 0)
    expect_equal(run(FALSE)$bloomindex, 0)
  }
})

test_that("both population kernels retain each bud's first bloom", {
  for (fn in list(PhenoFlex_pop, PhenoFlex_pop_slim)) {
    args <- list(temp = rep(8, 400), times = seq_len(400),
                 yc = c(0.1, 0.2, 0.3), zc = c(1, 2, 1e12))
    run <- function(stop, basic = FALSE)
      do.call(fn, c(args, list(stopatzc = stop, basic_output = basic)))
    early <- run(TRUE)
    full <- run(FALSE)
    expect_equal(full$bloomindex, early$bloomindex)
    expect_equal(run(FALSE, TRUE)$bloomindex, early$bloomindex)
    expect_equal(run(TRUE, TRUE)$bloomindex, early$bloomindex)
    for (bud in 1:2) {
      first <- which(full$z[, bud] >= args$zc[bud])[1]
      expect_true(first > 1 && first < 400)
      expect_equal(full$bloomindex[bud], first)
      expect_equal(full$z[seq_len(first), bud], early$z[seq_len(first), bud])
      expect_gt(full$z[400, bud], full$z[first, bud])
    }
    expect_equal(full$bloomindex[3], 0)
  }
})

test_that("direct population calls finish trajectories for cuts after bloom", {
  for (fn in list(PhenoFlex_pop, PhenoFlex_pop_slim)) {
    args <- list(temp = rep(8, 400), times = seq_len(400),
                 yc = c(0.1, 0.2), zc = c(1, 2), i_cut = c(200, 399))
    early <- do.call(fn, c(args, list(stopatzc = TRUE, basic_output = FALSE)))
    full <- do.call(fn, c(args, list(stopatzc = FALSE, basic_output = FALSE)))
    expect_true(all(early$bloomindex > 0 & early$bloomindex < 200))
    expect_equal(early, full)
    expect_true(all(early$z[400, ] > args$zc))
    expect_equal(as.vector(early$exp), rep(0, 4))
    basic <- do.call(fn, c(args, list(stopatzc = TRUE, basic_output = TRUE)))
    expect_equal(basic$bloomindex, full$bloomindex)
    expect_named(basic, "bloomindex")
  }
})
