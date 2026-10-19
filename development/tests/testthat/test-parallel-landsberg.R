test_that("Landsberg coupling uses exponential effectiveness and keeps every pool updating", {
  chill <- cbind(x=c(0,2,3,4,5,6), xs=10:15, final=c(0,0,1,2,3,4))
  saved <- chill
  heat <- rep(10,5)
  out <- apply_parallel_landsberg_structure(chill,heat,y0=1,yc=3,zc=10,
                                           stopatzc=FALSE,basic_output=FALSE)
  expect_equal(out$z,c(0,0,cumsum(10*(1-exp(-c(1,2,3,4))))))
  expect_equal(out$bloomindex,4)
  expect_equal(out$chill[6,],chill[6,])
  expect_identical(out$chill[1:5,],chill[1:5,])
  expect_identical(chill,saved)
  expect_identical(heat,rep(10,5))
  expect_equal(tail(diff(out$z),1),10*(1-exp(-4)))
  early <- apply_parallel_landsberg_structure(chill,heat,y0=1,yc=3,zc=10,
                                             basic_output=FALSE)
  expect_equal(early$bloomindex,4)
  expect_equal(early$z[5:6],c(0,0))
  expect_true(all(early$chill[5:6,]==0))
  expect_identical(apply_parallel_landsberg_structure(chill,heat,y0=1,yc=3,zc=10),
                   list(bloomindex=4L))
})

test_that("Landsberg handles initial chill, threshold crossing and numerical extremes", {
  run <- function(y, y0=1, yc=3, zc=100, heat=rep(1,length(y)-1)) {
    apply_parallel_landsberg_structure(matrix(y,ncol=1),heat,y0=y0,yc=yc,
                                       zc=zc,stopatzc=FALSE,basic_output=FALSE)
  }
  initial <- run(c(4,5,6))
  expect_equal(initial$chill[,1],c(4,5,6))
  expect_equal(initial$z,c(0,cumsum(1-exp(-c(5,6)))))
  overshoot <- run(c(0,4,5))
  expect_identical(overshoot,run(c(0,4,5),yc=100))
  expect_equal(overshoot$chill[,1],c(0,4,5))
  expect_equal(overshoot$z,c(0,cumsum(1-exp(-c(4,5)))))
  expect_equal(run(c(0,0,0))$z,c(0,0,0))
  expect_equal(run(c(0,1,2))$bloomindex,0)
  expect_equal(run(c(0,1,2),zc=1)$bloomindex,3)
  expect_equal(run(c(0,1e-20))$z[2],1e-20,tolerance=1e-30)
  expect_equal(run(c(0,1e300),y0=1e-300)$z[2],1)
  expect_error(run(c(0,1,2),heat=rep(.Machine$double.xmax,2),
                   y0=1e-300),"overflowed")
})

test_that("Landsberg validates scale and structure inputs", {
  chill <- matrix(c(0,1,2),ncol=1)
  for (bad in c(0,-1,Inf,NA_real_)) {
    expect_error(apply_parallel_landsberg_structure(chill,c(1,1),y0=bad),"y0")
    expect_error(apply_parallel_landsberg_structure(chill,c(1,1),y0=1,yc=bad),"yc")
    expect_error(apply_parallel_landsberg_structure(chill,c(1,1),y0=1,zc=bad),"zc")
  }
  expect_error(apply_parallel_landsberg_structure(chill,1,y0=1))
  expect_error(apply_parallel_landsberg_structure(chill,c(-1,1),y0=1))
  expect_error(apply_parallel_landsberg_structure(matrix(c(0,-1),ncol=1),1,y0=1),
               "non-negative")
  expect_error(apply_parallel_landsberg_structure(matrix(c(0,NA),ncol=1),1,y0=1))
})
