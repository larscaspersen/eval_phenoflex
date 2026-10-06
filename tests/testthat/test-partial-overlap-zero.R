test_that("zero overlap completes chill and holds Ca at zero while pools continue", {
  chill <- cbind(x=c(0,1,3,4,5),xs=10:14,y=c(0,.5,1.5,2,3))
  run <- function(...) apply_partial_overlap_structure(chill,c(9,0,3,4),yc=1,
                              b1=2,b2=4,b3=100,ol=0,...)
  out <- run(stopatzc=FALSE,basic_output=FALSE)
  expect_equal(out$chill,chill)
  expect_equal(out$z,c(0,0,0,3,7))
  expect_equal(out$bloomindex,5)
  expect_identical(run(),list(bloomindex=5L))
  # Lower b3 must have no effect because Ca stays zero despite overshoot.
  slow <- apply_partial_overlap_structure(chill,c(9,0,3,4),yc=1,b1=2,b2=4,
                      b3=0,ol=0,stopatzc=FALSE,basic_output=FALSE)
  expect_identical(out,slow)
  short <- apply_partial_overlap_structure(chill,c(9,0,6,4),yc=1,b1=2,b2=4,
                           b3=100,ol=0,basic_output=FALSE)
  expect_equal(short$bloomindex,4)
  expect_equal(short$z,c(0,0,0,6,0))
  expect_true(all(short$chill[5,]==0))
})

test_that("zero overlap handles initial completion and unmet chill", {
  initial <- apply_partial_overlap_structure(matrix(c(2,3,4),ncol=1),c(3,3),
                   yc=1,b1=2,b2=4,ol=0,basic_output=FALSE)
  expect_equal(initial$chill[,1],c(2,3,4))
  expect_equal(initial$z,c(0,3,6))
  expect_equal(initial$bloomindex,3)
  chill <- matrix(c(0,.2,.3),ncol=1)
  unmet <- apply_partial_overlap_structure(chill,c(100,100),yc=1,ol=0,
                                            basic_output=FALSE)
  expect_equal(unmet$chill,chill)
  expect_equal(unmet$z,c(0,0,0))
  expect_equal(unmet$bloomindex,0)
})

test_that("po_model zero overlap agrees with sequential bloom at b1+b2", {
  A1 <- exp(2/281)
  args <- list(temp=rep(8,20),times=0:19,yc=.5,E0=1,E1=2,
               A0=2*A1*exp(-1/281),A1=A1,Tf=8)
  out <- do.call(po_model,c(args,list(b1=2,b2=4,ol=0,basic_output=FALSE,
                                     stopatzc=FALSE)))
  expected <- do.call(seq_model,c(args,list(zc=6)))
  expect_gt(expected$bloomindex,0)
  expect_equal(out$bloomindex,expected$bloomindex)
  reached <- which(out$y>=.5)[1]
  expect_gt(tail(out$y,1),out$y[reached])
  full_chill <- do.call(calculate_chill_dynamic,args[names(args)!="yc"])
  expect_equal(out$x,unname(full_chill[,"x"]))
  expect_equal(out$y,unname(full_chill[,"y"]))
})
