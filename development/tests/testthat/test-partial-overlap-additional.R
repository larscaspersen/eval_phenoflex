test_that("partial overlap uses additional chill and excludes crossing overshoot", {
  # At completion (y=3, yc=2), Ca=0 and the requirement is 2+8=10.
  # One later chill unit lowers it to 2+8/2=6, reached in the final row.
  chill <- cbind(x=c(0,2,4,5), y=c(0,3,4,8))
  out <- apply_partial_overlap_structure(chill,c(4,2,1),yc=2,b1=2,b2=8,
                                         b3=log(2),ol=10,basic_output=FALSE)
  expect_equal(out$bloomindex,3)
  expect_equal(out$z,c(0,4,6,0))
  expect_equal(out$chill[1:3,],chill[1:3,])
  expect_true(all(out$chill[4,]==0))
  # Translate the entire chill trajectory and requirement: Ca is unchanged.
  shifted <- chill; shifted[,2] <- shifted[,2]+10
  translated <- apply_partial_overlap_structure(shifted,c(4,2,1),yc=12,b1=2,
                    b2=8,b3=log(2),ol=10,basic_output=FALSE)
  expect_equal(translated$z,out$z)
  expect_equal(translated$bloomindex,out$bloomindex)
  full <- apply_partial_overlap_structure(chill,c(4,2,1),yc=2,b1=2,b2=8,
                        b3=log(2),ol=10,stopatzc=FALSE,basic_output=FALSE)
  expect_equal(full$bloomindex,3)
  expect_equal(full$chill,chill)
})

test_that("zero additional chill requires b1+b2, after overlap ends while chill continues", {
  chill <- cbind(x=0:3,y=c(0,40,80,100))
  out <- apply_partial_overlap_structure(chill,c(30,40,50),yc=40,b1=20,b2=100,
                           b3=.1,ol=.5,stopatzc=FALSE,basic_output=FALSE)
  expect_equal(out$chill,chill)
  expect_equal(out$chill[,1],chill[,1])
  expect_equal(out$bloomindex,4)
  expect_equal(out$z,c(0,30,70,120))
  # Already-complete inputs use the initial row as their baseline.
  initial <- apply_partial_overlap_structure(matrix(c(40,41,50),ncol=1),c(5,1),
                          yc=20,b1=2,b2=8,b3=log(2),ol=10,basic_output=FALSE)
  expect_equal(initial$bloomindex,3)
  expect_equal(initial$z,c(0,5,6))
})

test_that("additional chill after overlap cannot change heat or bloom", {
  # Completion at y=1; overlap ends at y=2 and heat=4, giving Ca=1.
  # The fixed requirement is 4+8/2=8 despite later large chill increases.
  chill <- cbind(x=0:5, xs=10:15, y=c(0,1,2,20,30,40))
  run <- function(ch) apply_partial_overlap_structure(ch,c(2,2,1,3,2),
    yc=1,b1=4,b2=8,b3=log(2),ol=1,stopatzc=FALSE,basic_output=FALSE)
  out <- run(chill)
  expect_equal(out$chill,chill)
  expect_equal(out$z,c(0,2,4,5,8,10))
  expect_equal(out$bloomindex,5)
  changed <- chill
  changed[4:6,] <- 0 # Even a later declining submodel cannot alter coupling.
  other <- run(changed)
  expect_equal(other$chill,changed)
  expect_equal(other$z,out$z)
  expect_equal(other$bloomindex,out$bloomindex)
})
