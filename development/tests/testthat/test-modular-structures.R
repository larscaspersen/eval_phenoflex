test_that("parallel coupling weights updated chill and continues chill after the crossing row", {
  chill <- cbind(x=c(100,101,102,103), final=c(0,0.5,1,2))
  saved <- chill
  out <- apply_parallel_structure(chill, c(2,4,6), yc=1, zc=5, kmin=0.2,
                                  stopatzc=FALSE, basic_output=FALSE)
  expect_equal(out$z, c(0,1.2,5.2,11.2))
  expect_equal(out$bloomindex, 3)
  expect_equal(out$chill[,2], c(0,0.5,1,2))
  expect_equal(out$chill[,1], c(100,101,102,103))
  expect_identical(chill, saved)
  expect_identical(colnames(out$chill), colnames(chill))
  expect_equal(apply_parallel_structure(chill, c(2,4,6), yc=1, zc=5, kmin=0.2)$bloomindex, 3)
  early <- apply_parallel_structure(chill, c(2,4,6), yc=1, zc=5, kmin=0.2, basic_output=FALSE)
  expect_equal(early$z, c(0,1.2,5.2,0))
  expect_equal(early$chill[4,], setNames(c(0,0),colnames(chill)))
  for (k in c(0,1)) {
    out <- apply_parallel_structure(chill[,2,drop=FALSE], c(2,4,6), yc=1, zc=100,
                                    kmin=k, stopatzc=FALSE, basic_output=FALSE)
    expect_equal(out$z, if(k==0) c(0,1,5,11) else c(0,2,6,12))
    expect_equal(out$bloomindex, 0)
  }
})

test_that("partial overlap holds coupling chill while returning continuing pools", {
  chill <- cbind(x=100:104, final=c(0,0.5,1,2,3))
  saved <- chill
  run <- function(...) apply_partial_overlap_structure(chill, c(2,4,3,10), yc=1,
                                                        b1=5,b2=8,b3=log(2),ol=0.8,...)
  out <- run(stopatzc=FALSE,basic_output=FALSE)
  # Heat reaches b1*ol=4 at row 3. Hold its coupling value from row 4 onward.
  # Coupling y=1 and baseline=1 give requirement 5+8=13, reached at row 5, not row 4.
  expect_equal(out$chill,chill)
  expect_equal(out$chill[,1], chill[,1])
  expect_equal(out$z, c(0,0,4,7,17))
  expect_equal(out$bloomindex, 5)
  expect_equal(run()$bloomindex, 5)
  expect_identical(chill,saved)
  expect_identical(colnames(out$chill),colnames(chill))
  # Allow one more interval of chilling: y=2 gives additional chill=1 and lowers the requirement to 9.
  later <- apply_partial_overlap_structure(chill,c(2,4,3,10),yc=1,b1=5,b2=8,
                                            b3=log(2),ol=1.4,stopatzc=FALSE,basic_output=FALSE)
  expect_equal(later$chill,chill)
  expect_equal(later$bloomindex,5)
  early <- apply_partial_overlap_structure(chill,c(2,4,3,10),yc=1,b1=5,b2=8,
                                            b3=log(2),ol=1.4,basic_output=FALSE)
  expect_equal(early$bloomindex,5)
  expect_equal(early$z,c(0,0,4,7,17))
  zero <- apply_partial_overlap_structure(chill,c(2,4,3,10),yc=1,ol=0,basic_output=FALSE)
  expect_equal(zero$chill,chill)
  expect_equal(zero$z,c(0,0,4,7,17))
  expect_equal(zero$bloomindex,0)
  expect_equal(apply_partial_overlap_structure(chill[,2,drop=FALSE],c(2,4,3,10),
                yc=1,b1=5,b2=8,b3=log(2),ol=0.8)$bloomindex,5)
})

test_that("parallel and partial-overlap adapters compose modules with their respective chill rules", {
  A1 <- exp(2/281)
  args <- list(temp=rep(8,6),times=0:5,A0=10*A1*exp(-1/281),A1=A1,E0=1,E1=2,Tf=8)
  chill <- do.call(calculate_chill_dynamic,args)
  heat <- calculate_heat_gdh(args$temp,args$times)
  cases <- list(list(adapter=parallel_model, structure=apply_parallel_structure,
                     parameters=list(yc=1,zc=1e6,kmin=0.1)),
                list(adapter=po_model,structure=apply_partial_overlap_structure,
                     parameters=list(yc=1,b1=1,b2=1e6,b3=0,ol=1)))
  for(case in cases) {
    parameters <- c(case$parameters,list(stopatzc=FALSE,basic_output=FALSE))
    actual <- do.call(case$adapter,c(args,parameters))
    expected <- do.call(case$structure,c(list(chill=chill,heat=heat),parameters))
    expect_gt(chill[2,"x"],1)
    expect_equal(actual$x,unname(chill[,"x"]))
    expect_equal(actual$y,unname(chill[,"y"]))
    expect_equal(actual$x,unname(expected$chill[,"x"]))
    expect_equal(actual$y,unname(expected$chill[,"y"]))
    expect_equal(actual$z,expected$z)
    expect_equal(actual$bloomindex,expected$bloomindex)
    expect_equal(actual$xs,c(chill[-1,"xs"],0))
  }
})

test_that("new structure modules validate their contracts and parameter domains", {
  chill <- matrix(c(0,1,2),ncol=1)
  for(fn in list(apply_parallel_structure,apply_partial_overlap_structure)) {
    expect_error(fn(matrix(0,1,1),numeric()),"at least two")
    expect_error(fn(matrix(0,3,0),c(1,1)),"one column")
    expect_error(fn(chill,1),"nrow")
    expect_error(fn(chill,c(1,-1)),"non-negative")
    expect_error(fn(chill,c(1,Inf)),"finite")
    expect_error(fn(matrix(c(0,NA,2),ncol=1),c(1,1)),"finite")
    expect_error(fn(chill,c(1,1),yc=0),"yc")
  }
  expect_error(apply_parallel_structure(chill,c(1,1),kmin=1.1),"kmin")
  expect_error(apply_parallel_structure(chill,c(1,1),zc=0),"zc")
  expect_error(apply_parallel_structure(matrix(c(0,-1,2),ncol=1),c(1,1)),"non-negative accumulated")
  expect_error(apply_partial_overlap_structure(chill,c(1,1),b1=0),"b1")
  for(name in c("b2","b3","ol"))
    expect_error(do.call(apply_partial_overlap_structure,c(list(chill=chill,heat=c(1,1)),
                                                           setNames(list(-1),name))),name)
  expect_error(apply_partial_overlap_structure(chill,c(1,1),b1=1e308,ol=1e308),"b1\\*ol")
})
