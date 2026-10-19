test_that("unscaled GDH has the expected shape, units and interval alignment", {
  temp <- c(-5,4,14.5,25,30.5,36,40,999)
  heat <- calculate_heat_gdh_unscaled(temp,seq_along(temp))
  expect_equal(heat,c(0,0,.5,1,1-sqrt(2)/2,0,0),tolerance=1e-12)
  expect_equal(calculate_heat_gdh(temp,seq_along(temp)),21*heat)
  expect_equal(calculate_heat_gdh_unscaled(temp+273,seq_along(temp),
                 Tb=277,Tu=298,Tc=309,deg_celsius=FALSE),heat)
  expect_error(calculate_heat_gdh_unscaled(c(5,6),c(0,2)),"one-hour")
  expect_error(calculate_heat_gdh_unscaled(c(5,6),0:1,Tb=25),"Tb < Tu < Tc")
  expect_error(calculate_heat_gdh_unscaled(c(NA,6),0:1),"finite")
})

test_that("PhenoFlex weights start-of-interval chill and preserves pool outputs", {
  chill <- cbind(x=100:104,final=c(0,1,2,4,8))
  original <- chill
  heat <- c(100,4,6,8)
  expected <- c(0,0,4*plogis(-1),4*plogis(-1)+3,
                4*plogis(-1)+3+8*plogis(.5))
  out <- apply_phenoflex_structure(chill,heat,yc=2,zc=4,s1=.5,
                                    stopatzc=FALSE,basic_output=FALSE)
  expect_equal(out$z,expected)
  expect_equal(out$bloomindex,4)
  expect_identical(out$chill,chill)
  expect_identical(chill,original)
  early <- apply_phenoflex_structure(chill,heat,yc=2,zc=4,s1=.5,basic_output=FALSE)
  expect_equal(early$z,c(expected[1:4],0))
  expect_true(all(early$chill[5,]==0))
  expect_identical(apply_phenoflex_structure(chill,heat,yc=2,zc=4),list(bloomindex=4L))
  expect_equal(apply_phenoflex_structure(chill,heat,yc=2,zc=100)$bloomindex,0)
  # Explicitly cover the legacy sigmoid saturation limits.
  capped <- apply_phenoflex_structure(matrix(c(.1,100,101),ncol=1),c(1,1),
                       yc=40,s1=.8,basic_output=FALSE)
  expect_equal(capped$z,c(0,0,1))
})

test_that("modular PhenoFlex matches chillR with unscaled GDH", {
  for (temp in list(rep(8,2000),rep(c(2,8,15,25,30),400))) {
    times <- seq_along(temp)
    chill <- calculate_chill_dynamic(temp,times)
    heat <- calculate_heat_gdh_unscaled(temp,times)
    for (stop in c(TRUE,FALSE)) {
      actual <- apply_phenoflex_structure(chill,heat,yc=10,zc=10,s1=.5,
                                          stopatzc=stop,basic_output=FALSE)
      reference <- chillR::PhenoFlex(temp=temp,times=times,Tu=25,Tb=4,Tc=36,
                             yc=10,zc=10,s1=.5,stopatzc=stop,basic_output=FALSE)
      # chillR overwrites bloomindex when stopatzc=FALSE; this API retains
      # the first crossing. Derive that event from the reference trajectory.
      crossings <- which(reference$z >= 10)
      expected_bloom <- if (length(crossings)) crossings[1] else 0L
      expect_equal(actual$bloomindex,expected_bloom)
      if (all(temp==8)) expect_gt(actual$bloomindex,0)
      if (stop) expect_equal(actual$bloomindex,reference$bloomindex)
      expect_equal(actual$z,reference$z,tolerance=1e-10)
      expect_equal(unname(actual$chill[,"x"]),reference$x,tolerance=1e-10)
      expect_equal(unname(actual$chill[,"y"]),reference$y,tolerance=1e-10)
      scaled <- apply_phenoflex_structure(chill,calculate_heat_gdh(temp,times),
                     yc=10,zc=210,s1=.5,stopatzc=stop,basic_output=FALSE)
      expect_equal(scaled$bloomindex,actual$bloomindex)
      expect_equal(scaled$z,21*actual$z,tolerance=1e-10)
    }
  }
})

test_that("PhenoFlex structure validates inputs", {
  chill <- matrix(c(0,1,2),ncol=1)
  for (parameter in c("yc","zc","s1"))
    for (bad in c(0,-1,Inf,NA_real_))
      expect_error(do.call(apply_phenoflex_structure,
        c(list(chill=chill,heat=c(1,1)),setNames(list(bad),parameter))),parameter)
  expect_error(apply_phenoflex_structure(matrix(c(0,-1,2),ncol=1),c(1,1)),"non-negative accumulated")
  expect_error(apply_phenoflex_structure(chill,c(1,NA)),"finite")
  expect_error(apply_phenoflex_structure(chill,c(1,-1)),"non-negative")
  expect_error(apply_phenoflex_structure(chill,1),"nrow")
  expect_error(apply_phenoflex_structure(matrix(c(0,NA,2),ncol=1),c(1,1)),"finite")
})
