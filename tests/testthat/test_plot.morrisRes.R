context("Test of plot.morrisRes()")

test_that("plot.morrisRes() throws errors", {
  # An arbitrary Sobol'-result:
  arbit_res_sobol <- diag(7)
  class(arbit_res_sobol) <- "sobolRes"
  
  expect_error(plot.morrisRes(1:3))
  expect_error(plot.morrisRes("No character!"))
  expect_error(plot.morrisRes(diag(7)))
  expect_error(plot.morrisRes(arbit_res_sobol))
})
