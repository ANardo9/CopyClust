test_that("formatter works", {
  expect_equal(dim(CC_format(test_data_raw)), c(100, 478))
})
