test_that("clusterer works", {
  expect_equal(dim(CopyClust(data_input = CC_format(test_data_raw))), c(100, 1))
})
