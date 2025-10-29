context("test-utils.R")

devtools::document("/home/yyx/R/Project/R_code/SigBridgeR")

test_that("test-utils.R", {
    mat = matrix(1:100, nrow = 10)
    expect_equal(matrixStats::rowSums(mat), rowSums2(mat))
})
