library(testthat)

test_that("Input sanitation and warnings work", {
  # Decimals should be dropped with a warning
  expect_warning(res <- logStirling2(n = c(5, 5.5)), "n coerced to natural numbers.")
  expect_equal(rownames(res), c("5", "5"))

  # Zero and negatives should be dropped
  expect_warning(logStirling2(n = c(10, 0, -5)), "n < 1 dropped.")

  # All invalid inputs should trigger the stopifnot
  expect_error(suppressWarnings(logStirling2(n = c(-1, 0.5))), "All n are less than 1.")
  expect_error(suppressWarnings(logStirling2(n = 5, k = 0)), "All k are less than 1.")
})

test_that("Single value dispatch (stirling2direct bypass) is correct", {
  # n = 10, k = 5
  res_mat <- logStirling2(10, 5, as.matrix = TRUE)
  res_vec <- logStirling2(10, 5, as.matrix = FALSE)

  expect_true(is.matrix(res_mat))
  expect_false(is.matrix(res_vec))
  expect_equal(as.numeric(res_mat), res_vec)
})

test_that("Small n (n < 3) bypass logic is correct", {
  # Tests the nu0 logic and dimension matching
  res <- logStirling2(n = 1:5, k = 1:2, as.matrix = TRUE)
  expect_equal(dim(res), c(5, 2))

  # Log(S(1,1)) and Log(S(2,1)) are 0
  expect_equal(unname(res[1:2, 1]), c(0, 0))
})

test_that("Backend C++ Routing Dispatches Correctly", {
  # We assume standard correctness, but we want to ensure the R logic
  # doesn't crash when hitting the different branches.

  # 1. Row_C (Single row in a cache block)
  expect_no_error(logStirling2(n = 1050, k = 5:6))

  # 2. All_C (Dense sequence crossing a block boundary)
  # 1000:1005 should trigger All_C for block 1000
  expect_no_error(res_all <- logStirling2(n = 998:1002, as.matrix = TRUE))
  expect_equal(nrow(res_all), 5)

  # 3. Mult_C (Sparse sequence in a block)
  expect_no_error(res_mult <- logStirling2(n = c(1010, 1050, 1090), as.matrix = TRUE))
  expect_equal(nrow(res_mult), 3)
})

test_that("Formatting flags shape the output correctly", {
  # k = NULL tests
  res_null_mat <- logStirling2(n = c(4, 5), k = NULL, as.matrix = TRUE)
  expect_equal(dim(res_null_mat), c(2, 5)) # max(n) columns

  # ones = FALSE test
  res_no_ones <- logStirling2(n = 5, k = NULL, as.matrix = FALSE, ones = FALSE)
  # S(5,1) and S(5,5) should be removed. Length should be 3.
  expect_equal(length(res_no_ones), 3)
})
