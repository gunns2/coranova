test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("full coranova works", {
  expected_output = list(pB = 0.0004667585, pW = 8.749722e-60, pI = 0.8621334)
  expect_equal(perform_coranova_parametric(list(afr, eur), "pheno", c("pgs1", "pgs2", "pgs3")),expected_output,tolerance=1e-7)
})

#need to add more


