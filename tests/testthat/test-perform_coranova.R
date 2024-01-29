test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("full coranova works", {
  expected_output = list( pB = 0.0004667585, pW = 8.749722e-60, pI = 0.8621334)
  expect_equal(perform_coranova_parametric(list(afr, eur), "pheno", c("pgs1", "pgs2", "pgs3")),expected_output,tolerance=1e-7)
})


test_that("one pop, two pgs coranova works", {
  expected_output = list(pW = 3.43991e-20 ,diff = 0.1596351, se = 0.01734359, LCB = 0.1256417, UCB = 0.1936286)
  expect_equal(perform_coranova_parametric(list(afr), "pheno", c("pgs1", "pgs2")),expected_output,tolerance=1e-5, ignore_attr = TRUE)
})

test_that("two pop, one pgs coranova works", {
  expected_output = list(pB = 0.04006947 ,diff = -0.03692067, se = 0.01798349, LCB = -0.07216831, UCB = -0.00167304)
  expect_equal(perform_coranova_parametric(list(afr, eur), "pheno", c("pgs1")),expected_output,tolerance=1e-5, ignore_attr = TRUE)
})

test_that("populate R works", {
  expected_output <- c(0.2986534, 0.1390183, 0.1193595, 0.3355741, 0.1886227, 0.1582048 )
  names(expected_output) <- c("pgs1", "pgs2", "pgs3", "pgs1", "pgs2", "pgs3")
  expect_equal(populate_mu(list(cor(afr), cor(eur)), "pheno", c("pgs1", "pgs2", "pgs3")),expected_output,tolerance=1e-7)
})

test_that("populate V one group works", {
var1 <- (1 - 0.2986534^2)^2/5000
var2 <- (1 - 0.139018264^2)^2/5000
var3 <- (1 - 0.119359495^2)^2/5000

cov12 <- (0.5*(2*0.1811254 - 0.2986534*0.139018264)*(1 - 0.1811254^2 - 0.2986534^2 - 0.139018264^2) + 0.1811254^3)/5000
cov23 <- (0.5*(2*0.009823649 - 0.139018264*0.119359495)*(1 - 0.009823649^2 - 0.139018264^2 - 0.119359495^2) + 0.009823649^3)/5000
cov13 <- (0.5*(2*0.1505211 - 0.2986534*0.119359495)*(1 - 0.1505211^2 - 0.2986534^2 - 0.119359495^2) + 0.1505211^3)/5000

output <- matrix(c(var1, cov12, cov13, cov12, var2, cov23, cov13, cov23, var3), ncol = 3, byrow = T)
expect_equal(populate_V_one_group(cor(afr), 5000, "pheno", c("pgs1","pgs2", "pgs3")), output,tolerance=1e-7 )
})

test_that("populate V different sizes works", {
  output <- bdiag(populate_V_one_group(cor(eur), 5000, "pheno", c("pgs1","pgs2", "pgs3")), populate_V_one_group(cor(eur2000), 2000, "pheno", c("pgs1","pgs2", "pgs3")))
  expect_equal(populate_sigma(list(cor(eur), cor(eur2000)), "pheno", c("pgs1","pgs2", "pgs3"), list(5000, 2000)), output,tolerance=1e-7 )
})
