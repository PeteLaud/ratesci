library(ratesci)
context("Examples")

test_that("legacy & new methods match published examples", {

  # RISK DIFFERENCE

  # Newcombe RD example (d)
  # Miettinen-Nurminen
  expect_equal(
    unname(round(scoreci(x1 = 5, x2 = 0, n1 = 56, n2 = 29, skew = F)$estimates[, c(1, 3)], 4)),
    c(-0.0326, 0.1933)
  )
  # Mee
  expect_equal(
    unname(round(scoreci(x1 = 5, x2 = 0, n1 = 56, n2 = 29, skew = F, bcf = F)$estimates[, c(1, 3)], 4)),
    c(-0.0313, 0.1926)
  )
  # Newcombe/'Score'/Square&add
  expect_equal(
    unname(round(moverci(x1 = 5, x2 = 0, n1 = 56, n2 = 29, type = "wilson")[, c(1, 3)], 4)),
    c(-0.0381, 0.1926)
  )
  # MOVER-J
  expect_equal(
    unname(round(moverci(x1 = 5, x2 = 0, n1 = 56, n2 = 29, type = "jeff")[, c(1, 3)], 3)),
    c(-0.010, 0.177)
  )
  # Gart Nam 1990 for skewness corrected RD (unable to replicate other examples)
  expect_equal(
    unname(round(scoreci(x1 = 7, x2 = 2, n1 = 25, n2 = 25, bcf = F, skew = F)$estimates[, c(1, 3)], 3)),
    c(-0.014, 0.412)
  )
  expect_equal(
    unname(round(scoreci(x1 = 7, x2 = 2, n1 = 25, n2 = 25, bcf = F, skew = T)$estimates[, c(1, 3)], 3)),
    c(-0.014, 0.414)
  )

  # RELATIVE RISK

  # Gart Nam 1988 for skewness corrected RR
  expect_equal(
    unname(round(scoreci(x1 = c(8, 6), n1 = c(15, 10), x2 = c(4, 6), n2 = c(15, 20), contrast = "RR", RRtang = FALSE, bcf = F, skew = F)$estimates[, c(1, 3)], 3)),
    matrix(c(0.815, 5.336, 0.844, 4.594), byrow = T, nrow = 2)
  )
  expect_equal(
    unname(round(scoreci(x1 = c(8, 6), n1 = c(15, 10), x2 = c(4, 6), n2 = c(15, 20), contrast = "RR", RRtang = FALSE, bcf = F, skew = T)$estimates[, c(1, 3)], 3)),
    matrix(c(0.806, 6.148, 0.822, 4.954), byrow = T, nrow = 2)
  )
  # Li Tang Wong for Poisson RR (don't quite match!)
  expect_equal(
    unname(round(moverci(x1 = 15, x2 = 41, n1 = 19017, n2 = 28010, distrib = "poi", contrast = "RR", type = "jeff")[, c(1, 3)], 4)),
    c(0.2933, 0.9571)
  )
  expect_equal(
    unname(round(moverci(x1 = 60, x2 = 30, n2 = 54308.8, n1 = 51477.5, level = 0.9, distrib = "poi", contrast = "RR", type = "jeff")[, c(1, 3)], 4)),
    c(1.4668, 3.0586)
  )
  # Zou Donner for RR - they pre-round the Wilson interval to 3dps
  expect_equal(
    unname(round(moverci(x1 = 6, n1 = 22, x2 = 112, n2 = 842, distrib = "bin", contrast = "RR", type = "wilson")[, c(1, 3)], 4)),
    c(0.9735, 3.7295)
  )

  # Fagerland Newcombe for RR
  expect_equal(
    unname(round(moverci(x1 = c(24, 29, 7), n1 = c(73, 55, 18), x2 = c(53, 11, 1), n2 = c(65, 11, 18), contrast = "RR", type = "wilson")[, c(1, 3)], 3)),
    matrix(c(0.282, 0.563, 0.398, 0.761, 1.277, 40.751), byrow = T, nrow = 3)
  )
  expect_equal(
    unname(round(scoreci(x1 = c(24, 29, 7), n1 = c(73, 55, 18), x2 = c(53, 11, 1), n2 = c(65, 11, 18), contrast = "RR", RRtang = FALSE, skew = F)$estimates[, c(1, 3)], 3)),
    matrix(c(0.280, 0.560, 0.397, 0.742, 1.296, 42.544), byrow = T, nrow = 3)
  )

  # Tang alternative RR score examples, data from Mehrotra & Railkar
  expect_equal(
    unname(round(scoreci(
      x2 = c(26, 25, 11), n2 = c(453, 174, 70),
      x1 = c(21, 11, 8), n1 = c(464, 165, 69),
      contrast = "RR", RRtang = TRUE, skew = FALSE, bcf = FALSE,
      stratified = TRUE, weighting = "IVS"
    )$estimates[, c(1, 3)], 4)),
    c(0.4444, 0.9468)
  )
  expect_equal(
    unname(round(scoreci(
      x2 = c(26, 25, 11), n2 = c(453, 174, 70),
      x1 = c(21, 11, 8), n1 = c(464, 165, 69),
      contrast = "RR", RRtang = TRUE, skew = TRUE, bcf = FALSE,
      stratified = TRUE, weighting = "IVS"
    )$estimates[, c(1, 3)], 4)),
    c(0.4410, 0.9462)
  )

  # Tang OR score examples, data from Mehrotra & Railkar
  expect_equal(
    unname(round(scoreci(
      x2 = c(26, 25, 11), n2 = c(453, 174, 70),
      x1 = c(21, 11, 8), n1 = c(464, 165, 69),
      contrast = "OR", skew = FALSE, bcf = FALSE, ORbias = TRUE,
      stratified = TRUE
    )$estimates[, c(1, 3)], 4)),
    c(0.4155, 0.9466)
  )
  expect_equal(
    unname(round(scoreci(
      x2 = c(26, 25, 11), n2 = c(453, 174, 70),
      x1 = c(21, 11, 8), n1 = c(464, 165, 69),
      contrast = "OR", skew = TRUE, bcf = FALSE,
      stratified = TRUE
    )$estimates[, c(1, 3)], 4)),
    c(0.4124, 0.9461)
  )



  # Check CMH test equivalences stated by Tang
  x1 <- c(21, 11, 8)
  n1 <- c(464, 165, 69)
  x0 <- x2 <- c(26, 25, 11)
  n0 <- n2 <- c(453, 174, 70)
  event <- c(
    rep(1, x1[1]), rep(0, n1[1] - x1[1]), rep(1, x1[2]), rep(0, n1[2] - x1[2]), rep(1, x1[3]), rep(0, n1[3] - x1[3]),
    rep(1, x0[1]), rep(0, n0[1] - x0[1]), rep(1, x0[2]), rep(0, n0[2] - x0[2]), rep(1, x0[3]), rep(0, n0[3] - x0[3])
  )
  strat <- factor(c(rep(1, n1[1]), rep(2, n1[2]), rep(3, n1[3]), rep(1, n0[1]), rep(2, n0[2]), rep(3, n0[3])))
  trt <- factor(c(rep(1, sum(n1)), rep(0, sum(n0))))

  expect_equal(
    scoreci(x1, n1, x2, n2, stratified = T, contrast = "RD", skew = F, bcf = T, weighting = "MH", random = F)$pval[2],
    mantelhaen.test(x = event, y = trt, z = strat, correct = F)$p.value
  )
  expect_equal(
    scoreci(x1, n1, x2, n2, stratified = T, contrast = "RR", skew = F, bcf = T, weighting = "MH", random = F)$pval[2],
    mantelhaen.test(x = event, y = trt, z = strat, correct = F)$p.value
  )
  expect_equal(
    scoreci(x1, n1, x2, n2, stratified = T, contrast = "OR", skew = F, bcf = T, weighting = "INV", random = F)$pval[2],
    mantelhaen.test(x = event, y = trt, z = strat, correct = F)$p.value
  )



  # ODDS RATIO

  # Fagerland Newcombe for OR
  expect_equal(
    unname(round(moverci(
      x1 = c(24, 29, 7), n1 = c(73, 55, 18), x2 = c(53, 11, 1), n2 = c(65, 11, 18),
      contrast = "OR", type = "wilson"
    )[, c(1, 3)], 3)),
    matrix(c(0.050, 0.242, 0, 0.427, 1.395, 75.890), byrow = T, nrow = 3)
  )
  expect_equal(
    unname(round(scoreci(
      x1 = c(24, 29, 7), n1 = c(73, 55, 18), x2 = c(53, 11, 1), n2 = c(65, 11, 18),
      contrast = "OR", skew = FALSE, ORbias = FALSE
    )$estimates[, c(1, 3)], 3)),
    matrix(c(0.050, 0.245, 0, 0.417, 1.416, 76.428), byrow = T, nrow = 3)
  )

  # Gart 1985 Example 1
  # Gart
  expect_equal(
    unname(round(scoreci(x1 = 28, n1 = 99, x2 = 23, n2 = 182, contrast = "OR", bcf = F, skew = T, cc = T)$estimates[, c(1, 3)], 2)),
    c(1.40, 5.33)
  )
  # Cornfield
  expect_equal(
    unname(round(scoreci(
      x1 = 28, n1 = 99, x2 = 23, n2 = 182, contrast = "OR",
      ORbias = F, bcf = F, skew = F, cc = T
    )$estimates[, c(1, 3)], 2)),
    c(1.41, 5.30)
  )

  # Gart 1985 Example 2
  # Gart
  expect_equal(
    unname(round(scoreci(x1 = c(4, 2, 4, 1), n1 = c(16, 16, 18, 15), x2 = c(5, 3, 10, 3), n2 = c(79, 87, 90, 82), contrast = "OR", bcf = F, skew = T, cc = T, stratified = T)$estimates[, c(1, 3)], 2)),
    c(1.24, 7.12)
  )
  expect_equal(
    unname(round(scoreci(
      x1 = c(4, 2, 4, 1), n1 = c(16, 16, 18, 15), x2 = c(5, 3, 10, 3), n2 = c(79, 87, 90, 82),
      contrast = "OR", ORbias = F, bcf = F, skew = F, cc = T, stratified = T
    )$estimates[, c(1, 3)], 2)),
    c(1.29, 7.31)
  )
  # Cornfield examples - very close
  expect_equal(
    unname(round(scoreci(
      x1 = 3, n1 = 14, x2 = 60, n2 = 92, contrast = "OR", ORbias = F,
      bcf = F, skew = F, cc = T, level = pchisq(3.841, 1)
    )$estimates[, c(1, 3)], 4)),
    c(0.0296, 0.6228)
  )
  expect_equal(
    unname(round(scoreci(
      x1 = 3, n1 = 14, x2 = 60, n2 = 92, contrast = "OR", ORbias = F,
      bcf = F, skew = F, cc = T, level = 0.99
    )$estimates[, c(1, 3)], 4)),
    c(0.0209, 0.8797)
  )

  # Gart & Thomas 1982 Cornfield example
  expect_equal(
    unname(round(scoreci(
      x1 = 28, n1 = 51, x2 = 71, n2 = 230, contrast = "OR", ORbias = F,
      bcf = F, skew = F, cc = F
    )$estimates[, c(1, 3)], 3)),
    c(1.476, 5.035)
  )

  # Single proportion, Newcombe examples
  expect_equal(
    unname(round(exactci(x = c(15, 0, 1), n = c(148, 20, 29), distrib = "bin", midp = TRUE), 4)),
    matrix(c(0.0601, 0.1581, 0, 0.1391, 0.0017, 0.1585), byrow = T, nrow = 3)
  )

  expect_equal(
    unname(round(exactci(x = c(15, 0, 1), n = c(148, 20, 29), distrib = "bin", midp = FALSE), 4)),
    matrix(c(0.0578, 0.1617, 0, 0.1684, 0.0009, 0.1776), byrow = T, nrow = 3)
  )

  # Hartung & Knapp stratified RD example - Laud 2017 Table BII
  # SCAS
  expect_equal(
    unname(round(scoreci(
      x1 = c(15, 12, 29, 42, 14, 44, 14, 29, 10, 17, 38, 19, 21),
      x2 = c(9, 1, 18, 31, 6, 17, 7, 23, 3, 6, 12, 22, 19),
      n1 = c(16, 16, 34, 56, 22, 54, 17, 58, 14, 26, 44, 29, 38),
      n2 = c(16, 16, 34, 56, 22, 55, 15, 58, 15, 27, 45, 30, 38),
      stratified = TRUE, weighting = "IVS"
    )$estimates[, c(1, 3)], 4)),
    c(0.2441, 0.3698)
  )
  # TDAS
  expect_equal(
    unname(round(scoreci(
      x1 = c(15, 12, 29, 42, 14, 44, 14, 29, 10, 17, 38, 19, 21),
      x2 = c(9, 1, 18, 31, 6, 17, 7, 23, 3, 6, 12, 22, 19),
      n1 = c(16, 16, 34, 56, 22, 54, 17, 58, 14, 26, 44, 29, 38),
      n2 = c(16, 16, 34, 56, 22, 55, 15, 58, 15, 27, 45, 30, 38),
      stratified = TRUE, random = TRUE, skew = FALSE, weighting = "IVS"
    )$estimates[, c(1, 3)], 4)),
    c(0.1928, 0.4647)
  )

  # Tango
  expect_equal(
    unname(round(pairbinci(x = c(4, 9, 3, 16), contrast = "RD")$estimates[, c(1, 3)], 3)),
    c(-0.027, 0.390)
  )
  expect_equal(
    unname(round(pairbinci(x = c(43, 0, 1, 0), contrast = "RD", level = 0.90)$estimates[, c(1, 3)], 3)),
    c(-0.096, 0.037)
  )

  # Tang
  expect_equal(
    unname(round(pairbinci(x = c(446, 5, 16, 690), contrast = "RR", level = 0.9)$estimates[, c(1, 3)], 3)),
    c(0.958, 0.992) # nearly
  )
  expect_equal(
    unname(round(pairbinci(x = c(43, 0, 1, 0), contrast = "RR", level = 0.9)$estimates[, c(1, 3)], 3)),
    c(0.904, 1.039) # nearly
  )

  # paired methods examples from Agresti & Min 2005
  expect_equal(
    unname(round(pairbinci(x = c(53, 16, 8, 9), contrast = "RD")$estimates[, c(1, 3)], 3)),
    c(-0.020, 0.207)
  )
  expect_equal(
    unname(round(pairbinci(x = c(53, 16, 8, 9), contrast = "RD", method_RD = "TDAS")$estimates[, c(1, 3)], 4)),
    c(-0.0191, 0.2052)
  )
  expect_equal(
    unname(round(pairbinci(x = c(53, 16, 8, 9), contrast = "RD", method_RD = "Score")$estimates[, c(1, 3)], 3)),
    c(-0.020, 0.207)
  )
  expect_equal(
    unname(round(pairbinci(x = c(53, 16, 8, 9), contrast = "OR")$estimates[, c(1, 3)], 4)),
    c(0.8718, 4.8816)
  )

  # paired methods examples from Fagerland et al 2014
  expect_equal(
    unname(round(pairbinci(x = c(1, 1, 7, 12), contrast = "RD")$estimates[, c(1, 3)], 3)),
    c(-0.517, -0.026)
  )
  expect_equal(
    unname(round(pairbinci(x = c(1, 1, 7, 12), contrast = "RR")$estimates[, c(1, 3)], 3)),
    c(0.065, 0.907)
  )
  expect_equal(
    unname(round(pairbinci(x = c(1, 1, 7, 12), contrast = "OR", method_OR = "midp")$estimates, 3)),
    c(0.006, 0.924)
  )
})
