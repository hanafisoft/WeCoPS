# Reference values computed with the Python version of WeCoPS (same algorithm).

test_that("sim1: weights, eigenvalues, motifs, saliences", {
  b <- load_simulated(1)$blocks; a <- wecops(b, 5)
  expect_equal(block_weights(b, "rv"), c(1.74192226, 1.5731582, 1.03682773), tolerance = 1e-6)
  expect_equal(block_weights(b, "mfa"), c(2.09395304, 1.70731014, 1.03694352), tolerance = 1e-6)
  expect_equal(block_weights(b, "statis"), c(1.23112117, 1.11198201, 0.03036317), tolerance = 1e-6)
  expect_equal(block_weights(b, "nvar"), c(0.05, 0.05, 0.05), tolerance = 1e-6)
  expect_equal(a$eigenvalues[1:6], c(0.348809, 0.319415, 0.153054, 0.071993, 0.050192, 0.006184), tolerance = 1e-5)
  expect_identical(a$motif, c('common 1-2', 'common 1-2', 'proper 1', 'proper 2', 'proper 3'))
  expect_identical(a$motif_pc, c('common 1-2', 'common 1-2-3', 'common 1-2-3', 'proper 1', 'proper 2'))
  expect_equal(max(a$mixing$mu_star), 0.154, tolerance = 2e-3)
  expect_equal(unname(a$saliences), matrix(c(0.4525, 0.003, 0.2566, 0.2131, 0.0022, 0.5671, 0.0046, 0.1982, 0.0021, 0.1531, 0.0108, 0.9641, 0.0015, 0.001, 0.0005), nrow = 3, byrow = TRUE), tolerance = 2e-3)
})

test_that("wines: weights, eigenvalues, motifs, saliences", {
  b <- load_wines()$blocks; a <- wecops(b, 3)
  expect_equal(block_weights(b, "rv"), c(1.29413883, 1.34921713, 1.20887299, 1.5268729), tolerance = 1e-6)
  expect_equal(block_weights(b, "mfa"), c(1.33385238, 1.45217685, 1.2283083, 1.90666628), tolerance = 1e-6)
  expect_equal(block_weights(b, "statis"), c(0.67134425, 0.75532312, 0.66347789, 0.52062377), tolerance = 1e-6)
  expect_equal(block_weights(b, "nvar"), c(0.25, 0.33333333, 0.25, 0.33333333), tolerance = 1e-6)
  expect_equal(a$eigenvalues[1:4], c(0.583018, 0.20716, 0.13696, 0.037344), tolerance = 1e-5)
  expect_identical(a$motif, c('common 1-2-3-4', 'common 1-2-3-4', 'proper 4'))
  expect_identical(a$motif_pc, c('common 1-2-3-4', 'common 1-2-3-4', 'common 1-4'))
  expect_equal(max(a$mixing$mu_star), 0.0136, tolerance = 2e-3)
  expect_equal(unname(a$saliences), matrix(c(0.6497, 0.1936, 0.037, 0.6632, 0.2489, 0.0104, 0.7818, 0.1289, 0.0407, 0.2345, 0.2539, 0.466), nrow = 4, byrow = TRUE), tolerance = 2e-3)
})

test_that("butters: weights, eigenvalues, motifs, saliences", {
  b <- load_butters()$blocks; a <- wecops(b, 4)
  expect_equal(block_weights(b, "rv"), c(1.69319531, 1.06148598, 1.31321948), tolerance = 1e-6)
  expect_equal(block_weights(b, "mfa"), c(2.14915708, 1.06276952, 1.34582781), tolerance = 1e-6)
  expect_equal(block_weights(b, "statis"), c(0.74056519, 0.67244718, 0.83818372), tolerance = 1e-6)
  expect_equal(block_weights(b, "nvar"), c(0.005, 0.00117647, 0.00055066), tolerance = 1e-6)
  expect_equal(a$eigenvalues[1:5], c(0.57213, 0.175357, 0.087472, 0.075289, 0.034535), tolerance = 1e-5)
  expect_identical(a$motif, c('common 1-2-3', 'proper 1', 'proper 1', 'proper 3'))
  expect_identical(a$motif_pc, c('common 1-2-3', 'common 1-3', 'common 1-3', 'proper 1'))
  expect_equal(max(a$mixing$mu_star), 0.2858, tolerance = 2e-3)
  expect_equal(unname(a$saliences), matrix(c(0.1756, 0.46, 0.0361, 0.2085, 0.9209, 0.0118, 0.0087, 0.0073, 0.5679, 0.0238, 0.2886, 0.0214), nrow = 3, byrow = TRUE), tolerance = 2e-3)
})

test_that("the 28 verdicts of the article", {
  r <- compare_weighting(load_simulated(1)$blocks, 5, "rv")  # Simulated 1
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 1.0, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(1)$blocks, 5, "mfa")  # Simulated 1
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9999, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(1)$blocks, 5, "statis")  # Simulated 1
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9487, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(1)$blocks, 5, "nvar")  # Simulated 1
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 1.0, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_simulated(2)$blocks, 10, "rv")  # Simulated 2
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9995, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(2)$blocks, 10, "mfa")  # Simulated 2
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.998, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(2)$blocks, 10, "statis")  # Simulated 2
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9998, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(2)$blocks, 10, "nvar")  # Simulated 2
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9965, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(3)$blocks, 5, "rv")  # Simulated 3
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.0132, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(3)$blocks, 5, "mfa")  # Simulated 3
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.0076, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(3)$blocks, 5, "statis")  # Simulated 3
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9995, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_simulated(3)$blocks, 5, "nvar")  # Simulated 3
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 1.0, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_simulated(4)$blocks, 5, "rv")  # Simulated 4
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.7134, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(4)$blocks, 5, "mfa")  # Simulated 4
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.6112, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(4)$blocks, 5, "statis")  # Simulated 4
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.7611, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_simulated(4)$blocks, 5, "nvar")  # Simulated 4
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 1.0, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_wines()$blocks, 3, "rv")  # Wines
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9998, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_wines()$blocks, 3, "mfa")  # Wines
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9992, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_wines()$blocks, 3, "statis")  # Wines
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9996, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_wines()$blocks, 3, "nvar")  # Wines
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9988, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_butters()$blocks, 4, "rv")  # Butters R4
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.997, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_butters()$blocks, 4, "mfa")  # Butters R4
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9848, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_butters()$blocks, 4, "statis")  # Butters R4
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9991, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_butters()$blocks, 4, "nvar")  # Butters R4
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.3714, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_butters()$blocks, 3, "rv")  # Butters R3
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.5562, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_butters()$blocks, 3, "mfa")  # Butters R3
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.3975, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
  r <- compare_weighting(load_butters()$blocks, 3, "statis")  # Butters R3
  expect_identical(r$equivalent, TRUE); expect_equal(min(r$kappa), 0.9952, tolerance = 2e-3); expect_identical(r$pc_same, TRUE)
  r <- compare_weighting(load_butters()$blocks, 3, "nvar")  # Butters R3
  expect_identical(r$equivalent, FALSE); expect_equal(min(r$kappa), 0.263, tolerance = 2e-3); expect_identical(r$pc_same, FALSE)
})

test_that("results do not depend on the signs of the principal components", {
  b <- load_simulated(4)$blocks
  r1 <- compare_weighting(b, 5, "rv")
  U <- cpca(b, R = 5)$U
  cc1 <- ccp(b, U); cc2 <- ccp(b, sweep(U, 2, c(1, -1, 1, -1, 1), `*`))
  expect_equal(cc1$S, cc2$S, tolerance = 1e-8)
})

test_that("preprocessing gives unit inertia", {
  b <- prep_blocks(list(matrix(rnorm(20), 4), matrix(rnorm(12), 4)))
  expect_equal(vapply(b, function(X) sum(X^2), numeric(1)), c(1, 1))
})

test_that("d_kh^2 equals the salience", {
  a <- wecops(load_simulated(1)$blocks, 5)
  for (k in 1:3) expect_equal(diag(a$factorization$ccp[[k]]$D)^2, unname(a$saliences[k, ]))
})

test_that("Hungarian assignment is optimal", {
  set.seed(1)
  perms <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
  for (it in 1:20) {
    n <- sample(2:6, 1); a <- matrix(rnorm(n * n), n)
    best <- min(vapply(perms(1:n), function(p) sum(a[cbind(1:n, p)]), numeric(1)))
    h <- WeCoPS:::.hungarian(a)
    expect_equal(sum(a[cbind(1:n, h)]), best)
  }
})

test_that("reweighting threshold bound holds", {
  r <- compare_weighting(load_simulated(1)$blocks, 5, "mfa")
  expect_true(r$bound$guaranteed)
  expect_true(all(r$kappa_components >= r$bound$kappa_lower_bound - 1e-10))
})

test_that("labels are those of the oriented components", {
  b <- load_wines()$blocks
  r <- compare_weighting(b, 3, "rv", alpha = 0.01)
  expect_true(r$equivalent)
  expect_identical(r$T_default, block_labels(r$saliences_default, 0.01))
})
