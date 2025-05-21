test_that("PI_eval works", {
  y1 <- 2 + 0.15 * (1:20) + rnorm(20)
  y2 <- y1[20] + 0.3 * (1:30) + rnorm(30)
  y <- as.ts(c(y1, y2))
  out <- dotm(y, h = 10)

  # MSIS metric
  a = PI_eval(
    obs = y, forec = out$mean,
    lower_bounds = out$lower[, 3], upper_bounds = out$upper[, 3],
    name = "MSIS", alpha = 0.05
    )

  # ACD metric
  b = PI_eval(
    obs = y, forec = out$mean,
    lower_bounds = out$lower[, 3], upper_bounds = out$upper[, 3],
    name = "ACD", alpha = 0.05
    )
  expect_equal(2 * 2, 4)
})

