test_that("%||% falls back only on NULL", {
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% 2, 2)
  expect_equal(NA %||% 2, NA)
  expect_equal(list() %||% 2, list())
})

test_that("check_installed passes for present packages and errors otherwise", {
  expect_true(check_installed("stats"))
  expect_error(
    check_installed("a.package.that.does.not.exist"),
    "install.packages"
  )
  expect_error(
    check_installed("another.missing.package", bioc = TRUE),
    "BiocManager::install"
  )
})

test_that("hc_dist reproduces stats::dist for standard metrics", {
  set.seed(42)
  x <- matrix(stats::rnorm(120), nrow = 20)

  for (method in c("euclidean", "maximum", "manhattan", "canberra", "minkowski")) {
    expect_equal(
      as.vector(hc_dist(x, method = method)),
      as.vector(stats::dist(x, method = method)),
      info = method
    )
  }
})

test_that("hc_dist returns 1 - correlation for correlation metrics", {
  set.seed(42)
  x <- matrix(stats::rnorm(120), nrow = 20)

  for (method in c("pearson", "spearman", "kendall")) {
    expect_equal(
      as.vector(hc_dist(x, method = method)),
      as.vector(stats::as.dist(1 - stats::cor(t(x), method = method))),
      info = method
    )
  }
})

test_that("hc_dist rejects unknown metrics", {
  expect_error(hc_dist(matrix(1:6, nrow = 3), method = "nonsense"))
})

test_that("row_maxs matches apply", {
  set.seed(1)
  x <- matrix(stats::rnorm(50), nrow = 10)
  expect_equal(row_maxs(x), apply(x, 1, max))

  x[2, 3] <- NA
  expect_equal(row_maxs(x, na.rm = TRUE), apply(x, 1, max, na.rm = TRUE))
  expect_true(is.na(row_maxs(x)[2]))
})
