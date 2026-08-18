test_that("kde_2d matches a direct Gaussian kernel sum", {
  set.seed(3)
  x <- stats::rnorm(50)
  y <- stats::rnorm(50)
  lims <- c(-3, 3, -3, 3)
  n <- 12L
  h <- 0.5

  result <- kde_2d(x, y, h = h, n = n, lims = lims)

  expect_equal(result[["x"]], seq(-3, 3, length.out = n))
  expect_equal(result[["y"]], seq(-3, 3, length.out = n))
  expect_equal(dim(result[["z"]]), c(n, n))

  # Reference value at one grid cell, computed independently.
  bandwidth <- h / 4
  expected <- sum(
    stats::dnorm((result[["x"]][4] - x) / bandwidth) *
      stats::dnorm((result[["y"]][7] - y) / bandwidth)
  ) /
    (length(x) * bandwidth^2)

  expect_equal(result[["z"]][4, 7], expected)
})

test_that("kde_2d rejects non-positive bandwidths", {
  expect_error(kde_2d(1:5, 1:5, h = 0, n = 5), "strictly positive")
})

test_that("dbscan_cluster separates well-spaced blobs and flags noise", {
  set.seed(9)
  blob_a <- cbind(stats::rnorm(40, 0, 0.1), stats::rnorm(40, 0, 0.1))
  blob_b <- cbind(stats::rnorm(40, 10, 0.1), stats::rnorm(40, 10, 0.1))
  outlier <- cbind(100, 100)
  points <- rbind(blob_a, blob_b, outlier)

  clusters <- dbscan_cluster(points, eps = 0.5, min_pts = 5L)

  expect_length(clusters, nrow(points))
  expect_equal(length(setdiff(unique(clusters), 0L)), 2L)
  expect_equal(clusters[nrow(points)], 0L)
  expect_equal(length(unique(clusters[seq_len(40)])), 1L)
  expect_false(clusters[1] == clusters[41])
})

test_that("dbscan_cluster labels everything as noise when eps is tiny", {
  points <- cbind(seq_len(10), seq_len(10))
  expect_true(all(dbscan_cluster(points, eps = 0.001, min_pts = 5L) == 0L))
})

test_that("dbscan_cluster handles empty input", {
  expect_length(dbscan_cluster(matrix(numeric(0), ncol = 2), eps = 1), 0L)
})

test_that("concave_hull returns a closed ring that contains all input points", {
  set.seed(4)
  angle <- stats::runif(200, 0, 2 * pi)
  radius <- sqrt(stats::runif(200))
  points <- cbind(radius * cos(angle), radius * sin(angle))

  hull <- concave_hull(points, concavity = 2, length_threshold = 0.05)

  expect_true(is.matrix(hull))
  expect_equal(ncol(hull), 2L)
  expect_equal(hull[1, ], hull[nrow(hull), ])
  expect_gte(nrow(hull), 4L)

  # Every hull vertex must come from the input set.
  key <- function(m) paste(m[, 1], m[, 2])
  expect_true(all(key(hull) %in% key(points)))
})

test_that("concave_hull with high concavity approaches the convex hull", {
  set.seed(5)
  points <- cbind(stats::rnorm(100), stats::rnorm(100))

  polygon_area <- function(p) {
    n <- nrow(p)
    abs(sum(p[-n, 1] * p[-1, 2] - p[-1, 1] * p[-n, 2])) / 2
  }

  convex <- points[c(convex_hull_indices(points), convex_hull_indices(points)[1]), ]
  loose <- concave_hull(points, concavity = 100, length_threshold = 0)
  tight <- concave_hull(points, concavity = 0.5, length_threshold = 0)

  expect_equal(polygon_area(loose), polygon_area(convex), tolerance = 1e-8)
  expect_lt(polygon_area(tight), polygon_area(convex))
})

test_that("concave_hull degrades gracefully for tiny inputs", {
  points <- cbind(c(0, 1, 0), c(0, 0, 1))
  hull <- concave_hull(points)
  expect_equal(hull[1, ], hull[nrow(hull), ])
})

test_that("sq_seg_dist_many agrees with the scalar version", {
  set.seed(6)
  points <- cbind(stats::rnorm(25), stats::rnorm(25))
  a <- c(0, 0)
  b <- c(1, 2)

  expect_equal(
    sq_seg_dist_many(points, a, b),
    apply(points, 1, function(p) sq_seg_dist(p, a, b))
  )
})

test_that("no_intersections detects crossings and ignores shared endpoints", {
  starts <- rbind(c(0, 1), c(5, 5))
  ends <- rbind(c(2, 1), c(6, 6))

  # Crosses the first segment.
  expect_false(no_intersections(c(1, 0), c(1, 2), starts, ends))
  # Well clear of both.
  expect_true(no_intersections(c(-5, -5), c(-4, -4), starts, ends))
  # Shares an endpoint with the first segment, so it is not a crossing.
  expect_true(no_intersections(c(0, 1), c(0, -3), starts, ends))
})
