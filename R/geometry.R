#' Two-dimensional kernel density estimate on a regular grid
#'
#' Base-R replacement for `MASS::kde2d()` with a fixed, isotropic bandwidth.
#' The bandwidth is quartered and used as the standard deviation of a Gaussian
#' kernel, exactly as `MASS` does.
#'
#' @param x,y Numeric coordinate vectors of equal length.
#' @param h Bandwidth, recycled to both dimensions.
#' @param n Number of grid points per dimension.
#' @param lims Numeric vector `c(xmin, xmax, ymin, ymax)`.
#'
#' @returns A list with grid coordinates `x` and `y` and the density matrix `z`.
#' @keywords internal
#' @noRd
kde_2d <- function(x, y, h, n = 25, lims = c(range(x), range(y))) {
  n_points <- length(x)
  grid_x <- seq.int(lims[1], lims[2], length.out = n)
  grid_y <- seq.int(lims[3], lims[4], length.out = n)

  bandwidth <- rep(h, length.out = 2L) / 4
  if (any(bandwidth <= 0)) {
    stop("Bandwidths must be strictly positive.", call. = FALSE)
  }

  kernel_x <- matrix(
    stats::dnorm(outer(grid_x, x, "-") / bandwidth[1]),
    nrow = n
  )
  kernel_y <- matrix(
    stats::dnorm(outer(grid_y, y, "-") / bandwidth[2]),
    nrow = n
  )

  list(
    x = grid_x,
    y = grid_y,
    z = tcrossprod(kernel_x, kernel_y) /
      (n_points * bandwidth[1] * bandwidth[2])
  )
}

#' Density-based spatial clustering (DBSCAN)
#'
#' Base-R replacement for `fpc::dbscan()`. Points are visited in index order;
#' core points seed a cluster and their density-reachable neighbourhood is
#' expanded. Points that belong to no cluster are labelled `0`, matching the
#' `fpc` convention for noise.
#'
#' @param coordinates Numeric matrix or data frame with one row per point.
#' @param eps Neighbourhood radius.
#' @param min_pts Minimum number of points, including the point itself,
#'   required for a point to be a core point.
#'
#' @returns Integer vector of cluster assignments, `0` for noise.
#' @keywords internal
#' @noRd
dbscan_cluster <- function(coordinates, eps, min_pts = 5L) {
  coordinates <- as.matrix(coordinates)
  n_points <- nrow(coordinates)
  if (n_points == 0L) {
    return(integer(0))
  }

  distances <- as.matrix(stats::dist(coordinates))
  neighbours <- lapply(seq_len(n_points), function(i) which(distances[i, ] <= eps))
  is_core <- lengths(neighbours) >= min_pts

  cluster <- integer(n_points)
  current <- 0L

  for (point in seq_len(n_points)) {
    if (cluster[point] != 0L || !is_core[point]) {
      next
    }

    current <- current + 1L
    cluster[point] <- current

    # Grow the cluster through chains of core points.
    queue <- neighbours[[point]]
    while (length(queue) > 0L) {
      candidate <- queue[1L]
      queue <- queue[-1L]

      if (cluster[candidate] != 0L) {
        next
      }
      cluster[candidate] <- current

      if (is_core[candidate]) {
        queue <- c(queue, neighbours[[candidate]][
          cluster[neighbours[[candidate]]] == 0L
        ])
      }
    }
  }

  cluster
}

#' Convex hull of a set of points
#'
#' @param points Two-column numeric matrix.
#'
#' @returns Integer vector of row indices tracing the hull counter-clockwise.
#' @keywords internal
#' @noRd
convex_hull_indices <- function(points) {
  rev(grDevices::chull(points[, 1], points[, 2]))
}

#' Squared Euclidean distance between two points
#' @keywords internal
#' @noRd
sq_dist <- function(a, b) {
  sum((a - b)^2)
}

#' Squared distance from a point to a line segment
#' @keywords internal
#' @noRd
sq_seg_dist <- function(p, a, b) {
  segment <- b - a
  length_sq <- sum(segment^2)

  projection <- if (length_sq > 0) {
    max(0, min(1, sum((p - a) * segment) / length_sq))
  } else {
    0
  }

  sum((p - (a + projection * segment))^2)
}

#' Squared distances from many points to one line segment
#'
#' Vectorised form of [sq_seg_dist()].
#'
#' @param points Two-column numeric matrix.
#' @param a,b Segment endpoints as length-2 numeric vectors.
#'
#' @returns Numeric vector of squared distances.
#' @keywords internal
#' @noRd
sq_seg_dist_many <- function(points, a, b) {
  segment <- b - a
  length_sq <- sum(segment^2)

  dx <- points[, 1] - a[1]
  dy <- points[, 2] - a[2]

  projection <- if (length_sq > 0) {
    pmax(0, pmin(1, (dx * segment[1] + dy * segment[2]) / length_sq))
  } else {
    rep(0, nrow(points))
  }

  (dx - projection * segment[1])^2 + (dy - projection * segment[2])^2
}

#' Orientation of an ordered point triple
#'
#' @returns Positive when `c` lies left of the directed line `a -> b`.
#' @keywords internal
#' @noRd
orientation <- function(a, b, c) {
  (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
}

#' Test whether a candidate segment crosses any existing outline segment
#'
#' Segments that share an endpoint with the candidate are ignored, matching the
#' identity checks in the reference implementation.
#'
#' @param p1,p2 Endpoints of the candidate segment.
#' @param starts,ends Two-column matrices of existing segment endpoints.
#'
#' @returns `TRUE` when no proper intersection exists.
#' @keywords internal
#' @noRd
no_intersections <- function(p1, p2, starts, ends) {
  shares_endpoint <- (starts[, 1] == p1[1] & starts[, 2] == p1[2]) |
    (starts[, 1] == p2[1] & starts[, 2] == p2[2]) |
    (ends[, 1] == p1[1] & ends[, 2] == p1[2]) |
    (ends[, 1] == p2[1] & ends[, 2] == p2[2])

  candidates <- which(!shares_endpoint)
  if (length(candidates) == 0L) {
    return(TRUE)
  }

  cross <- function(index) {
    q1 <- starts[index, ]
    q2 <- ends[index, ]
    (orientation(p1, p2, q1) > 0) != (orientation(p1, p2, q2) > 0) &&
      (orientation(q1, q2, p1) > 0) != (orientation(q1, q2, p2) > 0)
  }

  !any(vapply(candidates, cross, logical(1)))
}

#' Concave hull of a set of points
#'
#' Base-R implementation of the `concaveman` algorithm (Park & Oh). Starting
#' from the convex hull, each edge longer than `length_threshold` is "dug in"
#' towards the closest remaining point, provided that point lies within
#' `edge_length^2 / concavity^2` of the edge, is closer to this edge than to
#' either neighbouring edge, and that the two replacement edges do not cross the
#' current outline. Both replacement edges are then reconsidered in turn.
#'
#' Replaces the `concaveman` package, which depends on `V8` and `sf` and
#' therefore on the GDAL, GEOS and PROJ system libraries.
#'
#' @param points Two-column numeric matrix of coordinates.
#' @param concavity A relative measure of concavity; larger values give
#'   simpler, more convex outlines.
#' @param length_threshold Edges shorter than this are never subdivided.
#'
#' @returns A two-column matrix of hull vertices, closed (first point repeated).
#' @keywords internal
#' @noRd
concave_hull <- function(points, concavity = 2, length_threshold = 0) {
  points <- unique(as.matrix(points))
  dimnames(points) <- NULL
  n_points <- nrow(points)

  if (n_points < 4L) {
    return(rbind(points, points[1L, , drop = FALSE]))
  }

  hull_index <- convex_hull_indices(points)
  if (length(hull_index) < 3L) {
    return(rbind(points, points[1L, , drop = FALSE]))
  }

  concavity <- max(0, concavity)
  sq_concavity <- max(concavity^2, .Machine[["double.eps"]])
  sq_length_threshold <- length_threshold^2

  # Circular doubly linked list over node slots; `point_of` maps slot to point.
  n_slots <- length(hull_index)
  point_of <- hull_index
  next_of <- c(seq_len(n_slots)[-1], 1L)
  prev_of <- c(n_slots, seq_len(n_slots - 1L))

  available <- setdiff(seq_len(n_points), hull_index)
  queue <- seq_len(n_slots)

  while (length(queue) > 0L && length(available) > 0L) {
    node <- queue[1L]
    queue <- queue[-1L]

    a_point <- point_of[node]
    b_point <- point_of[next_of[node]]
    a <- points[a_point, ]
    b <- points[b_point, ]

    sq_len <- sq_dist(a, b)
    if (sq_len < sq_length_threshold) {
      next
    }
    max_sq_len <- sq_len / sq_concavity

    previous <- points[point_of[prev_of[node]], ]
    after_next <- points[point_of[next_of[next_of[node]]], ]

    available_points <- points[available, , drop = FALSE]
    distances <- sq_seg_dist_many(available_points, a, b)

    within_reach <- which(distances <= max_sq_len)
    if (length(within_reach) == 0L) {
      next
    }
    within_reach <- within_reach[order(distances[within_reach])]

    # Current outline segments, for the self-intersection test.
    slots <- seq_along(point_of)
    starts <- points[point_of[slots], , drop = FALSE]
    ends <- points[point_of[next_of[slots]], , drop = FALSE]

    candidate <- NULL
    for (index in within_reach) {
      p <- available_points[index, ]
      distance_to_edge <- distances[index]

      # Prefer points that clearly belong to this edge rather than a neighbour.
      if (
        distance_to_edge >= sq_seg_dist(p, previous, a) ||
          distance_to_edge >= sq_seg_dist(p, b, after_next)
      ) {
        next
      }

      if (
        no_intersections(a, p, starts, ends) &&
          no_intersections(b, p, starts, ends)
      ) {
        candidate <- index
        break
      }
    }

    if (is.null(candidate)) {
      next
    }

    p <- available_points[candidate, ]
    if (min(sq_dist(p, a), sq_dist(p, b)) > max_sq_len) {
      next
    }

    # Splice the chosen point in between `node` and its successor.
    new_slot <- length(point_of) + 1L
    point_of[new_slot] <- available[candidate]
    next_of[new_slot] <- next_of[node]
    prev_of[new_slot] <- node
    prev_of[next_of[node]] <- new_slot
    next_of[node] <- new_slot

    available <- available[-candidate]
    queue <- c(queue, node, new_slot)
  }

  # Walk the linked list back into an ordered, closed ring.
  order_of_slots <- integer(0)
  slot <- 1L
  repeat {
    order_of_slots <- c(order_of_slots, slot)
    slot <- next_of[slot]
    if (slot == 1L) {
      break
    }
  }

  points[c(point_of[order_of_slots], point_of[1L]), , drop = FALSE]
}
