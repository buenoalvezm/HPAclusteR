# Estimate 2D density

Estimate 2D density

## Usage

``` r
get_density(x, y, h = 0.5, n = 100, lims = c(range(x), range(y)))
```

## Arguments

- x:

  x coordinates

- y:

  y coordinates

- h:

  bandwidth

- n:

  number of grid points

- lims:

  limits

## Value

Tidy tibble with x, y, z (density), x_coord, y_coord
