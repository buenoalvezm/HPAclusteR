# Scale a data matrix

Scale a data matrix

## Usage

``` r
hc_scale(data, method = "sample")
```

## Arguments

- data:

  Data to be scaled, with samples in rows and genes in columns.

- method:

  Scaling method: `"sample"` (default), `"gene"`, `"pareto"`,
  `"vector"`, `"min-max"`, `"max_sequencial"` or `"max"`.

## Value

Scaled data matrix with genes in rows and samples in columns.
