# Download annotation databases

Download annotation databases

## Usage

``` r
get_annot_dbs(db_loc = "databases", hpa_version = NULL, timeout = 1800)
```

## Arguments

- db_loc:

  Destination directory to save the downloaded databases (default:
  `"databases"`).

- hpa_version:

  Human Protein Atlas release to download. `NULL` (default) uses the
  current release; a number such as `24` pins an archived release.

- timeout:

  Download timeout in seconds (default: 1800). These files are large,
  and R's 60 second default is not enough.

## Value

A named list mapping database IDs to file paths.
