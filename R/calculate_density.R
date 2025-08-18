get_density <-
  function(x, y, h = 0.5, n = 100, lims = c(
             range(x),
             range(y)
           )) {
    g_density <-
      MASS::kde2d(x, y,
        h = h, n = n,
        lims = lims
      )


    g_density$z %>%
      as_tibble(rownames = "x") %>%
      gather(y, z, -x) %>%
      mutate(y = gsub("V", "", y)) %>%
      mutate_at(vars(x, y), as.integer) %>%
      left_join(enframe(g_density$x, "x", "x_coord"),
        by = "x"
      ) %>%
      left_join(enframe(g_density$y, "y", "y_coord"),
        by = "y"
      )
  }
