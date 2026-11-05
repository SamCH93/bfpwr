bfpwr_sim_grid <- function(name = c("short", "long")) {
    name <- match.arg(name)
    switch(name,
           short = seq(10, 500, by = 10),
           long = seq(10, 10000, by = 10))
}

bfpwr_sim_validate_look_grid <- function(look_grid) {
    stopifnot(
        is.numeric(look_grid),
        length(look_grid) >= 1,
        all(is.finite(look_grid)),
        all(look_grid >= 2),
        all(look_grid == round(look_grid)),
        all(diff(look_grid) > 0),
        all(look_grid %% 10 == 0)
    )
    invisible(look_grid)
}
