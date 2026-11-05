bfpwr_sim_fixture_specs <- function(families = NULL,
                                    modes = NULL,
                                    tags = NULL) {
    specs <- c(
        bfpwr_sim_fixture_specs_z(),
        bfpwr_sim_fixture_specs_t(),
        bfpwr_sim_fixture_specs_binomial()
    )
    keep <- rep(TRUE, length(specs))
    if (!is.null(families)) {
        families <- match.arg(families, c("z", "t", "binomial"),
                              several.ok = TRUE)
        keep <- keep & vapply(specs, function(x) x$family %in% families,
                              logical(1))
    }
    if (!is.null(modes)) {
        modes <- match.arg(modes, c("fixed", "sequential"), several.ok = TRUE)
        keep <- keep & vapply(specs, function(x) x$mode %in% modes,
                              logical(1))
    }
    if (!is.null(tags)) {
        tags <- as.character(tags)
        keep <- keep & vapply(specs, function(x) all(tags %in% x$tags),
                              logical(1))
    }
    specs[keep]
}

bfpwr_sim_find_fixture_spec <- function(fixture_set_id,
                                        specs = bfpwr_sim_fixture_specs()) {
    ids <- vapply(specs, function(x) x$fixture_set_id, character(1))
    hit <- which(ids == fixture_set_id)
    if (length(hit) != 1) {
        stop("fixture spec not found or not unique: ", fixture_set_id)
    }
    specs[[hit]]
}
