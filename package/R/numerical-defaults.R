## Immutable package defaults. Session options and explicit arguments may
## override these values; abs.tol follows the chosen rel.tol.
.bfpwr_defaults <- list(
    ngrid = 10000,
    tol = 1e-8,
    rel.tol = 1e-8,
    subdivisions = 1000,
    tail.eps = 1e-6,
    tail.nquad = 512
)

## Validate the complete update before changing any session options.
.bfpwr_validate_options <- function(settings) {
    if (length(names(settings)) != length(settings) ||
        any(!names(settings) %in% names(.bfpwr_defaults)) ||
        anyDuplicated(names(settings))) {
        stop("settings must have unique names from bfpwrOptions()")
    }
    for (name in names(settings)) {
        value <- settings[[name]]
        valid <- length(value) == 1 && is.numeric(value) &&
            is.finite(value) && value > 0
        if (valid && name %in% c("ngrid", "subdivisions", "tail.nquad")) {
            valid <- value == floor(value) &&
                value >= (if (name == "tail.nquad") 2 else 1) &&
                value <= .Machine$integer.max
        }
        if (valid && name == "tail.eps") valid <- value < 0.5
        if (!valid) stop("invalid numerical setting: ", name)
    }
    invisible(settings)
}

#' @title Numerical settings for Bayes factor designs
#'
#' @description Inspect or change the numerical settings used in subsequent
#'     calculations. Explicit function arguments override these session settings.
#'
#' @param ... Named settings, or a single named list of settings. With no
#'     arguments, return the effective settings without changing them.
#'
#' @details The settings are \code{ngrid} (sequential integration grid, default
#'     \code{10000}), \code{tol} (sample-size and t-boundary root tolerance,
#'     \code{1e-8}), \code{rel.tol} (t-Bayes-factor integration tolerance,
#'     \code{1e-8}), \code{subdivisions} (integration limit, \code{1000}),
#'     \code{tail.eps} (omitted t-predictive tail mass, \code{1e-6}), and
#'     \code{tail.nquad} (t-tail quadrature nodes per dimension, \code{512}).
#'     The absolute integration tolerance follows \code{rel.tol} unless
#'     \code{abs.tol} is supplied to a calculation. The internal z-boundary
#'     tolerance is fixed separately to preserve accuracy near point priors.
#'
#'     Values are stored in R options named \code{bfpwr.ngrid},
#'     \code{bfpwr.tol}, and so on. Changes last for the current R session.
#'     Saved design objects retain the controls needed to reproduce their plots.
#'
#' @return A named list of effective settings. When setting values, the previous
#'     settings are returned invisibly, so they can be restored with another call.
#'
#' @examples
#' bfpwrOptions()
#' old <- bfpwrOptions(ngrid = 20000, rel.tol = 1e-9)
#' bfpwrOptions(old)
#'
#' @export
bfpwrOptions <- function(...) {
    settings <- list(...)
    if (length(settings) == 1 && is.null(names(settings)) &&
        is.list(settings[[1]])) settings <- settings[[1]]
    .bfpwr_validate_options(settings)
    current <- lapply(names(.bfpwr_defaults), function(name) {
        getOption(paste0("bfpwr.", name), .bfpwr_defaults[[name]])
    })
    names(current) <- names(.bfpwr_defaults)
    .bfpwr_validate_options(current)
    if (!length(settings)) return(current)
    names(settings) <- paste0("bfpwr.", names(settings))
    options(settings)
    invisible(current)
}
