## Shared defaults for numerical accuracy. Keep explicit function arguments as
## the user-facing controls; simulation validation uses these same defaults by
## omitting numerical overrides. abs.tol follows the chosen rel.tol.
.bfpwr_defaults <- list(
    ngrid = 10000,
    tol = 1e-8,
    rel.tol = 1e-8,
    subdivisions = 1000,
    tail.eps = 1e-6,
    tail.nquad = 512
)

## Vectorize() copies defaults into a closure in base, while evaluating them
## in the original package function. Declare the binding for static checks.
utils::globalVariables(".bfpwr_defaults")
