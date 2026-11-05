library(tinytest)
library(bfpwr)

## Count actual setup calls in the installed namespace. Each full schedule
## prepares its predictive distribution, factor and grid once, then reuses
## their leading prefixes. Independent calls must prepare their own setup.
count_sequential_setup <- function(fun) {
    optionName <- "bfpwr.test.sequential_setup_counts"
    oldValue <- getOption(optionName)
    options(structure(list(c(predictive = 0L, factor = 0L, grid = 0L)),
                      names = optionName))
    targets <- list(
        list(name = "predpars", where = asNamespace("bfpwr"), key = "predictive"),
        list(name = "chol", where = baseenv(), key = "factor"),
        list(name = "ghalton", where = asNamespace("qrng"), key = "grid")
    )
    for (target in targets) {
        tracer <- substitute({
            counts <- getOption("bfpwr.test.sequential_setup_counts")
            counts[[key]] <- counts[[key]] + 1L
            options(bfpwr.test.sequential_setup_counts = counts)
        }, list(key = target$key))
        suppressMessages(trace(target$name, where = target$where,
                               tracer = tracer, print = FALSE))
    }
    on.exit({
        for (target in targets) {
            suppressMessages(untrace(target$name, where = target$where))
        }
        options(structure(list(oldValue), names = optionName))
    })
    result <- fun()
    list(result = result, counts = getOption(optionName))
}

zargs <- list(k1 = 0.1, k0 = 3, n = c(20, 40, 60, 80),
              se = 1/sqrt(c(20, 40, 60, 80)), pm = 0, psd = 1,
              dpm = 0.2, dpsd = 0.1, alternative = "greater", ngrid = 31)
z <- count_sequential_setup(function() do.call(pbf01seq, zargs))
expect_equal(z$counts, c(predictive = 1L, factor = 1L, grid = 1L),
             info = "a full z schedule prepares integration once")
targs <- list(k1 = 0.1, k0 = 3, n = c(20, 40, 60), dpm = 0.2,
              dpsd = 0.1, alternative = "greater", ngrid = 31)
t <- count_sequential_setup(function() do.call(ptbf01seq, targs))
expect_equal(t$counts, c(predictive = 1L, factor = 1L, grid = 1L),
             info = "a full t schedule reuses predictive parameters from tail setup")
again <- count_sequential_setup(function() {
    first <- do.call(pbf01seq, zargs)
    changed <- zargs
    changed$dpm <- -0.3
    changed$ngrid <- 19
    second <- do.call(pbf01seq, changed)
    list(first, second)
})
expect_equal(again$counts, c(predictive = 2L, factor = 2L, grid = 2L),
             info = "setup is scoped to each call and design")
expect_equal(again$result[[1]], z$result)
expect_false(isTRUE(all.equal(again$result[[1]]$cumpH1,
                              again$result[[2]]$cumpH1)))
single <- zargs
single$n <- single$n[1]
single$se <- single$se[1]
one <- count_sequential_setup(function() do.call(pbf01seq, single))
expect_equal(one$counts, c(predictive = 1L, factor = 0L, grid = 0L))

## Alternate integration controls bypass the deterministic grid and factor.
for (args in list(zargs, targs)) {
    args$method <- "pmvnorm"
    args$algorithm <- mvtnorm::GenzBretz(maxpts = 1234, abseps = 1e-7,
                                        releps = 1e-7)
    fun <- if (is.null(args$se)) ptbf01seq else pbf01seq
    alternate <- count_sequential_setup(function() do.call(fun, args))
    expect_equal(alternate$counts, c(predictive = 1L, factor = 0L, grid = 0L))
    expect_equal(alternate$result$integration$method, "pmvnorm")
    expect_equal(alternate$result$integration$algorithm, args$algorithm)
}

## A supplied stage callback belongs to adaptive/search callers and must not
## trigger the full-design integration setup.
local({
    callbacks <- integer()
    evalStage <- function(i, bounds) {
        callbacks <<- c(callbacks, i)
        expect_equal(length(bounds), i)
        list(pH1 = i/100, pH0 = i/200)
    }
    zbuild <- count_sequential_setup(function() {
        bfpwr:::.bfseq_build_z_design(
            k1 = 0.1, k0 = 3, se = c(0.3, 0.2, 0.1), n = c(10, 20, 30),
            pm = 0, psd = 1, dpm = 0.2, dpsd = 0.1, type = "normal",
            strict = TRUE, dots = list(ngrid = 31), evalStage = evalStage)
    })
    expect_equal(callbacks, 1:3)
    expect_equal(zbuild$counts, c(predictive = 0L, factor = 0L, grid = 0L))
    expect_equal(zbuild$result$cumpH1, cumsum((1:3)/100))
    callbacks <- integer()
    bounds <- lapply(1:3, function(i) {
        list(se = 1/sqrt(10*i), n1 = 10*i, n2 = 10*i,
             zk0 = -1, zk1 = 2)
    })
    tbuild <- count_sequential_setup(function() {
        bfpwr:::.bfseq_build_t_design(
            k1 = 0.1, k0 = 3, bounds = bounds, dpm = 0.2, dpsd = 0.1,
            plocation = 0, pscale = 1, pdf = 1, alternative = "greater",
            type = "one.sample", trange = "adaptive", strict = TRUE,
            tail.eps = 1e-6, tail.nquad = 512, dots = list(ngrid = 31),
            evalStage = evalStage)
    })
    expect_equal(callbacks, 1:3)
    expect_equal(tbuild$counts, c(predictive = 0L, factor = 0L, grid = 0L))
    expect_equal(tbuild$result$cumpH0, cumsum((1:3)/200))
})
