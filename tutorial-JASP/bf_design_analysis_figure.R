## Bayes factor design-analysis figure for the JASP tutorial.
##
## The bfpwr package functions are BF01-oriented. This script plots BF10 =
## 1 / BF01, so upward crossings support the one-sided effect hypothesis and
## downward crossings support the point-null hypothesis.

required_packages <- c("ggplot2", "patchwork")

check_packages <- function(packages) {
    missing_packages <- packages[!vapply(packages, requireNamespace,
                                         logical(1), quietly = TRUE)]
    if (length(missing_packages) > 0) {
        stop("Install required package(s): ",
             paste(missing_packages, collapse = ", "))
    }
}

script_directory <- function() {
    script_file <- commandArgs(trailingOnly = FALSE)
    script_file <- script_file[grepl("^--file=", script_file)]
    if (length(script_file) == 1) {
        return(dirname(normalizePath(sub("^--file=", "", script_file),
                                     winslash = "/")))
    }
    normalizePath(getwd(), winslash = "/")
}

load_bfpwr_functions <- function(root) {
    package_r_dir <- file.path(root, "package", "R")
    if (dir.exists(package_r_dir)) {
        r_files <- sort(list.files(package_r_dir, pattern = "[.]R$",
                                   full.names = TRUE))
        invisible(lapply(r_files, source))
    } else {
        library(bfpwr)
    }
}

check_packages(required_packages)

script_dir <- script_directory()
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/",
                           mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "package", "R")) &&
        dir.exists(file.path(script_dir, "package", "R"))) {
    repo_root <- script_dir
}
output_dir <- script_dir
load_bfpwr_functions(repo_root)

## Design and plotting settings ----------------------------------------------

set.seed(20260523)

n_trajectories <- 100
n_final <- 100
look_n <- seq(5, n_final, by = 5)
sample_size_breaks <- seq(0, n_final, by = 20)
x_padding <- 3
x_plot_limits <- c(0, n_final + x_padding)
plot_text_scale <- 1.5
axis_text_size <- 8.8 * plot_text_scale
trajectory_base_size <- 11 * plot_text_scale
margin_base_size <- 9 * plot_text_scale
strip_axis_text_size <- axis_text_size
strip_axis_title_size <- axis_text_size
trajectory_alpha <- 0.5
probability_label_size <- 2.45 * plot_text_scale
strip_probability_label_size <- 2.35 * plot_text_scale
condition_title_size <- 11 * plot_text_scale * 0.8
endpoint_jitter_width <- 0.75
sequential_panel_margin <- 10
sequential_undecided_bins <- 40
effects <- data.frame(
    effect = c("d = 0.5", "d = 0"),
    delta = c(0.5, 0),
    stringsAsFactors = FALSE
)

bf10_upper <- 10
bf10_lower <- 1 / 10
posterior_plot_limits <- c(0, 1)
bf10_breaks <- c(1 / 10, 1 / 3, 1, 3, 10)
bf10_labels <- c("1/10", "1/3", "1", "3", "10")

bf10_to_posterior <- function(bf10) {
    out <- bf10 / (1 + bf10)
    out[is.infinite(bf10) & bf10 > 0] <- 1
    pmin(pmax(out, 0), 1)
}

posterior_to_bf10 <- function(posterior) {
    posterior / (1 - posterior)
}

analysis_args <- list(
    plocation = 0,
    pscale = 1 / sqrt(2),
    pdf = 1,
    type = "one.sample",
    alternative = "greater"
)

analysis_design_grid <- function(delta, nodes = 21) {
    if (delta == 0) {
        return(data.frame(d = 0, weight = 1))
    }
    p <- (seq_len(nodes) - 0.5) / nodes
    data.frame(
        d = analysis_args$plocation +
            analysis_args$pscale * stats::qt((1 + p) / 2,
                                             df = analysis_args$pdf),
        weight = rep(1 / nodes, nodes)
    )
}

evidence_levels <- c("Evidence for effect", "Undecided", "Evidence for null")
evidence_cols <- c(
    "Evidence for effect" = "#2563EB",
    "Undecided" = "#8A8F98",
    "Evidence for null" = "#F97316"
)

## Evidence classification ----------------------------------------------------

evidence_zone <- function(bf10) {
    out <- ifelse(bf10 >= bf10_upper, "Evidence for effect",
                  ifelse(bf10 <= bf10_lower, "Evidence for null", "Undecided"))
    factor(out, levels = evidence_levels)
}

status_label <- function(evidence, delta) {
    evidence <- as.character(evidence)
    delta <- rep(delta, length.out = length(evidence))
    null_is_true <- delta == 0
    out <- ifelse(
        null_is_true,
        ifelse(evidence == "Evidence for null", "Conclusive",
               ifelse(evidence == "Evidence for effect", "Misleading",
                      "Undecided")),
        ifelse(evidence == "Evidence for effect", "Conclusive",
               ifelse(evidence == "Evidence for null", "Misleading",
                      "Undecided"))
    )
    factor(out, levels = c("Conclusive", "Undecided", "Misleading"))
}

## Bayes factor calculations --------------------------------------------------

bf10_one_sample <- function(x, n) {
    t_stat <- sqrt(n) * mean(x[seq_len(n)]) / stats::sd(x[seq_len(n)])
    1 / do.call(tbf01, c(list(t = t_stat, n = n), analysis_args))
}

draw_design_delta <- function(delta) {
    design <- analysis_design_grid(delta)
    sample(design$d, size = 1, prob = design$weight)
}

simulate_one_trajectory <- function(delta, effect, id) {
    true_delta <- draw_design_delta(delta)
    y <- stats::rnorm(n_final, mean = true_delta, sd = 1)
    bf10 <- vapply(look_n, function(n) bf10_one_sample(y, n), numeric(1))
    data.frame(
        effect = effect,
        delta = delta,
        true_delta = true_delta,
        id = id,
        n = c(0, look_n),
        bf10 = c(1, bf10)
    )
}

find_stop <- function(path) {
    crossed <- which(path$n > 0 & (path$bf10 >= bf10_upper |
                                      path$bf10 <= bf10_lower))
    if (length(crossed) == 0) {
        stop_i <- nrow(path)
    } else {
        stop_i <- crossed[1]
    }
    path[stop_i, c("effect", "delta", "id", "n", "bf10")]
}

## Analytical design-analysis probabilities ----------------------------------

fixed_probabilities <- function(delta) {
    design <- analysis_design_grid(delta)
    p_effect <- sum(design$weight * vapply(design$d, function(d) {
        suppressWarnings(do.call(ptbf01, c(
            list(k = 1 / bf10_upper, n = n_final, dpm = d, dpsd = 0),
            analysis_args
        )))
    }, numeric(1)))
    p_null <- sum(design$weight * vapply(design$d, function(d) {
        suppressWarnings(do.call(ptbf01, c(
            list(k = 1 / bf10_lower, n = n_final, dpm = d, dpsd = 0,
                 lower.tail = FALSE),
            analysis_args
        )))
    }, numeric(1)))
    data.frame(
        evidence = factor(evidence_levels, levels = evidence_levels),
        probability = c(p_effect, 1 - p_effect - p_null, p_null)
    )
}

sequential_stop_curve <- function(delta) {
    design <- analysis_design_grid(delta)
    effect_cumulative <- null_cumulative <- numeric(length(look_n))
    for (i in seq_len(nrow(design))) {
        seq_design <- suppressWarnings(do.call(ptbf01seq, c(
            list(k1 = 1 / bf10_upper, k0 = 1 / bf10_lower, n = look_n,
                 dpm = design$d[i], dpsd = 0),
            analysis_args
        )))
        effect_cumulative <- effect_cumulative +
            design$weight[i] * seq_design$cumpH1
        null_cumulative <- null_cumulative +
            design$weight[i] * seq_design$cumpH0
    }
    data.frame(
        n = look_n,
        effect_cumulative = effect_cumulative,
        null_cumulative = null_cumulative,
        effect_mass = diff(c(0, effect_cumulative)),
        null_mass = diff(c(0, null_cumulative))
    )
}

fixed_density <- function(delta, grid_length = 120) {
    design <- analysis_design_grid(delta)
    posterior_grid <- seq(0.001, 0.999, length.out = grid_length)
    bf10_grid <- posterior_to_bf10(posterior_grid)
    cdf_grid <- vapply(bf10_grid, function(bf10_value) {
        sum(design$weight * vapply(design$d, function(d) {
            suppressWarnings(do.call(ptbf01, c(
                list(k = 1 / bf10_value, n = n_final, dpm = d, dpsd = 0,
                     lower.tail = FALSE),
                analysis_args
            )))
        }, numeric(1)))
    }, numeric(1))
    cdf_grid <- pmin(pmax(cdf_grid, 0), 1)
    cdf_grid <- cummax(cdf_grid)
    posterior_mid <- (posterior_grid[-1] +
        posterior_grid[-length(posterior_grid)]) / 2
    density <- diff(cdf_grid) / diff(posterior_grid)
    density[!is.finite(density)] <- 0
    density <- pmax(density, 0)
    data.frame(
        delta = delta,
        posterior = posterior_mid,
        bf10 = posterior_to_bf10(posterior_mid),
        density = density,
        evidence = evidence_zone(posterior_to_bf10(posterior_mid))
    )
}

## Labels, counts, and plotting constants -------------------------------------

format_probability <- function(x) {
    ifelse(x < 0.005, "<1%", paste0(round(100 * x), "%"))
}

probability_counts <- function(probabilities, n = n_trajectories) {
    probabilities <- pmax(probabilities, 0)
    probabilities <- probabilities / sum(probabilities)
    raw_counts <- n * probabilities
    counts <- floor(raw_counts)
    remainder <- n - sum(counts)
    if (remainder > 0) {
        add_to <- order(raw_counts - counts, decreasing = TRUE)[seq_len(remainder)]
        counts[add_to] <- counts[add_to] + 1
    }
    counts
}

probability_status_label <- function(probability, evidence, delta) {
    paste(as.character(status_label(evidence, delta)),
          paste0("(", format_probability(probability), ")"),
          sep = "\n")
}

probability_status_label_inline <- function(probability, evidence, delta) {
    paste(as.character(status_label(evidence, delta)),
          paste0("(", format_probability(probability), ")"))
}

zone_midpoints <- data.frame(
    evidence = factor(evidence_levels, levels = evidence_levels),
    ymin = c(bf10_to_posterior(bf10_upper),
             bf10_to_posterior(bf10_lower),
             posterior_plot_limits[1]),
    ymax = c(posterior_plot_limits[2],
             bf10_to_posterior(bf10_upper),
             bf10_to_posterior(bf10_lower))
)
zone_midpoints$y <- (zone_midpoints$ymin + zone_midpoints$ymax) / 2

## Plotting helpers -----------------------------------------------------------

transparent_background_theme <- function() {
    ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "transparent",
                                                colour = NA),
        panel.background = ggplot2::element_rect(fill = "transparent",
                                                 colour = NA),
        legend.background = ggplot2::element_rect(fill = "transparent",
                                                  colour = NA),
        legend.box.background = ggplot2::element_rect(fill = "transparent",
                                                      colour = NA)
    )
}

make_condition_title <- function(hypothesis) {
    label <- if (hypothesis == 1) "Under~H[1]" else "Under~H[0]"
    ggplot2::ggplot() +
        ggplot2::annotate(
            "text",
            x = 0,
            y = 0.5,
            label = label,
            parse = TRUE,
            hjust = 0,
            vjust = 0.5,
            size = condition_title_size / 2.845276,
            fontface = "bold",
            colour = "#111827"
        ) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1),
                                 clip = "off") +
        ggplot2::theme_void(base_size = condition_title_size) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(
                0, sequential_panel_margin, 0, sequential_panel_margin
            )
        ) +
        transparent_background_theme()
}

deterministic_unit_interval <- function(x) {
    x <- sin(x * 12.9898) * 43758.5453
    x - floor(x)
}

add_endpoint_jitter <- function(points, width = endpoint_jitter_width) {
    points$x_endpoint <- points$n
    if (nrow(points) <= 1) {
        return(points)
    }

    evidence_index <- match(as.character(points$evidence), evidence_levels)
    key <- points$id * 1009 + points$n * 917 + evidence_index * 101
    jitter_order <- deterministic_unit_interval(key)
    point_groups <- split(seq_len(nrow(points)), points$n)

    for (group in point_groups) {
        group_size <- length(group)
        if (group_size <= 1) {
            next
        }
        ranks <- rank(jitter_order[group], ties.method = "first")
        offsets <- ((ranks - 0.5) / group_size - 0.5) * 2 * width
        points$x_endpoint[group] <- points$n[group] + offsets
    }

    points
}

make_trajectory_plot <- function(paths, final_points, sequential = FALSE,
                                 show_y = TRUE, show_x = TRUE) {
    final_points <- add_endpoint_jitter(final_points)

    p <- ggplot2::ggplot() +
        ggplot2::geom_hline(yintercept = bf10_to_posterior(c(bf10_lower,
                                                              bf10_upper)),
                            linewidth = 0.35, linetype = "22",
                            colour = "#343A40") +
        ggplot2::geom_hline(yintercept = bf10_to_posterior(1),
                            linewidth = 0.25,
                            colour = "#6B7280") +
        ggplot2::scale_y_continuous(
            limits = posterior_plot_limits,
            breaks = bf10_to_posterior(bf10_breaks),
            labels = bf10_labels
        ) +
        ggplot2::scale_x_continuous(
            breaks = sample_size_breaks,
            expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(xlim = x_plot_limits,
                                 ylim = posterior_plot_limits,
                                 clip = "on") +
        ggplot2::scale_fill_manual(values = evidence_cols,
                                    limits = evidence_levels,
                                    drop = FALSE,
                                    guide = "none") +
        ggplot2::scale_colour_manual(values = evidence_cols,
                                      limits = evidence_levels,
                                      drop = FALSE,
                                      guide = "none") +
        ggplot2::labs(
            x = if (show_x) "Sample Size" else NULL,
            y = if (show_y) expression(BF[10]) else NULL
        ) +
        ggplot2::theme_minimal(base_size = trajectory_base_size) +
        ggplot2::theme(
            axis.text = ggplot2::element_text(size = axis_text_size),
            axis.title = ggplot2::element_text(size = axis_text_size),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            legend.position = "none",
            plot.margin = ggplot2::margin(6, 6, 6, 6)
        ) +
        transparent_background_theme()

    if (sequential) {
        p <- p +
            ggplot2::geom_line(
                data = paths,
                ggplot2::aes(x = n, y = bf10_to_posterior(bf10),
                             group = id, colour = evidence),
                linewidth = 0.275,
                alpha = trajectory_alpha,
                lineend = "round"
            ) +
            ggplot2::geom_point(
                data = final_points,
                ggplot2::aes(x = x_endpoint, y = bf10_to_posterior(bf10),
                             fill = evidence),
                shape = 21,
                colour = "white",
                stroke = 0.35,
                size = 2.3
            )
    } else {
        p <- p +
            ggplot2::geom_line(
                data = paths,
                ggplot2::aes(x = n, y = bf10_to_posterior(bf10), group = id),
                linewidth = 0.225,
                alpha = trajectory_alpha,
                colour = "#4B5563",
                lineend = "round"
            ) +
            ggplot2::geom_point(
                data = final_points,
                ggplot2::aes(x = x_endpoint, y = bf10_to_posterior(bf10),
                             fill = evidence),
                shape = 21,
                colour = "white",
                stroke = 0.35,
                size = 2.3
            )
    }

    if (!show_y) {
        p <- p + ggplot2::theme(
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        )
    }
    if (!show_x) {
        p <- p + ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank()
        )
    }
    p
}

make_fixed_margin <- function(delta, probabilities) {
    density <- fixed_density(delta)
    label_data <- merge(zone_midpoints, probabilities, by = "evidence",
                        all.x = TRUE, sort = FALSE)
    label_data$status <- status_label(label_data$evidence, delta)
    label_data$label <- probability_status_label(
        probability = label_data$probability,
        evidence = label_data$evidence,
        delta = delta
    )
    upper_boundary <- bf10_to_posterior(bf10_upper)
    lower_boundary <- bf10_to_posterior(bf10_lower)
    label_data$y[label_data$status == "Conclusive" &
                     label_data$evidence == "Evidence for effect"] <-
        upper_boundary + 0.32 * (posterior_plot_limits[2] - upper_boundary)
    label_data$y[label_data$status == "Conclusive" &
                     label_data$evidence == "Evidence for null"] <-
        lower_boundary * 0.72
    upper_evidence_labels <- label_data$status %in% c("Conclusive",
                                                      "Misleading") &
        label_data$evidence == "Evidence for effect"
    lower_evidence_labels <- label_data$status %in% c("Conclusive",
                                                      "Misleading") &
        label_data$evidence == "Evidence for null"
    label_data$y[upper_evidence_labels] <- pmax(
        upper_boundary,
        label_data$y[upper_evidence_labels] -
            0.20 * (posterior_plot_limits[2] - upper_boundary)
    )
    label_data$y[lower_evidence_labels] <- pmin(
        lower_boundary,
        label_data$y[lower_evidence_labels] +
            0.20 * (lower_boundary - posterior_plot_limits[1])
    )
    max_density <- max(density$density, na.rm = TRUE)
    if (!is.finite(max_density) || max_density == 0) max_density <- 1
    label_data$x <- max_density * 2.55

    ggplot2::ggplot(density) +
        ggplot2::geom_segment(
            ggplot2::aes(x = 0, xend = density,
                         y = posterior, yend = posterior,
                         colour = evidence),
            linewidth = 0.8,
            alpha = 0.75
        ) +
        ggplot2::geom_text(
            data = label_data,
            ggplot2::aes(x = x, y = y, label = label,
                         colour = evidence),
            hjust = 1,
            size = probability_label_size,
            fontface = "bold",
            show.legend = FALSE
        ) +
        ggplot2::scale_y_continuous(
            limits = posterior_plot_limits,
            breaks = bf10_to_posterior(bf10_breaks),
            labels = NULL
        ) +
        ggplot2::scale_x_continuous(limits = c(0, max_density * 2.7),
                                    expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = posterior_plot_limits, clip = "off") +
        ggplot2::scale_colour_manual(values = evidence_cols,
                                      limits = evidence_levels,
                                      drop = FALSE,
                                      guide = "none") +
        ggplot2::theme_void(base_size = margin_base_size) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(6, 2, 6, 2)
        ) +
        transparent_background_theme()
}

integrate_final_regions <- function(final_regions, mean, sigma,
                                    method = "lpmvnorm", ...) {
    stopifnot(
        is.list(final_regions),
        is.numeric(mean),
        is.matrix(sigma),
        length(mean) == nrow(sigma),
        nrow(sigma) == ncol(sigma)
    )

    m <- length(mean)
    if (m == 1) {
        return(vapply(final_regions, function(region) {
            if (any(is.nan(region))) {
                return(0)
            }
            diff(stats::pnorm(q = c(region[1, ], region[2, ]),
                              mean = mean[1],
                              sd = sqrt(sigma[1, 1])))
        }, numeric(1)))
    }

    if (method == "lpmvnorm") {
        C <- t(chol(sigma))
        Ct <- mvtnorm::ltMatrices(C[lower.tri(C, diag = TRUE)], diag = TRUE)
        ngrid <- 1000
        w <- withr::with_seed(seed = 42, code = {
            t(qrng::ghalton(n = ngrid, d = m - 1))
        })
    }

    vapply(final_regions, function(region) {
        if (any(is.nan(region))) {
            return(0)
        }
        if (method == "lpmvnorm") {
            exp(mvtnorm::lpmvnorm(lower = region[1, ],
                                  upper = region[2, ],
                                  mean = mean,
                                  chol = Ct,
                                  M = ngrid,
                                  w = w,
                                  ...))
        } else {
            mvtnorm::pmvnorm(lower = region[1, ],
                             upper = region[2, ],
                             mean = mean,
                             sigma = sigma,
                             seed = 42,
                             keepAttr = FALSE,
                             ...)
        }
    }, numeric(1))
}

sequential_undecided_regions <- function(breaks) {
    m <- length(look_n)
    suppressWarnings({
        critical_null <- sapply(seq_along(look_n), function(i) {
            tcrit(k = 1 / bf10_lower, n1 = look_n[i], n2 = look_n[i],
                  plocation = analysis_args$plocation,
                  pscale = analysis_args$pscale, pdf = analysis_args$pdf,
                  alternative = analysis_args$alternative,
                  type = analysis_args$type)
        })
        critical_effect <- sapply(seq_along(look_n), function(i) {
            tcrit(k = 1 / bf10_upper, n1 = look_n[i], n2 = look_n[i],
                  plocation = analysis_args$plocation,
                  pscale = analysis_args$pscale, pdf = analysis_args$pdf,
                  alternative = analysis_args$alternative,
                  type = analysis_args$type)
        })
    })

    finite <- !is.nan(critical_null) & !is.nan(critical_effect)
    if (any(finite) && all(critical_effect[finite] >= critical_null[finite])) {
        direction <- "positive"
    } else if (any(finite) &&
               all(critical_effect[finite] < critical_null[finite])) {
        direction <- "negative"
    } else {
        stop("Inconsistent critical values: direction cannot be inferred.")
    }

    lower_boundary <- bf10_to_posterior(bf10_lower)
    upper_boundary <- bf10_to_posterior(bf10_upper)
    final_regions <- vector("list", length(breaks) - 1)

    for (i in seq_along(final_regions)) {
        bin_lower <- max(breaks[i], lower_boundary)
        bin_upper <- min(breaks[i + 1], upper_boundary)
        region <- matrix(NaN, nrow = 2, ncol = m)

        if (bin_lower < bin_upper) {
            if (direction == "positive") {
                for (j in seq_len(m - 1)) {
                    region[, j] <- c(
                        if (is.nan(critical_null[j])) -Inf else critical_null[j],
                        critical_effect[j]
                    )
                }
                region[, m] <- c(
                    tcrit(k = 1 / posterior_to_bf10(bin_lower),
                          n1 = n_final, n2 = n_final,
                          plocation = analysis_args$plocation,
                          pscale = analysis_args$pscale,
                          pdf = analysis_args$pdf,
                          alternative = analysis_args$alternative,
                          type = analysis_args$type),
                    tcrit(k = 1 / posterior_to_bf10(bin_upper),
                          n1 = n_final, n2 = n_final,
                          plocation = analysis_args$plocation,
                          pscale = analysis_args$pscale,
                          pdf = analysis_args$pdf,
                          alternative = analysis_args$alternative,
                          type = analysis_args$type)
                )
            } else {
                for (j in seq_len(m - 1)) {
                    region[, j] <- c(
                        critical_effect[j],
                        if (is.nan(critical_null[j])) Inf else critical_null[j]
                    )
                }
                region[, m] <- c(
                    tcrit(k = 1 / posterior_to_bf10(bin_upper),
                          n1 = n_final, n2 = n_final,
                          plocation = analysis_args$plocation,
                          pscale = analysis_args$pscale,
                          pdf = analysis_args$pdf,
                          alternative = analysis_args$alternative,
                          type = analysis_args$type),
                    tcrit(k = 1 / posterior_to_bf10(bin_lower),
                          n1 = n_final, n2 = n_final,
                          plocation = analysis_args$plocation,
                          pscale = analysis_args$pscale,
                          pdf = analysis_args$pdf,
                          alternative = analysis_args$alternative,
                          type = analysis_args$type)
                )
            }
        }
        final_regions[[i]] <- region
    }

    final_regions
}

sequential_undecided_histogram <- function(delta, probabilities,
                                           bins = sequential_undecided_bins) {
    evidence <- "Undecided"
    probability <- probabilities$probability[
        as.character(probabilities$evidence) == evidence
    ]
    if (length(probability) == 0) probability <- 0

    breaks <- seq(posterior_plot_limits[1], posterior_plot_limits[2],
                  length.out = bins + 1)
    final_regions <- sequential_undecided_regions(breaks)
    design <- analysis_design_grid(delta)
    bin_mass <- numeric(length(final_regions))
    for (i in seq_len(nrow(design))) {
        neff <- look_n
        se <- 1 / sqrt(neff)
        pars <- predpars(se = se, dpm = design$d[i], dpsd = 0)
        bin_mass <- bin_mass + design$weight[i] * integrate_final_regions(
            final_regions = final_regions,
            mean = pars$mean,
            sigma = pars$sigma
        )
    }
    total_mass <- sum(bin_mass, na.rm = TRUE)
    bin_mass <- if (total_mass > 0 && probability > 0) {
        probability * bin_mass / total_mass
    } else {
        numeric(length(final_regions))
    }

    data.frame(
        ymin = breaks[-length(breaks)],
        ymax = breaks[-1],
        y = (breaks[-length(breaks)] + breaks[-1]) / 2,
        density = bin_mass / diff(breaks),
        evidence = factor(evidence, levels = evidence_levels)
    )
}

make_sequential_undecided_margin <- function(histogram_data, probabilities,
                                             delta, x_limit) {
    label_data <- probabilities[
        as.character(probabilities$evidence) == "Undecided",
    ]
    label_data$label <- probability_status_label(
        label_data$probability,
        label_data$evidence,
        delta
    )
    label_data$x <- x_limit * 0.95
    label_data$y <- zone_midpoints$y[
        as.character(zone_midpoints$evidence) == "Undecided"
    ]

    ggplot2::ggplot(histogram_data) +
        ggplot2::geom_rect(
            ggplot2::aes(xmin = 0, xmax = density, ymin = ymin, ymax = ymax),
            fill = evidence_cols[["Undecided"]],
            alpha = 0.55
        ) +
        ggplot2::geom_text(
            data = label_data,
            ggplot2::aes(x = x, y = y, label = label, colour = evidence),
            hjust = 1,
            vjust = 0.5,
            size = probability_label_size,
            lineheight = 0.92,
            fontface = "bold",
            show.legend = FALSE
        ) +
        ggplot2::scale_y_continuous(
            limits = posterior_plot_limits,
            breaks = bf10_to_posterior(bf10_breaks),
            labels = NULL
        ) +
        ggplot2::scale_x_continuous(limits = c(0, x_limit), expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = posterior_plot_limits, clip = "off") +
        ggplot2::scale_colour_manual(values = evidence_cols,
                                      limits = evidence_levels,
                                      drop = FALSE,
                                      guide = "none") +
        ggplot2::theme_void(base_size = margin_base_size) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(
                0, sequential_panel_margin, 0, 2
            )
        ) +
        transparent_background_theme()
}

make_stop_strip <- function(stop_curve, boundary = c("effect", "null"),
                            show_x_axis = FALSE, probabilities = NULL,
                            delta = NULL, mass_limit = NULL) {
    boundary <- match.arg(boundary)
    if (boundary == "effect") {
        mass <- stop_curve$effect_mass
        colour <- evidence_cols[["Evidence for effect"]]
        signed_mass <- mass
    } else {
        mass <- stop_curve$null_mass
        colour <- evidence_cols[["Evidence for null"]]
        signed_mass <- -mass
    }
    if (is.null(mass_limit)) {
        mass_limit <- max(c(mass, 0.01), na.rm = TRUE) * 1.2
    }
    y_limits <- if (boundary == "effect") {
        c(0, mass_limit)
    } else {
        c(-mass_limit, 0)
    }
    strip_data <- data.frame(n = stop_curve$n,
                             mass = mass,
                             signed_mass = signed_mass)

    p <- ggplot2::ggplot(strip_data, ggplot2::aes(x = n)) +
        ggplot2::geom_col(ggplot2::aes(y = signed_mass),
                          width = min(diff(look_n)) * 0.8,
                          fill = colour,
                          alpha = 0.55) +
        ggplot2::scale_x_continuous(breaks = sample_size_breaks,
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(x = if (show_x_axis) "Sample Size" else NULL,
                      y = NULL) +
        ggplot2::coord_cartesian(xlim = x_plot_limits, ylim = y_limits,
                                 clip = "off") +
        ggplot2::theme_void(base_size = margin_base_size) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(
                0, sequential_panel_margin, 0, sequential_panel_margin
            )
        ) +
        transparent_background_theme()

    if (!is.null(probabilities) && !is.null(delta)) {
        label_evidence <- if (boundary == "effect") {
            "Evidence for effect"
        } else {
            "Evidence for null"
        }
        label_data <- probabilities[
            as.character(probabilities$evidence) == label_evidence,
        ]
        label_data$label <- probability_status_label_inline(
            label_data$probability,
            label_data$evidence,
            delta
        )
        label_data$x <- n_final + x_padding * 0.85
        label_data$y <- if (boundary == "effect") {
            y_limits[2] * 0.62
        } else {
            y_limits[1] * 0.52
        }

        p <- p +
            ggplot2::geom_text(
                data = label_data,
                ggplot2::aes(x = x, y = y, label = label,
                             colour = evidence),
                hjust = 1,
                vjust = 0.5,
                size = strip_probability_label_size,
                lineheight = 0.92,
                fontface = "bold",
                show.legend = FALSE
            ) +
            ggplot2::scale_colour_manual(values = evidence_cols,
                                          limits = evidence_levels,
                                          drop = FALSE,
                                          guide = "none")
    }

    if (show_x_axis) {
        p <- p +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    colour = "#4B5563",
                    size = strip_axis_text_size
                ),
                axis.ticks.x = ggplot2::element_line(
                    colour = "#6B7280",
                    linewidth = 0.25
                ),
                axis.ticks.length.x = grid::unit(2 * plot_text_scale, "pt"),
                axis.title.x = ggplot2::element_text(
                    colour = "#374151",
                    size = strip_axis_title_size,
                    margin = ggplot2::margin(t = 4)
                ),
                plot.margin = ggplot2::margin(
                    0, sequential_panel_margin, 3, sequential_panel_margin
                )
            )
    }

    p
}

## Outcome probabilities and sampling helpers ---------------------------------

fixed_probs <- lapply(effects$delta, fixed_probabilities)
seq_curves <- lapply(effects$delta, sequential_stop_curve)

sequential_probabilities_from_curve <- function(stop_curve) {
    p_effect <- tail(stop_curve$effect_cumulative, 1)
    p_null <- tail(stop_curve$null_cumulative, 1)
    data.frame(
        evidence = factor(evidence_levels, levels = evidence_levels),
        probability = c(p_effect, 1 - p_effect - p_null, p_null)
    )
}

fixed_outcome <- function(path) {
    as.character(evidence_zone(path$bf10[path$n == n_final]))
}

sequential_outcome <- function(path) {
    stop_point <- find_stop(path)
    as.character(evidence_zone(stop_point$bf10))
}

sample_paths_by_outcome <- function(delta, effect, target_probabilities,
                                    outcome_fun, batch_size = 50,
                                    max_batches = 200) {
    target_counts <- probability_counts(target_probabilities$probability)
    names(target_counts) <- as.character(target_probabilities$evidence)
    selected <- setNames(vector("list", length(evidence_levels)), evidence_levels)
    selected_counts <- setNames(integer(length(evidence_levels)), evidence_levels)
    kept <- 0
    candidates <- 0

    for (batch in seq_len(max_batches)) {
        for (j in seq_len(batch_size)) {
            if (all(selected_counts >= target_counts)) {
                break
            }
            candidates <- candidates + 1
            path <- simulate_one_trajectory(delta = delta, effect = effect,
                                            id = candidates)
            outcome <- outcome_fun(path)
            if (selected_counts[[outcome]] < target_counts[[outcome]]) {
                kept <- kept + 1
                path$id <- kept
                selected_counts[[outcome]] <- selected_counts[[outcome]] + 1
                selected[[outcome]][[selected_counts[[outcome]]]] <- path
            }
        }
        if (all(selected_counts >= target_counts)) {
            break
        }
    }

    if (!all(selected_counts >= target_counts)) {
        stop("Could not sample enough trajectories for: ",
             paste(names(target_counts)[selected_counts < target_counts],
                   collapse = ", "))
    }

    do.call(rbind, unlist(selected, recursive = FALSE))
}

## Generate simulated paths ----------------------------------------------------

fixed_trajectory_data <- do.call(
    rbind,
    lapply(seq_len(nrow(effects)), function(i) {
        sample_paths_by_outcome(
            delta = effects$delta[i],
            effect = effects$effect[i],
            target_probabilities = fixed_probs[[i]],
            outcome_fun = fixed_outcome
        )
    })
)

seq_probabilities <- lapply(seq_curves, sequential_probabilities_from_curve)
seq_trajectory_data <- do.call(
    rbind,
    lapply(seq_len(nrow(effects)), function(i) {
        sample_paths_by_outcome(
            delta = effects$delta[i],
            effect = effects$effect[i],
            target_probabilities = seq_probabilities[[i]],
            outcome_fun = sequential_outcome
        )
    })
)

## Assemble plots and save outputs --------------------------------------------

fixed_points <- fixed_trajectory_data[fixed_trajectory_data$n == n_final, ]
fixed_points$evidence <- evidence_zone(fixed_points$bf10)
fixed_points$status <- status_label(fixed_points$evidence, fixed_points$delta)

seq_stop_points <- do.call(
    rbind,
    lapply(split(seq_trajectory_data, interaction(seq_trajectory_data$effect,
                                                  seq_trajectory_data$id,
                                                  drop = TRUE)), find_stop)
)
seq_stop_points$evidence <- evidence_zone(seq_stop_points$bf10)
seq_stop_points$status <- status_label(seq_stop_points$evidence,
                                       seq_stop_points$delta)

seq_paths <- merge(seq_trajectory_data,
                   seq_stop_points[, c("effect", "id", "n", "evidence")],
                   by = c("effect", "id"), suffixes = c("", "_stop"))
seq_paths <- seq_paths[seq_paths$n <= seq_paths$n_stop, ]

sequential_stop_mass_limit <- max(
    c(0.01, unlist(lapply(seq_curves, function(curve) {
        c(curve$effect_mass, curve$null_mass)
    }))),
    na.rm = TRUE
) * 1.2
seq_undecided_histograms <- lapply(seq_len(nrow(effects)), function(i) {
    sequential_undecided_histogram(
        delta = effects$delta[i],
        probabilities = seq_probabilities[[i]]
    )
})
seq_undecided_x_limit <- max(
    c(0.01, vapply(seq_undecided_histograms, function(histogram_data) {
        max(histogram_data$density, na.rm = TRUE)
    }, numeric(1))),
    na.rm = TRUE
) * 2.7

fixed_plots <- vector("list", nrow(effects) * 2)
sequential_plots <- vector("list", nrow(effects) * 4)
for (i in seq_len(nrow(effects))) {
    effect_i <- effects$effect[i]
    row_paths <- fixed_trajectory_data[fixed_trajectory_data$effect == effect_i, ]
    row_fixed_points <- fixed_points[fixed_points$effect == effect_i, ]
    row_seq_paths <- seq_paths[seq_paths$effect == effect_i, ]
    row_seq_points <- seq_stop_points[seq_stop_points$effect == effect_i, ]

    fixed_plots[[2 * (i - 1) + 1]] <- make_trajectory_plot(
        paths = row_paths,
        final_points = row_fixed_points,
        sequential = FALSE,
        show_y = i == 1
    )
    fixed_plots[[2 * (i - 1) + 2]] <-
        make_fixed_margin(effects$delta[i], fixed_probs[[i]])

    sequential_plots[[4 * (i - 1) + 1]] <-
        make_stop_strip(seq_curves[[i]], "effect",
                        probabilities = seq_probabilities[[i]],
                        delta = effects$delta[i],
                        mass_limit = sequential_stop_mass_limit)
    sequential_plots[[4 * (i - 1) + 2]] <- make_trajectory_plot(
        paths = row_seq_paths,
        final_points = row_seq_points,
        sequential = TRUE,
        show_y = i == 1,
        show_x = FALSE
    ) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(
                0, sequential_panel_margin, 0, sequential_panel_margin
            )
        )
    sequential_plots[[4 * (i - 1) + 3]] <-
        make_sequential_undecided_margin(
            histogram_data = seq_undecided_histograms[[i]],
            probabilities = seq_probabilities[[i]],
            delta = effects$delta[i],
            x_limit = seq_undecided_x_limit
        )
    sequential_plots[[4 * (i - 1) + 4]] <-
        make_stop_strip(seq_curves[[i]], "null", show_x_axis = TRUE,
                        probabilities = seq_probabilities[[i]],
                        delta = effects$delta[i],
                        mass_limit = sequential_stop_mass_limit)
}

names(fixed_plots) <- LETTERS[seq_along(fixed_plots)]
names(sequential_plots) <- LETTERS[seq_along(sequential_plots)]

fixed_figure_plots <- c(
    list(make_condition_title(1), make_condition_title(0)),
    fixed_plots
)
names(fixed_figure_plots) <- LETTERS[seq_along(fixed_figure_plots)]

sequential_figure_plots <- c(
    list(make_condition_title(1), make_condition_title(0)),
    sequential_plots
)
names(sequential_figure_plots) <- LETTERS[seq_along(sequential_figure_plots)]

fixed_figure <- patchwork::wrap_plots(
    fixed_figure_plots,
    design = "AABB\nCDEF",
    widths = c(1, 0.38, 1, 0.38),
    heights = c(0.13, 1)
) &
    transparent_background_theme()

sequential_figure <- patchwork::wrap_plots(
    sequential_figure_plots,
    design = "AABB\nC#G#\nDEHI\nF#J#",
    widths = c(1, 0.32, 1, 0.32),
    heights = c(0.13, 0.18, 1, 0.18)
) &
    transparent_background_theme()

fixed_pdf <- file.path(output_dir, "bf_design_analysis_fixed.pdf")
fixed_png <- file.path(output_dir, "bf_design_analysis_fixed.png")
sequential_pdf <- file.path(output_dir, "bf_design_analysis_sequential.pdf")
sequential_png <- file.path(output_dir, "bf_design_analysis_sequential.png")

export_design_analysis <- !tolower(
    Sys.getenv("BFPWR_SUPPRESS_DESIGN_ANALYSIS_EXPORT")
) %in% c("1", "true", "yes")

if (export_design_analysis) {
    ggplot2::ggsave(fixed_pdf, fixed_figure, width = 10, height = 3.6,
                    device = grDevices::cairo_pdf, bg = "transparent")
    ggplot2::ggsave(fixed_png, fixed_figure, width = 10, height = 3.6,
                    dpi = 300, bg = "transparent")
    ggplot2::ggsave(sequential_pdf, sequential_figure, width = 10, height = 4.6,
                    device = grDevices::cairo_pdf, bg = "transparent")
    ggplot2::ggsave(sequential_png, sequential_figure, width = 10, height = 4.6,
                    dpi = 300, bg = "transparent")

    message("Wrote ", fixed_pdf)
    message("Wrote ", fixed_png)
    message("Wrote ", sequential_pdf)
    message("Wrote ", sequential_png)
}
