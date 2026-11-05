## Combined Bayes factor design-analysis figure for the JASP tutorial.
##
## This script reuses the fixed and sequential plot objects from
## bf_design_analysis_figure.R, stacks them, and adds row labels.

script_directory <- function() {
    script_file <- commandArgs(trailingOnly = FALSE)
    script_file <- script_file[grepl("^--file=", script_file)]
    if (length(script_file) == 1) {
        return(dirname(normalizePath(sub("^--file=", "", script_file),
                                     winslash = "/")))
    }
    normalizePath(getwd(), winslash = "/")
}

script_dir <- script_directory()
source_file <- file.path(script_dir, "bf_design_analysis_figure.R")
output_dir <- script_dir

previous_suppress <- Sys.getenv("BFPWR_SUPPRESS_DESIGN_ANALYSIS_EXPORT",
                                unset = NA)
Sys.setenv(BFPWR_SUPPRESS_DESIGN_ANALYSIS_EXPORT = "true")
source(source_file)
if (is.na(previous_suppress)) {
    Sys.unsetenv("BFPWR_SUPPRESS_DESIGN_ANALYSIS_EXPORT")
} else {
    Sys.setenv(BFPWR_SUPPRESS_DESIGN_ANALYSIS_EXPORT = previous_suppress)
}

make_row_label <- function(label) {
    ggplot2::ggplot() +
        ggplot2::annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = label,
            angle = 90,
            hjust = 0.5,
            vjust = 0.5,
            size = condition_title_size / 2.845276,
            fontface = "bold",
            colour = "#111827"
        ) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1),
                                 clip = "off") +
        ggplot2::theme_void(base_size = condition_title_size) +
        transparent_background_theme()
}

fixed_height <- 3.6
sequential_height <- 4.6
sequential_titled_row_heights <- c(0.13, 0.18, 1, 0.18)
sequential_untitled_row_heights <- c(0.18, 1, 0.18)
sequential_untitled_height <- sequential_height *
    sum(sequential_untitled_row_heights) / sum(sequential_titled_row_heights)

content_width <- 10
row_label_width <- 0.6
combined_width <- content_width + row_label_width
combined_height <- fixed_height + sequential_untitled_height

sequential_figure_without_titles <- patchwork::wrap_plots(
    sequential_plots,
    design = "A#E#\nBCFG\nD#H#",
    widths = c(1, 0.32, 1, 0.32),
    heights = sequential_untitled_row_heights
) &
    transparent_background_theme()

combined_figure <- patchwork::wrap_plots(
    list(
        make_row_label("Fixed Design"),
        fixed_figure,
        make_row_label("Sequential Design"),
        sequential_figure_without_titles
    ),
    design = "AB\nCD",
    widths = c(row_label_width, content_width),
    heights = c(fixed_height, sequential_untitled_height)
) &
    transparent_background_theme()

combined_pdf <- file.path(output_dir, "bf_design_analysis_combined.pdf")
combined_png <- file.path(output_dir, "bf_design_analysis_combined.png")

ggplot2::ggsave(combined_pdf, combined_figure,
                width = combined_width, height = combined_height,
                device = grDevices::cairo_pdf, bg = "transparent")
ggplot2::ggsave(combined_png, combined_figure,
                width = combined_width, height = combined_height,
                dpi = 300, bg = "transparent")

message("Wrote ", combined_pdf)
message("Wrote ", combined_png)
