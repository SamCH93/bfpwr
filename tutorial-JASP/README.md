# JASP Tutorial Figures

This folder contains the script and rendered figures for the Bayes factor
design-analysis illustration used in the JASP tutorial.

Run from the repository root with:

```sh
Rscript tutorial-JASP/bf_design_analysis_figure.R
```

To generate the stacked combined figure, run:

```sh
Rscript tutorial-JASP/bf_design_analysis_combined_figure.R
```

The script writes these files into this folder:

- `bf_design_analysis_fixed.pdf`
- `bf_design_analysis_fixed.png`
- `bf_design_analysis_sequential.pdf`
- `bf_design_analysis_sequential.png`
- `bf_design_analysis_combined.pdf`
- `bf_design_analysis_combined.png`

The calculations use the local `package/R` sources when the script is run from
this repository. If those sources are not available, the script falls back to
the installed `bfpwr` package.
