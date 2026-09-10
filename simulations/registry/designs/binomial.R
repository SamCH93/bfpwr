bfpwr_sim_design_cases_binomial <- function() {
    list(
        bfpwr_sim_design_binomial(
            id = "binom-dpoint-0p5-short",
            tier = "short",
            look_grid = "short",
            seed = 20263001,
            design_prior = list(family = "point", prob = 0.5),
            tags = c("point-null", "center-probability"),
            rationale = "Point-null binomial design at p = 0.5."
        ),
        bfpwr_sim_design_binomial(
            id = "binom-dpoint-0p6-short",
            tier = "short",
            look_grid = "short",
            seed = 20263002,
            design_prior = list(family = "point", prob = 0.6),
            tags = c("ordinary-alternative", "center-probability"),
            rationale = "Ordinary positive binomial alternative for analyses around p0 = 0.5."
        ),
        bfpwr_sim_design_binomial(
            id = "binom-dpoint-0p2-short",
            tier = "short",
            look_grid = "short",
            seed = 20263003,
            design_prior = list(family = "point", prob = 0.2),
            tags = c("low-probability", "point-design"),
            rationale = "Low-proportion binomial design for p0 = 0.2 examples and directional analyses."
        ),
        bfpwr_sim_design_binomial(
            id = "binom-dpoint-0p75-short",
            tier = "short",
            look_grid = "short",
            seed = 20263004,
            design_prior = list(family = "point", prob = 0.75),
            tags = c("high-probability", "point-design"),
            rationale = "High-proportion binomial design for p0 = 0.75 examples."
        ),
        bfpwr_sim_design_binomial(
            id = "binom-dbeta-60-40-short",
            tier = "short",
            look_grid = "short",
            seed = 20263005,
            design_prior = list(family = "beta", shape1 = 60, shape2 = 40,
                                lower = 0, upper = 1),
            tags = c("beta-design", "concentrated-design", "center-probability"),
            rationale = "Concentrated beta design prior centered near p = 0.6."
        ),
        bfpwr_sim_design_binomial(
            id = "binom-dbeta-1-1-h1gt0p5-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20263006,
            design_prior = list(family = "beta", shape1 = 1, shape2 = 1,
                                lower = 0.5, upper = 1),
            tags = c("truncated-beta-design", "directional-h1"),
            rationale = "Uniform truncated-beta design prior on the H1 side of p0 = 0.5."
        ),
        bfpwr_sim_design_binomial(
            id = "binom-dbeta-1-1-h0le0p5-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20263007,
            design_prior = list(family = "beta", shape1 = 1, shape2 = 1,
                                lower = 0, upper = 0.5),
            tags = c("truncated-beta-design", "directional-h0"),
            rationale = "Uniform truncated-beta design prior on the H0 side of p0 = 0.5."
        ),
        bfpwr_sim_design_binomial(
            id = "binom-dpoint-0p5-long",
            tier = "long",
            look_grid = "long",
            seed = 20263008,
            design_prior = list(family = "point", prob = 0.5),
            tags = c("point-null", "center-probability", "many-looks"),
            rationale = "Long-grid point-null binomial design at p = 0.5."
        ),
        bfpwr_sim_design_binomial(
            id = "binom-dpoint-0p55-long",
            tier = "long",
            look_grid = "long",
            seed = 20263009,
            design_prior = list(family = "point", prob = 0.55),
            tags = c("local-effect", "center-probability", "many-looks"),
            rationale = "Long-grid local binomial alternative for stepwise power and sample-size behavior."
        )
    )
}
