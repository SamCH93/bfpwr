bfpwr_sim_design_cases_z <- function() {
    list(
        bfpwr_sim_design_z(
            id = "z-dpoint-0-usd1-short",
            tier = "short",
            look_grid = "short",
            seed = 20261001,
            dpm = 0,
            dpsd = 0,
            usd = 1,
            tags = c("matched-usd1", "null"),
            rationale = "Point-null z design for the usd = 1 sampling model."
        ),
        bfpwr_sim_design_z(
            id = "z-dpoint-0p2-usd1-short",
            tier = "short",
            look_grid = "short",
            seed = 20261002,
            dpm = 0.2,
            dpsd = 0,
            usd = 1,
            tags = c("matched-usd1", "local-effect"),
            rationale = "Local positive-effect z design for the usd = 1 sampling model."
        ),
        bfpwr_sim_design_z(
            id = "z-dpoint-0p5-usd1-short",
            tier = "short",
            look_grid = "short",
            seed = 20261003,
            dpm = 0.5,
            dpsd = 0,
            usd = 1,
            tags = c("matched-usd1", "moderate-effect"),
            rationale = "Moderate positive-effect z design for the usd = 1 sampling model."
        ),
        bfpwr_sim_design_z(
            id = "z-dnorm-0p5-s0p1-usd1-short",
            tier = "short",
            look_grid = "short",
            seed = 20261004,
            dpm = 0.5,
            dpsd = 0.1,
            usd = 1,
            tags = c("matched-usd1", "uncertain-design", "moderate-effect"),
            rationale = "Aligned uncertain positive-effect z design with small design-prior spread."
        ),
        bfpwr_sim_design_z(
            id = "z-dnorm-0p2-s0p3-usd1-short",
            tier = "short",
            look_grid = "short",
            seed = 20261005,
            dpm = 0.2,
            dpsd = 0.3,
            usd = 1,
            tags = c("matched-usd1", "uncertain-design", "local-effect"),
            rationale = "Local uncertain z design with design-prior spread larger than its mean."
        ),
        bfpwr_sim_design_z(
            id = "z-dnorm-0p2-s0p5-usd1-short",
            tier = "short",
            look_grid = "short",
            seed = 20261006,
            dpm = 0.2,
            dpsd = 0.5,
            usd = 1,
            tags = c("matched-usd1", "uncertain-design", "diffuse-design"),
            rationale = "Diffuse local z design for checking behavior as design uncertainty increases."
        ),
        bfpwr_sim_design_z(
            id = "z-dnorm-0-s0p5-usd1-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20261007,
            dpm = 0,
            dpsd = 0.5,
            usd = 1,
            tags = c("matched-usd1", "mixed-signs", "uncertain-design"),
            rationale = "Mixed-sign z design prior centered on the null."
        ),
        bfpwr_sim_design_z(
            id = "z-dpoint-m0p3-usd1-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20261008,
            dpm = -0.3,
            dpsd = 0,
            usd = 1,
            tags = c("matched-usd1", "wrong-direction"),
            rationale = "Wrong-direction z design for directional analysis cases."
        ),
        bfpwr_sim_design_z(
            id = "z-dpoint-0-usd50-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20261009,
            dpm = 0,
            dpsd = 0,
            usd = 50,
            tags = c("extreme-noise", "null"),
            rationale = "Extreme-noise point-null z design for numerical stress testing."
        ),
        bfpwr_sim_design_z(
            id = "z-dpoint-10-usd50-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20261010,
            dpm = 10,
            dpsd = 0,
            usd = 50,
            tags = c("extreme-noise", "large-effect"),
            rationale = "Extreme-noise z design with a large fixed effect visible only with enough information."
        ),
        bfpwr_sim_design_z(
            id = "z-dnorm-10-s20-usd50-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20261011,
            dpm = 10,
            dpsd = 20,
            usd = 50,
            tags = c("extreme-noise", "uncertain-design", "diffuse-design"),
            rationale = "Extreme-noise z design with a very diffuse large-effect design prior."
        ),
        bfpwr_sim_design_z(
            id = "z-dpoint-0-usd1-long",
            tier = "long",
            look_grid = "long",
            seed = 20261012,
            dpm = 0,
            dpsd = 0,
            usd = 1,
            tags = c("matched-usd1", "null", "many-looks"),
            rationale = "Long-grid point-null z design for many-look sequential behavior."
        ),
        bfpwr_sim_design_z(
            id = "z-dpoint-0p2-usd1-long",
            tier = "long",
            look_grid = "long",
            seed = 20261013,
            dpm = 0.2,
            dpsd = 0,
            usd = 1,
            tags = c("matched-usd1", "local-effect", "many-looks"),
            rationale = "Long-grid local-effect z design for late stopping and sample-size checks."
        ),
        bfpwr_sim_design_z(
            id = "z-dnorm-0p2-s0p1-usd1-long",
            tier = "long",
            look_grid = "long",
            seed = 20261014,
            dpm = 0.2,
            dpsd = 0.1,
            usd = 1,
            tags = c("matched-usd1", "uncertain-design", "many-looks"),
            rationale = "Long-grid local uncertain z design with small design-prior spread."
        ),
        bfpwr_sim_design_z(
            id = "z-dnorm-0p2-s0p5-usd1-long",
            tier = "long",
            look_grid = "long",
            seed = 20261015,
            dpm = 0.2,
            dpsd = 0.5,
            usd = 1,
            tags = c("matched-usd1", "uncertain-design", "diffuse-design", "many-looks"),
            rationale = "Long-grid diffuse local z design for many-look behavior under design uncertainty."
        )
    )
}
