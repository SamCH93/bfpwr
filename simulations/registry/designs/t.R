bfpwr_sim_design_cases_t <- function() {
    list(
        bfpwr_sim_design_t(
            id = "t-two-dpoint-0-short",
            tier = "short",
            look_grid = "short",
            seed = 20262001,
            dpm = 0,
            dpsd = 0,
            type = "two.sample",
            tags = c("balanced", "null"),
            rationale = "Balanced two-sample point-null t design."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dpoint-0p5-short",
            tier = "short",
            look_grid = "short",
            seed = 20262002,
            dpm = 0.5,
            dpsd = 0,
            type = "two.sample",
            tags = c("balanced", "moderate-effect"),
            rationale = "Balanced two-sample moderate-effect t design."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dnorm-0p5-s0p1-short",
            tier = "short",
            look_grid = "short",
            seed = 20262003,
            dpm = 0.5,
            dpsd = 0.1,
            type = "two.sample",
            tags = c("balanced", "uncertain-design", "moderate-effect"),
            rationale = "Balanced two-sample uncertain t design with small design-prior spread."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dnorm-0p5-s0p5-short",
            tier = "short",
            look_grid = "short",
            seed = 20262004,
            dpm = 0.5,
            dpsd = 0.5,
            type = "two.sample",
            tags = c("balanced", "uncertain-design", "diffuse-design"),
            rationale = "Balanced two-sample uncertain t design with larger design-prior spread."
        ),
        bfpwr_sim_design_t(
            id = "t-one-dpoint-0-short",
            tier = "short",
            look_grid = "short",
            seed = 20262005,
            dpm = 0,
            dpsd = 0,
            type = "one.sample",
            tags = c("one-sample", "null"),
            rationale = "One-sample point-null t design matching the one-sample alternative geometry."
        ),
        bfpwr_sim_design_t(
            id = "t-one-dpoint-0p5-short",
            tier = "short",
            look_grid = "short",
            seed = 20262006,
            dpm = 0.5,
            dpsd = 0,
            type = "one.sample",
            tags = c("one-sample", "moderate-effect"),
            rationale = "One-sample moderate-effect t design."
        ),
        bfpwr_sim_design_t(
            id = "t-paired-dpoint-0-short",
            tier = "short",
            look_grid = "short",
            seed = 20262007,
            dpm = 0,
            dpsd = 0,
            type = "paired",
            tags = c("paired", "null"),
            rationale = "Paired point-null t design matching the paired alternative geometry."
        ),
        bfpwr_sim_design_t(
            id = "t-paired-dpoint-0p5-short",
            tier = "short",
            look_grid = "short",
            seed = 20262008,
            dpm = 0.5,
            dpsd = 0,
            type = "paired",
            tags = c("paired", "moderate-effect"),
            rationale = "Paired moderate-effect t design."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dpoint-m0p4-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20262009,
            dpm = -0.4,
            dpsd = 0,
            type = "two.sample",
            tags = c("balanced", "wrong-direction"),
            rationale = "Balanced two-sample wrong-direction t design for one-sided analysis cases."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dpoint-0-n2x1p1-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20262010,
            dpm = 0,
            dpsd = 0,
            type = "two.sample",
            n2_multiplier = 1.1,
            tags = c("unequal-groups", "null"),
            rationale = "Unequal two-sample point-null t design matching the unequal alternative geometry."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dpoint-0p35-n2x1p1-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20262011,
            dpm = 0.35,
            dpsd = 0,
            type = "two.sample",
            n2_multiplier = 1.1,
            tags = c("unequal-groups", "local-effect"),
            rationale = "Unequal two-sample positive-effect t design for checking n1/n2 handling."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dnorm-0-s0p5-short",
            tier = "adversarial",
            look_grid = "short",
            seed = 20262012,
            dpm = 0,
            dpsd = 0.5,
            type = "two.sample",
            tags = c("balanced", "mixed-signs", "uncertain-design"),
            rationale = "Balanced two-sample mixed-sign t design prior centered on the null."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dpoint-0-long",
            tier = "long",
            look_grid = "long",
            seed = 20262013,
            dpm = 0,
            dpsd = 0,
            type = "two.sample",
            tags = c("balanced", "null", "many-looks"),
            rationale = "Long-grid balanced two-sample point-null t design."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dpoint-0p2-long",
            tier = "long",
            look_grid = "long",
            seed = 20262014,
            dpm = 0.2,
            dpsd = 0,
            type = "two.sample",
            tags = c("balanced", "local-effect", "many-looks"),
            rationale = "Long-grid balanced two-sample local-effect t design."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dnorm-0p2-s0p1-long",
            tier = "long",
            look_grid = "long",
            seed = 20262015,
            dpm = 0.2,
            dpsd = 0.1,
            type = "two.sample",
            tags = c("balanced", "uncertain-design", "many-looks"),
            rationale = "Long-grid balanced two-sample local uncertain t design with small spread."
        ),
        bfpwr_sim_design_t(
            id = "t-two-dnorm-0p2-s0p5-long",
            tier = "long",
            look_grid = "long",
            seed = 20262016,
            dpm = 0.2,
            dpsd = 0.5,
            type = "two.sample",
            tags = c("balanced", "uncertain-design", "diffuse-design", "many-looks"),
            rationale = "Long-grid balanced two-sample diffuse local t design."
        )
    )
}
