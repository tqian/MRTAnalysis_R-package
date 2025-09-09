# MRTAnalysis

The **MRTAnalysis** package provides functions to conduct post-study analyses of Micro-Randomized Trials (MRTs), focusing on estimating **causal excursion effects**.

- **Proximal outcomes** (measured shortly after each treatment decision point):
  - `wcls()`: Primary analysis for continuous proximal outcomes. Implements weighted and centered least squares (the $k=1$ special case of Boruvka et al., 2018).
  - `emee()`: Primary analysis for binary proximal outcomes. Implements the estimator for marginal excursion effect (the $\Delta=1$ special case of Qian et al., 2021).
  - `emee2()`: Variant of `emee()`, centering treatment in the residual term. Basis for the sample size calculator in `MRTSampleSizeBinary`.

- **Distal outcomes** (measured once at end of study):
  - `dcee()`: Exploratory analysis for distal causal excursion effects in MRTs (Qian et al. 2025). Supports linear models and machine-learning learners (lm, gam, random forest, ranger, SuperLearner) with optional cross-fitting.

- **Mediated Effects Through Time-Varying Mediators to Distal Outcomes**
  -   - `mcee()`: Exploratory analysis for mediated causal excursion effects in MRTs, estimating **natural direct excursion effects (NDEE)** and **natural indirect excursion effects (NIEE)** through time-varying mediators. Supports GLM, GAM, random forest, ranger, and SuperLearner learners for fitting nuisance parameters.


## Installation

You can install the package from CRAN:

```r
install.packages("MRTAnalysis")
```

## Usage

See vignettes for detailed examples:

```r
library(MRTAnalysis)

# Proximal outcome analysis (continuous)
fit1 <- wcls(
  data = data_mimicHeartSteps,
  id = "userid", outcome = "logstep_30min",
  treatment = "intervention", rand_prob = 0.6,
  moderator_formula = ~1,
  control_formula = ~logstep_pre30min,
  availability = "avail"
)
summary(fit1)

# Distal outcome analysis
fit2 <- dcee(
  data = data_distal_continuous,
  id = "userid", outcome = "Y",
  treatment = "A", rand_prob = "prob_A",
  moderator_formula = ~1,
  control_formula = ~X,
  availability = "avail",
  control_reg_method = "lm"
)
summary(fit2)

# Mediation with distal outcome
fit3 <- mcee(
  data = data_time_varying_mediator_distal_outcome,
  id = "id", dp = "dp",
  outcome = "Y", treatment = "A", mediator = "M",
  availability = "I", rand_prob = "p_A",
  time_varying_effect_form = ~1,                # constant effects over time
  control_formula_with_mediator = ~ dp + M + X, # adjustment set
  control_reg_method = "glm"
)
summary(fit3)
```

## References

- Boruvka, A., Almirall, D., Witkiewitz, K., & Murphy, S. A. (2018). Assessing time-varying causal effect moderation in mobile health. *Journal of the American Statistical Association*, 113(523), 1112–1121. <doi:10.1080/01621459.2017.1305274>

- Qian, T., Yoo, H., Klasnja, P., Almirall, D., & Murphy, S. A. (2021). Estimating time-varying causal excursion effects in mobile health with binary outcomes. *Biometrika*, 108(3), 507–527. <doi:10.1093/biomet/asaa070>

- Qian, T. (2025). Distal Causal Excursion Effects: Modeling Long-Term Effects of Time-Varying Treatments in Micro-Randomized Trials. arXiv:2502.13500.

- Qian, T. (2025). Dynamic Causal Mediation Analysis for Intensive Longitudinal Data. arXiv:2506.20027.