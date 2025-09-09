## MRTAnalysis 0.3.0

- Added new functionality for mediated causal excursion effects in MRTs:
  -  Added `mcee()` function: streamlined workflow for estimating natural direct excursion effect (NDEE) and natural indirect excursion effect (NIEE) in micro-randomized trials (MRTs) with distal outcomes.
  -  Added two advanced wrappers:
    - `mcee_general()`: flexible configuration of nuisance models (p, q, eta, mu, nu) with support for multiple learners (glm, gam, lm, rf, ranger, sl).
  	-  `mcee_userfit_nuisance()`: allows users to inject externally fitted nuisance predictions.
  	-  Included config helper functions (`mcee_config_glm()`, `mcee_config_gam()`, `mcee_config_ranger()`, etc.) and `mcee_config_maker()` for building nuisance specifications to pass into `mcee_general()`.
  -  New dataset `data_time_varying_mediator_distal_outcome` included to illustrate usage.
  -  Added vignette "Time-Varying Causal Excursion Effect Mediation in MRT: Continuous Distal Outcomes" with detailed examples and best practices.


## MRTAnalysis 0.2.0

- Added new functionality for distal outcomes in MRTs:
  - Implemented `dcee()` for estimating distal causal excursion effects.
  - Supports flexible nuisance regression learners (`lm`, `gam`, `rf`, `ranger`, `SuperLearner`) with optional cross-fitting.
  - Provides small-sample t inference via `summary.dcee_fit()`, consistent with `wcls()` and `emee()`.
  - New synthetic dataset `data_distal_continuous` for examples and testing.
  - Added vignette: Exploratory Analysis for MRT: Distal Outcomes.
- Minor bug fixes and improvements to wcls() and emee() documentation.

## MRTAnalysis 0.1.2

- Fixed a bug in wcls when the randomization probability is
    time-varying.
- Now all variable inputs need to be in quotation marks; for example,
    from now on one should specify id = "userid" instead of id = userid.
    This is to allow dynamically specified column names.

## MRTAnalysis 0.1.1

- Updated vignette to improve clarify.

## MRTAnalysis 0.1.0

- Initial release.
