## MRTAnalysis 0.2.0

- Added new functionality for distal outcomes in MRTs:
	- Implemented `dcee()` for estimating distal causal excursion effects.
	- Supports flexible nuisance regression learners (`lm`, `gam`, `rf`, `ranger`, `SuperLearner`) with optional cross-fitting.
	- Provides small-sample t inference via `summary.dcee_fit()`, consistent with `wcls()` and `emee()`.
	- New synthetic dataset `data_distal_continuous` for examples and testing.
	- Added vignette: Exploratory Analysis for MRT: Distal Outcomes.
- Minor bug fixes and improvements to wcls() and emee() documentation.

## MRTAnalysis 0.1.2

-   Fixed a bug in wcls when the randomization probability is
    time-varying.
-   Now all variable inputs need to be in quotation marks; for example,
    from now on one should specify id = "userid" instead of id = userid.
    This is to allow dynamically specified column names.

## MRTAnalysis 0.1.1

-   Updated vignette to improve clarify.

## MRTAnalysis 0.1.0

-   Initial release.
