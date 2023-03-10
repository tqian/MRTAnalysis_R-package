# MRTanalaysis

This package provides functions to conduct primary analysis (estimating marginal causal excursion effect and moderated causal excursion effect) for Micro-randomized Trials with continuous or binary outcomes. 

Please read the vignettes for details on how to use this package.

The main functions are the following:

- wcls(): Primary analysis for continuous outcome. This is an implementation of the weighted and centered least squares method (the k=1 special case of Boruvka, A., Almirall, D., Witkiewitz, K., & Murphy, S. A. (2018). Assessing time-varying causal effect moderation in mobile health. Journal of the American Statistical Association, 113(523), 1112-1121.)
- emee(): Primary analysis for binary outcome. This is an implementation of the estimator for marginal excursion effect method (the Delta=1 special case of Qian, T., Yoo, H., Klasnja, P., Almirall, D., & Murphy, S. A. (2021). Estimating time-varying causal excursion effects in mobile health with binary outcomes. Biometrika, 108(3), 507-527.)


Additional functions that can be useful:

- emee2(): Primary analysis for binary outcome. This is a slightly altered version of emee(), where the treatment assignment indicator is also centered in the residual term. It would have similar (but not exactly the same) numerical output as emee(). This is the estimator based on which the sample size calculator for binary outcome MRT is developed. (See R package MRTSampleSizeBinary.)