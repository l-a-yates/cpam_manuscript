## Changepoint models for time series omics data with [`cpam`](https://l-a-yates.github.io/cpam/)
This repository contains the code to reproduce the simulations and results in the 
manuscript "Shape-constrained, changepoint additive models for time series omics data" 
by Yates et. al (2024).

For details about the [`cpam`](https://l-a-yates.github.io/cpam/) package, including installation, 
see the package [repository](https://github.com/l-a-yates/cpam) or the [website](https://github.com/l-a-yates/cpam).

### Package use and reproducing the case studies
An introduction to the package and code to reproduce the case studies in the 
manuscript can be found in the following tutorials
 - [Case Study 1](https://raw.githack.com/l-a-yates/cpam_manuscript/main/R/crisp.html)
 - [Case Study 2](https://raw.githack.com/l-a-yates/cpam_manuscript/main/R/torre.html)

There is also an introductory example using a small simulated data set to get started with the package
 - [Introductory Example](https://raw.githack.com/l-a-yates/cpam_manuscript/main/R/example.html)

### Reproducing the simulations
Code to reproduce the simulations in the manuscript can be found in the `R` directory.

 - The files `sim_calibration_co.R` and `sim_calibration_cc.R` contain the code to 
 reproduce the $p$-value calibration results for the case-only and case-control 
 models, respectively.
 - The files `sim_trends_co.R` and `sim_trends_cc.R` contain the code to reproduce
 
