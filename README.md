# FLicc

**FLicc: FLR tools for length-interval catch-curve models**

`FLicc` is an FLR-friendly TMB implementation of the Bayesian length-interval catch-curve model (*fishblicc*) by Medley (2025).

The package retains the biological core of the original equilibrium, length-based, multi-gear framework while providing a fast penalized-likelihood / MPD workflow in Template Model Builder (TMB) within FLR.

---

## Motivation

`FLicc` provides a practical and extensible implementation of a length-interval catch-curve model that is:

- fast to estimate,
- easy to diagnose,
- compatible with FLR workflows.

It is designed for replication, sensitivity analysis, and integration into advisory contexts.

---

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("YOUR-USER/FLicc")
```

---

## Multi-gear focus

`FLicc` supports multiple gears in a single equilibrium fit:

- joint fitting of multiple length-frequency datasets,
- estimation of gear-specific fishing mortality,
- estimation of gear-specific selectivity,
- use of relative catch information to inform estimation.

---

## Main assumptions

1. Approximate steady state  
2. von Bertalanffy growth  
3. Variability in asymptotic length  
4. Parametric selectivity  
5. Use of length-frequency and relative catch data  

---

## Typical workflow

A standard analysis consists of:

1. Preparing length-frequency data and life-history parameters  
2. Constructing the equilibrium stock representation  
3. Inspecting observed length distributions  
4. Fitting the model  
5. Diagnosing the fit  
6. Extracting outputs such as SPR and selectivity  

---

## Working with multiple years (FLR-style workflows)

Although `FLicc` is an equilibrium model, it integrates naturally into FLR workflows.

A common approach is to loop over years or iterations:

- Data structured as length × year × iteration  
- Fit the model separately for each year  
- Store outputs such as SPR, selectivity, and fishing mortality  

This enables:

- time-series of SPR  
- sensitivity analyses  
- embedding in simulation or MSE frameworks  

Note: each fit assumes equilibrium, so interpretation should focus on relative comparisons.

---

## Current implementation

- single-period equilibrium fits  
- multiple gears  
- population-at-length recursion  
- logistic and double-sided normal selectivity  
- SPR calculations  
- TMB-based estimation  

---

## Planned development

- expanded selectivity options  
- improved diagnostics  
- Beverton–Holt extensions  
- tighter FLR integration  

---

## References

Kell et al. (2007) FLR framework  
Kristensen et al. (2016) TMB  
Medley (2025) fishblicc model  
