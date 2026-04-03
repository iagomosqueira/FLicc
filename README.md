# FLicc
FLicc: FLR tools for length-interval catch-curve models
# FLicc

`FLicc` is an FLR-friendly TMB implementation of the **Bayesian length-interval catch-curve model** (**fishblicc**) by Medley (2025). 

The package aims to retain the biological core of the original equilibrium, length-based, **multi-gear** framework while providing a fast penalized-likelihood / MPD workflow in **Template Model Builder** (TMB; Kristensen et al. 2016) within **FLR** (Kell et al. 2007).
The original fishblicc framework is described by Medley (2025)
## Motivation

The purpose of `FLicc` is to provide a practical and extensible TMB implementation of a length-interval catch-curve model that is easy to fit, diagnose, and embed in FLR-style workflows. The package is intended as a bridge between the original fishblicc model and a faster estimation framework suitable for replication, sensitivity analysis, method development, and eventual FLR integration. :contentReference[oaicite:1]{index=1}

## Multi-gear focus

A central feature of `FLicc` is support for **multiple gears in a single equilibrium fit**. Gear-specific selectivity and fishing mortality can be estimated jointly while sharing the same underlying population and growth assumptions. This follows the motivation of the original fishblicc framework, which was developed specifically to analyse **length data from multi-gear fisheries**. 

In practical terms, this means `FLicc` is designed to:
- fit multiple observed length-frequency samples simultaneously,
- estimate gear-specific fishing mortality parameters,
- estimate gear-specific selectivity curves,
- and use relative catch information to help identify fishing mortality across gears. 

## Main assumptions

The current implementation follows the main equilibrium assumptions of the original fishblicc formulation:

1. **Approximate steady state**  
   The population is assumed to have been in approximate equilibrium for at least around a generation near the time of sampling, so the model is not a time-series model.
   
3. **von Bertalanffy mean growth**  
   Mean growth is assumed to follow the von Bertalanffy growth curve.

4. **Variability in asymptotic length**  
   The original fishblicc formulation represents individual variability in asymptotic size, and the TMB implementation follows this core idea through the fishblicc-style population-at-length recursion. 

5. **Parametric selectivity**  
   Selectivity is represented through parametric functions. The original implementation highlights logistic, normal, and double-sided normal forms, and the current `FLicc` TMB implementation focuses first on the forms needed for replication work, especially logistic and double-sided normal. 

6. **Length-frequency and relative catch information**  
   The model uses length-frequency data, together with relative catch information among gears, to infer gear-specific fishing mortality and stock status in equilibrium. 

## Current implementation

The current version of `FLicc` focuses on the core equilibrium fitting problem:

- single-period equilibrium fit,
- multiple gears,
- fishblicc-style population-at-length recursion,
- gear-specific fishing mortality,
- logistic and double-sided normal selectivity,
- SPR and related post hoc calculations,
- TMB-based penalized-likelihood / MPD estimation. :contentReference[oaicite:9]{index=9}

## Current simplifications

To keep the first TMB implementation transparent and robust, several simplifications have been made relative to the original fishblicc package:

1. **No full Bayesian MCMC workflow**  
   The current package is fitted in TMB as a penalized-likelihood / MPD model rather than through the full Bayesian `rstan` workflow used in the original fishblicc implementation.

2. **Single-period scope**  
   The present implementation targets the equilibrium single-period case first.

3. **Selective parameter fixing in examples**  
   In some replication examples, biological parameters such as `Linf` and `Mk` may be held fixed to isolate structural differences and stabilize early development.

4. **Reduced selectivity scope**  
   The initial TMB version focuses on the selectivity forms most important for early replication work.

5. **Post hoc extensions**  
   Some management-oriented calculations, such as Beverton–Holt steepness-based translations of SPR into depletion metrics, are handled post hoc rather than as part of the fitted objective function.

These simplifications are deliberate. The priority is to recover the main biological behavior of the original multi-gear length-based model in TMB, then progressively expand toward a fuller FLR-compatible toolkit. :contentReference[oaicite:11]{index=11}

## Planned development

Planned next steps include:
- closer alignment of the observation model with the original fishblicc implementation,
- expanded selectivity options,
- richer plotting and diagnostics,
- post hoc Beverton–Holt calculations using user-supplied steepness,
- and tighter FLR integration for data handling and output classes. :contentReference[oaicite:12]{index=12}

## References

Kell, L. T., Mosqueira, I., Grosjean, P., Fromentin, J.-M., Garcia, D., Hillary, R., Jardim, E., Mardle, S., Pastoors, M. A., Poos, J. J., Scott, F., and Scott, R. D. 2007. *FLR: an open-source framework for the evaluation and development of management strategies*. ICES Journal of Marine Science, 64(4): 640–646. :contentReference[oaicite:13]{index=13}

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., and Bell, B. M. 2016. *TMB: Automatic Differentiation and Laplace Approximation*. Journal of Statistical Software, 70(5): 1–21. :contentReference[oaicite:14]{index=14}

Medley, P. A. H. 2025. *A new Bayesian catch curve stock assessment model for the analysis of length data from multi-gear fisheries*. ICES Journal of Marine Science, 82(12): fsaf224. :contentReference[oaicite:15]{index=15}
