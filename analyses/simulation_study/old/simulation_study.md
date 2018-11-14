---
title: 'Simulation study'
output:
  html_document:
    toc: yes
    theme: united
    highlight: tango
    keep_md: yes
---



## Objectives
* how does DIABLO compares with existing integrative classification schemes (Concatenation and Ensemble)
* effect of the design (connection between datasets) on the types of variables selected (correlated vs. discriminatory) and error rate

\pagebreak

## Generate different types of variables and apply diablo to each type separately
* 3 datasets (effective sample size = 100; group1=100 observations, group2=100 observations)
* each dataset has four types of variables; lets explore them now
  + 30 variables that contribute to correlated & discriminatory components
  + 30 variables that contribute to correlated & nondiscriminatory components
  + 100 variables that contribute to uncorrelated & discriminatory components
  + 100 variables that contribute to uncorrelated & nondiscriminatory components

### Correlation structure for each set of simulated variables

#### corDis

<img src="simulation_study_files/figure-html/unnamed-chunk-1-1.png" width="100%" />

#### corNonDis

<img src="simulation_study_files/figure-html/unnamed-chunk-2-1.png" width="100%" />

#### unCorDis

<img src="simulation_study_files/figure-html/unnamed-chunk-3-1.png" width="100%" />

#### unCorNonDis

<img src="simulation_study_files/figure-html/unnamed-chunk-4-1.png" width="100%" />

## Simulation: vary noise and fold-change and compare with other schemes (concatenation/ensembles)
* concatenation_splsda, ensemble_splsda, diablo_full, diablo_null

<img src="Figures/simulationResults-1.png" width="100%" />

> **Simulation study: performance assessment and benchmarking.** Simulated datasets included different types of variables: correlated & discriminatory (corDis); uncorrelated & discriminatory (unCorDis); correlated & nondiscriminatory (corNonDis) and uncorrelated & nondiscriminatory (unCorNonDis) for different fold-changes between sample groups and different noise levels (see Supplementary Note). Integrative classifiers included DIABLO with either the full or null design, concatenation and ensemble-based sPLSDA classifiers and were all trained to select 90 variables across three multi-omics datasets. a) Classification error rates (10-fold cross-validation averaged over 50 simulations). Dashed line indicates a random performance (error rate = 50%). All methods perform similarly with the exception of DIABLO_full which has a higher error rate. b) Number of variables selected according to their type. DIABLO_full selected mainly variables that were correlated & discriminatory (corDis, red), whereas the other methods selected an equal number of correlated or uncorrelated discriminatory variables (corDis and unCorDis, red and blue).
