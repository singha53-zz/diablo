---
title: 'DIABLO – an integrative, multi-omics, multivariate method for multi-group classification'
author: 'Amrit Singh, Casey P. Shannon, Benoît Gautier, Florian Rohart, Michaël Vacher, Scott J. Tebbutt, and Kim-Anh Lê Cao'
output:
  html_document:
    toc: yes
    theme: united
    highlight: tango
    keep_md: yes
---



# Table of contents

## 1) [Multi-omics methods overview](https://github.com/singha53/diablo/tree/master/analyses/methods_overview)
* Describe types of multi-omic methods (unsupervised and supervised methods and variable selection (yes/no))
  
## 2) [Simulation Study](https://github.com/singha53/diablo/blob/master/analyses/simulation_study/simulation_study.md)
* how does DIABLO compares with existing integrative classification schemes (Concatenation and Ensemble)
* effect of the design (connection between datasets) on the types of variables selected (correlated vs. discriminatory) and error rate

## 3) [Benchmarking](https://github.com/singha53/diablo/blob/master/analyses/benchmarking/src/benchmarking_enrichmentConnectivity.md)
* apply various unsupervised and supervised approaches to real cancer datasets (colon, kidney, gbm and lung) and identify multi-omic panels of equal sizes and compare them with respect:
  + which variables were selected
  + overlap between variables
  + network characteristics
  + ability to discriminant between phenotypic groups
  + Numbers of significant gene sets using various gene set collections

## 4) [Case study 1: BRCA - Standard workflow of DIABLO](https://github.com/singha53/diablo/blob/master/analyses/casestudy1_brca/brca_analysis.md)
* demonstrate use of the tune function
* assess stability of panel features
* show ability of DIABLO to uncover features (as part of the multi-omic biomarker signature) already associated with breast cancer
* demonstrate possible visualizations with DIABLO models

## 5) [Case study 2: Asthma - Incoporation of other features (module-based analysis, cross-over study design) with the DIABLO framework](https://github.com/singha53/diablo/blob/master/analyses/casestudy2_asthma/src/asthma_analysis.md)
* describe module-based analysis and how to model cross-over study designs
* compare features selected using paired vs. unpaired approach
* explain biological insights

