---
title: "Methods overview"
author: "Amrit Singh"
date: '2018-03-20'
output:
  html_document:
    toc: yes
    theme: united
    highlight: tango
    keep_md: yes
---



![**Figure 1. Simulation study: performance assessment and benchmarking.** Simulated datasets included different types of variables: correlated & discriminatory (corDis); uncorrelated & discriminatory (unCorDis); correlated & nondiscriminatory (corNonDis) and uncorrelated & nondiscriminatory (unCorNonDis) for different fold-changes between sample groups and different noise levels (see Supplementary Note). Integrative classifiers included DIABLO with either the full or null design, concatenation and ensemble-based sPLSDA classifiers and were all trained to select 90 variables across three multi-omics datasets. a) Classification error rates (10-fold cross-validation averaged over 50 simulations). Dashed line indicates a random performance (error rate = 50%). All methods perform similarly with the exception of DIABLO_full which has a higher error rate. b) Number of variables selected according to their type. DIABLO_full selected mainly variables that were correlated & discriminatory (corDis, red), whereas the other methods selected an equal number of correlated or uncorrelated discriminatory variables (corDis and unCorDis, red and blue).](/Users/asingh/Dropbox/PROOF/Manuscript/mixOmics/diablo/analyses/methods_overview/Figure1_asEdit.png)





![**Supplementary Figure 1. Overview of approaches used for the integration of multiple high dimensional omics datasets using either unsupervised or supervised analyses.** Most integrative methods were developed for unsupervised analyses. Variable selection is an important feature of the methods to improve interpretation of these complex models. Various types of integrative methods are listed, ranging from Component-based that reduce the dimensionality of high-throughput omics datasets, Bayesian methods, Network-based and multi-step approaches which  include concatenation and ensemble approaches1. Concatenation-based approach combine multiple matrices and apply standard single omics analysis without taking into account the type of omics variable in the model. Ensemble-based approaches involve the development of independent models for each omics dataset, after which the outputs are combined using various voting schemes (e.g. majority vote, average vote). Methods name in courier font indicate the name of the R package. *Methods are coded in other languages are indicated below.  Abbreviations: JIVE: Joint and Individual Variation Explained6, *sMBPLS: sparse Multiblock Partial Least Squares (Matlab)7, SNMNMF: Sparse Network-regularized Multiple Non-negative Matrix Factorization (Matlab)8, MOFA: Multi-Omics Factor Analysis9, *CONEXIC: Copy Number and Expression In Cancer (Java)10, WGCNA: Weighted Gene Co-expression Network Analysis11, SNF: Similarity Network Fusion12, PANDA: Passing Attributes between Networks for Data Assimilation13, BCC: Bayesian Consensus Clustering14, *RIMBANET: Reconstructing Integrative Molecular Bayesian Networks (Perl)3; sPCA : sparse Principal Component Analysis15; sGCCA: sparse generalized canonical correlation analysis16; rGCCA : regularized generalized canonical correlation analysis17; NMF: Non-Negative Factorization (Matlab); MFA: Multiple Co-inertia Analysis (MCIA); Multiple Factor Analysis19; glmnet: Lasso and Elastic-Net Regularized Generalized Linear Models20; sPLSDA: sparse Partial Least Squares Discriminant Analysis21; stSVM Smoothed t-statistics Support Vector Machine22; GELnet: Generalized Elastic Net23; *ATHENA: Analysis Tool for Heritable and Environmental Network Associations (Perl)4; SVM: Support Vector Machine; RF: Random Forest24; GRridge: Adaptive group-regularized ridge regression25; *iBAG: integrative Bayesian Analysis of Genomics (R and Shiny)5](/Users/asingh/Dropbox/PROOF/Manuscript/mixOmics/diablo/analyses/methods_overview/DIABLO_workflow.png)

