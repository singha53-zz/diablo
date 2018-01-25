---
title: 'DIABLO – an integrative, multi-omics, multivariate method for multi-group classification'
author: 'Amrit Singh, Benoît Gautier, Casey P. Shannon, Michaël Vacher, Florian Rohart, Scott J. Tebbutt, and Kim-Anh Lê Cao'
bibliography: mybib.bib
thanks: "Replication files are available at https://github.com/singha53/diablo"
output:
  html_document:
    theme: united
    highlight: tango
    keep_md: yes
biblio-style: genome-biology.csl
---



# Table of contents

## 1) [Multi-omics methods overview]()
* Describe types of multi-omic methods (unsupervised and supervised methods and variable selection (yes/no))
  
## 2) [Simulation Study]()
* how does DIABLO compares with existing integrative classification schemes (Concatenation and Ensemble)
* effect of the design (connection between datasets) on the types of variables selected (correlated vs. discriminatory) and error rate

## 3) [Benchmarking]()
* apply various unsupervised and supervised approaches to real cancer datasets (colon, kidney, gbm and lung) and identify multi-omic panels of equal sizes and compare them with respect:
  + which variables were selected
  + overlap between variables
  + network characteristics
  + ability to discriminant between phenotypic groups
  + Numbers of significant gene sets using various gene set collections

## 4) [Case study 1: BRCA - Standard workflow of DIABLO]()
* demonstrate use of the tune function
* assess stability of panel features
* show ability of DIABLO to uncover features (as part of the multi-omic biomarker signature) already associated with breast cancer
* demonstrate possible visualizations with DIABLO models

## 5) [Case study 2: Asthma - Incoporation of other features (module-based analysis, cross-over study design) with the DIABLO framework]()
* describe module-based analysis and how to model cross-over study designs
* compare features selected using paired vs. unpaired approach
* explain biological insights

