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



![**Supplementary Figure 1. Overview of approaches used for the integration of multiple high dimensional omics datasets using either unsupervised or supervised analyses.** Most integrative methods were developed for unsupervised analyses. Variable selection is an important feature of the methods to improve interpretation of these complex models. Various types of integrative methods are listed, ranging from Component-based that reduce the dimensionality of high-throughput omics datasets, Bayesian methods, Network-based and multi-step approaches which  include concatenation and ensemble approaches1. Concatenation-based approach combine multiple matrices and apply standard single omics analysis without taking into account the type of omics variable in the model. Ensemble-based approaches involve the development of independent models for each omics dataset, after which the outputs are combined using various voting schemes (e.g. majority vote, average vote). Methods name in courier font indicate the name of the R package. *Methods are coded in other languages are indicated below.  Abbreviations: JIVE: Joint and Individual Variation Explained6, *sMBPLS: sparse Multiblock Partial Least Squares (Matlab)7, SNMNMF: Sparse Network-regularized Multiple Non-negative Matrix Factorization (Matlab)8, MOFA: Multi-Omics Factor Analysis9, *CONEXIC: Copy Number and Expression In Cancer (Java)10, WGCNA: Weighted Gene Co-expression Network Analysis11, SNF: Similarity Network Fusion12, PANDA: Passing Attributes between Networks for Data Assimilation13, BCC: Bayesian Consensus Clustering14, *RIMBANET: Reconstructing Integrative Molecular Bayesian Networks (Perl)3; sPCA : sparse Principal Component Analysis15; sGCCA: sparse generalized canonical correlation analysis16; rGCCA : regularized generalized canonical correlation analysis17; NMF: Non-Negative Factorization (Matlab); MFA: Multiple Co-inertia Analysis (MCIA); Multiple Factor Analysis19; glmnet: Lasso and Elastic-Net Regularized Generalized Linear Models20; sPLSDA: sparse Partial Least Squares Discriminant Analysis21; stSVM Smoothed t-statistics Support Vector Machine22; GELnet: Generalized Elastic Net23; *ATHENA: Analysis Tool for Heritable and Environmental Network Associations (Perl)4; SVM: Support Vector Machine; RF: Random Forest24; GRridge: Adaptive group-regularized ridge regression25; *iBAG: integrative Bayesian Analysis of Genomics (R and Shiny)5](https://github.com/singha53/diablo/blob/master/analyses/methods_overview/Figure1_asEdit.png)

References

1.	Huang, S., Chaudhary, K. & Garmire, L. X. More Is Better: Recent Progress in Multi-Omics Data Integration Methods. Front. Genet. 8, (2017).
2.	Akavia, U. D. et al. An Integrated Approach to Uncover Drivers of Cancer. Cell 143, 1005–1017 (2010).
3.	Zhu, J. et al. Stitching together multiple data dimensions reveals interacting metabolomic and transcriptomic networks that modulate cell regulation. PLoS Biol. 10, e1001301 (2012).
4.	Kim, D., Li, R., Dudek, S. M. & Ritchie, M. D. ATHENA: Identifying interactions between different levels of genomic data associated with cancer clinical outcomes using grammatical evolution neural network. BioData Min. 6, 23 (2013).
5.	Wang, W. et al. iBAG: integrative Bayesian analysis of high-dimensional multiplatform genomics data. Bioinformatics 29, 149–159 (2013).
6.	Lock, E. F., Hoadley, K. A., Marron, J. S. & Nobel, A. B. Joint and individual variation explained (JIVE) for integrated analysis of multiple data types. Ann. Appl. Stat. 7, 523–542 (2013).
7.	Zhang, S. et al. Discovery of multi-dimensional modules by integrative analysis of cancer genomic data. Nucleic Acids Res. 40, 9379–9391 (2012).
8.	Zhang, S., Li, Q., Liu, J. & Zhou, X. J. A novel computational framework for simultaneous integration of multiple types of genomic data to identify microRNA-gene regulatory modules. Bioinformatics 27, i401–i409 (2011).
9.	Argelaguet, R. et al. Multi-Omics factor analysis disentangles heterogeneity in blood cancer. bioRxiv 217554 (2017).
10.	An Integrated Approach to Uncover Drivers of Cancer: Cell. Available at: http://www.cell.com/abstract/S0092-8674(10)01293-6. (Accessed: 12th February 2018)
11.	Langfelder, P. & Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008).
12.	Wang, B. et al. Similarity network fusion for aggregating data types on a genomic scale. Nat. Methods 11, 333–337 (2014).
13.	Glass, K., Huttenhower, C., Quackenbush, J. & Yuan, G.-C. Passing messages between biological networks to refine predicted interactions. PLoS ONE 8, e64832 (2013).
14.	Lock, E. F. & Dunson, D. B. Bayesian consensus clustering. Bioinformatics 29, 2610–2616 (2013).
15.	Shen, H. & Huang, J. Sparse Principal Component Analysis via Regularized Low Rank Matrix Approximation. J. Multivar. Anal. 99, 1015–1034 (2007).
16.	Tenenhaus, A. et al. Variable selection for generalized canonical correlation analysis. Biostatistics 15, 569–583 (2014).
17.	González, I. et al. Highlighting relationships between heterogeneous biological data through graphical displays based on regularized canonical correlation analysis. J. Biol. Syst. 17, 173–199 (2009).
18.	van der Maaten, L. & Hinton, G. Visualizing Data using t-SNE. J. Mach. Learn. Res. 1, 1–48 (2008).
19.	Abdi, H., Williams, L. J. & Valentin, D. Multiple factor analysis: principal component analysis for multitable and multiblock data sets. Wiley Interdiscip. Rev. Comput. Stat. 5, 149–179 (2013).
20.	Zou, H. & Hastie, T. Regularization and variable selection via the elastic net. J. R. Stat. Soc. Ser. B Stat. Methodol. 67, 301–320 (2005).
21.	Lê Cao, K.-A., Boitard, S. & Besse, P. Sparse PLS discriminant analysis: biologically relevant feature selection and graphical displays for multiclass problems. BMC Bioinformatics 12, 253 (2011).
22.	Cun, Y. & Fröhlich, H. Network and data integration for biomarker signature discovery via network smoothed t-statistics. PLoS ONE 8, e73074 (2013).
23.	Sokolov, A., Carlin, D. E., Paull, E. O., Baertsch, R. & Stuart, J. M. Pathway-based genomics prediction using generalized elastic net. PLoS Comput Biol 12, e1004790 (2016).
24.	Breiman, L. Random forests. Mach. Learn. 45, 5–32 (2001).
25.	van de Wiel, M. A., Lien, T. G., Verlaat, W., van Wieringen, W. N. & Wilting, S. M. Better prediction by use of co-data: adaptive group-regularized ridge regression. Stat. Med. 35, 368–381 (2016).
26.	Rohart, F., Gautier, B., Singh, A. & Cao, K.-A. L. mixOmics: An R package for ‘omics feature selection and multiple data integration. PLOS Comput. Biol. 13, e1005752 (2017).

