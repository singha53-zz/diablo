---
title: "asthma_analysis"
bibliography: /Users/asingh/Dropbox/Manuscript/diablo/mybib.bib
thanks: "Replication files are available at https://github.com/singha53/diablo"
output:
  html_document:
    toc: yes
    theme: united
    highlight: tango
    keep_md: yes
biblio-style: /Users/asingh/Dropbox/Manuscript/diablo/genome-biology.csl
---



\pagebreak

# FEV1 profiles

![](asthma_analysis_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

# DIABLO

## tune keepX

![](asthma_analysis_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```
## $cells
## [1] 1 1
## 
## $gene.module
## [1] 8 3
## 
## $metabolite.module
## [1] 6 3
```

```
## # A tibble: 2 x 3
## # Groups:   comp [2]
##   keepX  comp     error
##   <chr> <chr>     <dbl>
## 1 2_3_5 comp1 0.1428571
## 2 1_3_3 comp2 0.1428571
```

## error rate of optimal keepX

```
## $max.dist
##                comp 1    comp 2
## post        0.2142857 0.1428571
## pre         0.2142857 0.1428571
## Overall.ER  0.2142857 0.1428571
## Overall.BER 0.2142857 0.1428571
## 
## $centroids.dist
##                comp 1    comp 2
## post        0.1428571 0.1428571
## pre         0.1428571 0.1428571
## Overall.ER  0.1428571 0.1428571
## Overall.BER 0.1428571 0.1428571
## 
## $mahalanobis.dist
##                comp 1    comp 2
## post        0.1428571 0.1428571
## pre         0.1428571 0.1428571
## Overall.ER  0.1428571 0.1428571
## Overall.BER 0.1428571 0.1428571
```

## AUC


```
## [1] TRUE
```

![](asthma_analysis_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
## Area under the curve: 96.43%
```

# Component plots

![](asthma_analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

# heatmap

![](asthma_analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

> the above heatmap shows cells as an individula feature name

![](asthma_analysis_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


# Circos plot

>Error in `row.names<-.data.frame`(`*tmp*`, value = value) : invalid 'row.names' length



## network



