---
title: "Benchmarking"
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



# SNF data description
  * The SNF datasets were part of the datasets used in the Nature Methods paper on Similarity Network Fusion (SNF); https://www.nature.com/articles/nmeth.2810
  * The cancer datasets include GBM (Brain), Colon, Kidney, Lung and Breast (the Breast cancer dataset was excluded in order to avoid confusion with the case study on Breast Cancer)
  * The datasets were obtained from: http://compbio.cs.toronto.edu/SNF/SNF/Software.html
  * Survival times were provided for each disease cohort. The median survival time was used to dictomize each response variables into low and high survival times.

## number of samples in each group


```r
addmargins(sapply(snf_group, table))
```

```
##      colon kidney gbm lung Sum
## high    33     61 105   53 252
## low     59     61 108   53 281
## Sum     92    122 213  106 533
```

## number of variables in each dataset
* mRNA transcripts or cpg probes that mapped to the same gene were averaged 


```r
sapply(snf_data, function(i) sapply(i, ncol))
```

```
##       colon kidney   gbm  lung
## mrna  17814  17665 12042 12042
## mirna   312    329   534   352
## cpg   13381  13680   750 13289
```

# Multi-omic biomarker panels

## Unsupervised

### JIVE


```r
jive_joints <- list(colon = t(do.call(rbind, jive_colon$joint)), kidney = t(do.call(rbind, 
    jive_kidney$joint)), gbm = t(do.call(rbind, jive_gbm$joint)), lung = t(do.call(rbind, 
    jive_lung$joint)))
colnames(jive_joints$colon) <- sapply(snf_data$colon, colnames) %>% unlist
colnames(jive_joints$kidney) <- sapply(snf_data$kidney, colnames) %>% unlist
colnames(jive_joints$gbm) <- sapply(snf_data$gbm, colnames) %>% unlist
colnames(jive_joints$lung) <- sapply(snf_data$lung, colnames) %>% unlist

snf_jive <- mapply(function(x, y) {
    # run sPCA
    pca_jive = spca(x, ncomp = 2, center = FALSE, scale = FALSE, keepX = rep(ncol(x), 
        2))
    
    ## scores
    pcs <- pca_jive$variates$X %>% as.data.frame %>% mutate(pheno = y, Method = "JIVE")
    
    ## features
    comp1 <- split(pca_jive$loadings$X[, 1], factor(sapply(strsplit(names(pca_jive$loadings$X[, 
        1]), "_"), function(i) i[1]), c("mrna", "mirna", "cpg")))
    comp2 <- split(pca_jive$loadings$X[, 2], factor(sapply(strsplit(names(pca_jive$loadings$X[, 
        2]), "_"), function(i) i[1]), c("mrna", "mirna", "cpg")))
    
    panel <- mapply(function(x, y) {
        c(names(x[order(abs(x), decreasing = TRUE)][1:30]), names(y[order(abs(y), 
            decreasing = TRUE)][1:30]))
    }, x = comp1, y = comp2, SIMPLIFY = FALSE)
    
    return(list(pcs = pcs, panel = panel))
}, x = jive_joints, y = snf_group, SIMPLIFY = FALSE) %>% zip_nPure()
```

### MOFA


```r
mofa_joints <- list(colon = MOFA_colon, kidney = MOFA_kidney, gbm = MOFA_gbm, 
    lung = MOFA_lung)

snf_mofa <- mapply(function(x, y, z) {
    ## scores
    pcs <- getFactors(x, factors = "all", as.data.frame = FALSE, include_intercept = TRUE) %>% 
        as.data.frame()
    colnames(pcs) <- c("intercept", "PC1", "PC2")
    pcs <- pcs[, -1] %>% mutate(pheno = y, Method = "MOFA")
    
    ## features
    features_mofa = getWeights(x, views = "all", factors = 1:2)
    features_mofa <- features_mofa[names(z)]
    panel <- mapply(function(x, y) {
        c(colnames(y)[order(abs(x[, 1]), decreasing = TRUE)[1:30]], colnames(y)[order(abs(x[, 
            2]), decreasing = TRUE)[1:30]])
    }, x = features_mofa, y = z, SIMPLIFY = FALSE)
    
    return(list(pcs = pcs, panel = panel))
}, x = mofa_joints, y = snf_group, z = snf_data, SIMPLIFY = FALSE) %>% zip_nPure()
```

### sGCCA


```r
## design matrix
design <- matrix(1, nrow = 3, ncol = 3)
rownames(design) <- colnames(design) <- names(snf_data$colon)
diag(design) <- 0

keepX = lapply(snf_mofa$panel, function(i) {
    lapply(i, function(i) {
        rep(length(i)/2, 2)
    })
})

snf_sgcca <- mapply(function(x, y, z) {
    result.unsupervised = wrapper.sgcca(X = x, ncomp = 2, keepX = z, design = design)
    
    ## scores
    pcs <- Reduce("+", result.unsupervised$variates)/2
    colnames(pcs) <- c("PC1", "PC2")
    pcs <- pcs %>% as.data.frame() %>% mutate(pheno = y, Method = "sGCCA")
    
    ## features
    feat1 <- lapply(selectVar(result.unsupervised, comp = 1), function(i) i["name"])
    feat2 <- lapply(selectVar(result.unsupervised, comp = 2), function(i) i["name"])
    
    panel <- list(mrna = c(feat1$mrna$name, feat2$mrna$name), mirna = c(feat1$mirna$name, 
        feat2$mirna$name), cpg = c(feat1$cpg$name, feat2$cpg$name))
    
    return(list(pcs = pcs, panel = panel))
}, x = snf_data, y = snf_group, z = keepX, SIMPLIFY = FALSE) %>% zip_nPure()
```

## Supervised

### Concatenation_sPLSDA


```r
snf_concat_splsda <- mapply(function(x, y) {
    X = do.call(cbind, x)
    concat <- splsda(X = do.call(cbind, x), Y = y, ncomp = 2, keepX = rep(ncol(X), 
        2))
    ## scores
    pcs <- concat$variates$X %>% as.data.frame
    colnames(pcs) <- c("PC1", "PC2")
    pcs <- pcs %>% mutate(pheno = y, Method = "Concatenation")
    
    ## features
    comp1 <- split(concat$loadings$X[, 1], factor(sapply(strsplit(names(concat$loadings$X[, 
        1]), "_"), function(i) i[1]), c("mrna", "mirna", "cpg")))
    comp2 <- split(concat$loadings$X[, 2], factor(sapply(strsplit(names(concat$loadings$X[, 
        2]), "_"), function(i) i[1]), c("mrna", "mirna", "cpg")))
    
    panel <- mapply(function(x, y) {
        c(names(x[order(abs(x), decreasing = TRUE)][1:30]), names(y[order(abs(y), 
            decreasing = TRUE)][1:30]))
    }, x = comp1, y = comp2, SIMPLIFY = FALSE)
    
    return(list(pcs = pcs, panel = panel))
}, x = snf_data, y = snf_group, SIMPLIFY = FALSE) %>% zip_nPure()
```

### Ensemble_spslda


```r
snf_ensemble_splsda <- mapply(function(x, y, z) {
    ensem <- mapply(function(a, b) {
        result <- splsda(X = a, Y = y, keepX = b, ncomp = 2)
        
        ## scores
        pcs <- result$variates$X %>% as.data.frame
        colnames(pcs) <- c("PC1", "PC2")
        
        ## panels
        feat <- c(selectVar(result, comp = 1)$name, selectVar(result, comp = 2)$name)
        
        return(list(pcs = pcs, feat = feat))
    }, a = x, b = z, SIMPLIFY = FALSE) %>% zip_nPure()
    
    ## scores
    pcs <- as.data.frame(Reduce("+", ensem$pcs)/3) %>% mutate(pheno = y, Method = "Ensemble")
    
    ## features
    panel <- ensem$feat
    return(list(pcs = pcs, panel = panel))
}, x = snf_data, y = snf_group, z = keepX, SIMPLIFY = FALSE) %>% zip_nPure()
```

### DIABLO_null


```r
## design matrix
design <- matrix(0, nrow = 3, ncol = 3)
rownames(design) <- colnames(design) <- names(snf_data$colon)
diag(design) <- 0

snf_diabloNull <- mapply(function(x, y, z) {
    result.supervised = block.splsda(X = x, Y = y, ncomp = 2, keepX = z, design = design)
    
    ## scores
    pcs <- Reduce("+", result.supervised$variates)/length(x)
    colnames(pcs) <- c("PC1", "PC2")
    pcs <- pcs %>% as.data.frame() %>% mutate(pheno = y, Method = "DIABLO_null")
    
    ## features
    feat1 <- lapply(selectVar(result.supervised, comp = 1), function(i) i["name"])
    feat2 <- lapply(selectVar(result.supervised, comp = 2), function(i) i["name"])
    
    panel <- list(mrna = c(feat1$mrna$name, comp2 = feat2$mrna$name), mirna = c(feat1$mirna$name, 
        comp2 = feat2$mirna$name), cpg = c(feat1$cpg$name, comp2 = feat2$cpg$name))
    
    return(list(pcs = pcs, panel = panel))
}, x = snf_data, y = snf_group, z = keepX, SIMPLIFY = FALSE) %>% zip_nPure()
```

### DIABLO_full


```r
## design matrix
design <- matrix(1, nrow = 3, ncol = 3)
rownames(design) <- colnames(design) <- names(snf_data$colon)
diag(design) <- 0

snf_diabloFull <- mapply(function(x, y, z) {
    result.supervised = block.splsda(X = x, Y = y, ncomp = 2, keepX = z, design = design)
    
    ## scores
    pcs <- Reduce("+", result.supervised$variates)/length(x)
    colnames(pcs) <- c("PC1", "PC2")
    pcs <- pcs %>% as.data.frame() %>% mutate(pheno = y, Method = "DIABLO_full")
    
    ## features
    feat1 <- lapply(selectVar(result.supervised, comp = 1), function(i) i["name"])
    feat2 <- lapply(selectVar(result.supervised, comp = 2), function(i) i["name"])
    
    panel <- list(mrna = c(feat1$mrna$name, comp2 = feat2$mrna$name), mirna = c(feat1$mirna$name, 
        comp2 = feat2$mirna$name), cpg = c(feat1$cpg$name, comp2 = feat2$cpg$name))
    
    return(list(pcs = pcs, panel = panel))
}, x = snf_data, y = snf_group, z = keepX, SIMPLIFY = FALSE) %>% zip_nPure()
```

\pagebreak

# Number of features per panel


```r
multiOmicPanels <- list(JIVE = snf_jive$panel, MOFA = snf_mofa$panel, sGCCA = snf_sgcca$panel, 
    Concatenation = snf_concat_splsda$panel, Ensemble = snf_ensemble_splsda$panel, 
    DIABLO_null = snf_diabloNull$panel, DIABLO_full = snf_diabloFull$panel)

allPanels <- list(JIVE = snf_jive$panel, MOFA = snf_mofa$panel, sGCCA = snf_sgcca$panel, 
    Concatenation = snf_concat_splsda$panel, Ensemble = snf_ensemble_splsda$panel, 
    DIABLO_null = snf_diabloNull$panel, DIABLO_full = snf_diabloFull$panel) %>% 
    rapply(., function(i) {
        sapply(strsplit(i, "_"), function(j) j[[2]])
    }, how = "list")

panels <- rapply(multiOmicPanels, function(i) {
    length(unique(i))
}, how = "list") %>% lapply(., function(j) {
    do.call(rbind, j) %>% as.data.frame() %>% mutate(disease = rownames(.))
}) %>% do.call(rbind, .) %>% as.data.frame() %>% mutate(method = sapply(strsplit(rownames(.), 
    "\\."), function(i) i[1])) %>% gather(omic, nFeat, -c(disease:method)) %>% 
    mutate(nFeat = as.numeric(nFeat)) %>% ggplot(aes(x = method, y = nFeat, 
    fill = omic)) + geom_bar(stat = "identity") + facet_grid(~disease) + customTheme(sizeStripFont = 10, 
    xAngle = 45, hjust = 1, vjust = 1, xSize = 10, ySize = 10, xAxisSize = 10, 
    yAxisSize = 10)
panels
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-10-1.png" width="100%" />

```r
## single omic panels
mrna <- lapply(allPanels, function(i) {
    lapply(i, function(j) {
        j[["mrna"]]
    })
})
mirna <- lapply(allPanels, function(i) {
    lapply(i, function(j) {
        j[["mirna"]]
    })
})
cpg <- lapply(allPanels, function(i) {
    lapply(i, function(j) {
        j[["cpg"]]
    })
})

mrna.cpg <- lapply(names(mrna), function(i) {
    x <- lapply(names(mrna$JIVE), function(j) {
        unique(unlist(strsplit(c(mrna[[i]][[j]], cpg[[i]][[j]]), ";")))
    })
    names(x) <- names(mrna$JIVE)
    x
})
names(mrna.cpg) <- names(mrna)

mrna.mirna.cpg <- lapply(names(mrna), function(i) {
    x <- lapply(names(mrna$JIVE), function(j) {
        unique(unlist(strsplit(c(mrna[[i]][[j]], mirna[[i]][[j]], cpg[[i]][[j]]), 
            ";")))
    })
    names(x) <- names(mrna$JIVE)
    x
})
names(mrna.mirna.cpg) <- names(mrna)
```

# Component plots


```r
allscores <- rbind(do.call(rbind, snf_jive$pcs), do.call(rbind, snf_mofa$pcs), 
    do.call(rbind, snf_sgcca$pcs), do.call(rbind, snf_concat_splsda$pcs), do.call(rbind, 
        snf_ensemble_splsda$pcs), do.call(rbind, snf_diabloNull$pcs), do.call(rbind, 
        snf_diabloFull$pcs)) %>% as.data.frame %>% mutate(Disease = rep(rep(names(snf_group), 
    sapply(snf_group, length)), 7))

allCompPlots <- ggplot(allscores, aes(x = PC1, y = PC2, group = pheno, color = pheno)) + 
    geom_point() + facet_wrap(Disease ~ Method, scales = "free", ncol = 7) + 
    stat_ellipse(level = 0.8) + customTheme(sizeStripFont = 10, xAngle = 0, 
    hjust = 0.5, vjust = 0.5, xSize = 10, ySize = 10, xAxisSize = 10, yAxisSize = 10) + 
    xlab("Component 1") + ylab("Component 2")
allCompPlots
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-11-1.png" width="100%" />

```r
pdf(paste0(WhereAmI, "results/Figures/multiOmicPanels_allcomponentplots.pdf"), 
    width = 13, height = 12)
allCompPlots
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Overlap in panels

## Colon

### Intersection plot


```r
colon_panels <- list(JIVE = unlist(snf_jive$panel$colon), MOFA = unlist(snf_mofa$panel$colon), 
    sGCCA = unlist(snf_sgcca$panel$colon), Concatenation = unlist(snf_concat_splsda$panel$colon), 
    Ensemble = unlist(snf_ensemble_splsda$panel$colon), DIABLO_null = unlist(snf_diabloNull$panel$colon), 
    DIABLO_full = unlist(snf_diabloFull$panel$colon))

colonInput <- fromList(colon_panels)
metadata <- data.frame(approaches = colnames(colonInput))
metadata$type <- "supervised"
metadata$type[metadata$approaches %in% c("JIVE", "sGCCA", "MOFA")] <- "unsupervised"

png(paste0(WhereAmI, "results/Figures/colon_overlap.png"), width = 7, height = 3.5, 
    units = "in", res = 300)
upset(colonInput, sets = colnames(colonInput), keep.order = TRUE, queries = list(list(query = intersects, 
    params = list("Concatenation", "DIABLO_null", "Ensemble"), active = TRUE, 
    color = "#56B4E9"), list(query = intersects, params = list("JIVE", "MOFA", 
    "sGCCA"), active = TRUE, color = "#F0E442")), set.metadata = list(data = metadata, 
    plots = list(list(type = "matrix_rows", column = "type", colors = c(supervised = "green", 
        unsupervised = "purple"), alpha = 0.5))))
grid.text("Colon", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
dev.off()
```

```
## quartz_off_screen 
##                 2
```

### Venn diagram


```r
venn(colon_panels, zcolor = "style", cexsn = 1, cexil = 1.3)
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-13-1.png" width="100%" />

## Kidney

### Intersection plot


```r
kidney_panels <- list(JIVE = unlist(snf_jive$panel$kidney), MOFA = unlist(snf_mofa$panel$kidney), 
    sGCCA = unlist(snf_sgcca$panel$kidney), Concatenation = unlist(snf_concat_splsda$panel$kidney), 
    Ensemble = unlist(snf_ensemble_splsda$panel$kidney), DIABLO_null = unlist(snf_diabloNull$panel$kidney), 
    DIABLO_full = unlist(snf_diabloFull$panel$kidney))

kidneyInput <- fromList(kidney_panels)
metadata <- data.frame(approaches = colnames(kidneyInput))
metadata$type <- "supervised"
metadata$type[metadata$approaches %in% c("JIVE", "sGCCA", "MOFA")] <- "unsupervised"
png(paste0(WhereAmI, "results/Figures/kidney_overlap.png"), width = 7, height = 3.5, 
    units = "in", res = 300)
upset(kidneyInput, sets = colnames(kidneyInput), keep.order = TRUE, queries = list(list(query = intersects, 
    params = list("Concatenation", "DIABLO_null", "Ensemble"), active = TRUE, 
    color = "#56B4E9"), list(query = intersects, params = list("JIVE", "MOFA", 
    "sGCCA"), active = TRUE, color = "#F0E442")), set.metadata = list(data = metadata, 
    plots = list(list(type = "matrix_rows", column = "type", colors = c(supervised = "green", 
        unsupervised = "purple"), alpha = 0.5))))
grid.text("Kidney", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
dev.off()
```

```
## quartz_off_screen 
##                 2
```

### Venn diagram


```r
venn(kidney_panels, zcolor = "style")
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-15-1.png" width="100%" />

## GBM

### Intersection plot


```r
gbm_panels <- list(JIVE = unlist(snf_jive$panel$gbm), MOFA = unlist(snf_mofa$panel$gbm), 
    sGCCA = unlist(snf_sgcca$panel$gbm), Concatenation = unlist(snf_concat_splsda$panel$gbm), 
    Ensemble = unlist(snf_ensemble_splsda$panel$gbm), DIABLO_null = unlist(snf_diabloNull$panel$gbm), 
    DIABLO_full = unlist(snf_diabloFull$panel$gbm))

gbmInput <- fromList(gbm_panels)
metadata <- data.frame(approaches = colnames(gbmInput))
metadata$type <- "supervised"
metadata$type[metadata$approaches %in% c("JIVE", "sGCCA", "MOFA")] <- "unsupervised"
png(paste0(WhereAmI, "results/Figures/gbm_overlap.png"), width = 7, height = 3.5, 
    units = "in", res = 300)
upset(gbmInput, sets = colnames(gbmInput), keep.order = TRUE, queries = list(list(query = intersects, 
    params = list("Concatenation", "DIABLO_null", "Ensemble"), active = TRUE, 
    color = "#56B4E9"), list(query = intersects, params = list("JIVE", "MOFA", 
    "sGCCA"), active = TRUE, color = "#F0E442")), set.metadata = list(data = metadata, 
    plots = list(list(type = "matrix_rows", column = "type", colors = c(supervised = "green", 
        unsupervised = "purple"), alpha = 0.5))))
grid.text("GBM", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
dev.off()
```

```
## quartz_off_screen 
##                 2
```

### Venn diagram


```r
venn(gbm_panels, zcolor = "style")
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-17-1.png" width="100%" />

## Lung

### Intersection plot


```r
lung_panels <- list(JIVE = unlist(snf_jive$panel$lung), MOFA = unlist(snf_mofa$panel$lung), 
    sGCCA = unlist(snf_sgcca$panel$lung), Concatenation = unlist(snf_concat_splsda$panel$lung), 
    Ensemble = unlist(snf_ensemble_splsda$panel$lung), DIABLO_null = unlist(snf_diabloNull$panel$lung), 
    DIABLO_full = unlist(snf_diabloFull$panel$lung))

lungInput <- fromList(lung_panels)
metadata <- data.frame(approaches = colnames(lungInput))
metadata$type <- "supervised"
metadata$type[metadata$approaches %in% c("JIVE", "sGCCA", "MOFA")] <- "unsupervised"
png(paste0(WhereAmI, "results/Figures/lung_overlap.png"), width = 7, height = 3.5, 
    units = "in", res = 300)
upset(lungInput, sets = colnames(lungInput), keep.order = TRUE, queries = list(list(query = intersects, 
    params = list("Concatenation", "DIABLO_null", "Ensemble"), active = TRUE, 
    color = "#56B4E9"), list(query = intersects, params = list("JIVE", "MOFA", 
    "sGCCA"), active = TRUE, color = "#F0E442")), set.metadata = list(data = metadata, 
    plots = list(list(type = "matrix_rows", column = "type", colors = c(supervised = "green", 
        unsupervised = "purple"), alpha = 0.5))))
grid.text("Lung", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
dev.off()
```

```
## quartz_off_screen 
##                 2
```

### Venn diagram


```r
venn(lung_panels, zcolor = "style")
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-19-1.png" width="100%" />

# Gene set enrichment analysis
  * We wished to assess the enrichment of the selected features across a variety of annotated gene sets in the MSigDB collection (http://software.broadinstitute.org/gsea/msigdb), in particular:
1. C1 - positional gene sets  for each human chromosome and cytogenetic band.
2. C2 – curated gene sets (Pathway Interaction DB [PID], Biocarta [BIOCARTA], Kyoto Encyclopedia of Genes and Genomes [KEGG], Reactome [REACTOME], and others)
3. C3 - 	motif gene sets  based on conserved cis-regulatory motifs from a comparative analysis of the human, mouse, rat, and dog genomes.
4.	C4 – computational gene sets (from the Cancer Gene Neighbourhoods [CGN] and Cancer Modules [CM] – citation available via: http://www.broadinstitute.org/gsea/msigdb/collections.jsp)
5. C5 - GO gene sets  consist of genes annotated by the same GO terms.
6.	C6 – ontologic gene sets (Gene sets represent signatures of cellular pathways which are often dis-regulated in cancer).
7. C7 - immunologic gene sets  defined directly from microarray gene expression data from immunologic studies.
8. H - hallmark gene sets  are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.
&
A. BTM - Blood Transcriptional Modules (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2727981/)
B. TISSUES - cell-specific expression from Benita et al. Blood 2008 http://www.bloodjournal.org/content/115/26/5376

  *  Significance of enrichment was determined using a hypergeometric test of the overlap between the selected features (mapped to official HUGO gene symbols or official miRNA symbols) and the various gene sets contained in the collections. Resulting p-values were corrected for multiple hypothesis using the Benjamini-Hochberg procedure applied across ALL genesets (10k+ tests – as pessimistic as possible). Adjusted p-values are reported in the fdr column.

## mRNA and CpGs


```r
enriched_sets <- rapply(mrna.cpg, function(i) {
    sear(i, "mrna") %>% dplyr::group_by(collection) %>% dplyr::summarise(sig = sum(fdr < 
        0.05)) %>% dplyr::filter(collection != "ARCHIVED")
}, how = "list") %>% lapply(., function(i) {
    do.call(rbind, i) %>% mutate(disease = rep(names(i), each = nrow(i[[1]])))
})
enrichedPathways <- do.call(rbind, enriched_sets) %>% mutate(method = rep(names(enriched_sets), 
    each = nrow(enriched_sets[[1]])))
enrichedPathways$type <- "supervised"
enrichedPathways$type[enrichedPathways$method %in% c("JIVE", "MOFA", "sGCCA")] <- "unsupervised"

enrichedPathways %>% spread(collection, sig)
```

```
## # A tibble: 28 x 13
##    disease        method         type   BTM    C1    C2    C3    C4    C5
##  *   <chr>         <chr>        <chr> <int> <int> <int> <int> <int> <int>
##  1   colon Concatenation   supervised     0     0     2    12     2     2
##  2   colon   DIABLO_full   supervised    16     1    43     0    37   118
##  3   colon   DIABLO_null   supervised     0     0     5     6     2     1
##  4   colon      Ensemble   supervised     0     0     3     2     3     1
##  5   colon          JIVE unsupervised     0     0     9    20     0    24
##  6   colon          MOFA unsupervised     9     0    44     0    13    68
##  7   colon         sGCCA unsupervised     0     0     3     4     0    27
##  8     gbm Concatenation   supervised     1     0   265    34    34   494
##  9     gbm   DIABLO_full   supervised    16     0   590    38    82   969
## 10     gbm   DIABLO_null   supervised     2     0   280    52    42   685
## # ... with 18 more rows, and 4 more variables: C6 <int>, C7 <int>,
## #   H <int>, TISSUES <int>
```

```r
enrichedPathways %>% spread(collection, sig) %>% write.csv(., paste0(WhereAmI, 
    "results/Tables/multiOmicPanels_biologicalEnrichment.csv"))
```

### which method is leads to the greatest number of signficant pathways?


```r
enrichedPathways %>% group_by(collection, disease) %>% dplyr::filter(sig == 
    max(sig)) %>% dplyr::select(disease, collection, method) %>% dplyr::group_by(disease, 
    collection) %>% dplyr::summarise(method = paste(as.character(method), collapse = " & ")) %>% 
    ungroup %>% group_by(method) %>% dplyr::summarise(n_sig = n()) %>% as.data.frame()
```

```
##                                                                        method
## 1                                                               Concatenation
## 2                                                                 DIABLO_full
## 3                                                                 DIABLO_null
## 4                                                   DIABLO_null & DIABLO_full
## 5                                                                    Ensemble
## 6                                                                        JIVE
## 7  JIVE & MOFA & sGCCA & Concatenation & Ensemble & DIABLO_null & DIABLO_full
## 8                                                                        MOFA
## 9                                          MOFA & Concatenation & DIABLO_null
## 10                                                               MOFA & sGCCA
## 11                                                                      sGCCA
##    n_sig
## 1      1
## 2     21
## 3      1
## 4      1
## 5      1
## 6      7
## 7      1
## 8      3
## 9      1
## 10     1
## 11     2
```

```r
enrichedPathways %>% group_by(type, disease, collection) %>% dplyr::filter(sig == 
    max(sig)) %>% dplyr::select(type, disease, collection, method) %>% dplyr::group_by(type, 
    disease, collection) %>% dplyr::summarise(method = paste(as.character(method), 
    collapse = " & ")) %>% ungroup %>% group_by(type, method) %>% dplyr::summarise(n_sig = n()) %>% 
    as.data.frame()
```

```
##            type                                               method n_sig
## 1    supervised                                        Concatenation     2
## 2    supervised                          Concatenation & DIABLO_null     1
## 3    supervised Concatenation & Ensemble & DIABLO_null & DIABLO_full     2
## 4    supervised                                          DIABLO_full    27
## 5    supervised                                          DIABLO_null     3
## 6    supervised                            DIABLO_null & DIABLO_full     1
## 7    supervised                                             Ensemble     4
## 8  unsupervised                                                 JIVE    19
## 9  unsupervised                                  JIVE & MOFA & sGCCA     5
## 10 unsupervised                                         JIVE & sGCCA     1
## 11 unsupervised                                                 MOFA     9
## 12 unsupervised                                         MOFA & sGCCA     1
## 13 unsupervised                                                sGCCA     5
```

# Connectivity


```r
saveRDS(multiOmicPanels, paste0(WhereAmI, "results/multiOmicPanels.rds"))

# also save various correlation matrices Pearson
cormat_pearson <- lapply(multiOmicPanels, function(i) {
    mapply(function(x, y) {
        do.call(cbind, y)[, unlist(x)] %>% cor(., method = "pearson")
    }, x = i, y = snf_data, SIMPLIFY = FALSE)
})
saveRDS(cormat_pearson, paste0(WhereAmI, "results/cormat_pearson.rds"))
```

## datasets of multi-omic panels


```r
data <- multiOmicPanels %>% purrr::map(~{
    # for every model for every experiment for every list of features and list
    # of datatypes
    purrr::map2(., snf_data, ~{
        purrr::map2(.x, .y, ~{
            .y[, .x]
        })
    })
})
names(data)
```

```
## [1] "JIVE"          "MOFA"          "sGCCA"         "Concatenation"
## [5] "Ensemble"      "DIABLO_null"   "DIABLO_full"
```

```r
# combine
data <- purrr::modify_depth(data, 2, purrr::reduce, cbind)
```

## adjacency matrices


```r
adj <- modify_depth(data, 2, cor)
```

## Number of connections


```r
cor_mat <- cormat_pearson
plots <- list()

nConnectionsDat <- cor_mat %>% modify_depth(2, ~tibble(cor = .[lower.tri(., 
    diag = F)])) %>% purrr::map(bind_rows, .id = "dataset") %>% bind_rows(.id = "model") %>% 
    mutate(cor = abs(cor)) %>% group_by(model, dataset) %>% mutate(cor = abs(cor)) %>% 
    summarise(cor_0.5 = sum(cor > 0.5), cor_0.65 = sum(cor > 0.55), cor_0.6 = sum(cor > 
        0.6), cor_0.65 = sum(cor > 0.65), cor_0.7 = sum(cor > 0.7), cor_0.75 = sum(cor > 
        0.75), cor_0.8 = sum(cor > 0.8), cor_0.85 = sum(cor > 0.85), cor_0.9 = sum(cor > 
        0.9), cor_0.95 = sum(cor > 0.95), cor_1 = sum(cor == 1)) %>% gather(cor, 
    nConnections, -c(model:dataset)) %>% mutate(cor = as.numeric(sapply(strsplit(cor, 
    "_"), function(i) i[2])))
nConnectionsDat$type <- "supervised"
nConnectionsDat$type[nConnectionsDat$model %in% c("JIVE", "MOFA", "sGCCA")] <- "unsupervised"

nconnections <- ggplot(nConnectionsDat, aes(x = cor, y = nConnections, fill = model, 
    color = model, linetype = type)) + geom_line(size = 1) + facet_wrap(~dataset, 
    nr = 1) + customTheme(sizeStripFont = 15, xAngle = 0, hjust = 0.5, vjust = 0.5, 
    xSize = 10, ySize = 10, xAxisSize = 10, yAxisSize = 10) + xlab("Absolute correlation coefficient cut-off") + 
    ylab("Number of edges")
nconnections
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-25-1.png" width="100%" />

```r
pdf(paste0(WhereAmI, "results/Figures/multiOmicsPanels_nConnections.pdf"), width = 14, 
    height = 4)
nconnections
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Network attributes


```r
df_list <- list()
corSeq <- seq(0.5, 0.95, 0.05)
for (i in 1:length(corSeq)) {
    cor_cutoff <- corSeq[i]
    df <- cor_mat %>% modify_depth(2, ~tibble(adj = I(list(.)))) %>% purrr::map(bind_rows, 
        .id = "data") %>% bind_rows(.id = "model")
    
    # graphs
    df$net <- df$adj %>% purrr::map(~{
        .[abs(.) < cor_cutoff] <- 0
        .[abs(.) >= cor_cutoff] <- 1
        .
    }) %>% purrr::map(igraph::graph_from_adjacency_matrix, mode = "lower", weighted = NULL, 
        diag = F)
    
    ## attributes
    df$clusters <- purrr::map(df$net, igraph::cluster_edge_betweenness)
    df$ncommunity <- purrr::map(df$clusters, length) %>% unlist()
    # df$modularity <- purrr::map(df$clusters, modularity) %>% unlist()
    # df$transitivity <- purrr::map(df$net, igraph::transitivity) %>% unlist()
    df$triads <- purrr::map(df$net, igraph::triad_census) %>% purrr::map(last) %>% 
        unlist()
    df$density <- purrr::map(df$net, igraph::edge_density) %>% purrr::map(last) %>% 
        unlist()
    df$cor_cutoff <- cor_cutoff
    df_list[[i]] <- df
    
}

attributesDat <- do.call(rbind, df_list) %>% gather(attributes, value, ncommunity:density)
attributesDat$type <- "supervised"
attributesDat$type[attributesDat$model %in% c("JIVE", "MOFA", "sGCCA")] <- "unsupervised"

netAttributes <- ggplot(attributesDat, aes(x = cor_cutoff, y = value, fill = model, 
    color = model, linetype = type)) + geom_point() + geom_line() + facet_grid(attributes ~ 
    data, scales = "free") + customTheme(sizeStripFont = 15, xAngle = 0, hjust = 0.5, 
    vjust = 0.5, xSize = 10, ySize = 10, xAxisSize = 10, yAxisSize = 10)
netAttributes
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-26-1.png" width="100%" />

```r
pdf(paste0(WhereAmI, "results/Figures/multiOmicsPanels_networkAttributes.pdf"), 
    width = 10, height = 10)
netAttributes
dev.off()
```

```
## quartz_off_screen 
##                 2
```

## all network plots


```r
df <- cor_mat %>%
modify_depth(2, ~ tibble(adj = I(list(.)))) %>%
purrr::map(bind_rows, .id = 'data') %>%
bind_rows(.id = 'model')

# graphs
df$net <- df$adj%>%
  purrr::map(~ {
    .[abs(.) < 0.25] <- 0
    .[abs(.) >= 0.25] <- 1
    .
  }) %>%
  purrr::map(igraph::graph_from_adjacency_matrix, mode = 'lower', weighted = NULL, diag = F)
df$clusters <- purrr::map(df$net, igraph::cluster_edge_betweenness)

df$ggraphs <- pmap(list(c = df$clusters, g = df$net) , function(c, g) {
  V(g)$community <- as.character(c$membership)
  V(g)$block <- gsub('^(.+)_.+$', '\\1', V(g)$name)
  ggraph(g, layout = 'igraph', algorithm = 'nicely') +
    # geom_edge_fan(colour = 'lightgrey', show.legend = FALSE) +
    geom_encircle(aes(x = x, y = y, group = community),
                  s_shape = 0.5, expand = 0.025, colour = 'lightgrey') +
    geom_node_point(aes(fill = block), shape = 21, colour = 'white', size = 4) +
    scale_x_continuous(expand = c(0.25, 0.25)) +
    scale_y_continuous(expand = c(0.25, 0.25)) +
    theme(aspect.ratio = 1.25,
          # legend.position = 'bottom',
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank())
})

# grab legend
legend <- get_legend(df$ggraphs[[1]])

# r1 <- plots$adjacency
nets <- lapply(unique(df$data), function(i){
  r3 <- filter(df, data == i)
  r3 <- map2(r3$ggraphs, r3$model, ~ .x + ggtitle(.y) + theme(legend.position = 'none'))
  r3 <- rev(r3)
  r3 <- plot_grid(
  r3[[1]],
  r3[[2]],
  r3[[3]],
  r3[[4]],
  r3[[5]],
  r3[[6]],
  r3[[7]],
  legend,
  nrow = 1, rel_widths = c(1, 1, 1, 1, 1, 1, 1, 0.4)
)
r3
})
names(nets) <- unique(df$data)

p <- plot_grid(nets$colon, nets$kidney, nets$gbm, nets$lung, nrow = 4, ncol = 1, labels = c('', ''))
p
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-27-1.png" width="100%" />

```r
pdf(paste0(WhereAmI, "results/Figures/multiOmicsPanels_networks.pdf"), width = 12, height = 10)
p
dev.off()
```

```
## quartz_off_screen 
##                 2
```

## network and component plot of the multi-omic panels derived using the colon cancer dataset


```r
nets$colon
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-28-1.png" width="100%" />

```r
pdf(paste0(WhereAmI, "results/Figures/multiOmicsPanels_networks_colon.pdf"), 
    width = 14, height = 6)
nets$colon
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
## add colon component plot
compPlot <- allscores %>% filter(Disease == "colon") %>% mutate(Method = factor(Method, 
    c("DIABLO_full", "DIABLO_null", "Ensemble", "Concatenation", "sGCCA", "MOFA", 
        "JIVE"))) %>% ggplot(aes(x = PC1, y = PC2, group = pheno, color = pheno)) + 
    geom_point() + facet_wrap(Disease ~ Method, scales = "free", ncol = 7) + 
    stat_ellipse(level = 0.8) + customTheme(sizeStripFont = 10, xAngle = 0, 
    hjust = 0.5, vjust = 0.5, xSize = 10, ySize = 10, xAxisSize = 10, yAxisSize = 10) + 
    xlab("Component 1") + ylab("Component 2")
compPlot
```

<img src="benchmarking_enrichmentConnectivity_files/figure-html/unnamed-chunk-28-2.png" width="100%" />

```r
pdf(paste0(WhereAmI, "results/Figures/multiOmicsPanels_compPlot_colon.pdf"), 
    width = 14, height = 3.5)
compPlot
dev.off()
```

```
## quartz_off_screen 
##                 2
```

