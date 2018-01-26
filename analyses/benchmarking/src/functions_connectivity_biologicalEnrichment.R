#####################################################
#
# functions_connectivity_biologicalEnrichment.R
#
#####################################################

library(org.Hs.eg.db)
library(igraph)
library(statnet)
library(RColorBrewer)

networkStats = function(adjMat, mode = "graph", maxlen = 6){
  
  nrelations <- network::network(adjMat,directed=FALSE)
  #gplot(nrelations, usearrows=FALSE)
  
  # Centrality
  ## Degree
  nodeDegree <- degree(nrelations)
  # Betweenness 
  bet <- betweenness(nrelations, gmode="graph") # Geographic betweenness
  # Closeness
  clo <- closeness(nrelations) # Geographic closeness
  # Closeness2
  closeness2 <- function(x){ # Create an alternate closeness function!
    geo <- 1/geodist(x)$gdist # Get the matrix of 1/geodesic distance
    diag(geo) <- 0 # Define self-ties as 0
    apply(geo, 1, sum) # Return sum(1/geodist) for each vertex
  }
  clo2 <- closeness2(nrelations) # Use our new function on contiguity data
  
  # Geodesic distance
  gdist <- geodist(nrelations)$gdist # matrix of geodesic distances
  
  # Eigenvector centrality score
  eigenCentrality <- evcent(nrelations)
  
  # Harary graph centrality    
  # The Harary graph centrality of a vertex v is equal to 1/(max_u d(v,u)), where d(v,u) 
  # is the geodesic distance from v to u. Vertices with low graph centrality scores are 
  # likely to be near the “edge” of a graph, while those with high scores are likely to 
  # be near the “middle.” Compare this with closeness, which is based on the reciprocal 
  # of the sum of distances to all other vertices (rather than simply the maximum).
  hararyCentrality <- graphcent(nrelations)
  
  # Prestige is the name collectively given to a range of centrality 
  # scores which focus on the extent to which one is nominated by 
  # others. default = domain: indegree within the reachability graph (Lin's unweighted measure)
  prestigeScore <- prestige(nrelations, cmode="domain")
  
  # Stress Centrality
  # The stress of a vertex, v, is given by
  # C_S(v) = sum( g_ivj, i,j: i!=j,i!=v,j!=v)
  # where g_ijk is the number of geodesics from i to k through j. 
  # Conceptually, high-stress vertices lie on a large number of 
  # shortest paths between other vertices; they can thus be thought 
  # of as “bridges” or “boundary spanners.” Compare this with 
  # betweenness, which weights shortest paths by the inverse of 
  # their redundancy.
  stressScore <- stresscent(nrelations)     #Compute stress scores
  
  # Graph-level indicies
  graphDensity <- gden(nrelations) # graph density
  ## Reciprocity
  reciprocityMeasure <- c("dyadic", "dyadic.nonnull", "edgewise", "edgewise.lrr", "correlation")
  recipScore <- unlist(lapply(reciprocityMeasure, function(i){
    grecip(nrelations, measure = i)
  }))
  names(recipScore) <- reciprocityMeasure
  
  transitivityMeasure <- c("weak", "strong", "weakcensus", "strongcensus", "rank", "correlation")
  transScore <- unlist(lapply(transitivityMeasure, function(i){
    gtrans(nrelations, mode = mode, measure = i)
  }))
  names(transScore) <- transitivityMeasure
  
  # Census
  dyadCensus <- dyad.census(nrelations) # M,A,N counts
  triadCensus <- triad.census(nrelations, mode = mode) # Directed triad census
  
  # Compute k-path or k-cycle census statistics
  #pathCensus <- kpath.census(nrelations, mode = mode, maxlen=maxlen, tabulate.by.vertex=FALSE) # Count paths of length <=6
  #cycleCensus <- kcycle.census(nrelations, mode = mode, maxlen=maxlen, tabulate.by.vertex=FALSE) # Count cycles of length <=6
  
  # Cycle Census information
  # A (maximal) clique is a maximal set of mutually adjacenct 
  # vertices. Cliques are important for their role as cohesive 
  # subgroups
  cliques <- clique.census(nrelations, mode = mode, tabulate.by.vertex=FALSE, enumerate=FALSE)$clique.count # Find maximal cliques
  
  # isolate: An isolated vertex is a vertex with 
  # degree zero; that is, a vertex that is not an 
  # endpoint of any edge
  isol <- isolates(nrelations) # Get the entire list of isolates
  isolatedVertices <- network.vertex.names(nrelations)[isol]
  
  # Compure k-Core
  kc <- kcores(nrelations, mode = mode, cmode="indegree")
  
  nodeIndices <- rbind(nodeDegree, bet, clo, clo2, eigenCentrality, hararyCentrality, prestigeScore, stressScore, kc)
  
  return(list(nodeIndices = nodeIndices, gdist = gdist, graphDensity = graphDensity, recipScore = recipScore, 
    transScore = transScore, dyadCensus = dyadCensus, triadCensus = triadCensus, cliques = cliques, isolatedVertices = isolatedVertices))
}

manyPanels = function(X.train, Y.train, concat_alpha, concat_lambda, single_alphaList, single_lambdaList){
  # Concatenation-Enet
  combinedDat <- do.call(cbind, X.train)
  result <- enet(X = combinedDat, Y = Y.train, alpha=concat_alpha, family="multinomial", lambda = concat_lambda, X.test = NULL, 
    Y.test = NULL, filter = "none", topranked = 50)   # run once to determine lambda then set it
  concat_enetPanel <- lapply(X.train, function(i){ intersect(colnames(i), result$enet.panel)})
  
  # Ensemble-Enet
  ensemble_enetPanel <- mapply(function(x, y, z){
    result <- enet(X = x, Y = Y.train, alpha=y, family="multinomial", lambda = z, X.test = NULL, 
      Y.test = NULL, filter = "none", topranked = 50)
    result$enet.panel
  }, x = X.train, y = single_alphaList, z = single_lambdaList)
  
  # DIABLO-Enet
  design <- matrix(1, nrow = length(X.train), ncol = length(X.train))
  diag(design) <- 0
  ncomp <- nlevels(Y.train) - 1
  list.keepX = lapply(ensemble_enetPanel, function(i){
    rep(round(length(i)/ncomp, 0), ncomp)
  })
  diabloMod = block.splsda(X = X.train, Y = Y.train, 
    ncomp = ncomp, keepX = list.keepX, design = design,
    scheme = "centroid")
  
  diabloPanel <- lapply(diabloMod$loadings[-(length(X.train)+1)], function(x)
    unique(as.character(as.matrix(apply(x, 2, function(i) names(i)[which(i != 0)])))))
  
  ## Single-omics (univariate tests)
  design <- model.matrix(~Y.train)
  singleFeat <- mapply(function(x, y){
    fit <- eBayes(lmFit(t(x), design, method = "ls"))
    top <- topTable(fit, coef = 1:ncol(design), adjust.method = "BH", n = nrow(fit))
    rownames(top)[1:y]
  }, x = X.train, y = lapply(ensemble_enetPanel, length))
  
  ## Panels
  panels = list(Concatenation = concat_enetPanel, Ensemble = ensemble_enetPanel,
    DIABLO = diabloPanel, Single = singleFeat)
  
  ## Panel length
  panelLength <- data.frame(Concatenation = unlist(lapply(concat_enetPanel, length)),
    Ensemble = unlist(lapply(ensemble_enetPanel, length)),
    DIABLO = unlist(lapply(diabloPanel, length)),
    Single = unlist(lapply(singleFeat, length)))
  
  return(list(panels = panels, panelLength = panelLength, concat_lambda = result$lambda))
}

## Turn panels into graph indices
graphIndices = function(panels = panels, X.train = X.train, cutoff = cutoff){
  
  ## Determine adjacency matrices
  adjMat <- lapply(panels, function(i){
    adjMat = cor(do.call(cbind, mapply(function(x, y){
      y[, x]
    }, x = i, y = X.train)))
    adjMat[abs(adjMat) < cutoff] <- 0
    adjMat[abs(adjMat) > cutoff] <- 1
    adjMat
  })
  
  # Estimate graph statistics
  graphs <- lapply(adjMat, function(x) networkStats(adjMat = x))
  
  ## graphDensity, recipScore, transScore, cliques
  gIndices <- lapply(graphs, function(i){
    data.frame(Value = c(i$graphDensity, i$recipScore["edgewise.lrr"], i$transScore["weakcensus"], length(i$cliques)),
      Statistic = c("Graph Density", "Edgewise.lrr", "Transitivity", "NumOfCliques"))
    
  })
  gIndicesDat <- do.call(rbind, gIndices) %>% mutate(Method = rep(names(gIndices), each = 4))
  
  # Dyads and Triads
  dyads <- do.call(rbind, lapply(graphs, function(x){x$dyadCensus}))
  dyads <- as.data.frame(dyads)[,-2] %>% mutate(Method = names(graphs)) %>% 
    gather(Type, Number, -Method)
  triads <- do.call(rbind, lapply(graphs, function(x){x$triadCensus}))
  triads <- as.data.frame(triads) %>% mutate(Method = names(graphs)) %>% 
    gather(Type, Number, -Method)
  
  # Number of isolates
  isolatedFeat <- lapply(graphs, function(x){ x$isolatedVertices})
  isolates <- do.call(rbind, lapply(isolatedFeat, function(x){
    unlist(lapply(X.train, function(y){
      length(intersect(x, colnames(y)))
    }))
  })) %>% as.data.frame %>% mutate(Method = names(graphs)) %>% 
    gather(Dataset, NumOfIsolates, -Method)
  
  return(list(adjMat = adjMat, gIndicesDat = gIndicesDat, dyads = dyads, triads = triads, isolates = isolates))
}


generate_network = function(data, panel, cutoff){
  combinedDat <- do.call(cbind, data)
  ## generate network object
  adjMat0 = cor(combinedDat[, panel])
  adjMat <- adjMat0
  adjMat[abs(adjMat) < cutoff] <- 0
  adjMat[abs(adjMat) > cutoff] <- 1
  netStat <- networkStats(adjMat = adjMat, mode = "graph", maxlen = 6)
  
  ## node colors
  dark = brewer.pal(n = length(data), name = "Set1")
  nodeColList = mapply(function(x, y){
    rep(y, length(intersect(colnames(x), panel)))
  }, x = data, y = names(data))
  nodeCol = dark[as.numeric(factor(unlist(nodeColList)))]
  
  ## plot network
  nrelations <- network::network(adjMat, directed=FALSE)
  gplot(nrelations, usearrows=FALSE, vertex.col = nodeCol)
  legend("topleft", names(data), col = dark, pch = 19, bty = "n")
}



## Bimap interface:
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

#################################################################################
## Calculate enrichment on multiple gene sets                                  ##
#################################################################################
gsa = function(weights, ID, bhfdr){
  ## pathways
  library("dcGSA"); library(org.Hs.eg.db)
  id <- "entrez"
  DATA_DIR <- "~/Dropbox/Methods/kleinstein-logminer-d0763790de98/Data"
  KEGG.file <- file.path(DATA_DIR, paste0("c2.cp.kegg.v4.0.", id, ".gmt"))
  Reactome.file <- file.path(DATA_DIR, paste0("c2.cp.reactome.v5.1.", id, ".gmt"))
  BTM.file <- file.path(DATA_DIR, paste0("BTM_for_GSEA_", id, "_20131008.gmt"))
  GO.file <- file.path(DATA_DIR, paste0("c5.all.v4.0.", id, ".gmt"))
  CELLS.file <- file.path(DATA_DIR, paste0("AbbasCellSubsets_edit20141007_", id, ".gmt"))
  
  geneSetLibraries <- c("KEGG", "Reactome", "GO", "BTM", "CELLS")
  results <- vector("list", length(geneSetLibraries))
  names(results) <- geneSetLibraries
  for(geneSetLibrary in geneSetLibraries){
    gsList <- do.call(readGMT, list(file=as.name(paste0(geneSetLibrary, ".file"))))
    gsList <- lapply(gsList, function(i){
      as.character(unlist(xx[i]))
    })
    genesets <- lapply(gsList, function(i){
      which(ID %in% i)
    })
    
    pval <- lapply(genesets, function(i){
      limma::geneSetTest(i, statistics = abs(weights),
        alternative="up", type="t",
        ranks.only = TRUE)
    })
    adjPval <- p.adjust(pval, "BH")
    enrichment_result <- data.frame(Pathway = names(pval), P.Value = unlist(pval), adj.P.Val = adjPval)
    results[[geneSetLibrary]] <- enrichment_result[order(enrichment_result$P.Value), ]
  }
  unlist(lapply(results, function(i){
    sum(i$adj.P.Val < bhfdr)
  }))
}


### sear
library(dplyr);
library(sear)

sear2 = function (input, type = c("mrna", "mirna", "both")){
  library(sear)
  data("collections", envir = environment())
  collections <- collections %>% 
    rowwise() %>% 
    mutate(members_mrna.mirna = list(unique(c(members_mrna, members_mirna))))
  
  type <- match.arg(type)
  tbl <- switch(type, mrna = dplyr::select(collections, collection:geneset, 
    members = members_mrna), mirna = dplyr::select(collections, 
      collection:geneset, members = members_mirna), both = dplyr::select(collections, 
        collection:geneset, members = members_mrna.mirna))
  uni <- tbl$members %>% unlist() %>% unique()
  recognized <- input[input %in% uni]
  if (length(recognized) < 10) {
    warning(sprintf("Submitted %s symbols, but only %s are recognized.", 
      length(input), length(recognized)))
  }
  
  input <- recognized
  tbl %>% dplyr::rowwise(.) %>% dplyr::mutate(n_input = length(input), 
    n_geneset = length(members), intersect = length(intersect(input, 
      members)), p_value = phyper(intersect - 1, n_geneset, 
        length(uni) - n_geneset, n_input, lower.tail = F)) %>% 
    dplyr::ungroup(.) %>% dplyr::group_by(collection) %>% 
    dplyr::mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
    dplyr::ungroup(.)
}

