#' @export gsoap_layout
#' @importFrom tsne tsne
#' @importFrom ProjectionBasedClustering Isomap SammonsMapping tSNE CCA KruskalStress
#' @importFrom philentropy distance
#' @importFrom packcircles circleRepelLayout
#' @importFrom WeightedCluster wcKMedRange wcKMedoids

create_association_matrix = function(l){
  # Extract instances and members
  instances = rep(names(l), times = sapply(l, length))
  members = unlist(l)
  # Create long format adjacency mat.
  am = data.frame(instances, members, adjacency = 1)
  # Reshape from long to wide format
  am = reshape(am, idvar = 'instances', timevar = 'members', direction = 'wide')
  am = am[,-1]
  # Set colnames and row names
  colnames(am) = unique(members)
  rownames(am) = unique(instances)
  # Replace NAs by zero
  am[is.na(am)] = 0
  # Get rid of the first column
  am = am[,-1]
  # Convert to matrix
  am = as.matrix(am)
  # Return adjacency matrix
  return(am)
}


calc_distance_matrix = function(m, distance.method = 'jaccard'){
  dist.mat = suppressMessages(philentropy::distance(m, distance.method))
  rownames(dist.mat) = rownames(m)
  colnames(dist.mat) = rownames(m)
  return(dist.mat)
}

non.diag = function(m){
  nd = upper.tri(m, diag = FALSE) | lower.tri(m, diag = FALSE)
  return(nd)
}

resolve.nondiag.zeros = function(M, k){
  n = nrow(M)
  m = ncol(M)
  nondiag.zeros = (M == 0 & non.diag(M))
  d = 1 - (k - 1)/k
  r = matrix(d, n, m)
  e = rnorm(n * m, 0., d/2.5)
  e = matrix(e, n, m)
  M[nondiag.zeros] = r[nondiag.zeros] + e[nondiag.zeros]
  return(M)
}

isomap_transformation = function(d, isomap.k = 3){
  res = ProjectionBasedClustering::Isomap(d,
                                          k = isomap.k,
                                          OutputDimension = 2,
                                          PlotIt = FALSE)
  return(res)
}


sammons_tranformation = function(d){
  res = ProjectionBasedClustering::SammonsMapping(d,
                                                  OutputDimension = 2,
                                                  PlotIt = FALSE)
  return(res)
}


tsne_transformation = function(d, tsne.perplexity = 30, tsne.iterations = 1e+3){
  res = ProjectionBasedClustering::tSNE(d,
                                        k = tsne.perplexity,
                                        OutputDimension = 2,
                                        Whitening = FALSE,
                                        PlotIt = FALSE,
                                        Iterations = tsne.iterations)
  return(res)
}


cca_transformation = function(d, cca.epochs = 10, cca.alpha0 = 0.5){
  res = ProjectionBasedClustering::CCA(d,
                                       Epochs = cca.epochs,
                                       OutputDimension = 2,
                                       alpha0 = cca.alpha0,
                                       PlotIt = FALSE)
  return(res)
}


min_max_scale = function(x)(x - min(x))/(max(x) - min(x))

create_layout = function(x, size, scale.factor = 1.0){

  # Normalize coordinates to 0-1 interval
  x = apply(x, 2, min_max_scale)

  # Rescale size by mean
  size.normed = (size / max(size)) / nrow(x) * scale.factor[1]

  # Calc raw radius
  radius = sqrt(size.normed / pi)

  # Create layout
  y = cbind(x, radius)

  return(y)
}


packing_simple = function(x, packing.maxiter = 1e+6){

  # Get range of values
  xrange = range(x[,1])
  yrange = range(x[,2])

  # Circle packing
  circles = packcircles::circleRepelLayout(x,
                                           xrange,
                                           yrange,
                                           xysizecols = 1:3,
                                           sizetype = "radius",
                                           maxiter = packing.maxiter,
                                           wrap = FALSE)

  # Take layout
  layout = circles$layout

  return(layout)
}


calc_closeness = function(dm, w){
  closeness = 1.
  if (class(dm) != 'matrix'){
    dm = as.matrix(dm)
  }
  if (all(dim(dm) > 1)){
    closeness = 1. / apply(dm , 1, weighted.mean, w)
  }
  return(closeness)
}


intracluster_closeness = function(cl, dm, w){
  icc = rep(1, length(cl))
  for (i in unique(cl)){
    idx = which(cl == i)
    icc[idx] = calc_closeness(dm[idx, idx], w[idx])
  }
  return(icc)
}


annotate_clusters = function(cl, icc, names){
  cn = names
  for (i in unique(cl)){
    idx = which(cl == i)
    n = names[idx]
    x = icc[idx]
    cn[idx] = n[which.max(x)]
  }
  return(cn)
}


select_clustering = function(cls, dm, w, stat = 'meta'){
  cq = lapply(cls, function(cl)WeightedCluster::wcClusterQuality(dm, cl, weights = w)$stats)
  cq = data.frame(do.call(rbind, cq))
  if (stat == 'meta'){
    cq = apply(cq, 2, rank)/ncol(cq)
    score = exp(rowMeans(log(cq)))
  } else {
    score = cq[,stat]
  }
  cl = cls[[which.max(score)]]
  return(cl)
}


hkclustering = function(dm, w, no.clusters = NULL, max.clusters = 5, hc.method = 'ward.D', cluster.stat = 'meta'){
  hc = hclust(dist(dm), method = hc.method, members = w)
  if (is.null(no.clusters)){
    max.clusters = min(max.clusters, nrow(dm))
    clustering_list = lapply(2:max.clusters, function(k)cutree(hc, k))
    clustering = select_clustering(clustering_list, dm, w, stat = cluster.stat)
  } else {
    no.clusters = min(no.clusters, nrow(dm))
    clustering = cutree(hc, no.clusters)
  }
  return(clustering)
}


#' A function to create a layout for GSOAP plot
#'
#' Some description
#'
#' More description
#'
#' @param x a data frame with the results of gene set over-representation analysis.
#' Must have a rownames indicating names of the genes sets and at least two columns,
#' one of which contains query gene members; and another one that contains respective
#' p-values (raw or adjusted for multiple testing).
#'
#' @param genes a character or integer, indicating name or index of the column containing genes members.
#' @param pvalues a character or integer, indicating name or index of the column containing p-values.
#' @param splitter a character to be used as a delimiter to parse the genes column.
#' @param distance a character indicating method used to calculate the distance/dissimilarity between instances.
#' Options include (but are not limited to) \emph{jaccard} (default), \emph{manhattan}, \emph{dice},
#' \emph{pearson}. For more details see \code{\link[philentropy]{distance}}.
#' @param projection a character indicating method used to project instances into 2-dimensional space based on their distance/dissimilarity..
#' Ooptions include \emph{iso} (isomap; default), \emph{mds} (multidimensional scaling), \emph{cca} (curvilinear component analysis), \emph{tsne} (t-distributed stochastic neighbor embedding),
#' @param scale.factor a positive real number to control dependence of the circle radius on the number of query gene members of the given gene set.
#' @param weighted a boolean indicating whether to use pathway \emph{significance}
#' (-log10(pvalue)) as a weight when closeness and clustering are calculated.
#' @param packing a boolean indicating whether to apply circle packing.
#' @param clustering a boolean indicating whether to apply clustering.
#' @param hc.method a character indicating method of hierarchical cluster to be used.
#' Options include: \emph{ward.D} (default), \emph{ward.D2}, \emph{single}, \emph{complete},
#' \emph{average}, \emph{mcquitty}, \emph{median} and \emph{centroid}.
#' @param no.clusters an integer indicating number of clusters, must be less than number of gene sets (rows).
#' @param max.clusters an integer indicating maximum number of clusters to consider, must be at least two.
#' @param cluster.stat an indicating statistic used to select optimal number of clusters.
#' Options are:
#' \itemize{
#'     \item \emph{meta} (default) is a combination of the methods listed below
#'     \item \emph{PBC} (point biserial correlation; default)
#'     \item \emph{HG} (Hubert's samma)
#'     \item \emph{HGSD} (Hubert’s samma - Somer's D)
#'     \item \emph{ASW} (average silhouette width)
#'     \item \emph{ASWw} (Average weighted silhouette width)
#'     \item \emph{CH} (Calinski-Harabasz index)
#'     \item \emph{R2} (R-squared)
#'     \item \emph{CHsq} (Calinski-Harabasz index using squared distances)
#'     \item \emph{R2sq} (R-squared using squared distances)
#'     \item \emph{HC} (Hubert’s C coefficient)
#' }
#'
#' @param isomap.k an integer indicating number of k nearest neighbors of
#' the \emph{isomap} projection.
#' @param tsne.perplexity an integer indicating \emph{tSNE} perplexity.
#' @param tsne.iterations an integer indicating maximum number of \emph{tSNE} iterations to perform.
#' @param cca.epochs an integer indicating \emph{CCA} training length.
#' @param cca.alpha0 a positive real number indicating \emph{CCA} initial step size.
#'
#' @return \code{layout} a data frame with x and y coordinates of
#'     the points representing the insttances, their size (radius) derived from
#'     the number of gene members; significance (-log10(p-value)), closeness,
#'     cluster membership and intracluster closeness.
#'
#' @author Tomas Tokar <tomastokar@gmail.com>
#'
#' @examples
#' data(pxgenes)
#'
#' l = gsoap_layout(pxgenes, 'Members', 'p.value')
#'
#' head(l$layout)
gsoap_layout = function(x,
                        genes,
                        pvalues,
                        splitter = '/',
                        distance = 'jaccard',
                        projection = 'iso',
                        scale.factor = 0.8,
                        weighted = TRUE,
                        packing = TRUE,
                        clustering = TRUE,
                        hc.method = 'ward.D',
                        no.clusters = NULL,
                        max.clusters = 8,
                        cluster.stat = 'meta',
                        isomap.k = 3,
                        tsne.perplexity = 30,
                        tsne.iterations = 1e+3,
                        cca.epochs = 10,
                        cca.alpha0 = 0.5){
  # -------------
  # Check inputs
  # -------------
  if (missing(x)){
    stop('Input is missing')
  }
  if (!is.data.frame(x)){
    stop('Input is not data frame')
  }
  if (!is.character(rownames(x))){
    stop('Input has missing or improper rownames')
  }
  if(!((genes %in% colnames(x))|(genes %in% 1:ncol(x)))){
    stop('Wrong `genes` value')
  }
  if(!((pvalues %in% colnames(x))|(pvalues %in% 1:ncol(x)))){
    stop('Wrong `pvalues` value')
  }
  if (!any(grepl(splitter, x[,genes]))){
    warning('Either `genes`, or `splitter` seem to be not correct.')
  }

  # --------------
  # Create layout
  # --------------
  # Extract query genes -- instances memberships
  memberships.list = setNames(strsplit(x[,genes], splitter), rownames(x))
  # Create association matrix
  asc.mat = create_association_matrix(memberships.list)
  # Get number of member genes
  no.members = rowSums(asc.mat)
  # Calculate distance matrix
  dist.mat = calc_distance_matrix(asc.mat, distance.method = distance)
  # --------------------------
  # Do projection to 2d space
  # --------------------------
  # Check for zeros appart of the main diagonal
  if (any(rowSums(dist.mat == 0.) > 1)){
    warning("Zero dissimilarity between non-identical entries.")
    k = ncol(asc.mat)
    dist.mat = resolve.nondiag.zeros(dist.mat, k)
  }
  if (projection == 'iso'){
    proj = suppressMessages(isomap_transformation(dist.mat,
                                                  isomap.k = isomap.k))
  }
  if (projection == 'mds'){
    #res = mds_transformation(d)
    proj = suppressMessages(sammons_tranformation(dist.mat))
  }
  if (projection == 'cca'){
    proj = suppressMessages(cca_transformation(dist.mat,
                                               cca.epochs,
                                               cca.alpha0))
  }
  if (projection == 'tsne'){
    proj = suppressMessages(tsne_transformation(dist.mat,
                                                tsne.perplexity,
                                                tsne.iterations))
  }
  # Do 2d projection
  xy = proj$ProjectedPoints
  # Calculate circle radius
  layout = create_layout(xy, no.members, scale.factor = scale.factor)
  # Circle packing
  if (packing){
    layout = packing_simple(layout)
  } else {
    layout = data.frame(layout)
  }
  # Set  colnames
  layout = setNames(layout, c('x', 'y', 'radius'))
  # Set rownames
  rownames(layout) = rownames(x)
  # Calculate number of members
  layout$size = no.members
  # Calculate significance
  layout$significance = -log10(x[,pvalues])
  # Set weights
  weights = rep(1, nrow(layout))
  if (weighted){
    weights = layout$significance
  }
  # Calculate closeness and add to layout
  layout$closeness = calc_closeness(dist.mat, weights)

  # ---------------------
  # Calculate distortion
  # ---------------------
  # Calculate euclidean distance between instances after projection
  dx = as.matrix(suppressMessages(philentropy::distance(layout[,1:2], method = 'euclidean')))
  # Calculate Kruskal stress after projection and packing
  stress = ProjectionBasedClustering::KruskalStress(dist.mat, dx)
  message(paste('Kruskall stress :', sprintf('%1.3f', stress)))
  # Calculate spearman correlation
  spcorr = cor(c(dist.mat), c(dx), method = 'spearman')
  message(paste('Rank correlation :', sprintf('%1.3f', spcorr)))

  # -----------------------
  # Extended functionality
  # -----------------------
  # Do clustering
  if (clustering){
    # Clustering
    layout$cluster = hkclustering(dist.mat,
                                  weights,
                                  no.clusters = no.clusters,
                                  max.clusters = max.clusters,
                                  hc.method = hc.method,
                                  cluster.stat = cluster.stat)
    # Calculate intracluster weights
    layout$intracluster_closeness = intracluster_closeness(layout$cluster,
                                                           dist.mat,
                                                           weights)
    # Add cluster names
    layout$cluster = annotate_clusters(layout$cluster,
                                       layout$intracluster_closeness,
                                       rownames(layout))
  }
  # Return
  return(layout)
}
