#' @export create_gsoap_layout
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


pamlustering = function(dm, w, no.clusters = NULL, max.clusters = 8, cluster.stat = 'PBC', boots = 100){
  if (is.null(no.clusters)){
    pam.stats = WeightedCluster::wcKMedRange(dm, 2:max.clusters, weights = w, R = boots)
    no.clusters = summary(pam.stats)[cluster.stat, 1]
  }
  clusters = WeightedCluster::wcKMedoids(dm, no.clusters, weights = w, cluster.only = T)
  clusters = factor(rownames(dm)[clusters])
  return(clusters)
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
#' @param distance.method a character indicating method used to calculate the distance/dissimilarity between instances.
#' Options include (but are not limited to) \emph{jaccard} (default), \emph{manhattan}, \emph{tanimoto},
#' \emph{intersection}. For more details see \code{\link[philentropy]{distance}}.
#' @param projection.method a character indicating method used to project instances into 2-dimensional space based on their distance/dissimilarity..
#' Ooptions include \emph{iso} (isomap; default), \emph{mds} (multidimensional scaling), \emph{cca} (curvilinear component analysis), \emph{tsne} (t-distributed stochastic neighbor embedding),
#' @param scale.factor a positive real number to control dependence of the circle radius on the number of query gene members of the given gene set.
#' @param weighted a boolean indicating whether to apply weights
#' (weight = \emph{1 / p-value}) when centrality and clustering are calculated.
#' @param log10.weights a boolean indicating whether weights should undergo
#' additional log10 tranformation
#' @param do.packing a boolean indicating whether to apply circle packing.
#' @param calc.centrality a boolean indicating whether to calculate centrality.
#' @param do.clustering a boolean indicating whether to apply clustering.
#' @param isomap.k an integer indicating number of k nearest neighbors of
#' the \emph{isomap} projection.
#' @param tsne.perplexity an integer indicating \emph{tSNE} perplexity.
#' @param tsne.iterations an integer indicating maximum number of \emph{tSNE} iterations to perform.
#' @param cca.epochs an integer indicating \emph{CCA} training length.
#' @param no.clusters an integer indicating number of clusters, must be less than number of gene sets (rows).
#' @param max.clusters an integer indicating maximum number of clusters to consider, must be at least two.
#' @param cluster.stat an indicating statistic used to select optimal number of clusters.
#' Options are:
#' \itemize{
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
#' @param pam.boots a positive integer indicating number of boostraps used to be used when selecting number of clusters.
#'
#' @return A \code{gsoap} object that is a list comprising following components
#' \itemize{
#'     \item \code{layout} data frame with x and y coordinates of
#'     the points representing the insttances, their size (radius) derived from
#'     the number of gene members; weight (-log10(p-value)), centrality
#'     and cluster membership.
#'     \item \code{stress} Kruskall stress caused by projection and circle packing.
#'     \item \code{spcorr} Spearman correlation between the original distances/dissimilarities
#'     between gene sets and those obtained after projection and circle packing.
#' }
#'
#' @author Tomas Tokar <tomastokar@gmail.com>
#'
#' @examples
#' data(pxgenes)
#'
#' l = create_gsoap_layout(pxgenes, 'Members', 'p.value')
#'
#' head(l$layout)
create_gsoap_layout = function(x,
                               genes,
                               pvalues,
                               splitter = '/',
                               distance.method = 'jaccard',
                               projection.method = 'iso',
                               scale.factor = 1.0,
                               weighted = TRUE,
                               log10.weights = TRUE,
                               do.packing = TRUE,
                               calc.centrality = TRUE,
                               do.clustering = TRUE,
                               isomap.k = 3,
                               tsne.perplexity = 30,
                               tsne.iterations = 1e+3,
                               cca.epochs = 10,
                               cca.alpha0 = 0.5,
                               no.clusters = NULL,
                               max.clusters = 8,
                               cluster.stat = 'PBC',
                               pam.boots = 100){ # add has genes
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
  dist.mat = calc_distance_matrix(asc.mat,
                                  distance.method = distance.method)
  # Check for zeros outside the diagonal
  if (any(rowSums(dist.mat == 0) > 1)){
    warning("Zero dissimilarity between non-identical entries.")
    if (projection.method == 'mds'){
      warning(paste('Projection method', dQuote(projection.method), 'cannot be used.'))
      projection.method = 'tsne'
      message(paste('Projection method set to', dQuote(projection.method)))
    }
    else if (projection.method == 'iso'){
      warning(paste('Projection method', dQuote(projection.method), 'cannot be used.'))
      projection.method = 'cca'
      message(paste('Projection method set to', dQuote(projection.method)))
    }

  }

  # --------------------------
  # Do projection to 2d space
  # --------------------------
  if (projection.method == 'iso'){
    proj = suppressMessages(isomap_transformation(dist.mat,
                                                  isomap.k = isomap.k))
  }
  if (projection.method == 'mds'){
    #res = mds_transformation(d)
    proj = suppressMessages(sammons_tranformation(dist.mat))
  }
  if (projection.method == 'cca'){
    proj = suppressMessages(cca_transformation(dist.mat,
                                               cca.epochs,
                                               cca.alpha0))
  }
  if (projection.method == 'tsne'){
    proj = suppressMessages(tsne_transformation(dist.mat,
                                                tsne.perplexity,
                                                tsne.iterations))
  }
  # Do 2d projection
  xy = proj$ProjectedPoints
  # Calculate circle radius
  layout = create_layout(xy, no.members, scale.factor = scale.factor)
  # Circle packing
  if (do.packing){
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

  # ---------------------
  # Calculate distortion
  # ---------------------
  # Calculate euclidean distance between instances after projection
  dx = as.matrix(suppressMessages(philentropy::distance(layout[,1:2], method = 'euclidean')))
  # Calculate Kruskal stress after projection and packing
  stress = ProjectionBasedClustering::KruskalStress(dist.mat, dx)
  # Calculate spearman correlation
  spcorr = cor(c(dist.mat), c(dx), method = 'spearman')

  # -----------------------
  # Extended functionality
  # -----------------------
  # Add p-values to layout
  # Set weights
  if (weighted){
    layout$Weight = 1. / x[,pvalues]
    if (log10.weights){
      layout$Weight = log10(layout$Weight)
    }
  } else {
    layout$Weight = rep(1, nrow(layout))
  }
  # Calculate centrality and add to layout
  if (calc.centrality){
    layout$Centrality = apply(dist.mat, 1, weighted.mean, layout$Weight)
    layout$Centrality = 1. - min_max_scale(layout$Centrality)
  }
  # Do clustering
  if (do.clustering){
    layout$Cluster = pamlustering(dist.mat,
                                  layout$Weight,
                                  no.clusters = no.clusters,
                                  max.clusters = max.clusters,
                                  cluster.stat = cluster.stat,
                                  boots = pam.boots)
  }
  # Wrap results together
  res = setNames(list(layout, stress, spcorr), c('layout','stress','spcorr'))
  # Return
  return(res)
}
