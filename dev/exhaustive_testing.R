library(gsoap)
data(pxgenes)

DISTANCE_METHODs = c('jaccard',
                     'cosine',
                     'dice',
                     'intersection')

PROJECTION_METHODs = c('iso',
                       'mds',
                       'cca',
                       'tsne')

PACKING_OPTIONs = c(TRUE, FALSE)
CENTRALITY_OPTIONs = c(TRUE, FALSE)
CLUSTERING_OPTIONs = c(TRUE, FALSE)

CLUSTER_STATS_OPTIONs = c('PBC',
                          'HGSD',
                          'ASW',
                          'ASWw',
                          'CH',
                          'R2',
                          'CHsq',
                          'R2sq',
                          'HC')

N = 100
all_args = list()
for (i in 1:N){
  args = list('distance.method' = sample(DISTANCE_METHODs, 1),
              'projection.method' = sample(PROJECTION_METHODs, 1),
              'weighted' = sample(c(TRUE, FALSE), 1),
              'log10.weights' = sample(c(TRUE, FALSE), 1),
              'do.packing' = sample(c(TRUE, FALSE), 1),
              'calc.centrality' = sample(c(TRUE, FALSE), 1),
              'do.clustering' = sample(c(TRUE, FALSE), 1),
              'cluster.stat' = sample(CLUSTER_STATS_OPTIONs, 1))
  all_args[[i]] = args
}

for (i in 1:100){
  args = all_args[[i]]
  print(i)
  print(args)
  l = create_gsoap_layout(pxgenes,
                          'Members',
                          'p.value',
                          distance.method = args$distance.method,
                          projection.method = args$projection.method,
                          weighted = args$weighted,
                          log10.weights = args$log10.weights,
                          do.packing = args$do.packing,
                          calc.centrality = args$calc.centrality,
                          do.clustering = args$do.clustering,
                          cluster.stat = args$cluster.stat)
}

