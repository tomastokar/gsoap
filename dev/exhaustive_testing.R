library(gsoap)
data(pxgenes)

DISTANCE_METHODs = c('jaccard',
                     'manhattan',
                     'dice',
                     'pearson')

PROJECTION_METHODs = c('iso',
                       'mds',
                       'cca',
                       'tsne')

PACKING_OPTIONs = c(TRUE, FALSE)
CLOSENESS_OPTIONs = c(TRUE, FALSE)
CLUSTERING_OPTIONs = c(TRUE, FALSE)

CLUSTER_STATS_OPTIONs = c('meta',
                          'PBC',
                          'HGSD',
                          'ASW',
                          'ASWw',
                          'CH',
                          'R2',
                          'CHsq',
                          'R2sq',
                          'HC')

sample_bool = function(){
  sample(c(TRUE, FALSE), 1)
}

N = 100
all_args = list()
for (i in 1:N){
  args = list('distance' = sample(DISTANCE_METHODs, 1),
              'projection' = sample(PROJECTION_METHODs, 1),
              'weighted' = sample_bool(),
              'log10.weights' = sample_bool(),
              'packing' = sample_bool(),
              'closeness' = sample_bool(),
              'clustering' = sample_bool(),
              'cluster.stat' = sample(CLUSTER_STATS_OPTIONs, 1, prob = c(10, rep(1, 9))))
  all_args[[i]] = args
}

for (i in 1:100){
  args = all_args[[i]]
  print(i)
  print(args)
  l = create_gsoap_layout(pxgenes,
                          'Members',
                          'p.value',
                          distance = args$distance,
                          projection = args$projection,
                          weighted = args$weighted,
                          log10.weights = args$log10.weights,
                          packing = args$packing,
                          closeness = args$closeness,
                          clustering = args$clustering,
                          cluster.stat = args$cluster.stat)
}

