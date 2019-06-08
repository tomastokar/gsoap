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

HC_METHODS = c('ward.D',
               'ward.D2',
               'single',
               'average',
               'median',
               'centroid')

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
              'hc.method' = sample(HC_METHODS, 1),
              'cluster.stat' = sample(CLUSTER_STATS_OPTIONs, 1, prob = c(10, rep(1, 9))))
  color = ifelse(args[['clustering']],
                 'Cluster',
                 ifelse(args[['closeness']],
                        'Closeness',
                        'p.value'))
  alpha = ifelse(args[['clustering']],
                 'Intracluster_closeness',
                 'p.value')
  viridis = sample(c('viridis', 'magma', 'plasma', 'inferno', 'cividis'), 1)
  args[['color']] = color
  args[['alpha']] = alpha
  args[['viridis']] = viridis
  args[['which.labels']] = sample(nrow(pxgenes), 5)
  all_args[[i]] = args
}

pdf('./img/exhtest.pdf')
for (i in 1:N){
  args = all_args[[i]]
  print(i)
  print(args)
  l = gsoap_layout(pxgenes,
                   'Members',
                   'p.value',
                   distance = args$distance,
                   projection = args$projection,
                   weighted = args$weighted,
                   log10.weights = args$log10.weights,
                   packing = args$packing,
                   closeness = args$closeness,
                   clustering = args$clustering,
                   hc.method = args$hc.method,
                   cluster.stat = args$cluster.stat)

  l$p.value = -log10(pxgenes$p.value)
  print(head(l))
  p = gsoap_plot(l,
                 as.color = args$color,
                 as.alpha = args$alpha,
                 show.size.guide = sample_bool(),
                 show.color.guide = sample_bool(),
                 show.alpha.guide = sample_bool(),
                 viridis.option = args$viridis,
                 which.labels = args$which.labels)
  plot(p)
}
dev.off()
