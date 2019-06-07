data("pxgenes")

OUTPUTS = c('layout', 'stress', 'spcorr')
COLNAMES = c("x", "y", "radius", "size", "Weight", "Closeness", "Cluster", "Intracluster_closeness")

test_that("Layout under default parametrization", {
  l = create_gsoap_layout(pxgenes, 'Members', 'p.value')
  expect_equal(class(l), 'data.frame')
  expect_equal(colnames(l), COLNAMES)
  expect_equal(nrow(l), nrow(pxgenes))
})

test_that("Plot layout", {
  l = create_gsoap_layout(pxgenes, 'Members', 'p.value')
  p = plot_gsoap(l, as.color = 'Cluster', as.alpha = 'Centrality')
  expect_equal(class(p), c('gg', 'ggplot'))
})
