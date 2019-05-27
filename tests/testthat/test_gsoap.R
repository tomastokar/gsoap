data("pxgenes")

OUTPUTS = c('layout', 'stress', 'spcorr')
COLNAMES = c("x", "y", "radius", "size", "Weight", "Centrality", "Cluster")

test_that("Layout under default parametrization", {
  l = create_gsoap_layout(pxgenes, 'Members', 'p.value')
  expect_equal(class(l), 'list')
  expect_equal(names(l), OUTPUTS)
  expect_equal(colnames(l$layout), COLNAMES)
  expect_equal(nrow(l$layout), nrow(pxgenes))
})

test_that("Plot layout", {
  l = create_gsoap_layout(pxgenes, 'Members', 'p.value')
  p = plot_gsoap(l$layout, as.color = 'Cluster', as.alpha = 'Centrality')
  expect_equal(class(p), c('gg', 'ggplot'))
})
