data("pxgenes")

OUTPUTS = c('layout', 'stress', 'spcorr')
COLNAMES = c("x", "y", "radius", "size", "significance", 'closeness', "cluster")

test_that("Layout under default parametrization", {
  l = gsoap_layout(pxgenes, 'Members', 'p.value')
  expect_equal(class(l), 'data.frame')
  expect_equal(colnames(l), COLNAMES)
  expect_equal(nrow(l), nrow(pxgenes))
})

test_that("Plot layout", {
  l = gsoap_layout(pxgenes, 'Members', 'p.value')
  p = gsoap_plot(l, as.color = 'cluster', as.alpha = 'closeness')
  expect_equal(class(p), c('gg', 'ggplot'))
})
