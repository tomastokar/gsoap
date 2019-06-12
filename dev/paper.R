library(gsoap)
library(ggpubr)

# # Load example data
# data("pxgenes")
#
# # Subset the example data
# N = 100
# pxgenes = head(pxgenes[order(pxgenes$FDR),], N)
#
# # -----------------
# # Generate layouts
# # -----------------
# projections = c('iso', 'mds', 'cca', 'tsne')
# layouts = list()
# for (projection in projections){
#   layouts[[projection]] = gsoap_layout(pxgenes,
#                                        'Members',
#                                        'p.value',
#                                        scale.factor = 0.9,
#                                        projection = projection,
#                                        packing = T,
#                                        no.clusters = 4)
# }
#
# # Add ANGPT1 logical index to layout
# # and order by sigificance
# for (projection in projections){
#   layouts[[projection]]$ANGPT1 = factor(grepl('ANGPT1', pxgenes$Members))
#   layouts[[projection]] = layouts[[projection]][order(layouts[[projection]]$significance, decreasing = F),]
# }
#
# save(layouts, file = './img/paper_layout.RData')

# -------------
# Create plots
# -------------
load('./img/paper_layout.RData')
f = function(p, a, L = 3){
  x = layouts[[p]]
  x = x[order(x['significance'], decreasing = T),]
  gsoap_plot(x,
             as.color = a,
             xlabel = paste(p, '1'),
             ylabel = paste(p, '2'),
             size.guide.loc = c(0., 1.),
             base.fontsize = 14,
             size.guide.fontsize = 9,
             size.guide.no = 3,
             viridis.direction = -1,
             which.labels = 1:L,
             label.fontsize = 12,
             label.alpha = 1.)
}

g = function(p, a1, a2, L = 3){
  x = layouts[[p]]
  x = x[order(x['significance'], decreasing = T),]
  gsoap_plot(x,
             as.color = a1,
             as.alpha = a2,
             xlabel = paste(p, '1'),
             ylabel = paste(p, '2'),
             size.guide.loc = c(0., 1.),
             base.fontsize = 14,
             size.guide.fontsize = 9,
             size.guide.no = 3,
             viridis.direction = -1,
             which.labels = 1:L,
             label.fontsize = 12,
             label.alpha = 1.)
}

aesthetics = c('significance', 'closeness', 'cluster', 'ANGPT1')
plots = list()
for (aesthetic in aesthetics){
  if (aesthetic %in% c('cluster', 'ANGPT1')){
    l = setNames(lapply(projections, g, aesthetic, 'significance'), projections)
  } else {
    l = setNames(lapply(projections, f, aesthetic), projections)
  }
  plots[[aesthetic]] = l
}

# --------------
# Arrange plots
# --------------
n = length(projections)
arranged = list()
for (aesthetic in aesthetics){
  arranged[[aesthetic]] = ggarrange(plotlist = plots[[aesthetic]],
                                    ncol = n,
                                    nrow = 1,
                                    common.legend = TRUE,
                                    legend = 'top')
}

final_plot = ggarrange(plotlist = arranged,
                       nrow = length(aesthetics),
                       ncol = 1,
                       align = 'hv')

pdf('./img/paper_plot.pdf', width = 15, height = 16)
plot(final_plot)
dev.off()
