# Load GSOAP package
library(gsoap)

# Load example dataset
data(pxgenes)

# Create GSOAP layout
layout = gsoap_layout(pxgenes, 'Members', 'p.value')

# Order instances by their significance
layout = layout[order(layout$significance, decreasing = TRUE),]

# Create GSOAP plot
p = gsoap_plot(layout, as.color = 'cluster', as.alpha = 'significance', which.label = 1:5)

pdf('./img/gsoap_example.pdf', width = 7, height = 5)
plot(p)
dev.off()

png('./img/gsoap_example.png', width = 7, height = 5, units = 'in', res = 900)
plot(p)
dev.off()
