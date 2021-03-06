# Load GSOAP package
library(gsoap)

# Load example dataset
data(pxgenes)

# Reduce to top 100 instances
pxgenes = head(pxgenes[order(pxgenes$FDR),], 100)

# Create GSOAP layout
layout = gsoap_layout(pxgenes, 'Members', 'p.value')

# Order instances by their significance
layout = layout[order(layout$significance, decreasing = TRUE),]

# Create GSOAP plot
p = gsoap_plot(layout, as.color = 'cluster', as.alpha = 'significance', which.label = 1:5)

# Plot to file
png('./img/gsoap_example.png', width = 7, height = 5, units = 'in', res = 900)
plot(p)
dev.off()
