# gsoap

A package for visualisation of Gene Set Over-representation Analysis results.

## Features


## Installation
```S
require(devtools)
install_github("tomastokar/gsoap", dependencies=T)
```

## Usage

### Load GSOAP package
```S
library(gsoap)
```

### Load example dataset 
```S
data("pxgenes")
```
Dataset if R data.frame with columns:
  * Source
  * Pathway
  * p.value
  * FDR
  * Members
  
### Create GSOAP layout
```S
gsoap_layout(pxgenes, 'Members', 'p.value')
```

### Plot GSOAP layout
```S
gsoap_plot(l, as.color = 'Cluster', as.alpha = 'Centrality')
```
