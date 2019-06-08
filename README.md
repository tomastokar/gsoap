# gsoap
=======

A package for visualisation of Gene Set Over-representation Analysis results.

## Features


## Installation
```S
require(devtools)
install_github("tomastokar/gsoap", dependencies=T)
```

## Usage

### Load GSOAP package
```
library(gsoap)
```

### Load example dataset 
```
data("pxgenes")
```
Dataset if R data.frame with columns:
  * Source
  * Pathway
  * p.value
  * FDR
  * Members
  
### Create GSOAP layout
```
gsoap_layout(pxgenes, 'Members', 'p.value')
```

### Plot GSOAP layout
```
gsoap_plot(l, as.color = 'Cluster', as.alpha = 'Centrality')
```
