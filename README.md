# gsoap

A package for visualisation of gene set over-representation enrichment analysis results.

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
The example dataset contains results of the over-representation analysis of 72 differentially expressed genes from [Tokar et al. 2018]. The analysis was performed using Pathway Data Integration Portal (pathDIP) [Rahmati et al., 2016]. It is an R data.frame containing following columns:
  * Source - original source of the pathway
  * Pathway - pathway name
  * p.value - statistical significance of the obtained enrichment
  * FDR - false discovery rate of the obtained enrichment
  * Members - list of query genes belonging to the given pathway (in the following format: ``GENE1/GENE2/..'')
  
### Create GSOAP layout
```S
layout = gsoap_layout(pxgenes, 'Members', 'p.value')
```
Obtained layout is datta.frame containing following columns:
  * x - x-coordinate of the circle representing given pathway
  * y - y-coordinate of the circle representing given pathway
  * radius - circle radius
  * size - number of pathway gene members (effect size)
  * importance - pathway importance defined as -log10(p.value)
  * closeness - pathway closeness
  * cluster - pathway cluster
  * intracluster_closeness - pathway intra-cluster closeness

### Create GSOAP plot
```S
gsoap_plot(layout, as.color = 'cluster', as.alpha = 'importance')
```

