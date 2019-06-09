# G

A package for visualisation of gene set over-representation enrichment analysis.

## Features

Per dafault, `gsoap_layout` will calculate Jaccard distance between instances (e.g. pathways, GO terms, etc.), i.e. will measure relative overlaps between their query genes. Multidimensional scaling (other options include tSNE, CCA, Isomap) is then applied to project instances into 2-dimensional space. Circle packing is then applied to increase visual clarity of the layout. Importance of the instance is calculated as -log10(pvalue) and is used as a instance weight to calculate instance closeness. Finally, cluster analysis is performed to identify clusters of instances.

Obtained layout is an R data.frame containing following columns:
  * x, y - coordinates of the circles representing pathways
  * radius - circle radius
  * size - number of pathway gene members (effect size)
  * importance - pathway importance defined as -log10(p.value)
  * closeness - pathway closeness
  * cluster - pathway cluster
  * intracluster_closeness - pathway intra-cluster closeness

Layout can be then visualized by `gsoap_plot`, that allows user to select columns to be used as color and alpha aesthetics, as well as the indices or names of the instances which should be annotated by labels.

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

### Create GSOAP plot
```S
gsoap_plot(layout, as.color = 'cluster', as.alpha = 'importance')
```

## References
 * Tokar, Tomas, et al. "Differentially expressed microRNAs in lung adenocarcinoma invert effects of copy number aberrations of prognostic genes." Oncotarget 9.10 (2018): 9137.
 * Rahmati, Sara, et al. "pathDIP: an annotated resource for known and predicted human gene-pathway associations and pathway enrichment analysis." Nucleic acids research 45.D1 (2016): D419-D426.
