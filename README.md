# gsoap

A package for visualisation of gene set over-representation enrichment analysis.

## Features
<p align="justify">
Per dafault, <code>gsoap_layout</code> will calculate Jaccard distance between instances (e.g. pathways, GO terms, etc.), i.e. will measure relative overlaps between their query genes. Multidimensional scaling (other options include tSNE, CCA, Isomap) is then applied to project instances into 2-dimensional space. Circle packing is then applied to increase visual clarity of the layout. Significance of the instance is calculated as *-log10(pvalue)* and is used as instance weight to calculate instance closeness. Finally, cluster analysis is performed to identify clusters of instances.
</p>

Obtained layout is an R data.frame that contains the following columns:
  * x, y - coordinates of the circles representing pathways
  * radius - circle radius
  * size - number of pathway gene members (effect size)
  * significance - pathway significance defined as -log10(p.value)
  * closeness - pathway closeness
  * cluster - pathway cluster membership

<p align="justify">
Layout can be then visualized by <code>gsoap_plot</code>, that allows user to select columns to be used as color and alpha aesthetics, as well as the indices or names of the instances which should be annotated by labels.
</p>

## Installation
```S
require(devtools)
install_github("tomastokar/gsoap", dependencies=T)
```

## Usage
### Example 1

<p align="justify">
Here We use example dataset provided by GSOAP. The example dataset contains results of the over-representation analysis of 72 differentially expressed genes from [Tokar et al. 2018]. The analysis was performed using Pathway Data Integration Portal (pathDIP) [Rahmati et al., 2016]. 
</p>

#### Load GSOAP package
```S
library(gsoap)
```

#### Load example dataset 
```S
data("pxgenes")
```
It is an R data.frame, whose rownames are pathway names and columns are:
  * Source - original source of the pathway
  * Pathway - pathway name
  * p.value - statistical significance of the obtained overlap
  * FDR - false discovery rate of the obtained enrichment
  * Members - list of query genes belonging to the given pathway (in the following format: ``GENE1/GENE2/..'')
  
#### Create GSOAP layout
Before creating the GSOAP layout, reduce the number of examples to 100. GSOAP works best with datasets containing up to 100 instances (rows). For bigger datasets, value of 'scale.factor' should be decreased.

```S
# Reduce to top 100 instances
pxgenes = head(pxgenes[order(pxgenes$FDR),], 100)

# Create layout using default parametrization
layout = gsoap_layout(pxgenes, 'Members', 'p.value')
```

#### Create GSOAP plot
Create GSOAP plot, using color to highlight the cluster membership and opacity to highlight significance. In addition, add labels of the 
5 most significant instances.
```S
# Order instances by their significance
layout = layout[order(layout$significance, decreasing = TRUE),]

# Create gsoap plot
gsoap_plot(layout, as.color = 'cluster', as.alpha = 'significance', which.label = 1:5)
```

![gsoap_example](https://user-images.githubusercontent.com/46754141/59848847-292e6680-9334-11e9-884e-1b5180fb9aa9.png)

### Example 2

Here we use the R package *clusterProfiler* [Yu et al., 2012] to perform GSOA on example set of genes. GSOAP is then applied on the obtained results.

#### Load libraries
```S
# Load clusterProfiler and annotation database
library(clusterProfiler)
library(org.Hs.eg.db)
```

#### Import example genes
```S
# Import example gene list
data(geneList)
```
Object `geneList` is an example vector of gene expression fold change, whose names are gene entrez IDs, provided by clusterProfiler. For more details see *clusterProfiler* project ![website]('https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html').

#### Perform GSOA
Reduce the genes to a differentially expressed subset; and perform GSOA on Gene Ontology (GO) biological processes (BP)
```S
# Reducte to differentially expressed genes
gene = names(geneList)[abs(geneList) > 2.0]

# Run enrichGO
x = enrichGO(gene = gene,
             ont  = "BP",
             OrgDb = org.Hs.eg.db,
             pvalueCutoff = 0.01,
             pAdjustMethod = "fdr",
             universe = names(geneList),
             minGSSize = 5,
             maxGSSize = 500,
             qvalueCutoff = 0.01,
             readable = FALSE)

# Convert to data frame and set BP description as rownames
x = as.data.frame(x, row.names = x$Description)
```
Obtained data.frame contains:
  * ID - GO biol. proc. ID
  * Description - a full name of the GO biol. proc.
  * GeneRatio - a fraction of query gene members of the given GO biol. proc. to the total number of query genes. 
  * BgRatio - a fraction of total number of genes of the given GO biol. to size of the ``gene universe''.
  * pvalue - a statistical significance of the obtained overlap
  * p.adjust - a statistical significance adjusted for multiple testing applying user-specified adjustment method.
  * qvalue - a statistical significance adjusted for multiple testing applying clusterProfiler custom adjustment method.
  * geneID -  list of query gene members of the given GO biol. proc. (in the following format: ``GENE1/GENE2/..'')
  * Count - total number of query gene members of the given GO biol. proc. (effect size).

#### Create GSOAP layout and plot
Create GSOAP layout using tSNE projection; decreate the scale factor to better accomodate the instances into the layout space. 

```S
l = gsoap_layout(x,
                 genes = 'geneID',
                 pvalues = 'p.adjust',
                 projection = 'tsne',
                 scale.factor = 0.8,
                 no.clusters = 5)
```

Take indices of all instances of the first cluster.
```S
idx = which(l$cluster == 'Cluster 1')
```

Plot GSOAP, whiile using color to highlight the cluster membership and opacity to highlight significance. In addition, add labels of instances from the first cluster.
```S
p = gsoap_plot(l,
               as.alpha = 'significance',
               as.color = 'cluster',
               which.labels = idx,
               viridis.option = 'plasma',
               viridis.direction = 1,
               viridis.range = c(.2, .8),
               size.guide.loc = c(1., 1.),
               label.fontsize = 10)
```

![gsoap_example_cluster_profiler](https://user-images.githubusercontent.com/46754141/59890165-a93ce680-939d-11e9-9d91-d244453c3e1f.png)


## References
 * Tokar, Tomas, et al. "Differentially expressed microRNAs in lung adenocarcinoma invert effects of copy number aberrations of prognostic genes." Oncotarget 9.10 (2018): 9137.
 * Rahmati, Sara, et al. "pathDIP: an annotated resource for known and predicted human gene-pathway associations and pathway enrichment analysis." Nucleic acids research 45.D1 (2016): D419-D426.
 * Yu, Guangchuang, et al. "clusterProfiler: an R package for comparing biological themes among gene clusters." Omics: a journal of integrative biology 16.5 (2012): 284-287.
