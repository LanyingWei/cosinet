# Cosinet Package

Cosinet is an R package that aims to bridge the gap between differential
co-expression network analysis and precision medicine by accurately
quantifying the degree of rewiring of co-expression networks at
the level of individual samples. Differential co-expression network analysis
is a valuable approach for understanding the underlying mechanisms of disease
and informing the development of personalized therapies, but current
approaches typically measure an average differential network across multiple
samples, making it difficult to apply in the field of precision medicine.
Cosinet utilizes gene expression data to determine the degree of similarity
between the gene co-expression patterns of a given sample and reference
conditions within the context of a function-specific differential
co-expression network. To obtain differential co-expression network, 
we provide R function for z-score-based differential co-expression analysis. 
Next, we calculate eigenvector centrality to get a ranked gene list, 
which represents the degree of hubness of each gene in the network.
Based on the ranking, one can further perform gene set enrichment analysis
(GSEA) to get enriched pathways, and use the core enrichment genes of a relevant 
pathway to form a function-specific differential co-expression sub-network. 
Finally, Cosinet scores can be calculated based on the specific 
sub-network.


## Installation

To install the cosinet package, you can use the `remotes` package:

``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("LanyingWei/cosinet")
```

## Vignette

For more details, please refer to the 
[online vignette](https://htmlpreview.github.io/?https://github.com/LanyingWei/cosinet/blob/main/vignettes/cosinet.html).
