#' @details
#'
#' Cosinet is an R package that aims to bridge the gap between differential
#' co-expression network analysis and precision medicine by accurately
#' quantifying the degree of rewiring of co-expression networks at
#' the level of individual samples. Differential co-expression network analysis
#' is a valuable approach for understanding the underlying mechanisms of disease
#' and informing the development of personalized therapies, but current
#' approaches typically measure an average differential network across multiple
#' samples, making it difficult to apply in the field of precision medicine.
#' Cosinet utilizes gene expression data to determine the degree of similarity
#' between the gene co-expression patterns of a given sample and reference
#' conditions within the context of a function-specific differential
#' co-expression network. The package employs a combination of techniques
#' such as differential co-expression analysis, network centrality calculation,
#' gene set enrichment analysis (GSEA), and a novel statistic that measures the
#' differences in statistical independence of a gene pair between two
#' conditions for a single sample. The general workflow is as follows:
#'
#'
#' @section Construct global differential co-expression network:
#' \itemize{
#' \item{\code{\link{matPreprocess}}}: Preprocess expression matrix for
#'     co-expression analysis. There should be two conditions of samples
#'     for comparison.
#' \item{\code{\link{DCMat}}}: Calculate differential co-expression between
#' the two conditions using a z-score based method.}
#' @section Calculate node centralities in differential co-expression network:
#' \itemize{
#' \item{\code{\link{eigenCentrality}}}: Calculate eigenvector centrality
#' using differential co-expression matrix.
#' \item{\code{\link{rankPlot}}}: Visualize the score and rank of gene
#' centrality.
#' }
#' @section Identify function-specific differential co-expression sub-networks:
#' Then the resulting gene centrality scores can be used as input for GSEA
#' using for example \code{GSEA} function from \code{clusterProfiler}
#' package. The core-enrichment genes of a relevant gene set can be used to
#' compute Cosinet scores.
#' @section Compute Cosinet scores based on differential co-expression sub-network:
#' \itemize{
#' \item{\code{\link{DCPlot}}}: Plot differential co-expression network for
#' a given set of genes.
#' \item{\code{\link{getCosi}}}: Calculate Cosinet score for a given set of
#' genes.
#' }
#' The final Cosinet score for each sample represents a measure of how
#' closely the gene co-expression patterns of that sample match the
#' patterns observed in the reference conditions within the differential
#' sub-network. A lower Cosinet score indicates that the gene co-expression
#' patterns of the given sample are more similar to those of the first
#' condition, while a higher Cosinet score indicates that they are more
#' similar to those of the second condition. Downstream analysis can then
#' be performed to further investigate the differences.
#'
#' @author Lanying Wei
#'
#' @seealso \code{\link{matPreprocess}}
#' @seealso \code{\link{DCMat}}
#' @seealso \code{\link{eigenCentrality}}
#' @seealso \code{\link{rankPlot}}
#' @seealso \code{\link{DCPlot}}
#' @seealso \code{\link{getCosi}}
#'
"_PACKAGE"
