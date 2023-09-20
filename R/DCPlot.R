#' Differential Co-expression Network Plotting
#'
#' Plot differential co-expression network for a given set of genes.
#'
#' This function plots correlation networks for a given set of genes separately
#' for two conditions. Edge weights represent the correlation coefficient
#' between two genes. The users can set various thresholds and colors to focus
#' on key changes in the network and understand the change in associations
#' between gene pairs under different conditions.
#'
#'
#' @param exprMat an expression matrix. Rows represent genes and
#' columns represent samples.
#' @param condition a vector that specifies the condition of each sample.
#' The vector should contain only two conditions. The order of the vector
#' should match the column order of the expression matrix.
#' @param genes a character vector of genes to plot.
#' @param dcMat a matrix of z-scores representing differential associations
#' of gene expression between two conditions. See {\link{DCMat}} for more
#' details. This matrix, along with posDCCut and negDCCut, are used to remove
#' edges between gene pairs that show low level of differential co-expression.
#' @param corMethod a character string indicating the correlation coefficient
#' to use. Can be either "spearman" (default) or "pearson".
#' @param posDCCut a positive numeric value that serves as positive
#' threshold for DCMat to hide unimportant edges. Gene pairs with DCMat values
#' below posDCCut and above negDCCut are hidden.
#' @param negDCCut a negative numeric value that serves as negative
#' threshold for DCMat to hide unimportant edges. Gene pairs with DCMat values
#' below posDCCut and above negDCCut are hidden.
#' @param minimum a numeric value between 0 and 1. Edges with absolute
#' weights less than this value are not shown.
#' @param maximum a numeric value between 0 and 1. The highest weight to scale
#' the edge widths to.
#' @param cut a numeric value between 0 and 1. Cut the scaling of edges in
#' width and color saturation. Edges with absolute weights over this value
#' have the strongest color intensity and become wider the stronger they are,
#' and edges with absolute weights under this value have the smallest width
#' and become vaguer the weaker the weight. If this is set to 0, no cutoff is
#' used and all edges vary in width and color.
#' @param posCol color of positive edges. Can be a vector of two to indicate
#' color of edges under 'cut' value and color of edges over 'cut' value.
#' @param negCol color of negative edges. Can be a vector of two to indicate
#' color of edges under 'cut' value and color of edges over 'cut' value.
#' @param layout the layout of the graph, can be "spring" or "circle". "spring"
#' gives a force embedded layout; "circle" places all nodes in a single
#' circle.
#' @param repulsion Scalar on the default repulse radius in the spring layout.
#' Defaults to 1. Setting this argument to lower values (e.g., 0.5) will cause
#' nodes in the spring layout to repulse each other less. This is especially
#' useful if a few unconnected nodes cause the giant component to visually
#' be clustered too much in the same place.
#' @param labelCex numeric scalar on the label size.
#' @param labelCol Character containing the color of the labels.
#' @param geneList1 a character vector of a list of genes to be colored.
#' Genes in this vector are labelled with color 'labelCol1'.
#' @param geneList2 a character vector of a list of genes to be colored.
#' Genes in this vector are labelled with color 'labelCol2'.
#' @param labelCol1 color of genes in geneList1.
#' @param labelCol2 color of genes in geneList2.
#' @param cond1Title plot title for the first condition.
#' @param cond2Title plot title for the second condition.
#' @param ... other parameters to be passed to \code{\link{qgraph}}.
#'
#' @return a list of two "qgraph" objects, each corresponds to one of the
#' correlation plots of the two input conditions. See \code{\link{qgraph}}
#' for more details.
#'
#' @examples
#'
#' # get expression data
#' data(exprBC)
#' # get sample information
#' data(samplesBC)
#' # get example genes
#' data(genesER)
#'
#' # calculate differential correlation matrix
#' dcMat <- DCMat(exprMat = exprBC,
#'                condition = samplesBC$er_status,
#'                corMethod = "spearman")
#'
#' # plot differential co-expression network using example genes
#' p <- DCPlot(exprMat = exprBC,
#'             condition = samplesBC$er_status,
#'             genes = genesER,
#'             dcMat = dcMat,
#'             posDCCut = 5,
#'             negDCCut = -5,
#'             corMethod = "spearman",
#'             minimum = 0.2)
#'
#' @author Lanying Wei
#'
#' @seealso \code{\link{DCMat}}
#' @seealso \code{\link{eigenCentrality}}
#'
#' @importFrom qgraph qgraph averageLayout
#' @importFrom graphics par
#'
#' @export


DCPlot <- function(exprMat,
                   condition,
                   genes,
                   dcMat,
                   corMethod=c("spearman", "pearson"),
                   posDCCut=3,
                   negDCCut=-posDCCut,
                   minimum=0.2,
                   maximum=1,
                   cut=0.8,
                   posCol=c("#4361ee", "#3f37c9"),
                   negCol=c("#ff5a5f", "#d90429"),
                   layout=c("spring", "circle"),
                   repulsion=1,
                   labelCex=1,
                   labelCol="#6c757d",
                   geneList1=NULL,
                   geneList2=NULL,
                   labelCol1="#e63946",
                   labelCol2="#072ac8",
                   cond1Title=NULL,
                   cond2Title=NULL,
                   ...) {

    stopifnot(!missing(exprMat),
              !missing(condition),
              !missing(genes),
              !missing(dcMat),
              length(condition) == ncol(exprMat),
              all(colnames(dcMat) == rownames(dcMat)))

    condition <- factor(condition)
    if (length(levels(condition)) != 2) {
        stop("There should be two levels of condition.")
    }
    if (!all(genes %in% rownames(exprMat))) {
        stop("Gene ",
             paste(genes[!genes %in% rownames(exprMat)], collapse = " "),
             " not in exprMat.")
    }
    if (!all(genes %in% rownames(dcMat))) {
        stop("Gene ",
             paste(genes[!genes %in% rownames(dcMat)], collapse = " "),
             " not in dcMat." )
    }

    corMethod <- match.arg(corMethod)
    layout <- match.arg(layout)

    corMat1 <-
        getCor(t(exprMat[genes, condition %in% levels(condition)[1]]),
               corMethod)
    corMat2 <-
        getCor(t(exprMat[genes, condition %in% levels(condition)[2]]),
               corMethod)

    corMat1[which(dcMat[genes, genes] > negDCCut &
                      dcMat[genes, genes] < posDCCut, arr.ind = TRUE)] <- 0
    corMat2[which(dcMat[genes, genes] > negDCCut &
                      dcMat[genes, genes] < posDCCut, arr.ind = TRUE)] <- 0


    corMat1[abs(corMat1) < minimum] <- 0
    corMat2[abs(corMat2) < minimum] <- 0
    genes2remove <- names(which(colSums(corMat1) == 0 & colSums(corMat2 == 0)))

    corMat1 <- corMat1[!rownames(corMat1) %in% genes2remove,
                       !colnames(corMat1) %in% genes2remove]

    corMat2 <- corMat2[!rownames(corMat2) %in% genes2remove,
                       !colnames(corMat2) %in% genes2remove]

    genes <- colnames(corMat1)

    corMat1 <-
        getCor(t(exprMat[genes, condition %in% levels(condition)[1]]),
               corMethod)
    corMat2 <-
        getCor(t(exprMat[genes, condition %in% levels(condition)[2]]),
               corMethod)

    corMat1[which(dcMat[genes, genes] > negDCCut &
                      dcMat[genes, genes] < posDCCut, arr.ind = TRUE)] <- 0
    corMat2[which(dcMat[genes, genes] > negDCCut &
                      dcMat[genes, genes] < posDCCut, arr.ind = TRUE)] <- 0


    L <- averageLayout(corMat1, corMat2,
                       repulsion = repulsion, layout = layout)

    labelColor <- rep(labelCol, ncol(corMat1))
    labelColor[colnames(corMat1) %in% geneList1] <- labelCol1
    labelColor[colnames(corMat1) %in% geneList2] <- labelCol2


    if (is.null(cond1Title)) {
        cond1Title <- levels(condition)[1]
    }

    if (is.null(cond2Title)) {
        cond2Title <- levels(condition)[2]
    }
    oldMfrow <- par()$mfrow
    par(mfrow = c(1, 2))

    p1 <- qgraph(corMat1, layout = L,
                 labels = colnames(corMat1),
                 label.cex = labelCex,
                 label.scale = FALSE,
                 minimum = minimum,
                 maximum = maximum,
                 cut = cut,
                 vsize = 0,
                 label.color = labelColor,
                 posCol = posCol,
                 negCol = negCol,
                 title = cond1Title,
                 ...)


    p2 <- qgraph(corMat2, layout = L,
                 labels = colnames(corMat1),
                 label.cex = labelCex,
                 label.scale = FALSE,
                 minimum = minimum,
                 maximum = maximum,
                 cut = cut,
                 vsize = 0,
                 label.color = labelColor,
                 posCol = posCol,
                 negCol = negCol,
                 title = cond2Title,
                 ...)

    par(mfrow = oldMfrow)

    res <- list(p1, p2)
    names(res) <- c(cond1Title, cond2Title)
    return(res)
}

