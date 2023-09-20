#' Rank Plot
#'
#' Draw the score and rank of genes on a scatter plot.
#'
#' @param geneScore a numeric vector of scores to rank, with gene as names.
#' @param genesToLabel a character vector specifying the genes to be labeled.
#' @param nLabel an integer specifying the number of top scored genes to be
#' labeled.
#' @param nIQR a positive numeric value indicating the fold of IQR to be used as
#' cutoff to color genes. Genes with scores > Q3 + nIQR * (Q3 - Q1) are colored
#' red.
#' @param labelSize a positive numeric value specifying the size of the label.
#' @param main a character string representing the plot title.
#'
#' @return A \code{ggplot} scatter plot object showing the ranked gene scores.
#'
#' @author Lanying Wei
#'
#' @examples
#' # get expression data
#' data(exprBC)
#' # get sample information
#' data(samplesBC)
#'
#' # calculate differential correlation matrix between ER status
#' dcMat <- DCMat(exprMat = exprBC,
#'                condition = samplesBC$er_status,
#'                corMethod = "spearman")
#'
#' # calculate gene centrality base on differential correlation matrix
#' geneCentrality <- eigenCentrality(mat = dcMat)
#'
#' # plot gene centrality rank plot
#' rankPlot(geneCentrality)
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats quantile
#'
#' @seealso \code{\link{eigenCentrality}}
#'
#'
#' @export

rankPlot <- function(geneScore,
                     genesToLabel=NULL,
                     nLabel=10,
                     nIQR=1.5,
                     labelSize=3,
                     main=NULL) {

    stopifnot(!missing(geneScore))

    Q1 <- quantile(geneScore, 0.25)
    Q3 <- quantile(geneScore, 0.75)
    cutoff <- Q3 + nIQR * (Q3 - Q1)

    data <- data.frame(Gene = names(geneScore),
                       Score = geneScore)
    data$Rank <- rank(-data$Score)
    data$TopGenes <- "No"
    data$TopGenes[data$Score > cutoff] <- "Yes"

    myColor <- c("No" = "#c0c0c0",
                 "Yes" = "#e63946")

    idx <- (data$Rank <= nLabel) | (data$Gene %in% genesToLabel)
    p <- ggplot(data, aes_string(x = "Rank",
                                 y = "Score",
                                 color = "TopGenes")) +
        geom_point(size = 0.5)
    if(sum(idx) > 0) {
        p <- p + ggrepel::geom_text_repel(aes_string(label = "Gene"),
                                          data = data[idx,],
                                          size = labelSize,
                                          max.overlaps = Inf)
    }

    p <- p + scale_color_manual(values = myColor) +
        scale_fill_manual(values = myColor) +
        labs(x = "Rank",
             y = "Score",
             title = main) +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "none")
    return(p)
}
