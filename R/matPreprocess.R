#' Expression Matrix Preprocessing
#'
#' Preprocess expression matrix for co-expression analysis.
#'
#' Preprocess expression matrix for co-expression analysis. Perform
#' log-transformation, remove rows with duplicate gene names
#' and filter genes expressed in too few samples.
#'
#' @param exprMat an expression matrix. Rows represent genes and
#' columns represent samples.
#' @param condition a vector that specifies the condition of each sample.
#' The vector should contain only two conditions. The order of the vector
#' should match the column order of the expression matrix.
#' @param minSample minimum number of samples required to have an expression
#' value over minExp for each condition. Genes that fail this criterion are
#' removed.
#' @param minExp minimum expression value for gene filtration.
#' @param logTrans logical indicator for performing log-transformation.
#' @param logBase the base with respect to which logarithms are computed.
#' @param pseudoCount pseudo count to add before log-transformation.
#'
#' @return a matrix of prepossessed gene expression values.
#'
#' @examples
#'
#' testMat <- matrix(data = 0:11, nrow = 3)
#' matPreprocess(exprMat = testMat,
#'               condition = c(1, 1, 2, 2),
#'               minSample = 2, minExp = 1, logTrans = TRUE,
#'               logBase = 2, pseudoCount = 0.1)
#'
#' @author Lanying Wei
#'
#' @export

matPreprocess <- function(exprMat,
                          condition,
                          minSample=20,
                          minExp=0,
                          logTrans=FALSE,
                          logBase=2,
                          pseudoCount=1) {

    stopifnot(!missing(exprMat),
              !missing(condition),
              length(condition) == ncol(exprMat))

    condition <- factor(condition)

    if (length(levels(condition)) != 2) {
        stop("There should be two levels of condition.")
    }

    exprMat <- matGeneFilter(exprMat = exprMat,
                             condition = condition,
                             minSample = minSample,
                             minExp = minExp)

    if (logTrans) {
        exprMat <- matLogTrans(exprMat = exprMat,
                               base = logBase,
                               pseudoCount = pseudoCount)
    }

    message("Done.")
    return(exprMat)
}


matLogTrans <- function(exprMat,
                        base=2,
                        pseudoCount=1) {

    exprMat <- log(exprMat + pseudoCount, base = base)
    message("Log", base,
            " transformed with pseudo count of ",
            pseudoCount, ".")
    return(exprMat)
}

matGeneFilter <- function(exprMat,
                          condition,
                          minSample=20,
                          minExp=0) {
    condition <- factor(condition)
    genesIdx2Keep <-
        rowSums(exprMat[, condition %in% levels(condition)[1]] > minExp,
                na.rm = TRUE) >= minSample &
        rowSums(exprMat[, condition %in% levels(condition)[2]] > minExp,
                na.rm = TRUE) >= minSample
    message(nrow(exprMat) - sum(genesIdx2Keep),
            " genes with expression value larger than ", minExp,
            " in less than ", minSample,
            " samples for any condition were removed.")
    exprMat <- exprMat[genesIdx2Keep, ]
    dupGenes <- rownames(exprMat)[duplicated(rownames(exprMat))]
    dupGeneIdx <- rownames(exprMat) %in% dupGenes
    if (length(dupGeneIdx) > 0) {
        exprMat <- exprMat[!dupGeneIdx, , drop = FALSE]
    }
    message(sum(dupGeneIdx),
            " rows with duplicate gene names were removed.")
    return(exprMat)
}


