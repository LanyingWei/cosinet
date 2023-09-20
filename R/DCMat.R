#' Differential Co-expression
#'
#' Calculate global differential co-expression network.
#'
#' This function calculates a global differential co-expression network between
#' two conditions using a z-score-based method. This method calculates z-scores
#' through the application of Fisher transformation to the correlation
#' coefficients and a subsequent z-test.
#'
#' @param exprMat an expression matrix. Rows represent genes and
#' columns represent samples.
#' @param condition a vector that specifies the condition of each sample.
#' The vector should contain only two conditions. The order of the vector
#' should match the column order of the expression matrix.
#' @param corMethod a character string indicating the correlation coefficient
#' to use. Can be either "spearman" (default) or "pearson".
#'
#' @return a matrix of z-scores representing the global differential
#' co-expression network between the two input conditions.
#'
#' @examples
#' # obtain expression data
#' data(exprBC)
#' # obtain sample information
#' data(samplesBC)
#'
#' # calculate differential correlation matrix
#' dcMat <- DCMat(exprMat = exprBC,
#'                condition = samplesBC$er_status,
#'                corMethod = "spearman")
#'
#' @author Lanying Wei
#'
#' @seealso \code{\link{DCPlot}}
#' @seealso \code{\link{eigenCentrality}}
#'
#' @export

DCMat <- function(exprMat,
                  condition,
                  corMethod=c("spearman", "pearson")) {

    stopifnot(!missing(exprMat),
              !missing(condition),
              length(condition) == ncol(exprMat))

    corMethod <- match.arg(corMethod)
    condition <- factor(condition)

    if (length(levels(condition)) != 2) {
        stop("There should be two levels of condition.")
    }

    message("Calculating differential coexpression matrix for ",
            levels(condition)[2],
            " (n=",
            sum(condition %in% levels(condition)[2]),
            ") ",
            "versus ", levels(condition)[1],
            " (n=",
            sum(condition %in% levels(condition)[1]),
            ") ...")

    exprMat1 <- exprMat[, condition %in% levels(condition)[1], drop = FALSE]
    exprMat2 <- exprMat[, condition %in% levels(condition)[2], drop = FALSE]

    const <- 1
    if (corMethod %in% 'spearman') {
        const <- 1.06
        message("Use Spearman correlation.")
    } else {
        message("Use Pearson correlation.")
    }

    message("Calculating correlation matrix for the first condition...")
    r1 <- getCor(t(exprMat1), corMethod)
    message("Calculating correlation matrix for the second condition...")
    r2 <- getCor(t(exprMat2), corMethod)
    message("Fisher transformation...")
    z1 <- atanh(r1)
    z2 <- atanh(r2)
    z <- (z2-z1)/sqrt(const/(sum(condition %in% levels(condition)[2])-3)+
                          const/(sum(condition %in% levels(condition)[1])-3))
    z[is.nan(z)] <- 0
    z[is.infinite(z)] <- 0
    diag(z) <- 0

    message("Done.")

    return(z)
}

getCor <- function(exprMat,
                   corMethod=c("spearman", "pearson")) {

    corMethod <- match.arg(corMethod)

    if (corMethod %in% "spearman") {
        exprMat <- apply(exprMat, 2, rank)
    }

    design <- matrix(data = 1,
                     nrow = nrow(exprMat),
                     ncol = 1)
    QR <- qr(design)
    E <- qr.qty(QR, exprMat)
    s2 <- colMeans(E[-1, ]^2)
    U <- t(t(E[-1, ])/sqrt(s2))
    c <- crossprod(U)/nrow(U)
    diag(c) <- 1

    return(c)
}


