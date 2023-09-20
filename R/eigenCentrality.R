#' Eigenvector Centrality
#'
#' Calculates eigenvector centrality using power iteration method.
#'
#' Perform power iteration to calculate eigenvector centrality. Eigenvector
#' centrality is used to measure the importance of a node by considering the
#' centrality of its neighbors.Nodes (genes) with high eigenvector centralities
#' are those that are highly connected to many other nodes which are highly
#' connected to many others (and so on). Using the differential co-expression
#' z-score matrix as input, the nodes with high eigenvector centralities
#' represent the genes that are located at the hubs of highly rewired
#' sub-networks.
#'
#' @param mat a symmetric matrix with zero values on the diagonal.
#' @param tolerance the error tolerance used to check convergence in power
#' iteration. The algorithm stops when the total absolute difference between
#' two consecutive iterations is less than the tolerance or when maxIter
#' iterations have been reached.
#' @param maxIter the maximum number of iterations allowed for the calculation.
#' If the tolerance is not reached before this limit, the algorithm stops.
#' @param sorted a logical indicator of whether the centrality values should
#' be returned in decreasing order (TRUE) or not (FALSE).
#'
#' @return a vector of eigenvector centrality values for each node.
#'
#' @examples
#' # get expression data
#' data(exprBC)
#' # get sample information
#' data(samplesBC)
#'
#' # calculate differential correlation matrix
#' dcMat <- DCMat(exprMat = exprBC,
#'                condition = samplesBC$er_status,
#'                corMethod = "spearman")
#'
#' # calculate gene centrality base on differential correlation matrix
#' geneCentrality <- eigenCentrality(mat = dcMat)
#'
#' @author Lanying Wei
#'
#' @seealso \code{\link{DCMat}}
#' @seealso \code{\link{rankPlot}}
#'
#' @export
#'

eigenCentrality <- function(mat,
                            tolerance=1e-6,
                            maxIter=500,
                            sorted=TRUE) {

    stopifnot(!missing(mat))

    if (!all(rownames(mat) == colnames(mat))) {
        stop("Symmetric matrix required.")
    }
    if (!all(diag(mat) == 0)) {
        stop("Diagonal values should all be zero.")
    }

    mat <- abs(mat)
    message("Use the absolute value of the matrix to proceed.")
    n <- ncol(mat)
    x0 <- rep(0, n)
    x1 <- rep(1 / n, n)
    iter <- 0

    while (sum(abs(x0 - x1)) > tolerance && iter < maxIter) {
        if (iter != 0 && iter %% 50 == 0) {
            message("Number of current iterations: ", iter)
        }
        x0 <- x1
        x1 <- x1 %*% mat
        m <- max(x1, na.rm = TRUE)
        x1 <- x1 / m
        iter <- iter + 1
    }

    message("Total number of iterations: ", iter)
    message("Total differences between the last two iterations: ",
            sum(abs(x0 - x1)))

    names(x1) <- colnames(mat)
    if (sorted) {
        x1 <- sort(x1, decreasing = TRUE)
    }

    return(x1)
}






