#' Cosinet Score
#'
#' Calculate Cosinet scores with a given set of genes.
#'
#' Calculate Cosinet scores with a given set of genes. Cosinet uses gene
#' expression data to determine the degree of similarity between the
#' gene co-expression patterns of a given sample and reference conditions
#' in the context of a function-specific differential co-expression network.
#' Cosinet scores quantifies the degree of rewiring of co-expression network at
#' the level of individual samples.
#'
#'
#' @param exprMat an expression matrix. Rows represent genes and
#' columns represent samples.
#' @param condition a vector that specifies the condition of each sample.
#' The vector should contain only two conditions. The order of the vector
#' should match the column order of the expression matrix.
#' @param genes a character vector of genes to be evaluated. These genes can be
#' the genes involved in a function-specific differential co-expression network.
#' Gene pairs generated from this gene list that pass the thresholds for
#' dcMat are used to calculate Cosinet scores.
#' @param sampleNames a character vector of sample names to calculate the
#' Cosinet scores for.
#' @param bin a numeric value between 0 and 1 that represents the fraction of
#' the bin size relative to the binRange length. Binned regions are used to
#' estimate gene pair dependency at a local region around an individual sample.
#' @param binRange a numeric value between 0 and 1 that represents the fraction
#' of the middle expression range relative to the full expression range. The
#' length of binRange is used as the base to generate the bin length. When
#' binRange is set to 1, the range between the maximum expression value and
#' the minimum expression value is used as the base. If users are concerned
#' about extreme data having a large impact on the bin size setting, a smaller
#' value of binRange, for example 0.95, can be set so that the expression range
#' between the 2.5% quantile and the 97.5% quantile is used as the base. A
#' binRange value less than 0.9 is not recommended, as this will exclude up to
#' 10% of the samples from the bin size calculation.
#' @param dcMat a matrix of z-scores representing differential associations
#' of gene expression between two conditions. See {\link{DCMat}} for more
#' details. This matrix, along with posDCCut and negDCCut, are used to 1. remove
#' gene pairs that show low level of differential co-expression, so that only
#' important changes in the network are considered in the Cosinet score
#' calculation; and 2. used as weights when calculating Cosinet scores.
#' @param posDCCut a positive numeric value that serves as positive
#' threshold for dcMat. Gene pairs with dcMat values
#' below posDCCut and above negDCCut are not used for the calculation.
#' @param negDCCut a negative numeric value that serves as negative
#' threshold for dcMat. Gene pairs with dcMat values
#' below posDCCut and above negDCCut are not used for the calculation.
#' @param weighted logical indicator for weighting each gene pair by the
#' absolute value of corresponding dcMat score during Cosinet score calculation.
#' @param threads number of threads used for calculation.
#' @return a matrix of the final Cosinet scores and the individual scores for
#' each gene pair before averaging. A lower Cosinet score indicates that the
#' gene co-expression patterns of the given sample are more similar to those of
#' the first condition, while a higher Cosinet score indicates that they are
#' more similar to those of the second condition.
#'
#'
#' @author Lanying Wei
#'
#' @seealso \code{\link{DCMat}}
#' @seealso \code{\link{eigenCentrality}}
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom stats p.adjust t.test
#' @importFrom utils combn
#'
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
#' # calculate Cosinet scores with example genes
#' scoreMat <- getCosi(exprMat = exprBC,
#'                     condition = samplesBC$er_status,
#'                     genes = genesER,
#'                     bin = 0.1,
#'                     dcMat = dcMat,
#'                     posDCCut = 3,
#'                     negDCCut = -3,
#'                     threads = 2)
#'
#'
#' @export

getCosi <- function(exprMat,
                    condition,
                    genes,
                    sampleNames=colnames(exprMat),
                    bin=0.1,
                    binRange=1,
                    dcMat=NULL,
                    posDCCut=3,
                    negDCCut=-posDCCut,
                    weighted=TRUE,
                    threads=2) {

    stopifnot(!missing(exprMat),
              !missing(condition),
              !missing(genes),
              length(condition) == ncol(exprMat),
              bin > 0 & bin < 1)

    condition <- factor(condition)

    if (length(levels(condition)) != 2) {
        stop("There should be two levels of condition.")
    }

    if (!all(genes %in% rownames(exprMat))) {
        stop("Gene ",
             paste(genes[!genes %in% rownames(exprMat)], collapse = " "),
             " not in exprMat.")
    }


    if (!all(sampleNames %in% colnames(exprMat))) {
        stop("Sample ",
             paste(sampleNames[!sampleNames %in% colnames(exprMat)],
                   collapse = " "),
             " not in exprMat.")
    }

    condition <- factor(condition)
    genePairs <- as.data.frame(t(combn(genes, 2)))
    weights <- apply(genePairs, 1, function(x) dcMat[x[1], x[2]])

    if (!is.null(dcMat) & !is.null(posDCCut) & !is.null(negDCCut)) {
        if (any(weights > posDCCut | weights < negDCCut)) {
            nGenePairsAll <- nrow(genePairs)
            genePairs <- genePairs[weights > posDCCut | weights < negDCCut, ,
                                   drop = FALSE]
            weights <- weights[weights > posDCCut | weights < negDCCut]
            message(nGenePairsAll - nrow(genePairs),
                    " gene pairs that failed to meet the dcMat thresholds",
                    " were removed from the Cosinet score calculation.")
        } else {
            stop("No gene pair meets the given dcMat thresholds.")
        }
    }

    deltaScores <- deltaScore(exprMat = exprMat,
                              genePairs = genePairs,
                              condition = condition,
                              sampleNames = sampleNames,
                              bin = bin,
                              binRange = binRange,
                              threads = threads)

    if (!is.null(dcMat) & weighted == TRUE) {
        message("Weight each gene pair by absolute dcMat value.")
        deltaScores <- deltaScores * abs(weights)
    }

    deltaScores <- t(deltaScores)
    Cosinet <- apply(deltaScores, 1, function(x) mean(x, na.rm = TRUE))

    finalRes <- cbind(Cosinet, deltaScores)
    finalRes <- as.data.frame(finalRes)

    message(ncol(deltaScores),
            " gene pairs contributed to the final Cosinet score calculation.")

    return(finalRes)
}

deltaScore <- function(exprMat,
                       genePairs,
                       condition,
                       sampleNames=colnames(exprMat),
                       bin=0.1,
                       binRange = binRange,
                       threads=2) {

    condition <- factor(condition)
    genes2Use <- unique(c(genePairs[, 1], genePairs[, 2]))
    exprMat1 <- exprMat[genes2Use, condition %in% levels(condition)[1]]
    exprMat2 <- exprMat[genes2Use, condition %in% levels(condition)[2]]

    cl <- makeCluster(threads)
    registerDoParallel(cl)

    message("Performing calculation using ", threads, " threads...")

    nSamples <- ceiling(length(sampleNames) / threads)
    nTimes <- min(threads, ceiling(length(sampleNames) / nSamples))

    i <- 1
    res <- foreach(i = seq_len(nTimes),
                   .combine = "c",
                   .inorder = TRUE,
                   .verbose = FALSE,
                   .export = c("rhoScore")
    ) %dopar% {
        start <- (i - 1) * nSamples + 1
        end <- min(i * nSamples, length(sampleNames))
        deltaList <- list()
        for (sampleName in sampleNames[start:end]) {
            if (!sampleName %in% colnames(exprMat1)) {
                m1 <- cbind(exprMat1, exprMat[genes2Use, sampleName, drop = FALSE])
                colnames(m1)[ncol(m1)] <- sampleName
            } else {
                m1 <- exprMat1
            }
            if (!sampleName %in% colnames(exprMat2)) {
                m2 <- cbind(exprMat2, exprMat[genes2Use, sampleName, drop = FALSE])
                colnames(m2)[ncol(m2)] <- sampleName
            } else {
                m2 <- exprMat2
            }
            rho1 <- rhoScore(m1, sampleName, genePairs,
                             bin = bin, binRange = binRange)
            rho2 <- rhoScore(m2, sampleName, genePairs,
                             bin = bin, binRange = binRange)
            rho1[rho1 < 0] <- 0
            rho2[rho2 < 0] <- 0
            delta <- rho2 - rho1

            deltaList[[sampleName]] <- delta
        }
        return(deltaList)
    }
    res <- as.data.frame(res)
    rownames(res) <- paste(genePairs[, 1], genePairs[, 2], sep = "_")

    stopCluster(cl)
    return(res)
}

rhoScore <- function(m,
                     sampleName,
                     genePairs,
                     bin=0.1,
                     binRange=1) {
    n <- ncol(m)
    sideRange <- (1 - binRange) / 2
    lowerExp <- apply(m[, !colnames(m) %in% sampleName], 1,
                      function(x) quantile(x = x,
                                           probs = sideRange,
                                           na.rm = TRUE,
                                           names = FALSE))
    upperExp <- apply(m[, !colnames(m) %in% sampleName], 1,
                      function(x) quantile(x = x,
                                           probs = 1 - sideRange,
                                           na.rm = TRUE,
                                           names = FALSE))
    binVal <- (upperExp - lowerExp) * bin * 0.5
    diffExp <- abs(m - m[, sampleName])
    diffExp[, sampleName] <- Inf
    binMat <- diffExp <= binVal
    binMat <- t(binMat * 1)
    rhoScores <- mapply(function(gene1, gene2) {
        x <- binMat[, gene1]
        y <- binMat[, gene2]
        (sum(x & y) / n - sum(x) / n * sum(y) / n) * 100
    }, genePairs[, 1], genePairs[, 2], SIMPLIFY = TRUE, USE.NAMES = FALSE)

    return(rhoScores)
}
