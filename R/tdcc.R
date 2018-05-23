#' @title
#' A Measure of Top-Down Correlation
#'
#' @description
#' With the use of Savage scores, the Top-Down Correlation Coefficient TDCC 
#' compares \code{b} rankings.
#'
#' @details
#' NOTE: As the implementation of the coefficient of concordance is still defective,
#' please use the Pearson coefficient!
#'
#' @param ranks [\code{matrix(nrow = b, ncol = k)}]\cr
#'   \code{(bxk)}-matrix of the ranks of the \code{k}
#'   variables for each of the \code{b} sensitivity analyses, ties are neglected,
#'   must be integers.
#' @param pearson [\code{logical(1)}]\cr
#'   Should the ordinary Pearson coefficient with
#'   Savage scores be computed (\code{b = 2})? Default is \code{FALSE}.
#' @param plot [\code{logical(1)}]\cr
#'   Should scatter plots showing rankings and Savage scores be created
#'   (\code{b = 2})? Default is \code{FALSE}.
#'
#' @return A named vector with components:
#' \itemize{
#'   \item{\code{kendall}}: Coefficient of concordance.
#'   \item{\code{pearson}}: Pearson coefficient (only if \code{pearson = TRUE}).
#' }
#'
#' @examples
#' ranking <- rbind(1:20,
#'                  c(1,3,2,4,16,10,19,12,18,17,
#'                    20,5,14,7,8,11,6,15,9,13))
#' tdcc(ranking, pearson = TRUE, plot = TRUE)
#'
#' @author Stefan Theers
#' @references R. L. Iman and W. J. Conover, 
#'   \emph{A Measure of Top-Down Correlation}, 
#'   Technometrics, Vol. 29, No. 3 (Aug., 1987), pp. 351--357.
#'
#' @export
#'
tdcc <- function(ranks, pearson = FALSE, plot = FALSE) {
  # input checks:
  assertMatrix(ranks, mode = "numeric", min.rows = 2, min.cols = 2)
  assertIntegerish(ranks)
  assertLogical(pearson, len = 1)
  assertLogical(plot, len = 1)
  if(pearson | plot) stopifnot(nrow(ranks) == 2)
  # number of sensitivity methods which should be compared:
  b <- nrow(ranks)
  # number of variables:
  k <- ncol(ranks)
  # Savage scores:
  Sfun <- function(i) sum(1 / i:k)
  S <- matrix(sapply(t(ranks), Sfun), nrow = b, byrow = TRUE)
           ## Warum nicht gleich?
           ## cscores(ranks[1, ], type="Savage")
  # Savage-Scores fuer i-tes Merkmal ueber alle b Rankings:
  S.sums <- colSums(S)
  # Kendall's Coefficient of Concordance:
  S1 <- sum(1 / 1:k)
  C.T <- 1 / (b^2 * (k - S1)) * sum(S.sums^2) - b^2 * k
  # scatter plots:
  if(plot) {
    old.pars <- graphics::par(no.readonly = TRUE)
    on.exit(do.call(graphics::par, old.pars))
    graphics::par(mfrow = c(1, 2))
    plot(ranks[1, ], ranks[2, ], main = "Scatterplot of Ranks",
         xlab = "Ranking by A", ylab = "Ranking by B")
    plot(S[1, ], S[2, ], main = "Scatterplot of Savage Scores",
         xlab = "Savage Score for A", ylab = "Savage Score for B")
  }
  # return value:
  if(!pearson) {
    retVal <- c(kendall = C.T)
  } else {
    r.T <- 1 / (k - S1) * (sum(apply(S, 2, prod)) - k)
    # Analogously to section 2, a permutation test could be implemented here:
    retVal <- c(kendall = C.T, pearson = r.T)
  }
  return(retVal)
}
