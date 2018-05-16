#' @title
#' A Measure of Top-Down Correlation
#'
#' @description
#' The Top-Down Correlation Coefficient TDCC compares \code{b} rankings
#' using Savage Scores.
#'
#' @details
#' NOTE: As the implementation of the coefficient of concordance is still defective,
#' please use the Pearson coefficient!
#'
#' @param ranks [\code{matrix(nrow = b, ncol = k)}]\cr
#'   \code{(bxk)}-matrix of the ranks of the \code{k}
#'   variables for each of the \code{b} SAs, ties are neglected,
#'   must be integers.
#' @param pearson [\code{logical(1)}]\cr
#'   logical, should the ordinary Pearson coefficient with
#'   Savage scores be computed (\code{b = 2})?
#' @param plot [\code{logical(1)}]\cr
#'   logical, scatter plots showing rankings and Savage scores
#'   (\code{b = 2})?
#'
#' @return A named vector with components:
#' \itemize{
#'   \item{\bold{kendall}}: Coefficient of concordance.
#'   \item{\bold{pearson}}: Pearson Coefficient (only if \code{pearson = TRUE}).
#' }
#'
#' @examples
#' ranking <- rbind(1:20,
#'                  c(1,3,2,4,16,10,19,12,18,17,
#'                    20,5,14,7,8,11,6,15,9,13))
#' tdcc(ranking, pearson = TRUE, plot = TRUE)
#'
#' @author Stefan Theers
#' @references Iman and Conover (1987): A Measure of Top-Down Correlation
#'
#' @export
#'
tdcc <- function(ranks, pearson = FALSE, plot = FALSE) {
  # Plausibilitaet:
  assertMatrix(ranks, mode = "numeric", min.rows = 2, min.cols = 2)
  assertIntegerish(ranks)
  assertLogical(pearson, len = 1)
  assertLogical(plot, len = 1)
  if(pearson | plot) stopifnot(nrow(ranks) == 2)
  # Anzahl zu vergleichender SA-Methoden:
  b <- nrow(ranks)
  # Anzahl Merkmale:
  k <- ncol(ranks)
  # Savage-Scores:
  Sfun <- function(i) sum(1 / i:k)
  S <- matrix(sapply(t(ranks), Sfun), nrow = b, byrow = TRUE)
           ## Warum nicht gleich?
           ## cscores(ranks[1, ], type="Savage")
  # Savage-Scores fuer i-tes Merkmal ueber alle b Rankings:
  S.sums <- colSums(S)
  # Kendall's Coefficient of Concordance:
  S1 <- sum(1 / 1:k)
  C.T <- 1 / (b^2 * (k - S1)) * sum(S.sums^2) - b^2 * k
  # Scatterplots:
  if(plot) {
    graphics::par(mfrow = c(1, 2))
    plot(ranks[1, ], ranks[2, ], main = "Scatterplot of Ranks",
         xlab = "Ranking by A", ylab = "Ranking by B")
    plot(S[1, ], S[2, ], main = "Scatterplot of Savage Scores",
         xlab = "Savage Score for A", ylab = "Savage Score for B")
    graphics::par(mfrow = c(1, 1))
  }
  if(!pearson)
    return(c(kendall = C.T))
  if(pearson) {
    r.T <- 1 / (k - S1) * (sum(apply(S, 2, prod)) - k)
    # An dieser Stelle kann analog zu Abschnitt 2 ein Permutationstest
    # implementiert werden.
    return(c(kendall = C.T, pearson = r.T))
  }
}
