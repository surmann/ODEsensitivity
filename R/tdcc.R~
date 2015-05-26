tdcc <-
function(ranks, pearson = FALSE, sp = FALSE) {
  # Plausibilitaet:
  assertMatrix(ranks, mode = "numeric", min.rows = 2, min.cols = 2)
  assertIntegerish(ranks)
  assertLogical(pearson, len = 1)
  assertLogical(sp, len = 1)
  if(pearson | sp) stopifnot(nrow(ranks) == 2)
  # Anzahl zu vergleichender SA-Methoden:
  b <- nrow(ranks)
  # Anzahl Merkmale:
  k <- ncol(ranks)
  # Savage-Scores:
  S.fun <- function(i) sum(1 / i:k)
  S <- matrix(sapply(t(ranks), S.fun), nrow = b, byrow = TRUE)
           ## Warum nicht gleich?
           ## cscores(ranks[1, ], type="Savage")
  # Savage-Scores fuer i-tes Merkmal ueber alle b Rankings:
  S.sums <- colSums(S)
  # Kendall's Coefficient of Concordance:
  S1 <- sum(1 / 1:k)
  C.T <- 1 / (b^2 * (k - S1)) * sum(S.sums^2) - b^2 * k
  # Scatterplots:
  if(sp) {
    par(mfrow = c(1, 2))
    plot(ranks[1, ], ranks[2, ], main = "Scatterplot of Ranks",
         xlab = "Ranking by A", ylab = "Ranking by B")
    plot(S[1, ], S[2, ], main = "Scatterplot of Savage Scores",
         xlab = "Savage Score for A", ylab = "Savage Score for B")
    par(mfrow = c(1, 1))
  }
  if(!pearson)
    return(C.T)
  if(pearson) {
    r.T <- 1 / (k - S1) * (sum(apply(S, 2, prod)) - k)
    # An dieser Stelle kann analog zu Abschnitt 2 ein Permutationstest
    # implementiert werden.
    return(list(kendall = C.T, pearson = r.T))
  }
}
