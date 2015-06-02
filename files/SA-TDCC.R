##--------------------------------------------------------------------##
## FOR1511: R-Code                                                    ##
##                                                                    ##
##   Top Down Correlation Coefficient                                 ##
##--------------------------------------------------------------------##

options(help_type = "html")
library("ODEnetwork")                          # Schwingungen
library("sensitivity")                         # fuer SA
library("boot")                                # fuer SA, Sobol
library("exactRankTests")                      # fuer Savage Scores
library("checkmate")                           # Plausibilitaet


##----------------------------------------------------------------------
## TDCC: Ergebnisvgl. verschiedener SA-Methoden
##----------------------------------------------------------------------

## tdcc - Funktion zur Bestimmung der TDCC fuer bel. b

## Eingabe
##   ranks   - (b x k)-Matrix der Raenge der k Merkmale fuer jede der
##             b SAs, ohne Bindungen, integers
##   pearson - logisch, soll der normale Pearsonkoeffizient mit
##             Savage Scores ausgegeben werden (mit b = 2)?
##   sp      - logisch, sollen Scatterplots fuer Rankings und Savage
##             Scores erzeugt werden (mit b = 2)?

## Ausgabe
##   TDCC

tdcc <- function(ranks, pearson = FALSE, sp = FALSE) {
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


##----------------------------------------------------------------------
## TDCC: Beispiel aus Iman & Canover (1987)
##----------------------------------------------------------------------

temp <- rbind(1:20,
              c(1,3,2,4,16,10,19,12,18,17,
                           20,5,14,7,8,11,6,15,9,13))
tdcc(temp, pearson = TRUE, sp = TRUE)

temp <- rbind(1:5, c(2,1,5,4,3))
tdcc(temp, TRUE)

temp <- rbind(1:20, 1:20)
tdcc(temp, TRUE)

temp <- rbind(1:20, 20:1)
tdcc(temp, TRUE)


##----------------------------------------------------------------------
## TDCC: Beispiel der Ishigami-Funktion
##----------------------------------------------------------------------

# SA mit Morris:
res.m <- morris(model = ishigami.fun, factors = 3, r = 4,
                design = list(type = "oat", levels = 5, grid.jump = 3))
mu.star <- apply(res.m$ee, 2, function(x) mean(abs(x)))
ranking.m <- rank(mu.star)                   # niedriger Rang = wichtig
# SA mit Sobol:
n <- 1000
X1 <- data.frame(matrix(runif(3 * n), nrow = n))
X2 <- data.frame(matrix(runif(3 * n), nrow = n))
res.s <- sobol2002(model = ishigami.fun, X1 = X1, X2 = X2, nboot = 100)
res.s$T     # totale S-Indizes
res.s$S     # Haupteffekt-S-Indizes
ranking.s <- rank(- res.s$T[, 1])
# TDCC bestimmen:
ranking <- rbind(as.integer(ranking.m), as.integer(ranking.s))
tdcc(ranking, pearson = TRUE, sp = TRUE)
# Es ex. tatsaechlich ein Unterschied in der Bewertung der Merkmale.
