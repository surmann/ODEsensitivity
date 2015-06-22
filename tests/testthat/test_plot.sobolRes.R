context("Testing plot.sobolRes()")

test_plot.sobolRes <- function() {
  set.seed(2015)
  kVec <- 2:20                        # Anzahl Inputfaktoren
  potsVec <- c(1:3, 5, 20, 100)       # Anzahl Zeitpunkte
  oneTest <- function(k, pots) {
    # Ein Sobol-Ergebnis erzeugen:
    res <- list(S = rbind(1:pots, matrix(runif(k * pots), ncol = pots)))
    res$T <- rbind(res$S[1, ], as.matrix(res$S[-1, ] * 2, ncol = pots))
    rownames(res$S) <- rownames(res$T) <-
      c("time", paste("input", 1:k, sep = ""))
    res <- setClasses(res, "sobolRes")

    # Testen:
    expect_true(plot(res))
    expect_true(plot(res, type = "p"))
    expect_true(plot(res, type = "l"))
    expect_true(plot(res, type = "b"))
    expect_true(plot(res, type = "c"))
    expect_true(plot(res, type = "n"))
    expect_true(plot(res, legendPos = "topleft"))
    expect_true(plot(res, legendPos = "bottomright"))
    expect_true(plot(res, legendPos = "top"))
    expect_error(plot(res, main = "Something."))
    expect_error(plot.sobolRes(1:3))
    expect_error(plot.sobolRes("no character!"))
    expect_error(plot.sobolRes(diag(7)))
    expect_error(plot.sobolRes(res <- setClasses(diag(7), "sobolRes")))
    expect_error(plot(res, type = 1))
    expect_error(plot(res, type = "No!"))
    expect_error(plot(res, legendPos = 1))
    expect_error(plot(res, legendPos = "No!"))
    expect_error(plot.sobolRes())
  }

  # Alle Anzahlen Zeitpunke und Inputfaktoren durrchlaufen:
  for(k in kVec) {
    for(pots in potsVec)
      oneTest(k, pots)
  }
}
