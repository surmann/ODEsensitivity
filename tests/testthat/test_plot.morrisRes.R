context("Testing plot.morrisRes()")

set.seed(2015)
kVec <- 2:20                        # Anzahl Inputfaktoren
potsVec <- c(1:3, 5, 20, 100)       # Anzahl Zeitpunkte

test_that("Plots are generated", {
  # Alle Anzahlen Zeitpunke und Inputfaktoren durchlaufen:
  for(k in kVec){
    for(pots in potsVec){
      # Ein Morris-Ergebnis erzeugen:
      res <- rbind(1:pots, matrix(runif((3*k) * pots), ncol = pots))
      rownames(res) <- c("time",
                         paste("mu", 1:k, sep = ""),
                         paste("mu.star", 1:k, sep = ""),
                         paste("sigma", 1:k, sep = ""))
      res <- setClasses(res, "morrisRes")
      
      # Testen:
      expect_true(plot(res))
      expect_true(plot(res, type = "trajec"))
      expect_true(plot(res, type = "sep"))
      expect_true(plot(res, main = "Something."))
    }
  }
})

test_that("Errors are thrown", {
  # Alle Anzahlen Zeitpunke und Inputfaktoren durchlaufen:
  for(k in kVec){
    for(pots in potsVec){
      # Ein Morris-Ergebnis erzeugen:
      res <- rbind(1:pots, matrix(runif((3*k) * pots), ncol = pots))
      rownames(res) <- c("time",
                         paste("mu", 1:k, sep = ""),
                         paste("mu.star", 1:k, sep = ""),
                         paste("sigma", 1:k, sep = ""))
      res <- setClasses(res, "morrisRes")
      
      # Testen:
      expect_error(plot.morrisRes(1:3))
      expect_error(plot.morrisRes("no character!"))
      expect_error(plot.morrisRes(diag(7)))
      expect_error(plot.morrisRes(res <- setClasses(diag(7), "sobolRes")))
    }
  }
})
