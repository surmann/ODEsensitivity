leftHandSide <-
function(X, Xinv) {
  norm <- function(x) sqrt(sum(x^2))

  y <- trunc(rnorm(ncol(Xinv), sd = 256L))
  b0 <- drop(Xinv %*% y)

  y0 <- drop(X %*% b0) ## \hat y0
  norm_y0 <- norm(y0)

  r0 <- y - y0
  norm_r0 <- norm(r0)

  i <- 19L
  while (norm_r0 < (norm_y0 * 2^i))
    i <- i - 2L
  if (i < -19L)
    i <- -21L
  ck <- 2^(c(-21L, -3L, 1L, 19L) - i)
  yk <- y0 + r0 %*% t(ck)
  list(yk = yk, beta = Xinv %*% yk, phi = atan(ck))
}
