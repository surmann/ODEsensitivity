\name{plot}
\alias{plot.intwo}
\alias{plot.trajectories}
\title{
Plotting the results of Morris SA
}
\description{
plot.intwo plots mu.star and sigma seperately versus time.
plot.trajectories plots mu.star versus sigma for every point of time.
}
\usage{
plot.intwo(res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
of type cf. example.
}
}
\details{
details.
}
\value{
none.
}
\author{
Stefan Theers
}
\seealso{
\link[sensitivity]{morris}
}
\examples{
# Modell aufstellen fuer die DGL y'(t) = alpha * y(t):
dgl.fun <- function(X, my.t) {
  exp(X[, 1] * my.t) + X[, 2] * my.t
  ## exp(X[, 1] * my.t)            # natuerlich die wahre Lsg.
}

# Fuehre fuer diesen einfachen Fall eine SA zu verschiedenen Zeitpunkten
# durch:
t.vec <- 0:10
x.fun <- function(my.t) {
  morris(model = dgl.fun, factors = 2, r = 100, my.t = my.t,
         design = list(type = "oat", levels = 100, grid.jump = 1))
}
one.run <- function(x.fun, my.t) {
  x <- x.fun(my.t)
  # analog zur Hilfeseite von morris()/ hoestgradig primitiv:
  k <- ncol(x$ee)
  mu <- mu.star <- sigma <- numeric(k)
  for(i in 1:k) {
    mu[i]      <- mean(x$ee[, i])
    mu.star[i] <- mean(abs(x$ee[, i]))
    sigma[i]   <- sd(x$ee[, i])
  }
  # Ergebnisse:
  res <- c(my.t, mu, mu.star, sigma)
  names(res) <- c("time",
                  paste("mu", 1:k, sep = ""),
                  paste("mu.star", 1:k, sep = ""),
                  paste("sigma", 1:k, sep = ""))
  return(res)
}
set.seed(2015)
res <- sapply(t.vec, one.run, x.fun = x.fun)

# Plots:
plot.intwo(res)
plot.trajectories(res)
}
