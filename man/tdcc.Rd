\name{tdcc}
\alias{tdcc}
\title{
A Measure of Top-Down Correlation
}
\description{
The Top-Down Correlation Coefficient TDCC compares \code{b} rankings
using Savage Scores.
}
\usage{
tdcc(ranks, pearson = FALSE, sp = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ranks}{
(bxk)-matrix of the ranks of the k variables for each of the b SAs,
ties are neglected, must be integers
}
  \item{pearson}{
logical, should the ordinary Pearson coefficient with Savage scores be
computed (b = 2)?
}
  \item{sp}{
logical, scatter plots showing rankings and Savage scoresn (b = 2)?
}
}
\details{
details.
}
\value{
A named list with components:
\item{kendall} TDCC computed using Kendall's Coefficient of Concordance
\item{pearson} TDCC computed using Pearsons's ordinary Correlation Coefficient
}
\references{
Iman and Conover (1987): A Measure of Top-Down Correlation
}
\author{
Stefan Theers
}
\seealso{
no.
}
\examples{
ranking <- rbind(1:20,
                 c(1,3,2,4,16,10,19,12,18,17,
                   20,5,14,7,8,11,6,15,9,13))
tdcc(ranking, pearson = TRUE, sp = TRUE)
}
