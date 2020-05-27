\name{rarexpected_fun}
\alias{rarexpected_fun}

\title{
Functional Dissimilarity Rarefaction Curves
}
\description{
The function calculates functional dissimilarity rarefaction curves using the method proposed by Ricotta et al. (2010).
}
\usage{
rarexpected_fun(comm,dist_f)
}

\arguments{
  \item{comm}{a community data frame or matrix with N plots as rows, S species as columns. Only presence/absence (1/0) of species in plots are allowed as entries.}
  \item{dist_f}{an object of class \code{dist} containing the pairwise functional distances among species. The distance matrix should be squared Euclidean or simply Euclidean.} 
}
\details{
The expected richness of species at an increasingly larger cumulative number of plots N is used for the calculation of Rao's quadratic entropy to obtain an expected functional dissimilarity in the N random plots. The expected richness is the sum of the probability for every species in the community  to find that species in the N random plots.
}
\value{
A vector containing the values of the expected functional dissimilarity for every sampling dimension is returned.
}
\references{
Botta-Dukat, Z. (2005) Rao's quadratic entropy as a measure of functional diversity
based on multiple traits. \emph{Journal of Vegetation Science}, \bold{16}, 533--540.

Rao, C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. \emph{Theoretical Population Biology}, \bold{21}, 24--43.

Ricotta, C., Burrascano, S., Blasi, C. (2010) Incorporating functional dissimilarities into sample-based rarefaction curves: from taxon resampling to functional resampling. \emph{Journal of Vegetation Science}, \bold{21}, 280--286.
}
\author{
Giovanni Bacaro \email{gbacaro@units.it}

Elisa Thouverai \email{elisa.th95@gmail.com}
}

\seealso{
\code{\link{ser_functional}}
}
\examples{
require(ade4)
require(vegan)
data(duneFVG)
data(duneFVG.tr8)
tr8_N<-duneFVG.tr8$traits.tot[,c(1,3,4)]
tr8_D<-data.frame(duneFVG.tr8$traits.tot[,2])
tr8_Q<-duneFVG.tr8$traits.tot[,5:15]
tr8dist<-dist.ktab(ktab.list.df(list(tr8_N,tr8_D,tr8_Q)),type=c('N','D','Q'))
d<-decostand(duneFVG$total, method='pa')
rare<-rarexpected_fun(d,tr8dist)
plot(rare, ylab="Rao QE", xlab="Number of plots", type="l")
}

\keyword{models}
