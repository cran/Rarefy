\name{rao_permuted}
\alias{rao_permuted}

\title{
Functional Rarefaction Curves on Standardized Species Number 
}
\description{
The function calculates functional rarefaction curves using Rao's quadratic entropy (Q) for a reference community. However, functional distance matrixes are calculated by selecting randomly species with their related functional traits from another larger community. The number of species randomly selected correspond to the number of species in the reference community matrix. The relative abundances of the standard community are used instead for the iterative calculation of the rarefacted index. 
}
\usage{
rao_permuted(comm_st,dist_f,random=99)
}

\arguments{
  \item{comm_st}{a dataframe or matrix of the standard community with N plots as rows, S species as columns. Only the abundances of species in plots are allowed as entries. Plot names should be provided as row names.}
  \item{dist_f}{an object of class \code{dist} containing the pairwise functional distances among the species of the larger community. The number of species must be larger than \code{comm_st}. The distance matrix should be squared Euclidean or simply Euclidean.}
  \item{random}{number of iterations  for the selection of species from the whole community. The functional rarefaction curve and associated Confidence Intervals are are calculated as the average of these iterations.}
}

\details{
Rao's quadratic entropy (Rao 1982) is commonly used as a measure of functional diversity. It incorporates both the relative abundance of species with their pairwise functional dissimilarities. It expresses the average difference between two randomly selected individuals with replacements.

The functional dissimilarity matrix is calculated using the function \code{dist.ktab} from the package \code{ade4}. A generalization of the Gower's general coefficient of distance is used to allow the treatment of various statistical types of variables when calculating the functional pairwise distances of the species.

To calculate the rarefaction curve, a standard number of species that correspond to the number of species of \code{comm_st} are randomly selected in \code{tr8} with their related functional traits, then a functional distance matrix is obtained and Q is calculated for 1,2,...,N plots using the mean relative abundances of the species in \code{comm_st}. This procedure is repeated \code{random} times and the average rarefaction curve is calculated.
}
\value{
An object of class \code{data.frame} with 3 columns is returned: 

  - Rarefaction: : mean of the values of the accumulation curves for all the sampling dimensions;
  
  - IC_up: upper confidence interval;
  
  - IC_low: lower confidence interval.
}
\references{
Botta-Dukat, Z. (2005) Rao's quadratic entropy as a measure of functional diversity
based on multiple traits. \emph{Journal of Vegetation Science}, \bold{16}, 533--540.

Pavoine S., Vallet, J., Dufour, A.-B., Gachet, S. and Daniel, H. (2009) On the challenge of treating various types of variables: Application for improving the measurement of functional diversity. \emph{Oikos}, \bold{118}, 391--402.

Rao, C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. \emph{Theoretical Population Biology}, \bold{21}, 24--43.

Ricotta, C., Pavoine, S., Bacaro, G.,  Acosta, A. (2012) Functional rarefaction for species abundance data. \emph{Methods in Ecology and Evolution}, \bold{3}, 519--525. 

Tordoni, E., Petruzzellis, F., Nardini, A., Savi, T., Bacaro, G. (2019) Make it simpler: alien species decrease functional diversity of coastal plant communities. \emph{Journal of Vegetation Science}, \bold{30}, 498–-509.

}
\author{
Elisa Thouverai \email{elisa.th95@gmail.com}

Enrico Tordoni \email{etordoni@units.it}
}

\seealso{
\code{\link{ser_functional}}
}
\examples{
require(ade4)
data(duneFVG)
data(duneFVG.tr8)
tr8_N<-duneFVG.tr8$traits.tot[,c(1,3,4)]
tr8_D<-data.frame(duneFVG.tr8$traits.tot[,2])
tr8_Q<-duneFVG.tr8$traits.tot[,5:15]
tr8dist<-dist.ktab(ktab.list.df(list(tr8_N,tr8_D,tr8_Q)),type=c('N','D','Q'))
## Rarefaction
rare<-rao_permuted(duneFVG$alien,tr8dist)
plot(rare$Rao, ylab="Rao QE", xlab="Number of plots", type="l", ylim=range(rare))
lines(rare[,2], lty=2)
lines(rare[,3], lty=2)
}
\keyword{models}
