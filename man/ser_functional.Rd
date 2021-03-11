\name{ser_functional}
\alias{ser_functional}

\title{
Functional Spatially-Explicit Rarefaction Curves
}
\description{
The function calculates classic and spatially-explicit rarefaction curves using different indexes of functional diversity.
}
\usage{
ser_functional(comm,dist_f=NULL,dist_xy,method=c("rao","chao","fun_div"), 
tau=NA,q=0,comparison=FALSE,resampling=99,fun_div=NULL,args=NULL,verbose=FALSE)
}

\arguments{
   \item{comm}{a community dataframe or matrix with N plots as rows, S species as columns. Only the abundances of species in plots are allowed as entries for the methods \code{rao} and \code{chao}. Plot names should be provided as row names.}
 \item{dist_f}{an object of class \code{dist} containing the pairwise functional distances among species. The distance matrix should be squared Euclidean or simply Euclidean.} 
  \item{dist_xy}{an object of class \code{dist} containing the pairwise geographic distances among the plots or an object of class, an object of class \code{vector} containing the order of the sampling units along a gradient or an object of class \code{matrix} or \code{data.frame} where each column represents a gradient along which the sampling units are ordered. The names of the labels must be the same as the rows of \code{comm}. If \code{dist_xy} is not of class dist, then a distance matrix is calculated using the Euclidean distance coefficient.}
  \item{method}{the diversity index for the calculation of the rarefaction curve, one among "rao", "chao" or "fun_div". See details.}
  \item{tau}{threshold of functional distinctiveness between any two species for \code{chao}. See details.}
  \item{q}{parameter that determines the sensitivity of the measure to the relative abundance of the species for \code{chao}.}
  \item{comparison}{if TRUE, both non-spatially explicit and spatially explicit rarefactions are calculated.}
  \item{resampling}{number of times plots (rows) are randomly resampled from comm to calculate the mean accumulation curve for the non-spatially-explicit rarefaction.}
  \item{fun_div}{a string with the name of the user-defined function for the diversity index in the rarefaction. The function must calculate the value of the chosen diversity index per plot and return a numeric vector with the values calculated.}
  \item{args}{a list with the arguments for fun_div. The value NA should be given in place of the community matrix in the list. The names of the elements must correspond to the names of the arguments of the function passed.}
  \item{verbose}{if TRUE, the arguments of \code{fun_div} are inserted interactively by the user.}
}
\details{
The available methods are:

\code{rao}: Rao functional diversity is based on the quadratic entropy of Rao (1982) and incorporates both the relative abundances of species and a measure of the pairwise functional distances between species. It expresses the average difference between two randomly selected individuals with replacements. \deqn{Q(\mathbf{p}_i,\mathbf{D})=\sum_{k=1}^S\sum_{l=1}^S p_{k}p_{l}d_{kl}}{Q(p_i,D)=sum_k,l p_k p_l d_kl} where \eqn{p=(p1,...,p_k,...,S)}{p=(p1,...,p_k,...,S)} is the vector of relative species abundances, \eqn{S} is  the  number  of  species, \eqn{\mathbf{D}=(d_{kl})}{D=(d_kl)} is the  matrix of functional dissimilarities  among  species, and \eqn{d_{kl}}{d_kl} is the  functional dissimilarity between species \eqn{k} and \eqn{l}.

\code{chao}: Chao et al. (2018) functional beta-diversity index FD quantifies the effective number of equally-distinct functional groups in the considered plots at the distinctiveness \code{tau} threshold. Any two species with functional distance greater than or equal to \code{tau}, are treated as functionally equally-distinct and as belonging to different functional groups with distance \code{tau}. For each pair of species with functional distance lower than \code{tau} but different from zero, only a proportion of individuals is considered functionally equally-distinct, the other proportion of individuals is considered functionally indistinct. If the pairwise distance is equal to zero, the two species are treated as belonging to the same functional group. After dividing the set of species to form functionally indistinct groups, the contribution of every species is quantified and then the FD of order \emph{q} is calculated using the Hill number of order \eqn{q}: \deqn{^{q}\textrm{FD}(\Delta(\tau))=\left ( \sum_{i=1}^{S} \nu_{i}(\tau)\left(\frac{a_i(\tau)}{n_{+}} \right )^{(1/1-q)} \right )}{FD(tau)= (sum_i v_i(tau)(a_i /(n+))^1/1-q))} where \eqn{a_{i}(\tau)}{a_i(tau)} is the combined abundance of all functionally-indistinct individuals from species \eqn{i}, \eqn{v_{i}(\tau)=n_{i}/a_{i}(\tau)}{v_i(tau)= n_i/a_i(tau)} represents the attribute contribution of species \eqn{i} for a threshold level \code{tau} (\eqn{n_{i}}{n_i} is the abundance of species \eqn{i}), \eqn{n_+}{n+} is the total number of individuals in the community and \eqn{q} is the parameter that determines the sensitivity of the measure to the relative abundance of the species.
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

Chao, A., Chiu, C.-H., Villeger, S., Sun, I.F., Thorn, S., Lin, Y., Chiang, J.-M., B. Sherwin, W. (2018) An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. \emph{Ecological Monographs}, \bold{89}, e01343.

Rao, C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. \emph{Theoretical Population Biology}, \bold{21}, 24--43.

Ricotta, C., Pavoine, S., Bacaro, G.,  Acosta, A. (2012) Functional rarefaction for species abundance data. \emph{Methods in Ecology and Evolution}, \bold{3}, 519--525. 

}
\author{
Elisa Thouverai \email{elisa.th95@gmail.com}

with contributions of Sandrine Pavoine.
}

\seealso{
\code{\link{rare_Rao}}

\code{\link{rare_alpha}} 
}
\examples{
\dontrun{
#Time consuming
require(ade4)
data(duneFVG)
data(duneFVG.xy)
data(duneFVG.tr8)
tr8_N<-duneFVG.tr8$traits.tot[,c(1,3,4)]
tr8_D<-data.frame(duneFVG.tr8$traits.tot[,2])
tr8_Q<-duneFVG.tr8$traits.tot[,5:15]
tr8dist<-dist.ktab(ktab.list.df(list(tr8_N,tr8_D,tr8_Q)),type=c('N','D','Q'))
tr8dist<-cailliez(tr8dist) 
dist_sp<-dist(duneFVG.xy$tot.xy,'euclidean')
## Rao index
## Take care that the following scripts may be time-consuming 
rare<-ser_functional(duneFVG$total,tr8dist,dist_sp,method='rao',comparison=TRUE)
plot(rare [,1], ylab="Rao QE", xlab="Number of plots", type="l", ylim=range(rare, na.rm=TRUE))
lines(rare[,2], lty=2)
lines(rare[,3], lty=2)
lines(rare[,4 ], col=2)
lines(rare[,5], lty=2, col=2)
lines(rare[,6], lty=2, col=2)
legend("bottomright",  legend=c("spatially-explicit rarefaction", "classic rarefaction"),
lty=1, col=1:2)
# dashed lines represent the confidence interval
## Chao index
rareC<-ser_functional(duneFVG$total,tr8dist,dist_sp,method='chao',tau=0.6,q=2,comparison
=TRUE)
plot(rareC [,1], ylab="Chao", xlab="Number of plots", type="l", ylim=range(rareC, na.rm
=TRUE))
lines(rareC[,2], lty=2)
lines(rareC[,3], lty=2)
lines(rareC[,4 ], col=2)
lines(rareC[,5], lty=2, col=2)
lines(rareC[,6], lty=2, col=2)
legend("bottomright",  legend=c("spatially-explicit rarefaction", "classic rarefaction"),
lty=1, col=1:2)
# dashed lines represent the confidence interval
}}
\keyword{models}
