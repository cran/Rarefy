\name{rare_alpha}
\alias{rare_alpha}

\title{Alpha Diversity Rarefaction Curve}

\description{
The function calculates spatial and non-spatially explicit rarefaction curves using different indexes of alpha diversity.
}

\usage{
rare_alpha(comm,dist_xy=NULL,method=c("HCDT","hill","fun_div"),q=0,random=99,fun_div=NULL
,args=NULL,verbose=FALSE,spatial=FALSE,mean=FALSE)
}

\arguments{
  \item{comm}{a community data frame or matrix with N plots as rows, S species as columns. Both the presence/absence (1/0) or the abundances of species in plots are allowed as entries. Plot names should be provided as row names.}
  \item{dist_xy}{an object of class \code{dist} containing the pairwise geographic distances among the plots. The names of the labels must be the same as the rows of \code{comm}. The argument could be empty (default=NULL) if non-spatial rarefaction curves are considered.}
  \item{method}{the diversity index for the calculation of the rarefaction curve, one among "HCDT","hill" or "fun_div". See details.}
  \item{q}{parameter that determines the sensitivity of the measure to the relative abundance of the species.}
   \item{random}{number of iterations to construct the non-spatially explicit rarefaction curve}
  \item{fun_div}{a string with the name of the user-defined function for the diversity index in the rarefaction. The function must calculate the value of the chosen diversity index per plot and return a numeric vector with the values calculated.}
  \item{args}{a list with the arguments for fun_div. The value NA should be given in place of the community matrix in the list. The names of the elements must correspond to the names of the arguments of the function passed.}
   \item{verbose}{if TRUE, the arguments of \code{fun_div} are inserted interactevely by the user.}
  \item{spatial}{if TRUE, the function calculates the spatially-explicit rarefaction.}
  \item{mean}{if TRUE, the mean of the relative abundances of the species in all the sampling dimensions is used for the calculation of the rarefaction instead of the total abundance, giving less importance to abundant species.}
}

\details{
The avaiable methods are:

\code{HCDT}: HCDT entropy is a generalization of the standard coefficient of entropy (Harvda and Charvat 1967; Daroczy 1970; Tsallis 1988): \deqn{HCDT=\frac{1-\left(\sum_{i}^{S}(p_i)^q \right )}{q-1}}{HCDT= 1-(sum_i (p_i)^q)/(q-1)} where where \eqn{S} is the number of species, \eqn{p_i} is the relative abundance of species \eqn{i} and \eqn{q} is the parameter that regulates the sensitivity to species abundance. When \eqn{q} tends to 1, the equation becomes the Shannon entropy: \deqn{HCDT=-\sum_{i=1}^{S}p_ilog(p_i)}{HCDT= -sum_i p_i log(p_i)}

\code{hill}: Hill numbers is a class of measures that obeys to the replication principle and integrates species richness and species abundances. The parameter \eqn{q}, called 'order', regulates the sensitivity of the index to the species abundance: with \eqn{q=0} the value of the index correspond to the species richness, with \eqn{q=1} the measure tends to the exponential of Shannon index, and with \eqn{q=2} it corresponds to the inverse of Simpson index: \deqn{^{q}\textrm{D}=\left( \sum_{i=1}^{S}p_i^q\right )^{1/(1-q)}}{D=(sum_i (p_i)^q)^(1/(1-q))} where \eqn{S} is the number of species, \eqn{p_i} is the relative abundance of species \eqn{i} and \emph{q} is the order.

\code{hill} doesn't behave like a non-decreasing function of the increasing number of plots \emph{N}. This behavior is typical of non-concave indexes. 
}
\value{
An object of class \code{data.frame} with 3 columns is returned: 

  - Rarefaction:: mean of the values of the accumulation curves for all the sampling dimensions
  
  - IC_up: upper confidence interval;
  
  - IC_low: lower confidence interval.
}

\references{
Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K., Ellison, A. M (2014). Rarefaction and extrapolation with Hill numbers: A framework for sampling and estimation in species diversity studies. \emph{Ecological Monographs}, \bold{84}, 45--67.

Daroczy, Z. (1970) Generalized information functions. \emph{Information and Control}, \bold{16}, 36--51.

Havrda, M., Charvat F. (1967) Quantification method of classification processes: concept of structural alpha-entropy. \emph{Kybernetik}, \bold{3}, 30--35.

Ricotta, C., Pavoine, S., Bacaro, G., Acosta, A. (2012) Functional rarefaction for species abundance data. \emph{Methods in Ecology and Evolution}, \bold{3}, 519--525.

Tsallis, C. (1988) Possible generalization of Boltzmann-Gibbs statistics. \emph{Journal of Statistical Physics}, \bold{52}, 480--487.

}
\author{
Elisa Thouverai \email{elisa.th95@gmail.com}

with contributions of Sandrine Pavoine.
}

\examples{
\donttest{
data(duneFVG)
data(duneFVG.xy)
dist_sp<-vegdist(duneFVG.xy$tot.xy,'euclidean')
## Spatial explicit rarefaction with HCDT
rare<-rare_alpha(duneFVG$total,dist_sp,q=2,random=99,spatial=TRUE)
## Non-spatially explicit rarefaction with Hill
rare1<-rare_alpha(duneFVG$total,dist_sp,method='hill',q=2,mean=TRUE)
##fun_div with the list of arguments passed 
#to set the argument, digit the correspondent number
#digit the number correspondent to "stop" when the process is complete
a<-list(NA,'shannon')
names(a)<-c('x','index')
rare2<-rare_alpha(duneFVG$total,method="fun_div",random=999,fun_div='diversity',args=a,mean
=TRUE)
rare3<-rare_alpha(duneFVG$total,dist_sp,method="fun_div",fun_div='diversity',args=a,spatial
=TRUE,mean=TRUE)
##fun_div without the list of arguments passed (verbose=TRUE)
rare4<-rare_alpha(duneFVG$total,method="fun_div",random=999,fun_div='diversity',mean=TRUE)
rare5<-rare_alpha(duneFVG$total,dist_sp,method="fun_div",fun_div='diversity',spatial=TRUE
,mean=TRUE)
}
}

\keyword{models}
