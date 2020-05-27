\name{rare_beta}
\alias{rare_beta}

\title{
Beta Diversity Rarefaction Curve
}
\description{
The function calculates spatial and non-spatially explicit rarefaction curves using different indexes of beta diversity.
}
\usage{
rare_beta(comm,dist_xy=NULL,method=c("whittaker","jaccard","bray","cody"
,"fun_div"),random=99,fun_div=NULL,args=NULL,verbose=FALSE,spatial=FALSE)
}

\arguments{
 \item{comm}{a community dataframe or matrix with N plots as rows, S species as columns. Both the presence/absence (1/0) or the abundances of species in plots are allowed as entries. Plot names should be provided as row names.}
  \item{dist_xy}{an object of class \code{dist} containing the pairwise geographic distances among the plots. The names of the labels must be the same as the rows of \code{comm}. The argument could be empty (default=NULL) if non-spatial rarefaction curves are considered.}
  \item{method}{the diversity index for the calculation of the rarefaction curve, one among "whittaker","jaccard","bray","cody" (only for spatially-explicit rarefaction) or "fun_div". See details.}
  \item{random}{number of iterations to construct the non-spatially explicit rarefaction curve}
  \item{fun_div}{a string with the name of the function you want to use to calculate the rarefaction. The function must calculate the value of the chosen index per plot and return a pairwise distance matrix with the values.}
  \item{args}{a list with the arguments for fun_div. The value \code{NA} should be given in place of the community matrix in the list. The names of the elements must correspond to the names of the arguments of the function passed.}
   \item{verbose}{if TRUE, the arguments of \code{fun_div} are inserted interactevely by the user.}
  \item{spatial}{if TRUE, the function calculates the spatially-explicit rarefaction.}
}
\details{
The avaiable methods are:

\code{whittaker}: Whittaker's species turnover calculates how many times there is a change in species composition among the plots: \deqn{\beta_w=\gamma/\alpha -1}{beta_w= gamma/alpha - 1} 

\code{jaccard}: Jaccard similarity coefficient is a pairwise measure of similarity between plots: \deqn{\beta_j=a/(\alpha_1+\alpha_2-a)}{beta_j= a/(alpha_1 + alpha_2 – a)} where \eqn{a} is the number of species in common between two plots, and \eqn{\alpha_1} and \eqn{\alpha_2} are the values of alpha diversity (species richness) of the 2 plots compared. The rarefaction curve is calculated with the mean of the pairwise diversities in N random plots.

\code{bray}: Bray-Curtis dissimilarity is a pairwise measure of similarity between plots weighted by the abundances of the species: \deqn{\beta_{bray}=\frac{\sum_{i} \left | x_i - x_j \right |}{\sum_{i} x_i+x_j}}{beta_bray= sum_i (abs(x_i - x_j))/sum_i(x_i + x_j)} where \eqn{x_i} and \eqn{x_j} are the abundances of the species \eqn{x} in the plots \eqn{i} and \eqn{j}.

\code{cody}: Cody index is defined as the rate at which species are being replaced in censuses at each point on the habitat gradient and is fixed for samples arranged along gradients of environmental change: \deqn{\beta_c=[g(H)+l(H)]/2}{beta_c= [g(H) + l(H)]/2} where \eqn{g(H)} is the number of species gained along the habitat gradient \eqn{H} and \eqn{l(H)} is the number of species lost.

In some cases, the rarefaction doesn't behave like a non-decreasing function of plot number \emph{N}. This behavior is typical of non-concave indexes. 
}
\value{
An object of class \code{data.frame} with 3 columns is returned: 

  - Rarefaction: : mean of the values of the accumulation curves for all the sampling dimensions;
  
  - IC_up: upper confidence interval;
  
  - IC_low: lower confidence interval.
}
\references{
Baselga, A. (2013) Multiple site dissimilarity quantifies compositional heterogeneity among several sites, while average pairwise dissimilarity may be misleading. \emph{Ecography}, \bold{36}, 124--128.

Koleff, P., Gaston, K.J., Lennon, J.J. (2003) Measuring beta diversity for presence-absence data. \emph{Journal of Animal Ecology}, \bold{72}, 367--382.

Wilson, M.V., Shmida, A. (1984) Measuring Beta Diversity with Presence-Absence Data. \emph{Journal of Ecology}, \bold{72}, 1055--1064.
}
\author{
Elisa Thouverai \email{elisa.th95@gmail.com}

Enrico Tordoni \email{etordoni@units.it}

with contributions of Sandrine Pavoine.
}

\examples{
\donttest{
data(duneFVG)
data(duneFVG.xy)
dist_sp<-vegdist(duneFVG.xy$tot.xy,'euclidean')
## Spatial explicit rarefaction
## Take care that the following scripts may be a bit time-consuming 
## depending on your computer (less than a minute to a few minutes)
rare<-rare_beta(duneFVG$total,dist_sp,method="bray",spatial=TRUE)
rare1<-rare_beta(duneFVG$total,method="bray",random=99)
plot(rare [,1], ylab="Rao QE", xlab="Number of plots", type="l", ylim=range(rare, na.rm
=TRUE))
lines(rare[,2], lty=2)
lines(rare[,3], lty=2)
plot(rare1 [,1], ylab="Rao QE", xlab="Number of plots", type="l", ylim=range(rare, na.rm
=TRUE))
lines(rare1[,2], lty=2)
lines(rare1[,3], lty=2)

}}

\keyword{models}

