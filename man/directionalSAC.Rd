\name{directionalSAC}
\alias{directionalSAC}

\title{
Spatial Explicit Rarefaction Curves
}
\description{
The function calculates directional and non-directional accumulation curves of species diversity as a function of sampling effort.
}
\usage{
directionalSAC(community,gradient)
}

\arguments{
  \item{community}{a community dataframe with N plots as rows, S species as columns. Both presence/absence and species abundances are allowed as entries. Plot names should be provided as row names on the dataframe.}
  \item{gradient}{a vector of numeric values, a matrix, a dataframe or an object of class \code{dist}.}
}
\details{
If \code{gradient} is a vector, then plots are ordered along a single spatial or environmental gradient. The length of the vector must be the number of rows in object community. Values in the vector must be in the same order as plots in object community. If names are given to each value of the vector, then the row names in \code{community} must be equal as those in \code{gradient}. For example, values in the vector may be the latitude of each plot or an environmental variable such as the temperature.

If \code{gradient} is a matrix or a data frame, then each column of the matrix must be a gradient along which plots are ordered. Only numeric values are allowed. For example, values in the matrix may be environmental data such as the temperature, the precipitation, the elevation. The result given by the function is an average over all gradients specified by the matrix. The number of rows in the matrix must be the same as the number of rows in \code{community}. If the matrix has row names, then they should be the same as in \code{community}. There must be as many columns in \code{gradient} as there are gradients of interest (for example, as many columns as there are environmental variables).

If \code{gradient} is of class \code{dist}, then the result given by the function is an average over all possible directional accumulation curves. In that case, \code{gradient} contains any pairwise dissimilarity/distance measure among plots. If names are given in \code{gradient}, they should be the same as row names in \code{community}. As specified in the main text, adjacent plots are combined step by step using the specified distance among plots as a constraining factor. In the simplest case, given a set of N plots, for each plot, the first, second, ..., k-th nearest neighbor are determined and a directional species accumulation curve is constructed using the resulting sequence of plots. This procedure is repeated for all plots, generating N directional accumulation curves from which a mean spatially explicit curve is calculated. The resulting curve is thus an intermediate solution between a non-directional accumulation curve and a pure directional curve in which all plots are ordered along a single spatial or environmental gradient.

}
\value{
An object of class \code{data.frame} is returned containing the following statistics:

  - N_SCR: Directional species accumulation curve.

  - N_Exact: Non directional species accumulation curve (classic accumulation curve).

  - Alpha_dir: directional mean number of species in the M plots (for details on the calculation, see Ricotta et al., 2019).

  - Beta_M_dir: Directional beta diversity as a function of sampling effort M ( for details on the calculation, Ricotta et al., 2019).

  - Beta_N_dir: Normalized directional beta diversity.

  - Beta_M: Non-directional beta diversity as a function of sampling effort M.

  - Beta_N: Normalized non-directional beta diversity as a function of sampling effort M.

  - Beta_Autocor: A normalized measure of autocorrelation for directional beta diversity calculated as the normalized difference between directional and non-directional beta.

}
\references{
Ricotta, C., Acosta, A., Bacaro, G., Carboni, M.,
Chiarucci, A., Rocchini, D., Pavoine, S. (2019) Rarefaction of beta diversity. \emph{Ecological Indicators}, \bold{107}, 105606. \doi{10.1016/j.ecolind.2019.105606}.

Chiarucci, A., Bacaro, G., Rocchini, D., Ricotta, C., Palmer, M. W.,
Scheiner, S. M. (2009) Spatially constrained rarefaction: incorporating
the autocorrelated structure of biological communities into sample-based
rarefaction. \emph{Community Ecology}, \bold{10(2)}, 209--214.
}
\author{
Sandrine Pavoine \email{sandrine.pavoine@mnhn.fr}

Giovanni Bacaro \email{gbacaro@units.it}
}

\examples{
require(vegan)
data(mite)
data(mite.xy)
comm_matrix <- mite 

# Spatially-explicit curves can be obtained as follows

spatialdist <- dist(mite.xy) # to calculate the geographic 
# distance between plots, i.e. the Euclidean distance 
# between the coordinates of the plots)

betas <- directionalSAC(comm_matrix, spatialdist) # to calculate directional 
# and non directional beta diversity

plot(1:70, betas$N_Exact, xlab="M", ylab="Species richness", ylim=range(c(betas$N_Exact, 
betas$N_SCR, betas$Alpha, mean(apply(comm_matrix, 1, function(x) length(x[x>0]))))))
points(1:70,rep( mean(apply(comm_matrix, 1, function(x) length(x[x>0]))), 70), pch=2)
points(1:70, betas$N_SCR, pch=3)
points(1:70, betas$Alpha_dir, pch=4)
legend("right", legend=c("Non-directional SAC", 
    "Non-directional alpha diversity", "Directional SAC", 
    "Directional alpha diversity"), pch=1:4) 
# M is the number of plots

plot(1:70, betas$Beta_M, xlab="M", ylab="Beta diversity", 
 ylim=range(c(betas$Beta_M_dir, betas$Beta_M)))
points(1:70, betas$Beta_M_dir, pch=2)
legend("right", legend=c("Non-directional beta", "Directional beta"), pch=1:2) 

plot(2:70, betas$Beta_N[2:70], xlab="M", ylab="Normalized beta diversity", 
    ylim=range(c(betas$Beta_N_dir[2:70], betas$Beta_N[2:70])))
points(2:70, betas$Beta_N_dir[2:70], pch=2)
legend("right", legend=c("Non-directional beta", "Directional beta"), pch=1:2) 

plot(2:70, betas$Beta_Autocor[2:70], xlab="M", 
    ylab="Normalized measure of autocorrelation")
}

\keyword{models}

