\name{duneFVG}
\alias{duneFVG}
\docType{data}
\title{
Abundance of Vascular Plants Data in Northen Adriatic Dunes
}
\description{
This data set is a list which contains information about the abundance of native and alien vascular plants sampled in 128 16 m^{2} squared plots in the Northern Adriatic Dunes. See Tordoni et al. (2019) for details.
}
\usage{
data("duneFVG")
}

\format{
  A list that contains 3 objects:
  \describe{
    \item{\code{total}}{a numeric dataframe with plots as rows and species as columns. It contains the abundances of all the species collected.}
    \item{\code{native}}{a numeric dataframe with plots as rows and species as columns. It contains the abundances of the native species collected.}
    \item{\code{alien}}{a numeric dataframe with plots as rows and species as columns. It contains the abundances of the alien species collected.}
  }
  There are two linked lists to \code{duneFVG}, one with the geographic coordinates of the plots (\code{\link{duneFVG.xy}}) and one with the functional traits of the species (\code{\link{duneFVG.tr8}}).

}
\source{
Tordoni, E., Petruzzellis, F., Nardini, A., Savi, T., Bacaro, G. (2019) Make it simpler: alien species decrease functional diversity of coastal plant communities. \emph{Journal of Vegetation Science}, \bold{30}, 498--509.
}

\examples{
data(duneFVG)
}
\keyword{datasets}
