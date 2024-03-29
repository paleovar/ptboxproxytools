% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Proxytools_tools.R
\name{paleodata_explained_variance}
\alias{paleodata_explained_variance}
\title{Compute explained variance of an extracted signal}
\usage{
paleodata_explained_variance(
  xin,
  signal_type = "pca",
  signal_components = 1,
  reference_signal = stats::prcomp(xin)$x[, 1]
)
}
\arguments{
\item{xin}{Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate}

\item{signal_type}{Method to extract signals;
Implemented methods:
'pca' Principal component analysis
'ca' Correspondance analysis
'prc' Principal curves (princurve package)
'custom_via_rda' Uses redundancy analysis to compute the explained variance of a specified signal (e.g. explained variance of AP signal for pollen assemblage record)}

\item{signal_components}{Components for which explained variance should be computed}

\item{reference_signal}{Specified signal, if signal_type='custom_via_rda' is selected}
}
\value{
Explained variance (given as a fraction) of the selected components (double between 0 and 1); if multiple components are specified,
the joint explained variance is returned; if input is Proxytibble, output is proxytibble where explained variance is given in the proxy_data column
}
\description{
Compute explained variance of an extracted signal
}
\examples{
# Load ice core example data
library(PTBoxProxydata)
mng <- ProxyDataManager()
icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
# Compute PCA from multivariate zoo's
icecoredata_expl_var <- paleodata_explained_variance(icecoredata, 'pca')
print(paste0("Explained variance by PC1 for EDC time series: ",icecoredata_expl_var$proxy_data[[1]]))
print(paste0("Explained variance by PC1 for EPICA DML time series: ",icecoredata_expl_var$proxy_data[[2]]))
# We can do the same by first computing the PCA and then extract explained variance of the signal via redundancy analysis
icecoredata_pca <- paleodata_signal_extraction(icecoredata, 'pca')
icecoredata_expl_var_rda <- paleodata_explained_variance(icecoredata$proxy_data[[1]], 'custom_via_rda', reference_signal = icecoredata_pca$proxy_data[[1]])
print(icecoredata_expl_var_rda)

}
\seealso{
\link{prcomp} (from `stats`) for principal component analysis

\link{cca} (from `vegan`) for correspondance analysis

\link{principal_curve} (from `princurve`) for principal curves

\link{rda} (from `vegan`) for redundancy analysis
}
