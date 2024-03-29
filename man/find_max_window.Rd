% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Proxytools_tools.R
\name{find_max_window}
\alias{find_max_window}
\title{Extract maximal window}
\usage{
find_max_window(
  xin,
  t_min = min(zoo::index(xin)),
  t_max = max(zoo::index(xin)),
  min_res = t_max - t_min,
  max_step = t_max - t_min,
  min_length = 0
)
}
\arguments{
\item{xin}{Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate}

\item{t_min}{Lower time limit}

\item{t_max}{Upper time limit}

\item{min_res}{Minimal mean inter-sample range}

\item{max_step}{Maximal time step}

\item{min_length}{Minimal length of signal}
}
\value{
Proxytibble with proxy data in `zoo::zoo` format or irregular time series object (`zoo::zoo`), with data restricted to the found values, or NA if no window fulfills criteria
}
\description{
Extract maximal window fulfilling properties: minimal/maximal dates, minimal resolution (mean time step), maximal time step (to avoid large gaps), minimal length
}
\examples{
# Load Monticchio example data from PTBoxProxydata
library(PTBoxProxydata)
mng <- PTBoxProxydata::ProxyDataManager()
monticchiodata <- PTBoxProxydata::load_set(mng,'monticchio_testset',zoo_format = 'zoo', force_file=TRUE)
# Find the longest period in the Monticchio record in interval 30-60ka BP with mean temporal resolution below 250yrs, a maximal step length of 500yrs, and at least 20 samples
monticchiodata_highres <- find_max_window(monticchiodata,30000,60000,250,500,20)
# Plot restricted data
plot(monticchiodata_highres$proxy_data[[1]])
# Or directly for the AP zoo:
plot(find_max_window(monticchiodata$proxy_data[[1]],30000,60000,250,500,20))

}
