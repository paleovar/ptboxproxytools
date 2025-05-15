# HELP FUNCTIONS
# bestNormalize

### Exported helpers ----

#' Reverse the time axis of `zoo::zoo`
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#'
#' @return zoo with returned time axis
#' @export
#'
#' @examples
#' # Load Monticchio example data from PTBoxProxydata
#' library(PTBoxProxydata)
#' mng <- PTBoxProxydata::ProxyDataManager()
#' monticchiodata <- PTBoxProxydata::load_set(mng,'monticchio_testset',zoo_format = 'zoo', force_file=TRUE)
#' # Change from yrs BP to yrs AP (i.e. forward instead of backward) time axis for proxytibble
#' monticchiodata_rev <- rev_time_axis(monticchiodata)
#' # Change from yrs BP to yrs AP (i.e. forward instead of backward) time axis for zoo
#' monticchiozoo_rev <- rev_time_axis(monticchiodata$proxy_data[[1]])
#' # Plot zoo data
#' plot(monticchiozoo_rev)
#' plot(monticchiodata_rev$proxy_data[[1]])
#'
rev_time_axis <- function(xin) UseMethod('rev_time_axis')

#' @export
rev_time_axis.zoo <- function(xin) {
    xin <- zoo::rev.zoo(xin)
    zoo::index(xin) <- -rev(zoo::index(xin))
    return(xin)
}

#' @export
rev_time_axis.Proxytibble <- function(xin) {
    if (all(class(xin[[Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`rev_time_axis` only implemented for `zoo_format == 'zoo'`")
    return(PTBoxProxydata::apply_proxy(xin, fun = rev_time_axis.zoo))
}


#' Remove NAs and aggregate non-unique timesteps in zoo objects
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param remove_na Flag if NAs should be removed
#' @param aggregation Flag if non-unique timesteps should be merged
#' @param aggregation_fun Function for merging non-unique timesteps
#'
#' @return Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`)
#' @export
#'
#' @examples
#' # Load Monticchio example data from PTBoxProxydata
#' library(PTBoxProxydata)
#' mng <- PTBoxProxydata::ProxyDataManager()
#' monticchiodata <- PTBoxProxydata::load_set(mng,'monticchio_testset',zoo_format = 'zoo', force_file=TRUE)
#' # Clean proxy data time series
#' monticchiodata_cleaned <- clean_timeseries(monticchiodata)
#' # Clean zoo time series
#' monticchiozoo_cleaned <- clean_timeseries(monticchiodata$proxy_data[[1]])
#' # Plot zoo data
#' plot(monticchiozoo_cleaned)
#' plot(monticchiodata_cleaned$proxy_data[[1]])
#'
#' @seealso
#' \link{aggregate} (from `stats`) for merging non-unique timesteps
#'
clean_timeseries <- function(xin,remove_na=TRUE,aggregation=TRUE,aggregation_fun=mean) UseMethod('clean_timeseries')

#' @export
clean_timeseries.zoo <-function(xin,remove_na=TRUE,aggregation=TRUE,aggregation_fun=mean) {
    if (remove_na==TRUE) {
        xin <- zoo::as.zoo(na.omit(xin))
    }
    if (aggregation==TRUE) {
        xin <- zoo:::aggregate.zoo(xin,zoo::index(xin),aggregation_fun)
    }
    return(xin)
}

#' @export
clean_timeseries.Proxytibble <- function(xin,remove_na=TRUE,aggregation=TRUE,aggregation_fun=mean) {
    if (all(class(xin[[Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`clean_timeseries` only implemented for `zoo_format == 'zoo'`")
    return(PTBoxProxydata::apply_proxy(xin, fun = clean_timeseries.zoo,remove_na=remove_na,aggregation=aggregation,aggregation_fun=aggregation_fun))
}


#' Center and rescale irregular time series (NA values are removed from computing means / standard deviations)
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param center logical center to zero mean?
#' @param scale logical normalize to unit standard deviation?
#'
#' @return Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @export
#'
#' @examples
#' # Load Monticchio example data from PTBoxProxydata
#' library(PTBoxProxydata)
#' mng <- PTBoxProxydata::ProxyDataManager()
#' monticchiodata <- PTBoxProxydata::load_set(mng,'monticchio_testset',zoo_format = 'zoo', force_file=TRUE)
#' # Normalize proxy data in proxytibble
#' monticchiodata_normalized <- normalize(monticchiodata)
#'
#' # Normalize single zoo
#' monticchiozoo_normalized <- normalize(monticchiodata$proxy_data[[1]])
#' # Plot zoo data
#' plot(monticchiozoo_normalized)
#' plot(monticchiodata_normalized$proxy_data[[1]])
#'
normalize <- function(xin,center=TRUE,scale=TRUE) UseMethod('normalize')

#' @export
normalize.zoo <- function(xin, center = TRUE, scale = TRUE) {
    if (!("matrix" %in% class(zoo::coredata(xin)))) {
        if (center == TRUE) {
            xin <- xin - mean(xin, na.rm = TRUE)
        }
        if (scale == TRUE) {
            xin <- xin / sd(xin)
        }
        return(xin)
    } else {
        if (center == TRUE) {
            xin <- PTBoxProxydata::zoo_applyfix(xin, function(xx)
                xx - mean(xx, na.rm = TRUE))
        }
        if (scale == TRUE) {
            xin <- PTBoxProxydata::zoo_applyfix(xin, function(xx)
                xx / sd(xx))
        }
        return(xin)
    }
}

#' @export
normalize.Proxyzoo <- function(xin, center = TRUE, scale = TRUE) {
    if (center == TRUE) {
        xin <- PTBoxProxydata::zoo_applyfix.Proxyzoo(xin, function(xx)
            xx - mean(xx, na.rm = TRUE))
    }
    if (scale == TRUE) {
        xin <- PTBoxProxydata::zoo_applyfix.Proxyzoo(xin, function(xx)
            xx / sd(xx))
    }
    return(xin)
}

#' @export
normalize.Proxytibble <- function(xin, center = TRUE, scale = TRUE) {
    if (!any(class(xin[[Proxytibble_colnames_proxy_data()]][[1]]) %in% c('zoo', 'Proxyzoo')))
        stop("`normalize` only implemented for `zoo_format == 'zoo'`")
    return(PTBoxProxydata::apply_proxy(
        xin,
        fun = normalize, #.zoo,
        center = center,
        scale = scale
    ))
}


#' Remove empty rows (i.e. consisting only of zeros) in a matrix
#'
#' @param xin Matrix
#'
#' @return Matrix
#' @export
#'
#' @examples
#' testmatrix <- matrix(c(1,0,1,0),ncol=2)
#' print(testmatrix)
#' print(remove_empty_rows(testmatrix))
#'
remove_empty_rows <- function(xin) {
    nonempty_rows <- which(rowSums(xin)!=0)
    return(list(xin=xin[nonempty_rows,],rows=nonempty_rows))
}

#' Convert factors to numerics (non-numeric values are converted to NA)
#'
#' @param x factor vector
#'
#' @return numerics vector
#' @export
#'
as.numeric.factor <- function(x) {
    as.numeric(levels(x))[x]
}

#' Calculate spatial mean of a climate data field
#'
#' @param lon Vector of longitudes
#' @param lat Vector of latitudes
#' @param clim_field Climate data field (with dimension #lon x #lat)
#'
#' @return Number (spatial mean)
#' @export
#'
#' @examples
#' # Create dummy field
#' lon <- seq(0,90,by=10)
#' lat <- seq(0,70,by=10)
#' clim_field <- matrix(1:80,nrow=10)
#' # Compute spatial mean
#' clim_mean <- spatial_means(lon,lat,clim_field)
#' # Print spatial mean
#' print(clim_mean)
#'
spatial_means <- function(lon,lat,clim_field) {
    spatial_weights <- array(0,dim=c(length(lon),length(lat)))
    for (i in 1:length(lat)) {spatial_weights[,i] <- (cos((lat[i]+mean(diff(lat))/2)*pi/180)+cos((lat[i]-mean(diff(lat))/2)*pi/180))/2}
    spatial_weights[which(is.na(clim_field[1:(length(lon)*length(lat))]))] <- NA
    spatial_weights <- spatial_weights/sum(spatial_weights,na.rm=T)
    return(sum(spatial_weights*clim_field,na.rm=T))
}


#' Replicate data vector of zoo object (time axis does not change)
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param times Number of replications
#'
#' @return Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @export
#'
#' @examples
#' # Load Monticchio example data from PTBoxProxydata
#' library(PTBoxProxydata)
#' mng <- PTBoxProxydata::ProxyDataManager()
#' monticchiodata <- PTBoxProxydata::load_set(mng,'monticchio_testset',zoo_format = 'zoo', force_file=TRUE)
#' # Clean proxy data time series
#' monticchiodata_rep <- rep_zoo(monticchiodata,5)
#' # Clean zoo time series
#' monticchiozoo_rep <- rep_zoo(monticchiodata$proxy_data[[1]],5)
#' # Plot zoo data
#' plot(monticchiodata_rep$proxy_data[[1]])
#' plot(monticchiozoo_rep)
#'
rep_zoo <- function(xin, times=1) UseMethod('rep_zoo')

#' @export
rep_zoo.zoo <- function(xin,times=1) {
    return(zoo::zoo(array(rep(zoo::coredata(xin),times=times),dim=c(length(zoo::coredata(xin)),times)),order.by=zoo::index(xin)))
}

#' @export
rep_zoo.Proxytibble <- function(xin,times = 1) {
    if (all(class(xin[[Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`rep_zoo` only implemented for `zoo_format == 'zoo'`")
    return(PTBoxProxydata::apply_proxy(xin, fun = rep_zoo.zoo, times = times))
}


#' Calculate credible intervals for samples (e.g. output from Monte Carlo or Markov Chain Monte Carlo methods)
#'
#' @param samples Vector of samples
#' @param level Credible level
#' @param simple If TRUE computation through quantiles, otherwise smallest interval containing the respective probability mass is computed
#'
#' @return Vector with upper and lower credible interval
#' @export
#'
#' @examples
#' testsamples <- rnorm(1000,0,1)
#' level <- 0.95
#' simple <- TRUE
#' testci <- calculate_ci(testsamples, level,simple)
#' print(testci)
#'
calculate_ci <- function(samples,level,simple=TRUE) {
    if (simple==TRUE) {
        ci_min <- quantile(samples,(1-level)/2,na.rm=TRUE)
        ci_max <- quantile(samples,(1+level)/2,na.rm=TRUE)
    } else {
        num_samples <- length(samples)
        samples <- sort(samples)
        tmp <- c(samples[1],samples[ceiling(level*num_samples)])
        j <- 1
        while ((ceiling(level*num_samples)+j) < num_samples) {
            if ((samples[ceiling(level*num_samples)+j]-samples[j+1]) < (tmp[2]-tmp[1])) {
                tmp <- c(samples[j],samples[ceiling(level*num_samples)+j])
            }
            j <- j+1
        }
        ci_min <- tmp[1]
        ci_max <- tmp[2]
    }
    return(c(ci_min,ci_max))
}

#' Binning
#'
#' @param xin Input zoo
#' @param xout Output time axis
#' @param bin_width Width of bins if interpolation method is "binning". Defaults to the mean sample resolution (no variable bin sizes are supported at the moment)
#' @param binning_function How should values within one bin be averaged? Default is "mean"
#'
#' @returns Zoo with binned values at the provided interpolation dates
#' @export
binning <- function(xin, xout, bin_width = mean(diff(xout)), binning_function = mean) {
    xout <- zoo::zoo(sapply(xout, function(x) binning_function(paleodata_windowing(xin,x-bin_width/2,x+bin_width/2),na.rm=TRUE)),order.by=xout)
    return(xout)
}


#' Anti-aliasing linear interpolation method from Baudouin et al. 2025 ("bwr25")
#'
# "bwr25" combines linear interpolation to high resoluton to reduce aliasing and binning the high-resolution output to the desired resolution
# At the moment, it only works with an equally spaced xout
# It uses a fixed 10 times xout resolution for anti-aliasing
#'
#' @param xin Input zoo
#' @param xout Output time axis (only equally spaced xout are supported at the moment)
#'
#' @returns Zoo with interpolated values
#' @export
interp_bwr25 <- function(xin,xout) {
    n <- length(xout)
    xout_highres <- seq(xout[1]*1.45 - 0.45*xout[2], xout[n]*1.45 - 0.45*xout[n-1], (xout[2]-xout[1])/10)
    xout_highres <- zoo::zoo(approx(x = zoo::index(xin), y = zoo::coredata(xin), xout = xout_highres)$y, order.by=xout_highres)
    sel <- seq(1, length(xout_highres), 10)
    y <- rep(0, length(xout))
    for (j in 1:10) {
        y <- y + zoo::coredata(xout_highres)[sel+j-1]
    }
    return(zoo::zoo(y/10, order.by=xout))
}

#' Gaussian kernel interpolation and smoothing
#'
#' @param xin Input too
#' @param xout Output time axis
#' @param smooth_scale Smoothing scale of the Gaussian kernel
#' @param pass Gain at the smoothing scale
#'
#' @returns Zoo with smoothed and interpolated values at the provided interpolation dates
#' @export
#'
#' @seealso
#' \link{ksmooth} (from `stats`) is used for the smoothing and interpolation, but with modified definition of the kernel width to align better with timescale point of view
#'
gkinterp <- function(xin, xout, smooth_scale, pass = 0.5) {
    xout <- zoo::zoo(ksmooth(zoo::index(xin), zoo::coredata(xin), bandwidth = smooth_scale/pi*sqrt(log(1/pass)/2)*4*qnorm(0.75, 0, 1), x.points = xout, kernel = "normal")$y,order.by=xout)
    # sqrt(log(1/pass)/2)/pi comes from the fourier transform of the gaussian kernel: the gain is "pass" when apply at the smoothing scale (cutoff period)
    # 4*qnorm(0.75, 0, 1) comes from the definition of ksmooth
    return(xout)
}


#' Remove samples in interpolated zoo that are far away from original samples
#'
#' @param xin_raw Original proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param xin_interp Zoo with interpolated values (all timeseries from xin_raw need to be interpolated to the same time axis), xin_interp can have dimension 1 (time), 2 (time x records or sites x time), or 3 (sites x time x records)
#' @param max_dist Cutoff distance from nearest raw sample beyond which values in xin_interp are set to NA
#'
#' @return zoo of same dimension as xin_interp
#' @export
#'
remove_extrapolated_samples <- function(xin_raw, xin_interp,max_dist=5*mean(diff(zoo::index(xin_interp)))) UseMethod('remove_extrapolated_samples')
#' @export
remove_extrapolated_samples.zoo <- function(xin_raw, xin_interp,max_dist=5*mean(diff(zoo::index(xin_interp)))) {
    for (j in 1:length(zoo::index(xin_interp))) {
        if (min(abs(zoo::index(xin_raw) - zoo::index(xin_interp)[j])) > max_dist) {
            if (!("matrix" %in% class(zoo::coredata(xin_interp)))) {
                xin_interp[j] <- NA
            } else {
                xin_interp[j,] <- NA
            }
        }
    }
    return(xin_interp)
}
#' @export
remove_extrapolated_samples.Proxytibble <- function(xin_raw, xin_interp,max_dist=5*mean(diff(zoo::index(xin_interp)))) {
    for (i in 1:dim(proxytibble)[1]) {
        for (j in 1:length(zoo::index(xin_interp))) {
            if (min(abs(zoo::index(xin_raw$proxy_data[[i]]) - zoo::index(xin_interp)[j])) > max_dist) {
                if (length(dim(xin_interp)) == 2) {
                    xin_interp[i,j] <- NA
                } else {
                    xin_interp[i,j,] <- NA
                }
            }
        }
    }
    return(xin_interp)
}

#' Efficient resampling of elements for uncertainty quantification with bootstrapping
#'
#' @param x Input data
#' @param size Number of output elements, defaults to length(x)
#' @param replace Replace entries?, defaults to TRUE
#' @param ... Other parameters for sample.int
#'
#' @returns Shuffled version of input data
#' @export
resample <- function(x, size=length(x), replace=TRUE, ...) {
    return(x[sample.int(length(x), size=size, replace=replace, ...)])
}

### Private helpers ----

# Extract the non-negative part of a vector
nonneg <- function(x) {
    if (class(x) != "matrix") {
        return(sapply(x,function(y) max(y,0)))
    } else {
        return(apply(x,c(1,2),function(y) max(y,0)))
    }
}


# Extract the non-positive part of x
nonpos <- function(x) {
    if (class(x) != "matrix") {
        retrun(sapply(x,function(y) min(y,0)))
    } else {
        return(apply(x,c(1,2),function(y) min(y,0)))
    }
}


# Quantile transformation
qtrafo <- function(x,weighted=FALSE) {
    # Quantile transformation
    trafo <- bestNormalize::orderNorm(x)$x.t
    # Potentially rescale by standard deviation to preserve variance
    if (weighted == TRUE) {
        trafo <- trafo * sd(x)
    }
    return(trafo)
}

