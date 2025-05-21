### (Statistical) tools for paleodata processing
# princurve
# VGAM
# vegan
# PaleoSpec
# nest
# bestNormalize

#' Extract all samples in a given time period from an irregular time series
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param start_date Lower end of time window (included in the output `zoo::zoo`)
#' @param end_date Upper end of time window (included in the output `zoo::zoo`)
#'
#' @return Proxytibble with proxy data in `zoo::zoo` format or irregular time series object (`zoo::zoo`) containing only the samples in the defined window
#' @export
#'
#' @examples
#' # Load Monticchio example data from PTBoxProxydata
#' library(PTBoxProxydata)
#' mng <- PTBoxProxydata::ProxyDataManager()
#' monticchiodata <- PTBoxProxydata::load_set(mng,'monticchio_testset',zoo_format = 'zoo', force_file=TRUE)
#' # Reduce to samples between 30ka BP and 60ka BP
#' monticchiodata <- paleodata_windowing(monticchiodata,30000,60000)
#' # Extract zoo with samples from 30ka BP to 40ka BP
#' monticchiozoo <- paleodata_windowing(monticchiodata$proxy_data[[1]],30000,40000)
#' # Plot zoo data
#' plot(monticchiozoo)
#'
paleodata_windowing <- function(xin,start_date,end_date) UseMethod('paleodata_windowing')

#' @export
paleodata_windowing.zoo <- function(xin, start_date, end_date) {
    sel <-
        which(zoo::index(xin) >= start_date & zoo::index(xin) <= end_date)
    return(zoo::zoo(xin[sel, drop = FALSE],
                    order.by = zoo::index(xin)[sel]))
}

#' @export
paleodata_windowing.Proxytibble <-
    function(xin, start_date, end_date) {
        if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
            stop("`paleodata_windowing` only implemented for `zoo_format == 'zoo'`")
        return(
            PTBoxProxydata::apply_proxy(
                xin,
                fun = paleodata_windowing.zoo,
                start_date = start_date,
                end_date = end_date
            )
        )
    }


#' Interpolation functions for irregular time series
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param xout Dates of the interpolated time series
#' @param method Type of interpolation, either 'linear' (piecewise linear interpolation), 'nn' (nearest neighbor), spline' (spline interpolation), 'lh14' (using the interpolation method from Laepple and Huybers 2014, Rehfeld et al. 2018, implemented in 'PaleoSpec'), 'binning' (aggregation in equally-sized bins), 'loess' (loess interpolation), 'bwr25' (binning with initial linear interpolation ot high resolution to mitigate aliasing, from Baudouin et al. 2025), 'gk' (Gaussian kernel smoothing and interpolation)
#' @param remove_na Flag if NAs should be removed after the interpolation
#' @param aggregation Flag if non-unique timesteps should be merged after the interpolation
#' @param aggregation_fun Function for merging non-unique timesteps
#' @param lh14_lowpass Scaling factor for the lowpass frequency. Default is 1.2 (same as in PaleoSpec::MakeEquidistant)
#' @param lh14_length Scaling factor for the Length of the filter. Default is 5 (same as in PaleoSpec::MakeEquidistant)
#' @param bin_width Width of bins if interpolation method is "binning". Defaults to the mean sample resolution (no variable bin sizes are supported at the moment)
#' @param binning_function How should values within one bin be averaged? Default is "mean"
#' @param loess_span "span" parameter in loess fitting, controls degree of smoothing
#' @param gk_antialiasing Should linear interpolation to higher resolution be applied prior to smoothing to avoid aliasing. Default is "TRUE"
#' @param gk_smooth_scale Smoothing scale of the Gaussian kernel
#' @param gk_pass Gain at the smoothing scale. Default is 0.5
#'
#' @return Proxytibble with interpolated proxy data in `zoo::zoo` format or interpolated irregular time series object (`zoo::zoo`)
#' @export
#'
#' @examples
#' # Load ice core example data
#' library(PTBoxProxydata)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # Interpolate to a regular 100yr resolution between 30ka BP and 60ka BP using spline interpolation
#' icecoredata_interpolated <- paleodata_interpolation(icecoredata,"spline",seq(30000,60000,by=100))
#' # Plot datasets
#' plot(icecoredata_interpolated$proxy_data[[1]])
#' plot(icecoredata_interpolated$proxy_data[[2]])
#' # Interpolate to a regular 1000yr resolution between 30ka BP and 60ka BP using interpolation method from PaleoSpec package (named "lh14" after Laepple and Huybers 2014, this method is optimized for computation of spectral densities)
#' edc_data_interpolated <- paleodata_interpolation(icecoredata$proxy_data[[1]],"lh14",seq(30000,60000,by=100))
#' # Plot zoo data
#' plot(edc_data_interpolated)
#'
#' @seealso
#' \link{optim} (from `stats`) for linear and nn interpolation
#'
#' \link{spline} (from `stats`) for spline interpolation
#'
#' \link{MakeEquidistant} (from `PaleoSpec`) for 'lh14' interpolation (optimized for computation of spectral densities from irregular time series)
#'
#' \link{loess} (from `stats`) for loess interpolation
#'
#' \link{ksmooth} (from `stats`) for Gaussian kernel interpolation
#'
paleodata_interpolation <- function(xin,
                                    xout,
                                    method,
                                    remove_na = TRUE,
                                    aggregation = TRUE,
                                    aggregation_fun = mean,
                                    remove_extrapolated_values = FALSE,
                                    max_dist = 5*mean(diff(zoo::index(xout))),
                                    lh14_lowpass = 1.2,
                                    lh14_length = 5,
                                    bin_width = mean(diff(xout)),
                                    binning_function = mean,
                                    loess_span = 0.25,
                                    gk_antialiasing = TRUE,
                                    gk_smooth_scale = NULL,
                                    gk_pass = 0.5)
    UseMethod('paleodata_interpolation')

#' @export
paleodata_interpolation.zoo <-
    function(xin,
             xout,
             method,
             remove_na = TRUE,
             aggregation = TRUE,
             aggregation_fun = mean,
             remove_extrapolated_values = FALSE,
             max_dist = 5*mean(diff(zoo::index(xout))),
             lh14_lowpass = 1.2,
             lh14_length = 5,
             bin_width = mean(diff(xout)),
             binning_function = mean,
             loess_span = 0.25,
             gk_antialiasing = TRUE,
             gk_smooth_scale = NULL,
             gk_pass = 0.5) {
        if (!method %in% c("linear","nn","spline","lh14","binning","loess","bwr25","gk")) {
            stop("`method` not supported")
        }
        if (method == "linear") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                cln <- colnames(xin)
                xout <- zoo::zoo(
                    approx(zoo::index(xin), zoo::coredata(xin), xout=xout, method = "linear")$y,
                    order.by = xout
                )
                colnames(xout) <- cln
            } else {
                xout <- PTBoxProxydata::zoo_apply(xin,
                                                  function(xx) {
                                                      xo <- zoo::zoo(
                                                          approx(zoo::index(xx), zoo::coredata(xx), xout=xout, method = "linear")$y,
                                                          order.by = xout)
                                                      return(xo)
                                                  },
                                                  out_index = xout)
            }
        }
        if (method == "nn") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                cln <- colnames(xin)
                xout <- zoo::zoo(
                    approx(zoo::index(xin), zoo::coredata(xin), xout=xout, method = "constant")$y,
                    order.by = xout
                )
                colnames(xout) <- cln
            } else {
                xout <- PTBoxProxydata::zoo_apply(xin,
                                                  function(xx) {
                                                      xo <- zoo::zoo(
                                                          approx(zoo::index(xx), zoo::coredata(xx), xout=xout, method = "constant")$y,
                                                          order.by = xout)
                                                      return(xo)
                                                  },
                                                  out_index = xout)
            }
        }
        if (method == "lh14") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                cln <- colnames(xin)
                xout <- zoo::zoo(
                    PaleoSpec::MakeEquidistant(zoo::index(xin),
                                               xin,
                                               time.target = xout,
                                               k=lh14_length,
                                               kf=lh14_lowpass),
                    order.by = xout
                )
                colnames(xout) <- cln
            } else {
                xout <- PTBoxProxydata::zoo_apply(xin,
                                 function(xx) {
                                     xo <- zoo::zoo(
                                         PaleoSpec::MakeEquidistant(zoo::index(xx),
                                                                    xx,
                                                                    time.target = xout,
                                                                    k=lh14_length,
                                                                    kf=lh14_lowpass),
                                         order.by = xout)
                                     return(xo)
                                 },
                                 out_index = xout)
            }
        }
        if (method == "spline") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                xout <- zoo::zoo(
                    spline(zoo::index(xin),
                           xin,
                           xout = xout)$y,
                    order.by = xout
                )
            } else {
                xout <- PTBoxProxydata::zoo_apply(xin,
                                                  function(xx)
                                                      spline(zoo::index(xx),
                                                             xx,
                                                             xout = xout)$y,
                                                  out_index = xout)
            }
        }
        if (method == "binning") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                xout <- binning(xin,
                            xout = xout,
                            bin_width = bin_width,
                            binning_function = binning_function)
            } else {
                xout <- PTBoxProxydata::zoo_apply(xin,
                                                  function(xx)
                                                      binning(xx,
                                                              xout = xout,
                                                              bin_width = bin_width,
                                                              binning_function = binning_function),
                                                  out_index = xout)
            }
        }
        if (method == "loess") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                xout <- zoo::zoo(predict(loess(y ~ x,
                                  data.frame(x = zoo::index(xin), y = zoo::coredata(xin)),
                                  span = loess_span, degree = 2),
                            data.frame(x=xout)),
                    order.by = xout)
            } else {
                xout <- PTBoxProxydata::zoo_apply(xin,function(xx)
                                                             predict(loess(y ~ x,
                                                                            data.frame(x = zoo::index(xx), y = zoo::coredata(xx)),
                                                                            span = loess_span, degree = 2),
                                                                            data.frame(x=xout)),
                                                             out_index = xout)
            }
        }
        if (method == "bwr25") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                xout <- interp_bwr25(xin,xout)
            } else {
                xout <- PTBoxProxydata::zoo_apply(xin,function(xx)
                                                            interp_bwr25(xx,xout),
                                                            out_index = xout)
            }
        }
        if (method == "gk") {
            if (gk_antialiasing == TRUE) {
                xin <- paleodata_interpolation(xin,
                                               seq(min(xout)-2*gk_smooth_scale, max(xout)+2*gk_smooth_scale, by=min(diff(xout))/10),
                                               method="linear",remove_na==TRUE)
            }
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                xout <- gkinterp(xin, xout, smooth_scale = gk_smooth_scale, pass = gk_pass)
            } else {
                xout <- PTBoxProxydata::zoo_apply(xin,function(xx)
                                                            gkinterp(xx, xout = xout, smooth_scale = gk_smooth_scale, pass = gk_pass),
                                                            out_index = xout)
            }
        }
        if (remove_extrapolated_values == TRUE) {
            xout <- remove_extrapolated_samples(xin,xout,max_dist=max_dist)
        }
        xout <- clean_timeseries(xout, remove_na=remove_na, aggregation=aggregation, aggregation_fun=aggregation_fun)
        return(xout)
    }

#' @export
paleodata_interpolation.Proxytibble <-
    function(xin,
             xout,
             method,
             remove_na = TRUE,
             aggregation = TRUE,
             aggregation_fun = mean,
             remove_extrapolated_values = FALSE,
             max_dist = 5*mean(diff(zoo::index(xout))),
             lh14_lowpass = 1.2,
             lh14_length = 5,
             bin_width = mean(diff(xout)),
             binning_function = mean,
             loess_span = 0.25,
             gk_antialiasing = TRUE,
             gk_smooth_scale = NULL,
             gk_pass = 0.5) {
        if (!method %in% c("linear","nn","spline","lh14","binning","loess","bwr25","gk")) {
            stop("`method` not supported")
        }
        if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
            stop("`paleodata_interpolation` only implemented for `zoo_format == 'zoo'`")
        return(
            PTBoxProxydata::apply_proxy(
                xin,
                fun = paleodata_interpolation.zoo,
                method = method,
                xout = xout,
                remove_na = remove_na,
                aggregation = aggregation,
                aggregation_fun = aggregation_fun,
                remove_extrapolated_values = remove_extrapolated_values,
                max_dist = max_dist,
                lh14_lowpass = lh14_lowpass,
                lh14_length = lh14_length,
                bin_width = bin_width,
                binning_function = binning_function,
                loess_span = loess_span,
                gk_antialiasing = gk_antialiasing,
                gk_smooth_scale = gk_smooth_scale,
                gk_pass = gk_pass
            )
        )
    }


#' Gaussian filtering of irregular time series (using functions from `nest` package)
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param filter_type Type of filter, either 'detrend' (high pass), 'smooth' (low pass), or 'bandpass' (high and low pass)
#' @param filter_scales Upper and lower cut-off periods for bandpass filtering
#' @param detr_scale Cut-off period for detrending
#' @param smooth_scale Cut-off period for smoothing
#'
#' @return Proxytibble with filtered proxy data in `zoo::zoo` format, or filtered irregular time series object (`zoo::zoo`)
#' @export
#'
#' @examples
#' # Load ice core example data
#' library(PTBoxProxydata)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # Detrend the data with 10kyr cutoff timescale
#' icecoredata_detrended <- paleodata_filtering(icecoredata, 'detrend', detr_scale=10000)
#' # Smooth the data with 10kyr cutoff timescale
#' icecoredata_smoothed <- paleodata_filtering(icecoredata, 'smooth', smooth_scale=10000)
#' # Apply bandpass filter for timescales from 1kyr to 10kyr
#' icecoredata_filtered <- paleodata_filtering(icecoredata, 'bandpass', filter_scales=data.frame(lower=1000,upper=10000))
#' # Plot results
#' plot(icecoredata_detrended$proxy_data[[1]])
#' plot(icecoredata_smoothed$proxy_data[[1]])
#' plot(icecoredata_filtered$proxy_data[[1]])
#'
#' @seealso
#' \link{gaussbandpass} (from `nest`) for specifics of the Gaussian smoothing / detrending / bandpass filtering
#'
paleodata_filtering <- function(xin,filter_type,filter_scales=NULL,detr_scale=NULL,smooth_scale=NULL) UseMethod('paleodata_filtering')

#' @export
paleodata_filtering.zoo <-
    function(xin,
             filter_type,
             filter_scales = NULL,
             detr_scale = NULL,
             smooth_scale = NULL) {
        if (!filter_type %in% c("smooth","detrend","bandpass")) {
            stop("`filter_type` not supported")
        }
        if (filter_type == "bandpass") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                return(
                    nest::gaussbandpass(xin, per1 = filter_scales$lower, per2 = filter_scales$upper)$filt
                )
            } else {
                return(PTBoxProxydata::zoo_apply(xin, function(xx)
                    nest::gaussbandpass(xx, per1 = filter_scales$lower, per2 = filter_scales$upper)$filt))
            }
        }
        if (filter_type == "detrend") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                return(nest::gaussdetr(xin, tsc.in = detr_scale)$detr)
            } else {
                return(PTBoxProxydata::zoo_apply(xin, function(xx)
                    nest::gaussdetr(xx, tsc.in = detr_scale)$detr))
            }
        }
        if (filter_type == "smooth") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                return(nest::gaussdetr(xin, tsc.in = smooth_scale)$Xsmooth)
            } else {
                return(PTBoxProxydata::zoo_apply(xin, function(xx)
                    nest::gaussdetr(xx, tsc.in = smooth_scale)$Xsmooth))
            }
        }
    }

#' @export
paleodata_filtering.Proxytibble <-
    function(xin,
             filter_type,
             filter_scales = NULL,
             detr_scale = NULL,
             smooth_scale = NULL) {
        if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
            stop("`paleodata_filtering` only implemented for `zoo_format == 'zoo'`")
        return(
            PTBoxProxydata::apply_proxy(
                xin,
                fun = paleodata_filtering.zoo,
                filter_type = filter_type,
                filter_scales = filter_scales,
                detr_scale = detr_scale,
                smooth_scale = smooth_scale
            )
        )
    }



#' Transform data in an irregular time series object
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param transformation_type Type of transformation;
#' Implemented methods:
#' 'logit' (logit transformation),
#' invlogit' (inverse logit transformation),
#' 'probit' (probit transformation),
#' 'invprobit' (inverse probit transformation),
#' 'sqrt' (square-root transformation)
#' 'quadratic' (quadratic transformation)
#' 'log' (logarithmic transformation)
#' 'exp' (exponential transformation)
#' 'quantile' (quantile transformation to transform the time series into a standard Gaussian distribution)
#' 'weighted quantile' (quantile transformation with re-scaling to preserve the variance in the original time series, i.e. results in a Gaussian distribution with mean zero and standard deviation according to the orginal standard deviation)
#' 'nonneg' (Set negative values to zero)
#' 'normalize' (center and standardize the time series)
#' 'standardize' (standardize the time series, i.e. divide by its standard deviation)
#' 'center' (center the time series, i.e. subtract its mean)
#'
#' @return Proxytibble with transformed proxy data in `zoo::zoo` format, or transformed irregular time series object (`zoo::zoo`)
#' @export
#'
#' @examples
#' # Load Monticchio example data from PTBoxProxydata
#' library(PTBoxProxydata)
#' mng <- PTBoxProxydata::ProxyDataManager()
#' monticchiodata <- PTBoxProxydata::load_set(mng,'monticchio_testset',zoo_format = 'zoo', force_file=TRUE)
#' # Transform the arboreal pollen signal using logit, square root, and quantile transformations
#' monticchio_logit <- paleodata_transformation(monticchiodata$proxy_data[[1]]/100,transformation_type="logit")
#' monticchio_sqrt <- paleodata_transformation(monticchiodata$proxy_data[[1]],transformation_type="sqrt")
#' monticchio_quantile <- paleodata_transformation(monticchiodata$proxy_data[[1]],transformation_type="quantile")
#' # Plot zoo data
#' plot(monticchio_logit)
#' plot(monticchio_sqrt)
#' plot(monticchio_quantile)
#'
#' @seealso
#' \link{logitlink} (from `VGAM`) for 'logit' and 'invlogit' transformations
#'
#' \link{probitlink} (from `VGAM`) for 'probit' and 'invprobit' transformations
#'
#' \link{orderNorm} (from `bestNormalize`) for 'quantile' transformation
#'
#' \link{normalize} (from `PTBoxProxytools`) for 'normalize', 'standardize', and 'center' transformations
#'
paleodata_transformation <- function(xin, transformation_type="logit") UseMethod('paleodata_transformation')

#' @export
paleodata_transformation.zoo <-
    function(xin, transformation_type = "logit") {
        if (!transformation_type %in% c("logit","invlogit","probit","invprobit","sqrt","quadratic","log","exp","quantile","weighted_quantile","nonneg","normalize","standardize","center")) {
            stop("`method` not supported")
        }
        if (transformation_type == "logit") {
            data_trafo <- VGAM::logitlink(xin, inverse = FALSE)
        }
        if (transformation_type == "invlogit") {
            data_trafo <- VGAM::logitlink(xin, inverse = TRUE)
        }
        if (transformation_type == "probit") {
            data_trafo <- VGAM::probitlink(xin, inverse = FALSE)
        }
        if (transformation_type == "invprobit") {
            data_trafo <- VGAM::probitlink(xin, inverse = TRUE)
        }
        if (transformation_type == "sqrt") {
            data_trafo <- sqrt(xin)
        }
        if (transformation_type == "quadratic") {
            data_trafo <- xin ^ 2
        }
        if (transformation_type == "log") {
            data_trafo <- log(xin)
        }
        if (transformation_type == "exp") {
            data_trafo <- exp(xin)
        }
        if (transformation_type == "quantile") {
            if (!requireNamespace("bestNormalize", quietly = TRUE)) stop("suggested package `bestNormalize` required for `transformation_type = 'quantile'`")
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                data_trafo <- qtrafo(as.numeric(xin), weighted = FALSE)
            } else {
                data_trafo <-
                    apply(xin, 2, function(xx)
                        qtrafo(as.numeric(xx), weighted = FALSE))
            }
        }
        if (transformation_type == "weighted_quantile") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                data_trafo <- qtrafo(as.numeric(xin), weighted = TRUE)
            } else {
                data_trafo <-
                    apply(xin, 2, function(xx)
                        qtrafo(as.numeric(xx), weighted = TRUE))
            }
        }
        if (transformation_type == "nonneg") {
            data_trafo <- nonneg(zoo::coredata(xin))
        }
        if (transformation_type == "normalize") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                data_trafo <- normalize(xin)
            } else {
                data_trafo <- PTBoxProxydata::zoo_apply(xin, normalize)
            }
        }
        if (transformation_type == "standardize") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                data_trafo <- normalize(xin, center = FALSE, scale = TRUE)
            } else {
                data_trafo <-
                    PTBoxProxydata::zoo_apply(xin, normalize, center = FALSE, scale = TRUE)
            }
        }
        if (transformation_type == "center") {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                data_trafo <- normalize(xin, center = TRUE, scale = FALSE)
            } else {
                data_trafo <-
                    PTBoxProxydata::zoo_apply(xin, normalize, center = TRUE, scale = FALSE)
            }
        }
        return(clean_timeseries(zoo::zoo(data_trafo, order.by = zoo::index(xin))))
    }

#' @export
paleodata_transformation.Proxytibble <-
    function(xin, transformation_type = "logit") {
        if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
            stop("`paleodata_transformation` only implemented for `zoo_format == 'zoo'`")
        return(
            PTBoxProxydata::apply_proxy(xin,
                        fun = paleodata_transformation.zoo,
                        transformation_type = transformation_type)
        )
    }



#' Extract low-dimensional signals from multivariate datasets (e.g. pollen assemblages)
#'
#' @param xin Proxytibble with multivariate proxy data in `zoo::zoo` format, or a multivariate irregular time series object (`zoo::zoo`)
#' @param signal_type Method to extract signals;
#' Implemented methods:
#' 'pca' Principal component analysis
#' 'dca' Detrended correspondance analysis
#' 'ca' Correspondance analysis
#' 'prc' Principal curves
#' @param signal_components Components that should be returned (e.g. c(2,3) for 'pca' returns the second and third principal components)
#'
#' @return Proxytibble with extracted signals in proxy data, or irregular time series object (`zoo::zoo`) containing the extracted signals
#' @export
#'
#' @examples
#' # Load ice core example data
#' library(PTBoxProxydata)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # Compute PCA from multivariate zoo's
#' icecoredata_pca <- paleodata_signal_extraction(icecoredata, 'pca')
#' plot(icecoredata_pca$proxy_data[[1]])
#' plot(icecoredata_pca$proxy_data[[2]])
#' # Compute principal curve from multivariate zoo's
#' icecoredata_prc <- paleodata_signal_extraction(icecoredata, 'prc')
#' plot(icecoredata_prc$proxy_data[[1]])
#' plot(icecoredata_prc$proxy_data[[2]])
#'
#' @seealso
#' \link{prcomp} (from `stats`) for principal component analysis
#'
#' \link{decorana} (from `vegan`) for detrended correspondance analysis
#'
#' \link{cca} (from `vegan`) for correspondance analysis
#'
#' \link{principal_curve} (from `princurve`) for principal curves
#'
paleodata_signal_extraction <- function(xin,signal_type,signal_components=1) UseMethod('paleodata_signal_extraction')

#' @export
paleodata_signal_extraction.zoo <- function(xin,signal_type,signal_components=1) {
    if (!signal_type %in% c("pca","dca","ca","prc")) {
        stop("`signal_type` not supported")
    }
    if (signal_type == "pca") {
        signal <- stats::prcomp(xin)$x[,signal_components]
    }
    if (signal_type == "dca") {
        if (!requireNamespace("vegan", quietly = TRUE)) stop("suggested package `vegan` required for `signal_type = 'dca'`")
        signal <- vegan::decorana(xin)$rproj[,signal_components]
    }
    if (signal_type == "ca") {
        if (!requireNamespace("vegan", quietly = TRUE)) stop("suggested package `vegan` required for `signal_type = 'ca'`")
        signal <- vegan::cca(xin)$CA$u[,signal_components]
    }
    if (signal_type == "prc") {
        # Principal curve algorithm from princurve (original)
        if (!requireNamespace("princurve", quietly = TRUE)) stop("suggested package `princurve` required for `signal_type = 'princurve'`")
        signal <- princurve::principal_curve(as.matrix(xin))$lambda
    }
    return(zoo::zoo(data.frame(signal),order.by=zoo::index(xin)))
}

#' @export
paleodata_signal_extraction.Proxytibble <- function(xin,signal_type,signal_components=1) {
    if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`paleodata_signal_extraction` only implemented for `zoo_format == 'zoo'`")
    return(
        PTBoxProxydata::apply_proxy(xin,
                    fun = paleodata_signal_extraction.zoo,
                    signal_type = signal_type,
                    signal_components = signal_components)
    )
}

#' Variance computation by integration of the spectrum, based on PaleoSpec::GetVarFromSpectra()
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param freq.start Vector containing the start frequencies of the fitting interval(s)
#' @param freq.end Vector containing the start frequencies of the fitting interval(s)
#' @param target denotes  "raw" or "logsmooth", default "raw"
#' @param interpolation Logical: should interpolation be applied?
#' @param method Parameters for interpolation (see paleodata_interpolation)
#' @param xout Parameters for interpolation (see paleodata_interpolation). If NULL the interpolation is set to the mean temporal resolution of each proxy time series.
#' @param transformation Logical: should transformation be applied?
#' @param transformation_type Parameters for transformation (see paleodata_transformation)
#' @param detrend Logical: should time series be linearly detrend before computing the spectrum
#' @param df_log Width of the log-smoother in log units
#' @param bLog TRUE: average in the log space of the power, FALSE: arithmetic average
#'
#' @return List with estimated variance
#' @export
#'
#' @examples
#' # Load ice core example data
#' library(PTBoxProxydata)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # compute variance
#' testout <- icecoredata %>% paleodata_varfromspec(.,  target="raw", freq.start = 5e-5, freq.end = 5e-4)
#' print(testout)
#'
#' @seealso
#' \link{GetVarFromSpectra} (from `PaleoSpec`) for regular time series version
#'
#' \link{paleodata_spectrum} (from `PTBoxProxytools`) for computation of spectra
#'
paleodata_varfromspec <-
    function(xin,
             freq.start = NULL,
             freq.end = NULL,
             target = "raw",
             interpolation = TRUE,
             method = "lh14",
             xout = NULL,
             transformation = FALSE,
             transformation_type = "normalize",
             detrend = TRUE,
             df_log = 0.05,
             bLog = FALSE)
        UseMethod('paleodata_varfromspec')


#' @export
paleodata_varfromspec.zoo <-
    function(xin,
             freq.start = NULL,
             freq.end = NULL,
             target = "raw",
             interpolation = TRUE,
             method = "lh14",
             xout = NULL,
             transformation = FALSE,
             transformation_type = "normalize",
             detrend = TRUE,
             df_log = 0.05,
             bLog = FALSE) {

        GetVar <- function (spec, f, dfreq = NULL)
        {
            if (f[1] >= f[2])
                stop("f1 must be < f2")
            freqVector <- spec$freq
            if (f[1] < PaleoSpec::FirstElement(freqVector)) {
                warning("f[1] is smaller than the lowest frequency in the spectrum, set to the lowest frequency")
                f[1] <- PaleoSpec::FirstElement(freqVector)
            }
            if (f[2] > PaleoSpec::LastElement(freqVector)) {
                warning("f[2] is larger than the highest frequency in the spectrum, set to the highest frequency")
                f[2] <- PaleoSpec::LastElement(freqVector)
            }
            Intp <- function (freqRef, spec)
            {
                result <- list()
                result$freq <- freqRef
                result$spec <- approx(spec$freq, spec$spec, freqRef)$y
                class(result) <- "spec"
                return(result)
            }
            if (is.null(dfreq))
                dfreq <- min(diff(spec$freq)[1]/5, (f[2] - f[1])/100)
            newFreq <- seq(from = f[1], to = f[2], by = dfreq)
            vars <- mean(Intp(newFreq, spec)$spec) * (freq.end - freq.start) * 2
            return(vars)
        }

        xspec <- paleodata_spectrum(xin,
                                  interpolation,
                                  method,
                                  xout,
                                  transformation,
                                  transformation_type,
                                  detrend = detrend,
                                  df_log = df_log,
                                  bLog = bLog)
        # If xin is multivariate and therefore a list of spectra is compute use lapply, otherwise not
        if  (dim(as.matrix(xin))[2] > 1) {
            variance = lapply(xspec, function(xx) {GetVar(xx[[target]], f=c(freq.start, freq.end))})
        } else {
            variance = GetVar(xspec[[target]], f=c(freq.start, freq.end))
        }

        return(variance)
    }

#' @export
paleodata_varfromspec.Proxytibble <- function(xin,
                                          freq.start = NULL,
                                          freq.end = NULL,
                                          target = "raw",
                                          interpolation = TRUE,
                                          method = "lh14",
                                          xout = NULL,
                                          transformation = FALSE,
                                          transformation_type = "normalize",
                                          detrend = TRUE,
                                          df_log = 0.05,
                                          bLog = FALSE) {
    if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`paleodata_varfromspec` only implemented for `zoo_format == 'zoo'`")
    return(
        PTBoxProxydata::apply_proxy(
            xin,
            freq.start = freq.start,
            freq.end = freq.end,
            target = target,
            fun = paleodata_varfromspec.zoo,
            interpolation = interpolation,
            xout = xout,
            method = method,
            transformation = transformation,
            transformation_type = transformation_type,
            detrend = detrend,
            df_log = df_log,
            bLog = bLog
        )
    )
}


#' Function to compute scaling of the spectrum following the multi-taper methodology applied in Laepple and Huybers 2014, Rehfeld et al. 2018
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param freq.start Vector containing the start frequencies of the fitting interval(s)
#' @param freq.end Vector containing the start frequencies of the fitting interval(s)
#' @param target either "raw" or "logsmooth", default "raw"
#' @param interpolation Logical: should interpolation be applied?
#' @param method Parameters for interpolation (see paleodata_interpolation)
#' @param xout Parameters for interpolation (see paleodata_interpolation). If NULL the interpolation is set to the mean temporal resolution of each proxy time series.
#' @param transformation Logical: should transformation be applied?
#' @param transformation_type Parameters for transformation (see paleodata_transformation)
#' @param detrend Logical: should time series be linearly detrend before computing the spectrum
#' @param df_log Width of the log-smoother in log units
#' @param bLog TRUE: average in the log space of the power, FALSE: arithmetic average
#'
#' @return List with scaling coefficient, standard deviation and spectra (object description in SpecMTM and LogSmooth functions of PaleoSpec package)
#' @export
#'
#' @examples
#' # Load ice core example data
#' library(PTBoxProxydata)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # compute scaling
#' testout <- icecoredata %>% paleodata_scaling(., target="raw")
#' print(testout)
#'
#' @seealso
#' \link{paleodata_spectrum} (from `PTBoxProxytools`) for computation of spectra
#'
#' \link{SlopeFit} (from `PaleoSpec`) for computation of spectral slope
#'
paleodata_scaling <-
    function(xin,
             freq.start = NULL,
             freq.end = NULL,
             target = "raw",
             interpolation = TRUE,
             method = "lh14",
             xout = NULL,
             transformation = FALSE,
             transformation_type = "normalize",
             detrend = TRUE,
             df_log = 0.05,
             bLog = FALSE)
        UseMethod('paleodata_scaling')


#' @export
paleodata_scaling.zoo <-
    function(xin,
             freq.start = NULL,
             freq.end = NULL,
             target = "raw",
             interpolation = TRUE,
             method = "lh14",
             xout = NULL,
             transformation = FALSE,
             transformation_type = "normalize",
             detrend = TRUE,
             df_log = 0.05,
             bLog = FALSE) {

        xin <- paleodata_spectrum(xin,
                                 interpolation,
                                 method,
                                 xout,
                                 transformation,
                                 transformation_type,
                                 detrend = detrend,
                                 df_log = df_log,
                                 bLog = bLog)

        scaling = lapply(xin, function(xx) {PaleoSpec::SlopeFit(xx[[target]], freq.start, freq.end)})

        return(scaling)
}

#' @export
paleodata_scaling.Proxytibble <- function(xin,
                                          freq.start = NULL,
                                          freq.end = NULL,
                                          target = "raw",
                                          interpolation = TRUE,
                                          method = "lh14",
                                          xout = NULL,
                                          transformation = FALSE,
                                          transformation_type = "normalize",
                                          detrend = TRUE,
                                          df_log = 0.05,
                                          bLog = FALSE) {
    if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`paleodata_spectrum` only implemented for `zoo_format == 'zoo'`")
    return(
        PTBoxProxydata::apply_proxy(
            xin,
            freq.start = freq.start,
            freq.end = freq.end,
            target = target,
            fun = paleodata_scaling.zoo,
            interpolation = interpolation,
            xout = xout,
            method = method,
            transformation = transformation,
            transformation_type = transformation_type,
            detrend = detrend,
            df_log = df_log,
            bLog = bLog
        )
    )
}


#' Function to compute spectra (raw and smoothed using log-smooth) following the multi-taper methodology applied in Laepple and Huybers 2014, Rehfeld et al. 2018
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param interpolation Logical: should interpolation be applied?
#' @param method Parameters for interpolation (see paleodata_interpolation)
#' @param xout Parameters for interpolation (see paleodata_interpolation). If NULL the interpolation is set to the mean temporal resolution of each proxy time series.
#' @param transformation Logical: should transformation be applied?
#' @param transformation_type Parameters for transformation (see paleodata_transformation)
#' @param detrend Logical: should time series be linearly detrend before computing the spectrum
#' @param df_log Width of the log-smoother in log units
#' @param bLog TRUE: average in the log space of the power, FALSE: arithmetic average
#'
#' @return List with raw and log-smoothed spectra (object description in SpecMTM and LogSmooth functions of PaleoSpec package)
#' @export
#'
#' @examples
#' # Load ice core example data
#' library(PTBoxProxydata)
#' library(PaleoSpec)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # First compute spectra
#' testout <- icecoredata$proxy_data[[1]] %>% paleodata_spectrum(.)
#' # Second plot spectra (output is list with raw and log-smoothed spectra that can be plotted using PaleoSpec())
#' PaleoSpec::LPlot(testout$raw)
#' PaleoSpec::LLines(testout$logsmooth)
#'
#' @seealso
#' \link{paleodata_interpolation} (from `PTBoxProxytools`) for interpolation to regular time axis
#'
#' \link{paleodata_transformation} (from `PTBoxProxytools`) for time series transformation
#'
#' \link{SpecMTM} (from `PaleoSpec`) for computation of spectrum with multi-taper methodology
#'
#' \link{LogSmooth} (from `PaleoSpec`) for log smoothing of spectrum
#'
paleodata_spectrum <-
    function(xin,
             interpolation = TRUE,
             method = "lh14",
             xout = NULL,
             transformation = FALSE,
             transformation_type = "normalize",
             detrend = TRUE,
             df_log = 0.05,
             bLog = FALSE)
        UseMethod('paleodata_spectrum')


#' @export
paleodata_spectrum.zoo <-
    function(xin,
             interpolation = TRUE,
             interpolation_method = "lh14",
             interpolation_xout = NULL,
             transformation = FALSE,
             transformation_type = "normalize",
             detrend = TRUE,
             df_log = 0.05,
             bLog = FALSE) {
        if(is.null(interpolation_xout)){
            interpolation_xout = seq(from = zoo::index(xin)[1], zoo::index(xin)[length(zoo::index(xin))],
                                  by = mean(diff(zoo::index(xin))))
            }
        if (interpolation == TRUE) {
            xin <- paleodata_interpolation(xin, xout = interpolation_xout, method = interpolation_method)
        }
        if (transformation == TRUE) {
            xin <- paleodata_transformation(xin, transformation_type)
        }
        if (! ("matrix" %in% class(zoo::coredata(xin)))) {
            spectrum <- list()
            spectrum$raw <- PaleoSpec::SpecMTM(stats::as.ts(xin), detrend = detrend)
            spectrum$logsmooth <- PaleoSpec::LogSmooth(spectrum$raw, df_log = df_log, bLog = bLog)
            return(spectrum)
        } else {
            spectra <-
                apply(xin, 2, function(xx) {
                    spectrum <- list()
                    spectrum$raw <-
                        PaleoSpec::SpecMTM(stats::as.ts(zoo::zoo(xx, order.by = zoo::index(xin))), detrend = detrend)
                    spectrum$logsmooth <-
                        PaleoSpec::LogSmooth(spectrum$raw)
                    return(spectrum)
                })
            return(spectra)
        }
    }

#' @export
paleodata_spectrum.Proxytibble <- function(xin,
                                            interpolation = TRUE,
                                            interpolation_method = "lh14",
                                            interpolation_xout = NULL,
                                            transformation = FALSE,
                                            transformation_type = "normalize",
                                            detrend = TRUE,
                                            df_log = 0.05,
                                            bLog = FALSE) {
    if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`paleodata_spectrum` only implemented for `zoo_format == 'zoo'`")
    return(
        PTBoxProxydata::apply_proxy(
            xin,
            fun = paleodata_spectrum.zoo,
            interpolation = interpolation,
            interpolation_xout = interpolation_xout,
            interpolation_method = interpolation_method,
            transformation = transformation,
            transformation_type = transformation_type,
            detrend = detrend,
            df_log = df_log,
            bLog = bLog
        )
    )
}


#' Compute explained variance of an extracted signal
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param signal_type Method to extract signals;
#' Implemented methods:
#' 'pca' Principal component analysis
#' 'ca' Correspondance analysis
#' 'prc' Principal curves (princurve package)
#' 'custom_via_rda' Uses redundancy analysis to compute the explained variance of a specified signal (e.g. explained variance of AP signal for pollen assemblage record)
#' @param signal_components Components for which explained variance should be computed
#' @param reference_signal Specified signal, if signal_type='custom_via_rda' is selected
#'
#' @return Explained variance (given as a fraction) of the selected components (double between 0 and 1); if multiple components are specified,
#' the joint explained variance is returned; if input is Proxytibble, output is proxytibble where explained variance is given in the proxy_data column
#' @export
#'
#' @examples
#' # Load ice core example data
#' library(PTBoxProxydata)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # Compute PCA from multivariate zoo's
#' icecoredata_expl_var <- paleodata_explained_variance(icecoredata, 'pca')
#' print(paste0("Explained variance by PC1 for EDC time series: ",icecoredata_expl_var$proxy_data[[1]]))
#' print(paste0("Explained variance by PC1 for EPICA DML time series: ",icecoredata_expl_var$proxy_data[[2]]))
#' # We can do the same by first computing the PCA and then extract explained variance of the signal via redundancy analysis
#' icecoredata_pca <- paleodata_signal_extraction(icecoredata, 'pca')
#' icecoredata_expl_var_rda <- paleodata_explained_variance(icecoredata$proxy_data[[1]], 'custom_via_rda', reference_signal = icecoredata_pca$proxy_data[[1]])
#' print(icecoredata_expl_var_rda)
#'
#' @seealso
#' \link{prcomp} (from `stats`) for principal component analysis
#'
#' \link{cca} (from `vegan`) for correspondance analysis
#'
#' \link{principal_curve} (from `princurve`) for principal curves
#'
#' \link{rda} (from `vegan`) for redundancy analysis
#'
paleodata_explained_variance <- function(xin,signal_type="pca",signal_components=1,reference_signal=stats::prcomp(xin)$x[,1]) UseMethod('paleodata_explained_variance')

#' @export
paleodata_explained_variance.zoo <- function(xin,signal_type="pca",signal_components=1,reference_signal=stats::prcomp(xin)$x[,1]) {
    if (!signal_type %in% c("pca","ca","prc","custom_via_rda")) {
        stop("`signal_type` not supported")
    }
    if (signal_type == "pca") {
        signal <- stats::prcomp(xin)
        expl_var <- sum(signal$sdev[signal_components]^2)/sum(signal$sdev^2)
    }
    if (signal_type == "ca") {
        if (!requireNamespace("vegan", quietly = TRUE)) stop("suggested package `vegan` required for `signal_type = 'ca'`")
        signal <- vegan::cca(xin)
        expl_var <- sum(signal$CA$eig[signal_components])/sum(signal$CA$eig)
    }
    if (signal_type == "prc") {
        # Principal curve algorithm from princurve (original)
        if (!requireNamespace("princurve", quietly = TRUE)) stop("suggested package `princurve` required for `signal_type = 'prc'`")
        signal <- princurve::principal_curve(as.matrix(xin))
        expl_var <- 1-signal$dist / sum((as.matrix(xin)-array(rep(colMeans(xin),each=dim(xin)[1]),dim=dim(xin)))^2)
    }
    if (signal_type == "custom_via_rda") {
        # Compute the explained variance of a given signal based on redundancy analysis
        # (i.e. explained variance from optimized linear regression of data and reference signal)
        if (!requireNamespace("vegan", quietly = TRUE)) stop("suggested package `vegan` required for `signal_type = 'custom_via_rda'`")
        signal <- vegan::rda(X=xin,Y=reference_signal)
        expl_var <- sum(signal$CCA$eig)/sum(signal$CCA$eig,signal$CA$eig)
    }
    return(expl_var)
}

#' @export
paleodata_explained_variance.Proxytibble <- function(xin,signal_type="pca",signal_components=1,reference_signal=stats::prcomp(xin)$x[,1]) {
    if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`paleodata_explained_variance` only implemented for `zoo_format == 'zoo'`")
    return(
        PTBoxProxydata::apply_proxy(
            xin,
            fun = paleodata_explained_variance.zoo,
            signal_type=signal_type,
            signal_components=signal_components,
            reference_signal=reference_signal
        )
    )
}


#' Extract maximal window
#'
#' Extract maximal window fulfilling properties: minimal/maximal dates, minimal resolution (mean time step), maximal time step (to avoid large gaps), minimal length
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param t_min Lower time limit
#' @param t_max Upper time limit
#' @param min_res Minimal mean inter-sample range
#' @param max_step Maximal time step
#' @param min_length Minimal length of signal
#'
#' @return Proxytibble with proxy data in `zoo::zoo` format or irregular time series object (`zoo::zoo`), with data restricted to the found values, or NA if no window fulfills criteria
#' @export
#'
#' @examples
#' # Load Monticchio example data from PTBoxProxydata
#' library(PTBoxProxydata)
#' mng <- PTBoxProxydata::ProxyDataManager()
#' monticchiodata <- PTBoxProxydata::load_set(mng,'monticchio_testset',zoo_format = 'zoo', force_file=TRUE)
#' # Find the longest period in the Monticchio record in interval 30-60ka BP with mean temporal resolution below 250yrs, a maximal step length of 500yrs, and at least 20 samples
#' monticchiodata_highres <- find_max_window(monticchiodata,30000,60000,250,500,20)
#' # Plot restricted data
#' plot(monticchiodata_highres$proxy_data[[1]])
#' # Or directly for the AP zoo:
#' plot(find_max_window(monticchiodata$proxy_data[[1]],30000,60000,250,500,20))
#'
find_max_window <- function(xin,t_min=min(zoo::index(xin)),t_max=max(zoo::index(xin)),min_res=t_max-t_min,max_step=t_max-t_min,min_length=0) UseMethod('find_max_window')

#' @export
find_max_window.zoo <-
    function(xin,
             t_min = min(zoo::index(xin)),
             t_max = max(zoo::index(xin)),
             min_res = t_max - t_min,
             max_step = t_max - t_min,
             min_length = 0) {
        indices <- which(zoo::index(xin) >= t_min & zoo::index(xin) <= t_max)
        sample_ages <- zoo::index(xin)[indices]
        max_window <- c(1, 1)
        if (length(sample_ages) > 1) {
            for (i in 1:(length(sample_ages) - 1)) {
                for (j in (i + 1):length(sample_ages)) {
                    if ((sample_ages[j] - sample_ages[i]) >= min_length) {
                        if ((mean(diff(sample_ages[i:j])) <= min_res) &
                            (max(diff(sample_ages[i:j])) <= max_step) &
                            ((sample_ages[j] - sample_ages[i]) > sample_ages[max_window[2]] - sample_ages[max_window[1]])) {
                            max_window <- c(i, j)
                        }
                    }
                }
            }
        }
        if (diff(max_window) > 0) {
            if (! ("matrix" %in% class(zoo::coredata(xin)))) {
                return(xin[(min(indices) - 1) + max_window[1]:max_window[2]])
            } else {
                return(xin[(min(indices) - 1) + max_window[1]:max_window[2], ])
            }
        } else {
            return(NA)
        }
    }

#' @export
find_max_window.Proxytibble <- function(xin,t_min=min(zoo::index(xin)),t_max=max(zoo::index(xin)),min_res=t_max-t_min,max_step=t_max-t_min,min_length=0) {
    if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`find_max_window` only implemented for `zoo_format == 'zoo'`")
    return(
        PTBoxProxydata::apply_proxy(
            xin,
            fun = find_max_window.zoo,
            t_min=t_min,
            t_max=t_max,
            min_res=min_res,
            max_step=max_step,
            min_length=min_length
        )
    )
}



#' Constructs stacks of records from a proxytibble (list of proxy records)
#'
#' @param xin Input proxy data
#' @param nr_samples Number of samples
#' @param lon_min Lon min
#' @param lon_max Lon max
#' @param lat_min Lat min
#' @param lat_max Lat max
#' @param time Time
#' @param period_start Start of period
#' @param period_end End of period
#' @param sample_interval Sample interval
#' @param use_complete_records_only Should only records with values at all timesteps be used?
#' @param restrict_to_window Should maximal window fulfilling quality criteria be used?
#' @param window_tmin Start of period for window
#' @param window_tmax End of period for window
#' @param window_minres Min resolution
#' @param window_maxstep Maximal time step
#' @param window_minlength Minimum length of window
#' @param xin_interp Interpolated values
#' @param age_uncertainty Method for including age uncertainty
#' @param age_uncertainty_std Std. for age uncertainty
#' @param value_uncertainty Method for including proxy uncertainty
#' @param value_uncertainty_std Std. for proxy uncertainty
#' @param interpolation_method Interpolation width
#' @param bin_width Binning width
#' @param binning_function Binning function
#' @param site_normalization Should site values be normalized?
#' @param bootstrapping Use bootstrapping over sites?
#' @param stacking_method Stacking method
#' @param gridbox_size Gridbox size for grid-based interpolation methods
#' @param land_area_only Use land areas only for latitudinal weighting
#' @param dist_exp Exponent in computing weights with avgdist
#' @param within_group_method Averaging method within groups
#' @param group_weights Group weight
#' @param group_maps List of maps to compute group weights
#' @param anomalies Should anomalies from time mean be reported for each output sample
#' @param rescale Should each output sample be rescaled to unit variance
#'
#' @returns List with ensemble mean, ensemble quantiles, and ensemble members of stack
#' @export
site_mean <- function(xin,
                      nr_samples=1000,
                      lon_min=-180,
                      lon_max=180,
                      lat_min=-90,
                      lat_max=90,
                      time = NULL,
                      period_start = min(time),
                      period_end = max(time),
                      sample_interval = mean(diff(time)),
                      use_complete_records_only=FALSE,
                      restrict_to_window=FALSE,
                      window_tmin = period_start,
                      window_tmax = period_end,
                      window_minres = 2500,
                      window_maxstep = 5000,
                      window_minlength = (period_end - period_start)/2,
                      window_maxdistfromsample = 3000,
                      xin_interp = NULL,
                      age_uncertainty="white_noise",
                      age_uncertainty_std=1000,
                      value_uncertainty="white_noise",
                      value_uncertainty_std=5,
                      interpolation_method="binning",
                      bin_width=sample_interval,
                      binning_function=mean,
                      site_normalization=FALSE,
                      bootstrapping=TRUE,
                      stacking_method="site_mean",
                      gridbox_size=c(20,10),
                      land_area_only=TRUE,
                      dist_exp=1,
                      within_group_method="avgdist",
                      group_weights=NULL,
                      group_maps=NULL,
                      anomalies=TRUE,
                      rescale=FALSE) {

    # Prep: set time axis
    if (is.null(time)) {
        time <- seq(period_start,period_end,sample_interval)
    }
    # Proxytibble "xin" needs lon, lat, and proxydata (assumed to be either list of zoo's [add proxyzoo later for age ensembles] or array of interpolated values [so far only case of ensembles of interpolated values included, add single member case later])
    # Steps:
    # 0) Record selection based on regional limits and record quality
    tmp_ind <- which(xin$lon >= lon_min & xin$lon <= lon_max & xin$lat >= lat_min & xin$lat <= lat_max)
    xin <- xin[tmp_ind,]
    if (!is.null(xin_interp)) {
        xin_interp <- xin_interp[tmp_ind,,]
    }
    if (restrict_to_window == TRUE) {
        windowed_proxydata <- lapply(xin$proxy_data, find_max_window, t_min = window_tmin,t_max = window_tmax, min_res = window_minres, max_step = window_maxstep, min_length = window_minlength)
        xin <- xin[which(!is.na(windowed_proxydata)),]
        if (!is.null(xin_interp)) {
            xin_interp <- xin_interp[which(!is.na(windowed_proxydata)),,]
        }
        windowed_proxydata <- windowed_proxydata[which(!is.na(windowed_proxydata))]
    }
    cat(dim(xin)[1]," used for computing curve \n",sep="")

    # 1) Sample age and proxy uncertainties (none, whitenoise, ensemble)
    # TBD

    # 2) Interpolate (and smooth) to common axis, if no xin_interp is provided
    # CASE 1: Input data is proxytibble with (raw) proxy records (i.e. not interpolated to common time axis, no prior treatment like smoothing or interpolation is assume)
    if (is.null(xin_interp)) {
        # Interpolations methods: binning --> OK, lh14 --> looks good to me except for a few very low res records which might have to be removed, spline --> DON'T USE, PRODUCES TO MANY STRANGE VALUES, High-res linear + Gaussian kernel --> OK, loess --> OK, GAMLS --> !
        var_at_sites <- array(NA,dim=c(dim(xin)[1],length(time)))
        for (i in 1:dim(xin)[1]) {
            # ATTENTION: NEED TO UPDATE THE PALEODATA_INTERPOLATIO FUNCTION WITH extrapolate_dist (set values to NA which are further away from range of record than extrapolate_dist) AND max_dist_from_input_sample (set values to NA which are further away from nearest input time point than max_dist_from_input_sample) PARAMETERS
            # var_at_sites[i,] <- paleodata_interpolation(xin$proxy_data[[i]],method=interpolation_method,xout=time,extrapolate_dist=sample_interval/2,max_dist_from_input_sample=Inf,remove_na = FALSE, aggregation = FALSE, bin_width = bin_width, binning_function = binning_function)
            var_at_sites[i,] <- paleodata_interpolation(xin$proxy_data[[i]],
                                                        xout=time,
                                                        method=interpolation_method,
                                                        remove_na = FALSE,
                                                        aggregation = FALSE,
                                                        bin_width = bin_width,
                                                        binning_function = binning_function)
        }
    } else {
        # CASE 2: xin_interp (ensembles of interpolated values) is provided from pre-computed interpolation, e.g. Bayesian gamls models which take long to compute for large ensembles of records
        if (is.null(nr_samples)) {
            nr_samples <- dim(xin_interp)[3]
        }
        if (nr_samples == dim(xin_interp)[3]) {
            var_at_sites <- xin_interp
        } else {
            if (nr_samples < dim(xin_interp)[3]) {
                var_at_sites <- xin_interp[,,resample(1:dim(xin_interp)[3],size=nr_samples,replace=FALSE)]
            } else {
                var_at_sites <- xin_interp[,,resample(1:dim(xin_interp)[3],size=nr_samples,replace=TRUE)]
            }
        }
    }

    if (restrict_to_window == TRUE) {
        if (length(dim(var_at_sites)) == 2) {
            var_at_sites <- t(remove_extrapolated_samples(windowed_proxydata,zoo::zoo(t(var_at_sites),order.by=time),max_dist=window_maxdistfromsample))
        } else {
            var_at_sites <- remove_extrapolated_samples(windowed_proxydata,var_at_sites,time,max_dist=window_maxdistfromsample)
        }
    }
    if (use_complete_records_only == TRUE) {
        window_time <- time[which(time >= window_tmin & time <= window_tmax)]
        if (length(dim(var_at_sites)) == 2) {
            ind_complete_records <- which(apply(var_at_sites[,which(time >= window_tmin & time <= window_tmax)], 1, function(x) (length(which(!is.na(x))) == length(window_time))))
            var_at_sites <- var_at_sites[ind_complete_records,]
        } else {
            ind_complete_records <- which(apply(apply(var_at_sites[,which(time >= window_tmin & time <= window_tmax),],1:2,mean), 1, function(x) (length(which(!is.na(x))) == length(window_time))))
            var_at_sites <- var_at_sites[ind_complete_records,,]
        }
        xin <- xin[ind_complete_records,]
    }

    # 4) Normalize records --> unclear if that is useful, values of records are not directly comparable, but also normalization will lead to major changes even for records without any changes
    if (site_normalization == TRUE) {
        if (length(dim(var_at_sites)) == 2) {
            var_at_sites <- t(apply(var_at_sites,1,normalize))
        } else {
            var_at_sites <- aperm(apply(var_at_sites,c(1,3),normalize),c(2,3,1))
        }
    }

    # 5) bootstrap over records
    if (bootstrapping == TRUE) {
        bootstrap_samples <- sapply(1:nr_samples,function(i) resample(x=1:dim(var_at_sites)[1],replace=TRUE))
    } else {
        bootstrap_samples <- sapply(1:nr_samples,function(i) 1:dim(var_at_sites)[1])
    }

    # 6) Spatial mean (site_mean, avgdist_weighted_mean [based on average distance to other records], gridbox_weighted_mean [10x5°boxes + mean over all boxes with data], gridbox_and_lat_weighted_mean [10x5°boxes + zonal mean + weighted mean over all latitudinal bands])
    if (stacking_method == "groupweighted" & is.null(group_weights)) {
        group_weights <- compute_group_weights_from_maps(group_maps = group_maps)
        group_weights <- group_weights[sort(names(group_weights),index=TRUE)$ix]
    }
    if (length(dim(var_at_sites)) == 2) {
        stack_samples <- sapply(1:nr_samples, function(k) stack_records(site_data = list(var=var_at_sites[bootstrap_samples[,k],], lon=xin$lon[bootstrap_samples[,k]], lat=xin$lat[bootstrap_samples[,k]], group=xin$group[bootstrap_samples[,k]]),
                                                                        stacking_method = stacking_method,
                                                                        lon_min = lon_min,
                                                                        lon_max = lon_max,
                                                                        lat_min = lat_min,
                                                                        lat_max = lat_max,
                                                                        gridbox_size = gridbox_size,
                                                                        land_area_only=land_area_only,
                                                                        dist_exp = dist_exp,
                                                                        within_group_method=within_group_method,
                                                                        group_weights=group_weights,
                                                                        group_maps=group_maps))
    } else {
        stack_samples <- sapply(1:nr_samples, function(k) stack_records(site_data = list(var=var_at_sites[bootstrap_samples[,k],,k], lon=xin$lon[bootstrap_samples[,k]], lat=xin$lat[bootstrap_samples[,k]], group=xin$group[bootstrap_samples[,k]]),
                                                                        stacking_method = stacking_method,
                                                                        lon_min = lon_min,
                                                                        lon_max = lon_max,
                                                                        lat_min = lat_min,
                                                                        lat_max = lat_max,
                                                                        gridbox_size = gridbox_size,
                                                                        land_area_only=land_area_only,
                                                                        dist_exp = dist_exp,
                                                                        within_group_method=within_group_method,
                                                                        group_weights=group_weights,
                                                                        group_maps=group_maps))
    }

    # 7) Compute summary statistics
    if (anomalies == TRUE) {
        stack_samples <- apply(stack_samples,2,normalize,scale=FALSE)
    }
    if (rescale == TRUE) {
        stack_samples <- apply(stack_samples,2,normalize)
    }

    # 8) Compute summary statistics
    stack_quantiles <- apply(stack_samples,1,quantile,c(0.05,0.1,0.25,0.5,0.75,0.9,0.95),na.rm=TRUE)
    stack_means <- apply(stack_samples,1,mean,na.rm=TRUE)
    stack_included_records <- apply(var_at_sites,2,function(x) length(which(!is.na(x))))

    # Return time, means, quantiles, and samples as list
    return(list(time=time, means=stack_means, samples=stack_samples, quantiles = stack_quantiles, included_records = stack_included_records))
}

#' Wrapper function for processing of datasets using various processing methods
#'
#' Input data --> filtering --> interpolation --> time restriction (windowing) --> transformation --> signal extraction --> output data
#'
#' @param xin Proxytibble with proxy data in `zoo::zoo` format, or irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param filtering Logical: should filtering be applied
#' @param filter_type Parameters for filtering (see paleodata_filtering)
#' @param filter_scales Parameters for filtering (see paleodata_filtering)
#' @param detr_scale Parameters for filtering (see paleodata_filtering)
#' @param smooth_scale Parameters for filtering (see paleodata_filtering)
#' @param interpolation Logical: should interpolation be applied?
#' @param method Parameters for interpolation (see paleodata_interpolation)
#' @param xout Parameters for interpolation (see paleodata_interpolation)
#' @param windowing Logical: should restriction to time window be applied?
#' @param start_date Parameters for windowing (see paleodata_windowing)
#' @param end_date Parameters for windowing (see paleodata_windowing)
#' @param transformation Logical: should transformation be applied?
#' @param transformation_type Parameters for transformation (see paleodata_transformation)
#' @param signal_extraction Logical: should signal extraction be applied?
#' @param signal_type Parameters for signal extraction (see paleodata_signal_extraction)
#' @param signal_components Parameters for signal extraction (see paleodata_signal_extraction)
#'
#' @return Proxytibble with proxy data in `zoo::zoo` format or irregular time series object (`zoo::zoo`) containing the processed signals
#' @export
#'
#' @examples
#' #' # Load ice core example data
#' library(PTBoxProxydata)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # Step-by-step processing of ice core data: 1) Bandpass filter, 2) Interpolate to common time axis, 3) Restrict to 30-60ka BP, 4) Normalize, 5) Compute 1. PC
#' icecoredata_processed <- paleodata_processing(icecoredata, filtering = TRUE, filter_type = "bandpass", filter_scales = data.frame(lower=1000, upper=10000),
#'                                                            interpolation = TRUE, method = "lh14", xout = seq(20000,70000,by=100),
#'                                                            windowing = TRUE, start_date = 30000, end_date = 60000,
#'                                                            transformation = TRUE, transformation_type = "normalize",
#'                                                            signal_extraction = TRUE, signal_type = "pca", signal_components = 1)
#' # Plotting
#' plot(icecoredata_processed[[1]]$proxy_data[[1]])
#' plot(icecoredata_processed[[2]]$proxy_data[[1]])
#' plot(icecoredata_processed[[3]]$proxy_data[[1]])
#'
#' @seealso
#' \link{paleodata_filtering} (from `PTBoxProxytools`) for filtering
#'
#' \link{paleodata_interpolation} (from `PTBoxProxytools`) for interpolation
#'
#' \link{paleodata_windowing} (from `PTBoxProxytools`) for time period restriction
#'
#' \link{paleodata_transformation} (from `PTBoxProxytools`) for transformation
#'
#' \link{paleodata_signal_extraction} (from `PTBoxProxytools`) for signal extraction
#'
paleodata_processing <-
    function(xin,
             filtering = FALSE,
             filter_type = NULL,
             filter_scales = NULL,
             detr_scale = NULL,
             smooth_scale = NULL,
             interpolation = FALSE,
             interpolation_method = NULL,
             interpolation_xout = NULL,
             interpolation_gk_smooth_scale=NULL,
             windowing = FALSE,
             start_date = NULL,
             end_date = NULL,
             transformation = FALSE,
             transformation_type = NULL,
             signal_extraction = FALSE,
             signal_type = NULL,
             signal_components = NA) UseMethod('paleodata_processing')

#' @export
paleodata_processing.zoo <- function(xin,
                                 filtering=FALSE,filter_type=NULL,filter_scales=NULL,detr_scale=NULL,smooth_scale=NULL,interpolation=FALSE,
                                 interpolation_method=NULL,interpolation_xout=NULL,interpolation_gk_smooth_scale=NULL,windowing=FALSE,start_date=NULL,end_date=NULL,
                                 transformation=FALSE,transformation_type=NULL,signal_extraction=FALSE,signal_type=NULL,signal_components=NA) {
    # Filtering
    if (filtering == TRUE) {
        xin <- paleodata_filtering(xin,filter_type,filter_scales,detr_scale,smooth_scale)
    }
    # Interpolation
    if (interpolation == TRUE) {
        xin <- paleodata_interpolation(xin,xout=interpolation_xout,method=interpolation_method,gk_smooth_scale=interpolation_gk_smooth_scale)
    }
    # Time restriction (windowing)
    if (windowing == TRUE) {
        xin <- paleodata_windowing(xin,start_date,end_date)
    }
    # Transformation
    if (transformation == TRUE) {
        xin <- paleodata_transformation(xin,transformation_type)
    }
    # Signal Extraction
    if (signal_extraction == TRUE) {
        xin <- paleodata_signal_extraction(xin,signal_type,signal_components)
    }
    return(xin)
}


#' @export
paleodata_processing.Proxytibble <- function(xin,
                                             filtering = FALSE,
                                             filter_type = NULL,
                                             filter_scales = NULL,
                                             detr_scale = NULL,
                                             smooth_scale = NULL,
                                             interpolation = FALSE,
                                             interpolation_method = NULL,
                                             interpolation_xout = NULL,
                                             interpolation_gk_smooth_scale=NULL,
                                             windowing = FALSE,
                                             start_date = NULL,
                                             end_date = NULL,
                                             transformation = FALSE,
                                             transformation_type = NULL,
                                             signal_extraction = FALSE,
                                             signal_type = NULL,
                                             signal_components = NA) {
    if (all(class(xin[[PTBoxProxydata::Proxytibble_colnames_proxy_data()]][[1]]) != 'zoo'))
        stop("`paleodata_processing` only implemented for `zoo_format == 'zoo'`")
    return(
        PTBoxProxydata::apply_proxy(
            xin,
            fun = paleodata_processing.zoo,
            filtering = filtering,
            filter_type = filter_type,
            filter_scales = filter_scales,
            detr_scale = detr_scale,
            smooth_scale = smooth_scale,
            interpolation = interpolation,
            interpolation_method = interpolation_method,
            interpolation_xout = interpolation_xout,
            interpolation_gk_smooth_scale=interpolation_gk_smooth_scale,
            windowing = windowing,
            start_date = start_date,
            end_date = end_date,
            transformation = transformation,
            transformation_type = transformation_type,
            signal_extraction = signal_extraction,
            signal_type = signal_type,
            signal_components = signal_components
        )
    )
}

#' Wrapper function to apply different processing chains to an irregular time series object
#'
#' Multiple versions of the form: input data --> filtering --> interpolation --> time restriction (windowing) --> transformation --> signal extraction --> output data
#'
#' @param xin Irregular time series object (`zoo::zoo`), xin can be multivariate
#' @param processing_name Names of the different processing applications
#' @param filtering Vector of logicals: should filtering be applied
#' @param filter_type Vector of parameters for filtering (see paleodata_filtering)
#' @param filter_scales Vector of parameters for filtering (see paleodata_filtering)
#' @param detr_scale Vector of parameters for filtering (see paleodata_filtering)
#' @param smooth_scale Vector of parameters for filtering (see paleodata_filtering)
#' @param interpolation Vector of logicals: should interpolation be applied?
#' @param method Vector of parameters for interpolation (see paleodata_interpolation)
#' @param xout Vector of parameters for interpolation (see paleodata_interpolation)
#' @param windowing Vector of logicals: should restriction to time window be applied?
#' @param start_date Vector of parameters for windowing (see paleodata_windowing)
#' @param end_date Vector of parameters for windowing (see paleodata_windowing)
#' @param transformation Vector of logicals: should transformation be applied?
#' @param transformation_type Vector of parameters for transformation (see paleodata_transformation)
#' @param signal_extraction Vector of logicals: should signal extraction be applied?
#' @param signal_type Vector of parameters for signal extraction (see paleodata_signal_extraction)
#' @param signal_components Vector of parameters for signal extraction (see paleodata_signal_extraction)
#'
#' @return List of processed signals, provided as list of irregular time series objects (`zoo::zoo`). Fields `$data` contain the processed signals.
#' @export
#'
#' @examples
#' #' # Load ice core example data
#' library(PTBoxProxydata)
#' mng <- ProxyDataManager()
#' icecoredata <- load_set(mng,'icecore_testset',zoo_format = 'zoo')
#' # Step-by-step processing of ice core data in multiple different forms (here, smooting, detrending, and bandpass filtering is applied at the same time): 1) Smooth / Detrend / Bandpass filter, 2) Interpolate to common time axis, 3) Restrict to 30-60ka BP, 4) Normalize, 5) Compute 1. PC
#' icecoredata_processed <- paleodata_multiprocessing(xin = icecoredata$proxy_data[[1]], processing_name = c("Smoothing","Detrending","Bandpass filtering"),
#'                                                            filtering = rep(TRUE,times=3), filter_type = c("smooth","detrend","bandpass"), detr_scale = rep(1000,times=3), smooth_scale = rep(10000,times=3), filter_scales = data.frame(lower=rep(1000,times=3), upper=rep(10000,times=3)),
#'                                                            interpolation = rep(TRUE,times=3), method = rep("lh14",times=3), xout = list(seq(20000,70000,by=100),seq(20000,70000,by=100),seq(20000,70000,by=100)),
#'                                                            windowing = rep(TRUE,times=3), start_date = rep(30000,times=3), end_date = rep(60000,times=3),
#'                                                            transformation = rep(TRUE,times=3), transformation_type = rep("normalize",times=3),
#'                                                            signal_extraction = rep(TRUE,times=3), signal_type = rep("pca",times=3), signal_components = rep(1,times=3))
#' # Plotting
#' plot(icecoredata_processed[[1]]$data)
#' plot(icecoredata_processed[[2]]$data)
#' plot(icecoredata_processed[[3]]$data)
#'
#' @seealso
#' \link{paleodata_processing} (from `PTBoxProxytools`) for wrapper of individual processing chains
#'
#' \link{paleodata_filtering} (from `PTBoxProxytools`) for filtering
#'
#' \link{paleodata_interpolation} (from `PTBoxProxytools`) for interpolation
#'
#' \link{paleodata_windowing} (from `PTBoxProxytools`) for time period restriction
#'
#' \link{paleodata_transformation} (from `PTBoxProxytools`) for transformation
#'
#' \link{paleodata_signal_extraction} (from `PTBoxProxytools`) for signal extraction
#'
paleodata_multiprocessing <-
    function(xin,
             processing_name,
             filtering = rep(FALSE, times = length(processing_name)),
             filter_type = NULL,
             filter_scales = NULL,
             detr_scale = NULL,
             smooth_scale = NULL,
             interpolation = rep(FALSE, times = length(processing_name)),
             interpolation_method = NULL,
             interpolation_xout = NULL,
             interpolation_gk_smooth_scale = NULL,
             windowing = rep(FALSE, times = length(processing_name)),
             start_date = NULL,
             end_date = NULL,
             transformation = rep(FALSE, times = length(processing_name)),
             transformation_type = NULL,
             signal_extraction = rep(FALSE, times = length(processing_name)),
             signal_type = NULL,
             signal_components = NA) UseMethod('paleodata_multiprocessing')

#' @export
paleodata_multiprocessing.zoo <- function(xin,processing_name,filtering=rep(FALSE,times=length(processing_name)),
                                      filter_type=NULL,filter_scales=NULL,detr_scale=NULL,smooth_scale=NULL,
                                      interpolation=rep(FALSE,times=length(processing_name)),interpolation_method=NULL,
                                      interpolation_xout=NULL,interpolation_gk_smooth_scale=NULL,windowing=rep(FALSE,times=length(processing_name)),start_date=NULL,
                                      end_date=NULL,transformation=rep(FALSE,times=length(processing_name)),transformation_type=NULL,
                                      signal_extraction=rep(FALSE,times=length(processing_name)),signal_type=NULL,signal_components=NA) {
    data_processings <- list()
    for (i in 1:length(processing_name)) {
        data_processings[[i]] <- list(name=processing_name[i])
        data_processings[[i]]$data <- paleodata_processing.zoo(xin,filtering=filtering[i],filter_type=filter_type[i],
                                                           filter_scales=filter_scales[i,],detr_scale=detr_scale[i],
                                                           smooth_scale=smooth_scale[i],interpolation=interpolation[i],
                                                           interpolation_method=interpolation_method[i],interpolation_xout=interpolation_xout[[i]],
                                                           interpolation_gk_smooth_scale=interpolation_gk_smooth_scale[i],
                                                           windowing=windowing[i],start_date=start_date[i],end_date=end_date[i],
                                                           transformation=transformation[i],transformation_type=transformation_type[i],
                                                           signal_extraction=signal_extraction[i],signal_type=signal_type[i],
                                                           signal_components=signal_components[[i]])
    }
    return(data_processings)
}
