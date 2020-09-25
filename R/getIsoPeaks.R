#' Isotopic peaks fitting for specified m/z
#'
#' Petforms peak fitting at provided 'm/z' value, and also at m/z values
#' corresponding to the 13C, 15N, 2H and 18O isotopes. Returns data_frame
#' object containing information about fitting results.
#'
#'
#' @param ss Spectrum formatted as numeric matrix with 2 columns (m/z and i),
#' sorted by mz
#' @param mz m/z value corresponding to the monoisotopic mass of singly charged
#' ion
#' @param width m/z interval for searching the isotopic peaks. Default value is
#' 0.002 which is good for resolution 50-60K (at m/z=200).
#' @param npoint Minimum number of data points in the peak. If peak have less
#' points it will be rejected.
#' @param tol m/z tolerance for searching monoisotopic peak.
#' @param fixSigma If True, then peak width is only fitted for monoisotopic
#' peak, and then it is fixed for the rest isotopes.
#' @param plot If True, fitted peaks are plotted using ggplot.
#' @param peaks Character vector defining which peaks to fit (13C, 15N, 2H or 18O).
#'
#' @return Dataframe with peak parameters
#' @export
#'
#' @examples
get_isopeaks <- function(ss, mz, width = 0.002, npoint = 10, tol = 0.01, fixSigma = T, plot = F, peaks=c("13C", "15N", "2H", "18O")) {
    # Checking input data
    stopifnot(length(dim(ss)) == 2)
    # here we'll put our result
    res <- data.frame()

    # Fitting monoisotopic peak
    p_int <- findInterval(c(mz - tol, mz + tol), ss[, 1])
    if (diff(p_int) < npoint) {
        warning(sprintf("No peak found at mz=%.4f with tol=%.4f", mz, tol))
        return(data.frame())
    }
    sub_ss <- ss[p_int[1]:p_int[2], ]

    istart = iend = which.max(sub_ss[, 2])
    while (sub_ss[istart, 2] > 0 & istart > 1) istart <- istart - 1
    while (sub_ss[iend, 2] > 0 & iend < nrow(sub_ss)) iend <- iend + 1

    if ((iend - istart) < npoint) {
        warning(sprintf("Less then %d datapoints in peak at mz=%.4f", npoint, mz))
        return(data.frame())
    }

    sub_ss <- sub_ss[istart:iend, ]
    # adjusting initial parameters for gaussian fit
    startPars <- c(mu = sub_ss[which.max(sub_ss[, 2]), 1], sigma = sd(sub_ss[, 1]), I = sub_ss[which.max(sub_ss[, 2]), 2] * sd(sub_ss[, 1])/2)
    model = fitPeak(sub_ss[, 1], sub_ss[, 2], startPars = startPars, maxAttempts = 10)
    r2 <- rSquared(model, sub_ss[, 2])
    if (plot) {
        plots <- list()
        plots[["mono"]] <- plot_peakfit(sub_ss, model) + ggtitle("mono")
    }
    if (class(model) != "try-error" & r2 > 0.9) {
        v <- coef(model)
        res <- tidy(model) %>% mutate(peak = "0") %>% select(peak, term, estimate) %>% spread(term, estimate) %>% mutate(R2 = r2, isoshift = 0, theorshift = 0,
            masserror = 0, isoratio = 1)
        for (isotope_ in c("13C", "15N", "2H", "18O")) {
            element_ = sub("\\d+", "", isotope_)
            isoshift_ = diff((isotopes %>% filter(element == element_) %>% filter(abundance > 0.9 | isotope == isotope_))$mass)
            isomz <- v["mu"] + isoshift_
            p_int <- findInterval(c(isomz - width, isomz + width), ss[, 1])
            if (diff(p_int) > npoint) {
                sub_ss <- ss[p_int[1]:p_int[2], ]
                istart = iend = which.max(sub_ss[, 2])
                while (sub_ss[istart, 2] > 0 & istart > 1) istart <- istart - 1
                while (sub_ss[iend, 2] > 0 & iend < nrow(sub_ss)) iend <- iend + 1
                if ((iend - istart + 1) > npoint) {
                  sub_ss <- sub_ss[istart:iend, ]
                  startPars <- c(isomz, v["sigma"], v["I"] * (isotopes %>% filter(isotope == isotope_ & element == element_))$abundance * 4)
                  model_ = fitPeak(sub_ss[, 1], sub_ss[, 2], startPars = startPars, maxAttempts = 10, sigma = ifelse(fixSigma, v["sigma"], 0))
                  r2_ = rSquared(model_, sub_ss[, 2])

                  if (class(model_) != "try-error" & r2_ > 0.5) {
                    if (plot) {
                      plots[[isotope_]] <- plot_peakfit(sub_ss, model_) + ggtitle(isotope_)
                    }
                    v_ = coef(model_)
                    res <- res %>% bind_rows(tidy(model_) %>% mutate(peak = isotope_) %>% select(peak, term, estimate) %>% spread(term, estimate) %>% mutate(R2 = r2_,
                      isoshift = v_["mu"] - v["mu"], theorshift = isoshift_, masserror = v_["mu"] - v["mu"] - isoshift_, isoratio = v_["I"]/v["I"]))
                  }
                }
            }
        }
    }
    if (plot && nrow(res) > 0) {
        do.call(grid.arrange, plots)
    }
    return(res)
}


#' Isotopic peaks fitting for specified m/z
#'
#' Petforms peak fitting at provided 'm/z' value, and also fits fine isotopic
#' structure of +1 isotopic peak corresponding to the 13C, 15N and 2H isotopes.
#' Unlike \code{get_isopeaks}, isotopic peaks are not fitted separately, but
#' rather as a single envelope with fixed m/z values for peaks.
#' Returns data_frame object containing information about fitting results.
#'
#'
#' @param ss Spectrum formatted as numeric matrix with 2 columns (m/z and i),
#' sorted by mz
#' @param mz m/z value corresponding to the monoisotopic mass of singly charged
#' ion
#' @param width m/z interval for searching the isotopic peaks. Default value is
#' 0.002 which is good for resolution 50-60K (at m/z=200).
#' @param npoint Minimum number of data points in the peak. If peak have less
#' points it will be rejected.
#' @param tol m/z tolerance for searching monoisotopic peak.
#' @param fixSigma If True, then peak width is only fitted for monoisotopic
#' peak, and then it is fixed for the rest isotopes.
#' @param plot If True, fitted peaks are plotted using ggplot.
#'
#' @return Dataframe with peak parameters
#' @export
#'
#' @examples
get_isopeaks2 <- function(ss, mz, width = 0.002, npoint = 10, tol = 0.01, fixSigma = T, plot = F) {
    # Checking input data
    stopifnot(length(dim(ss)) == 2)
    # here we'll put our result
    res <- data.frame()

    # Fitting monoisotopic peak
    p_int <- findInterval(c(mz - tol, mz + tol), ss[, 1])
    if (diff(p_int) < npoint) {
        warning(sprintf("No peak found at mz=%.4f with tol=%.4f", mz, tol))
        return(data.frame())
    }
    sub_ss <- ss[p_int[1]:p_int[2], ]

    istart = iend = which.max(sub_ss[, 2])
    while (sub_ss[istart, 2] > 0 & istart > 1) istart <- istart - 1
    while (sub_ss[iend, 2] > 0 & iend < nrow(sub_ss)) iend <- iend + 1

    if ((iend - istart) < npoint) {
        # warning(sprintf('Less then %d datapoints in peak at mz=%.4f', npoint, mz))
        return(data.frame())
    }

    sub_ss <- sub_ss[istart:iend, ]
    # adjusting initial parameters for gaussian fit
    startPars <- c(mu = sub_ss[which.max(sub_ss[, 2]), 1], sigma = sd(sub_ss[, 1]), I = sub_ss[which.max(sub_ss[, 2]), 2] * sd(sub_ss[, 1])/2)
    model = fitPeak(sub_ss[, 1], sub_ss[, 2], startPars = startPars, maxAttempts = 10)
    r2 <- rSquared(model, sub_ss[, 2])
    spec1 = matrix(ncol = 2, nrow = 0)
    if (plot & class(model) != "try-error") {
        plots <- list()
        plots[["mono"]] <- plot_peakfit(sub_ss, model) + ggtitle("mono")
    }

    if (class(model) != "try-error" & r2 > 0.9) {
        v <- coef(model)
        for (isotope_ in c("13C", "15N", "2H")) {
            element_ = sub("\\d+", "", isotope_)
            isoshift_ = diff((isotopes %>% filter(element == element_) %>% filter(abundance > 0.9 | isotope == isotope_))$mass)
            isomz <- v["mu"] + isoshift_
            sub_int <- findInterval(c(isomz - 4 * v["sigma"], isomz + 4 * v["sigma"]), ss[, 1])
            sub_ss <- ss[sub_int[1]:sub_int[2], ]
            spec1 <- rbind(spec1, sub_ss)
        }

    }
    startPars <- c(rc = 0.05, rn = 0.004, rh = 0.001, dm = 0)
    model_ <- fit_isopeaks(x = spec1[, 1], y = spec1[, 2], startPars = startPars, I = v["I"], mu = v["mu"], sigma = v["sigma"], maxAttempts = 25)
    r2_ = rSquared(model_, spec1[, 2])
    if (class(model_) != "try-error" && r2_ > 0.7) {
        if (plot)
            plots[["1"]] <- plot_peakfit(spec1, model_) + ggtitle("+1")
        tm <- try(tidy(model_), silent = T)
        if (class(tm) != "try-error")
            res <- res %>% bind_rows(tm %>% mutate(R2 = r2_, I = v["I"], sigma = v["sigma"], mu = v["mu"]))
    }
    if (plot && nrow(res) > 0) {
        do.call(grid.arrange, plots)
    }
    return(res)
}











get_isopeaks_ <- function(ss, mz, width = 0.002, npoint = 10, tol = 0.01, fixSigma = T, plot = F, peaks=c("13C", "15N", "2H", "18O")) {
    # Checking input data
    stopifnot(length(dim(ss)) == 2)
    # here we'll put our result
    res <- data.frame()

    # Fitting monoisotopic peak
    p_int <- findInterval(c(mz - tol, mz + tol), ss[, 1])
    if (diff(p_int) < npoint) {
        warning(sprintf("No peak found at mz=%.4f with tol=%.4f", mz, tol))
        return(data.frame())
    }
    sub_ss <- ss[p_int[1]:p_int[2], ]

    istart = iend = which.max(sub_ss[, 2])
    istart = istart-1
    dv_sub <- diff(sub_ss)
    while (dv_sub[istart, 2] > 0 & istart > 1) istart <- istart - 
        1
    while (dv_sub[iend, 2] < 0 & iend < nrow(sub_ss)) iend <- iend + 
        1
    if ((iend - istart) < npoint) {
        warning(sprintf("Less then %d datapoints in peak at mz=%.4f", npoint, mz))
        return(data.frame())
    }
    istart <-  istart + 1
    sub_ss <- sub_ss[istart:iend, ]
    # adjusting initial parameters for gaussian fit
    startPars <- c(mu = sub_ss[which.max(sub_ss[, 2]), 1], sigma = sd(sub_ss[, 1]), I = sub_ss[which.max(sub_ss[, 2]), 2] * diff(range(sub_ss[, 1]))/2)
    model = fitPeak(sub_ss[, 1], sub_ss[, 2], startPars = startPars, maxAttempts = 10)
    r2 <- rSquared(model, sub_ss[, 2])
    if (plot) {
        plots <- list()
        plots[["mono"]] <- plot_peakfit(sub_ss, model) + ggtitle("mono")
    }
    if (class(model) != "try-error" & r2 > 0.9) {
        v <- coef(model)
        res <- tidy(model) %>% mutate(peak = "0") %>% select(peak, term, estimate) %>% spread(term, estimate) %>% mutate(R2 = r2, isoshift = 0, theorshift = 0,
            masserror = 0, isoratio = 1)
        for (isotope_ in c("13C", "15N", "2H", "18O")) {
            element_ = sub("\\d+", "", isotope_)
            isoshift_ = diff((isotopes %>% filter(element == element_) %>% filter(abundance > 0.9 | isotope == isotope_))$mass)
            isomz <- v["mu"] + isoshift_
            p_int <- findInterval(c(isomz - width, isomz + width), ss[, 1])
            if (diff(p_int) > npoint) {
                sub_ss <- ss[p_int[1]:p_int[2], ]
                istart = iend = which.max(sub_ss[, 2])
                istart = istart-1
                dv_sub <- diff(sub_ss)
                while (dv_sub[istart, 2] > 0 & istart > 1) istart <- istart - 1
                while (dv_sub[iend, 2] < 0 & iend < nrow(sub_ss)) iend <- iend + 1
                if ((iend - istart + 1) > npoint) {
                  istart = istart + 1
                  sub_ss <- sub_ss[istart:iend, ]
                  startPars <- c(mu = sub_ss[which.max(sub_ss[, 2]), 1], sigma = sd(sub_ss[, 1]), I = sub_ss[which.max(sub_ss[, 2]), 2] * diff(range(sub_ss[, 1]))/2)
                  model_ = fitPeak(sub_ss[, 1], sub_ss[, 2], startPars = startPars, maxAttempts = 10)
                  r2_ = rSquared(model_, sub_ss[, 2])

                  if (class(model_) != "try-error" & r2_ > 0.5) {
                    if (plot) {
                      plots[[isotope_]] <- plot_peakfit(sub_ss, model_) + ggtitle(isotope_)
                    }
                    v_ = coef(model_)
                    res <- res %>% bind_rows(tidy(model_) %>% mutate(peak = isotope_) %>% select(peak, term, estimate) %>% spread(term, estimate) %>% mutate(R2 = r2_,
                      isoshift = v_["mu"] - v["mu"], theorshift = isoshift_, masserror = v_["mu"] - v["mu"] - isoshift_, isoratio = v_["I"]/v["I"]))
                  }
                }
            }
        }
    }
    if (plot && nrow(res) > 0) {
        do.call(grid.arrange, plots)
    }
    return(res)
}
