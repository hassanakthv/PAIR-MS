#' fit_isopeaks
#'
#' Fits gaussian peak for the data
#' \code{ y = I/sigma * exp(-1/2*(x-mu)^2 / sigma^2) }
#' using nls function with initial parameters specified by named vector
#' \code{startPars}. If \code{sigma > 0} orivuded, then only two parameters are
#' fitted using nls - intensity \code{I} and position \code{mu}.
#'
#' @param x vector of x-values;
#' @param y vector of y-values;
#' @param startPars names vector of initial parameters
#' \code{mu}, \code{I} and \code{sigma};
#' @param maxAttempts maximum number of attempts, see below;
#' @param sigma if \code{sigma>0} provided, two-parameter fitting is performed,
#' fixing \code{sigma} at provided value;
#'
#' @return Returns an object of class \code{nls} if the fitting was successful.
#'
#' @export
#'
#' @examples
#' fitPeak(x=x, y=y, startPars=c('mu' = 70.065, 'I' = 5e3, 'sigma' = 3e-4))
#' fitPeak(x=x, y=y, startPars=c('mu' = 70.065, 'I' = 5e3, 'sigma' = 3e-4), sigma=3e-4)
fit_isopeaks <- function(x, y, startPars = c(rc = 0.02, rn = 0.003, rh = 4e-04, dm = 0), I = 10000, sigma = 0.001, mu = 70, maxAttempts = 1000) {
    dc = diff((isotopes %>% filter(element == "C") %>% filter(abundance > 0.9 | isotope == "13C"))$mass)
    dn = diff((isotopes %>% filter(element == "N") %>% filter(abundance > 0.9 | isotope == "15N"))$mass)
    dh = diff((isotopes %>% filter(element == "H") %>% filter(abundance > 0.9 | isotope == "2H"))$mass)
    fctStr <- function(x, rc, rn, rh, dm) {
        r <- abs(rc) * I/sigma * exp(-1/2 * (x - mu - dc - dm)^2/sigma^2)
        r <- r + abs(rn) * I/sigma * exp(-1/2 * (x - mu - dn - dm)^2/sigma^2)
        r <- r + abs(rh) * I/sigma * exp(-1/2 * (x - mu - dh - dm)^2/sigma^2)
        return(r)
    }
    fitFct <- as.formula("y ~ fctStr(x, rc, rn, rh, dm)")
    
    xVec <- x
    yVec <- y
    
    varyPars <- 0
    attempts <- 0
    repeatLoop <- TRUE
    ## Check if number of non-missing values is sufficient (NLS can only handle data with at least three non-missing values)
    validValues <- !is.na(yVec)
    if (sum(validValues) <= 2) {
        m <- NA
        class(m) <- "try-error"
    } else {
        ## Perform fit
        while (repeatLoop & attempts < maxAttempts) {
            parTmp <- startPars * (1 + varyPars * (runif(1, -0.1, 0.1)))
            m <- try(nls(formula = fitFct, start = parTmp, data = list(x = xVec, y = yVec), na.action = na.exclude, algorithm = "port"), silent = TRUE)
            attempts <- attempts + 1
            varyPars <- 1
            if (class(m) != "try-error") {
                mcoef <- try(tidy(m), silent = TRUE)
                if (class(mcoef) != "try-error") 
                  repeatLoop <- FALSE
            }
        }
    }
    return(m)
}
