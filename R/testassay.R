#' Hypothesis testing procedure for assay validation
#'
#' Specify the assay results, procedure parameters, and the model assumptions,
#' and this function calculates the relevant statistics
#'
#' @param x The vector of assay values
#' @param m The vector of values indicating sample membership
#' @param n The vector of values indicating replicate membership
#' @param q The confidence level, typically 0.8 or 0.9
#' @param model String specifying the distribution for the assay values. Valid
#'   values are "normal" or "lognormal"
#' @param constant String specifying whether the variance is assumed to be
#'   constant over the levels ("variance") or the coefficient of variation
#'   ("CV").
#' @param data Data frame or environment in which to look for x
#'
#' @export
#' @aliases testassay
#' @return An object of class "assaytest", which is a list of components including a data frame of the relevant statistics calculated on x. Print, summary, predict, and plot methods are available.
#' @details The list has the following components
##' \itemize{
##'  \item{sumtab}{Table summarizing the experiment, includes mean values, SD or CV estimates, and upper confidence limits on those. }
##'  \item{Umax}{The maximum of the upper limits on the SD or CV, used in the effective SD interval calculation}
##'  \item{n}{The number of samples per level}
##'  \item{m}{The number of levels}
##'  \item{q}{The confidence level}
##'  \item{model}{The assumed model}
##'  \item{constant}{The parameter assumed to be constant (either SD or CV)}
##'  \item{alpha}{The alpha level, calculated as (1 - q)^(1 / m)}
##'  \item{x}{The data vector supplied by the user}
##' }
#'

testassay <- function(x, m, n, q = .9, model = "normal", constant = "variance", data = NULL) {

  stopifnot(q > 0 & q < 1)
  stopifnot(model %in% c("normal", "lognormal"))
  stopifnot(constant %in% c("var", "variance", "sd", "SD", "CV", "cv"))

  if(!is.null(data)){

    x <- eval(substitute(x), envir = data)
    m <- eval(substitute(m), envir = data)
    n <- eval(substitute(n), envir = data)

  }

  ## need same number of levels per sample

  tabm <- table(m)
  if(!all(tabm == max(tabm))) {

    stop("Needs the same number of levels per sample. ")
  }

  inn <- length(unique(n))
  inm <- length(tabm)

  alpha <- (1 - q)^(1 / inm)

  meanAA <- tapply(x, m, mean)
  varAA <- tapply(x, m, var)

  if(constant %in% c("var", "variance", "sd", "SD")) {
    ## mean and variance by m
    constant <- "SD"

    Uj <- sqrt(
      (inn - 1) * varAA / qchisq(alpha, inn - 1)
    )

    Umax <- max(Uj)
    sumtab <- data.frame(sample = names(meanAA), mean = meanAA,
                         sd = sqrt(varAA), U.SD = Uj)
    sumtab <- sumtab[order(sumtab$mean), ]
    rownames(sumtab) <- NULL

  }

  if(constant %in% c("cv", "CV")){
    constant <- "CV"

    cvAA <- sqrt(varAA) / meanAA
    ## construct upper confidence limit for cv

    rootfunc<-function(cvin, T){
      cvin - sqrt(inn)/qt(alpha, inn - 1, sqrt(inn) / T,lower.tail=TRUE)
    }

    Uj <- sapply(cvAA, function(cvout) uniroot(function(x) rootfunc(cvout, x), c(1e-4, 1e4))$root)
    Umax <- max(Uj)
    sumtab <- data.frame(sample = names(meanAA), mean = meanAA,
                         cv = cvAA, U.CV = Uj)
    sumtab <- sumtab[order(sumtab$mean), ]
    rownames(sumtab) <- NULL

  }

  ## assemble the object and give it class "assaytest"

  structure(list(sumtab = sumtab, Umax = Umax, n = inn, m = inm, q = q,
                 model = model, constant = constant, alpha = alpha, x = x), class = "assaytest")



}


#' Print the results of an assaytest
#'
#' @export
#'
#' @param x An object of class "assaytest", as created by the function \link{testassay}
#' @param ... Additional arguments, currently unused
#'
print.assaytest <- function(x, ...){

  cat("Results of a", paste(x$m, x$n, round(x$q * 100), sep = ":"), "% procedure. \n")
  cat("Assuming a", x$model, "model with constant", x$constant, "\n \n")

  print(x$sumtab, ...)


  cat("\n", paste0("U.", x$constant),  "is the one-sided upper confidence limit at alpha =", 1 - x$alpha)
  cat("\nBetween assay values of", min(x$sumtab$mean), "to", max(x$sumtab$mean),
      "the assay passes the procedure with bound", x$Umax, "\n")
  cat("Use predict to obtain confidence intervals for future observations")

}



#' Construct effective standard deviation intervals for observed assay values
#'
#' Computes effective standard deviation intervals for observed assay results.
#' These intervals have greater than 68.27% coverage.
#'
#' @export
#'
#' @param object An object of class "assaytest"
#' @param newdata A vector of observed values
#' @param ... additional arguments
#'
#' @return A data frame with the observed values, lower, and upper confidence
#'   limits
#'

predict.assaytest <- function(object, newdata, ...){

  if(missing(newdata)){

    newdata <- object$x

  }

  if(any(min(object$sumtab$mean) > newdata) | any(max(object$sumtab$mean) < newdata)) {

    warning("New observations outside the range of the assay validation procedure!")

  }

  if(object$model == "normal" & object$constant == "SD"){

    data.frame(obs = newdata, lower = newdata - object$Umax, upper = newdata + object$Umax)

  } else if(object$model == "normal" & object$constant == "CV"){

    normConstCVCI(newdata, object$Umax, ...)

  } else if(object$model == "lognormal" & object$constant == "CV"){

    lognormConstCVCI(newdata, object$Umax, ...)

  }

}


#' Normal constant CV model
#'
#' @importFrom stats pnorm qchisq qnorm qt uniroot var
#' @param y Observed value
#' @param tau Upper limit of CI for the CV
#' @param conf.level Confidence level
#' @param eps Limit

normConstCVCI <- function(y, tau, conf.level = .6827, eps=.Machine$double.eps^.25){

  alpha <- 1 - conf.level
  if (tau * qnorm(1-alpha)>=1){
    Upper <- Inf
    Lower <- y / (1 + tau * qnorm(conf.level))
  } else {
    ##  make log-centered confidence interval
    ## want CI of form: log(y) +/-  r
    ## So we want
    ## rL= rU
    ## rL = log(1+tau*qnorm(1-alphaL))
    ## rU = - log(1-tau*qnorm(1-alphaU))
    ##   where alphaL+alphaU = alpha
    ## Make a rootfunc on log(alphaL)
    rootfunc<-function(logaL){
      # Note: qnorm(10^-18) = 1-qnorm(1-10^-18)
      #   computationally, left side is more accurate
      #   right side goes to -Inf too quick
      -log(1 - tau*qnorm(exp(logaL))) - log(1 + tau*qnorm(alpha-exp(logaL)))

    }
    ## find limits, no real lower limt
    ## set to -700 so that aL approx 2 x 10^-9
    logaL.lower<- log(.Machine$double.xmin)+10
    ## For upper limits, we need
    ## 1.  aL <= alpha
    ## 2. 1+ tau*qnorm(1-aL) > 0    =>  aL < pnorm(1/tau)
    ## 3. 1 - tau*qnorm(1-alpha+aL) > 0 => aL < -1 + alpha + pnorm(1/tau)
    logaL.upper<- log( min(alpha-eps, pnorm(1/tau), -1+alpha+pnorm(1/tau)) - eps )
    if (rootfunc(logaL.lower)>0){
      Lower<- 0
      Upper<- y/(1+tau*qnorm(alpha))
    } else {
      logaL<- uniroot( rootfunc, c(logaL.lower, logaL.upper) )$root
      aL<- exp(logaL)
      aU<- alpha - aL
      ## check that rL=rU
      ## rL<- log(1-tau*qnorm(aL))
      ## rU<- - log(1+tau*qnorm(aU))
      Upper <- y / (1+tau*qnorm(aU))
      Lower <- y / (1-tau*qnorm(aL))
    }
  }

  #c(loglower=log(Lower),logy=log(y),logupper=log(Upper))
  data.frame(obs = y, lower = Lower, upper = Upper)
}


#' log-normal constant CV model
#'
#' @param y Observed value
#' @param tau Upper limit of CI for the CV
#' @param conf.level Confidence level

lognormConstCVCI<- function(y, tau, conf.level=.6827){
  # get log-centered confidence interval
  alpha<- 1 - conf.level
  eta<-log(tau ^ 2 + 1)
  rootfunc<-function(fL){
    -sqrt(eta)* qnorm(1-alpha*fL)+eta/2 +
      sqrt(eta)* qnorm(1-alpha*(1-fL))+eta/2
  }
  fL<- uniroot(rootfunc,c(.5,0))$root
  rL<- qnorm(1-alpha*fL)*sqrt(eta) - eta/2
  rU<- qnorm(1-alpha*(1-fL))*sqrt(eta) + eta/2
  r<- mean(c(rL,rU))
  #r10<- log10( exp(r) )
  #c(check=rL-rU,r=r,lower=y*exp(-r),
  # upper=y*exp(r),loglower=log(y)-r,logy=log(y),logupper=log(y)+r,
  # r10=r10,
  # log10lower=log10(y) - r10,
  # log10y=log10(y),
  # log10upper=log10(y) + r10)

  data.frame(obs = y, lower=y*exp(-r), upper=y*exp(r) )
}


#' Growth Inhibition Assay
#'
#' Data from a growth inhibition assay experiment. Samples were run repeatedly
#' on different assays, for two different types of parasites. gia is the value
#' of interest, and the meanAAgia is the sample level mean, which is the best
#' estimate of the "true" gia level for that sample. varAAgia is the sample level
#' variance.
#'
#' @format A data frame with six variables: \code{parasite}, \code{assay},
#'   \code{plate}, \code{elisa}, \code{gia}, \code{sample}, \code{meanAAgia},
#'   and \code{varAAgia}
#' @aliases gia

"gia"
