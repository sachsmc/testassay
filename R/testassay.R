#' Hypothesis testing procedure for assay validation for precision
#'
#' Does an m:n:q procedure for assay validation for precision.
#' Returns an object of class 'assaytest'. There is are \code{\link[=predict.assaytest]{predict}}
#' and \code{\link[=print.assaytest]{print}} methods for that class.
#'
#' @details
#' The m:n:q procedure uses m different samples that have different levels of the true value
#' with n replicates for each sample. The output is a 100q percent upper limit of the bound on the
#' precision parameter when the true values within the range of values for the m samples all follow either a
#' a constant coefficient of variation model or a  constant
#' standard deviation model (same as a constant variance model) (see \code{constant} argument).
#'
#' For example, if the 4:4:90 percent procedure using a normal model with a constant variance model
#' returns a bound on the standard deviation  (the Umax element of the assaytest class) of 7.9
#' then under the assumptions we have 90 percent confidence that the true SD is less than 7.9.
#'
#' The \code{\link[=predict.assaytest]{predict}} method gives effective standard deviation intervals (i.e., 68.27 pct CIs)
#' for the expected response from subsequent observed values from the assay.
#'
#'
#' @param x The vector of assay values
#' @param m The vector of values indicating sample membership
#' @param n The vector of values indicating replicate membership
#' @param q The confidence level, typically 0.8 or 0.9
#' @param model String specifying the distribution for the assay values. Valid
#'   values are "normal" or "lognormal"
#' @param constant String specifying whether the standard deviation is assumed to be
#'   constant over the levels ("SD") or the coefficient of variation
#'   is assumed constant over the levels ("CV"). The values "sd", "var", or "variance" may be used for "SD", and "cv" may be used
#'    for "CV".
#' @param data Data frame or environment in which to look for x
#'
#' @export
#' @aliases testassay
#' @return An object of class "assaytest", which is a list of components including a data frame of the relevant statistics calculated on x. Print, summary, predict, and plot methods are available.
#' The list has the following components
##' \itemize{
##'  \item{sumtab}{Table summarizing the experiment, includes mean values, SD or CV estimates, and upper confidence limits on those. }
##'  \item{Umax}{The maximum of the upper limits on the SD or CV, used in the effective SD interval calculation}
##'  \item{n}{The number of samples per level}
##'  \item{m}{The number of levels}
##'  \item{q}{The confidence level}
##'  \item{model}{The assumed model}
##'  \item{constant}{The parameter assumed to be constant (either 'SD' or 'CV').  }
##'  \item{alpha}{The alpha level, calculated as (1 - q)^(1 / m)}
##'  \item{x}{The data vector supplied by the user}
##' }
#'
#' @references
#' Fay, MP, Sachs, MC, and Miura, K (2016). A Hypothesis Testing Framework for Validating and Assay for Precision
#' (unpublished manuscript).
#'
#' @examples
#' # reproduce Table 3 of Fay, Sachs and Miura
#' I<- gia$parasite=="3D7" & gia$meanAAgia<80
#' treD7.test<-testassay(x=gia, m=sample, n=assay, q=.9,
#'   data=subset(gia, parasite=="3D7" & meanAAgia<80))
#' treD7.test
#' # get estimated effective standard deviation intervals (68.27 percent CIs)
#' # for observed values 21.4 and 65.9
#' # using results from testassay
#' predict(treD7.test,c(21.4,65.9))
#'
testassay <- function(x, m, n, q = .9, model = "normal", constant = "SD", data = NULL) {

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
#' These intervals have at least 68.27 percent coverage.
#'
#' @details
#' Takes the \code{Umax} element from the \code{assaytest} object and treats it as the known precision
#' parameter. For the constant SD model, the effective standard deviation interval for observed data
#' value y is (y-Umax, y+Umax). For the constant CV models the effective SD interval uses either
#' \code{\link{normConstCVCI}} (for the "normal" model) or  \code{\link{lognormConstCVCI}} (for the "lognormal" model).
#'
#' Although \code{Umax} is an upper bound (not an estimate) of the precision parameter, simulations have shown
#' that treating   \code{Umax} as the true precision parameter
#' gives effective SD intervals with coverage of at least 68.27 percent (see Fay, Sachs, and Miura, 2016).
#'
#' @references
#' Fay, MP, Sachs, MC, and Miura, K (2016). A Hypothesis Testing Framework for Validating and Assay for Precision
#' (unpublished manuscript).
#'
#' @export
#'
#' @param object An object of class "assaytest"
#' @param newdata A vector of observed values. If missing, uses object$x.
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


#' Log-centered confidence intervals from a Normal constant coeffficient of variation model
#'
#' Assume Y is normal with mean mu>0 and coefficient of variation theta, then Y/mu ~ N(1, theta^2).
#' Get log-centered confidence intervals (when possible), meaning intervals such that log(y) +/- r(theta), where
#' r(theta) is a constant function of theta.
#'
#' @importFrom stats pnorm qchisq qnorm qt uniroot var
#' @param y vector of observed values, should be positive
#' @param theta coefficient of variation (assumed known)
#' @param conf.level Confidence level
#' @param eps a small number used in the algorithm (look at code before changing)
#'
#'
#' @return A list with the following components
##' \itemize{
##'  \item{obs}{ y }
##'  \item{lower}{ lower confidence limit on mu=E(Y)}
##'  \item{upper}{ upper confidence limit on mu=E(Y)}
##' }
#'
#' @examples
#' # defaults to 68.27 percent confidence level, same level as Normal plus or minus 1 std dev.
#' normConstCVCI(3.4,.6)
#' # symmetric on log scale
#' log(normConstCVCI(3.4,.6))
#'
#' @export
normConstCVCI <- function(y, theta, conf.level = .6827, eps=.Machine$double.eps^.25){

  alpha <- 1 - conf.level
  if (theta * qnorm(1-alpha)>=1){
    Upper <- Inf
    Lower <- y / (1 + theta * qnorm(conf.level))
  } else {
    ##  make log-centered confidence interval
    ## want CI of form: log(y) +/-  r
    ## So we want
    ## rL= rU
    ## rL = log(1+theta*qnorm(1-alphaL))
    ## rU = - log(1-theta*qnorm(1-alphaU))
    ##   where alphaL+alphaU = alpha
    ## Make a rootfunc on log(alphaL)
    rootfunc<-function(logaL){
      # Note: qnorm(10^-18) = 1-qnorm(1-10^-18)
      #   computationally, left side is more accurate
      #   right side goes to -Inf too quick
      -log(1 - theta*qnorm(exp(logaL))) - log(1 + theta*qnorm(alpha-exp(logaL)))

    }
    ## find limits, no real lower limt
    ## set to -700 so that aL approx 2 x 10^-9
    logaL.lower<- log(.Machine$double.xmin)+10
    ## For upper limits, we need
    ## 1.  aL <= alpha
    ## 2. 1+ theta*qnorm(1-aL) > 0    =>  aL < pnorm(1/theta)
    ## 3. 1 - theta*qnorm(1-alpha+aL) > 0 => aL < -1 + alpha + pnorm(1/theta)
    logaL.upper<- log( min(alpha-eps, pnorm(1/theta), -1+alpha+pnorm(1/theta)) - eps )
    if (rootfunc(logaL.lower)>0){
      Lower<- 0
      Upper<- y/(1+theta*qnorm(alpha))
    } else {
      logaL<- uniroot( rootfunc, c(logaL.lower, logaL.upper) )$root
      aL<- exp(logaL)
      aU<- alpha - aL
      ## check that rL=rU
      ## rL<- log(1-theta*qnorm(aL))
      ## rU<- - log(1+theta*qnorm(aU))
      Upper <- y / (1+theta*qnorm(aU))
      Lower <- y / (1-theta*qnorm(aL))
    }
  }

  #c(loglower=log(Lower),logy=log(y),logupper=log(Upper))
  data.frame(obs = y, lower = Lower, upper = Upper)
}


#' log-normal constant CV model
#'
#' This function gets confidence intervals on mu=E(Y) assuming Y is lognormal and the coefficient of variation is known.
#'
#' @param y Observed value
#' @param theta coefficient of variation (assumed known)
#' @param conf.level Confidence level
#'
#' @details
#' Let Y be lognormal, so that log(Y) is normal with mean xi and variance eta.
#' Then E(Y) =mu = exp(xi+eta/2) and
#' Var(Y)=sigma^2 = mu^2 (exp(eta)-1),
#' so that the coefficient of variation is sigma/mu = sqrt( exp(eta)-1).
#' We want to get log-centered confidence intervals on mu, meaning intervals such that
#' log(y) +/- r(theta), where
#' r(theta) is a constant function of theta.
#'
#' @return A list with the following components
##' \itemize{
##'  \item{obs}{ y }
##'  \item{lower}{ lower confidence limit on mu=E(Y)}
##'  \item{upper}{ upper confidence limit on mu=E(Y)}
##' }
#'
#' @examples
#' # defaults to 68.27 percent confidence level, same level as Normal plus or minus 1 std dev.
#' lognormConstCVCI(3.4,.6)
#' # compare to normal constant CV model result
#' normConstCVCI(3.4,.6)
#'
#' @export
lognormConstCVCI<- function(y, theta, conf.level=.6827){
  # get log-centered confidence interval
  alpha<- 1 - conf.level
  eta<-log(theta ^ 2 + 1)
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
#' on different assays, for two different strains of parasites (3d7 and FVO).
#' \code{elisa} is a measure of the amount of antibody and is measured once for each sample.
#' \code{sample} is a unique name for each sample and is defined as \code{paste(gia$parasite,gia$elisa,sep=".")}.
#' \code{gia} is the value
#' of interest, and the \code{meanAAgia} is the sample level mean, which is the best
#' estimate of the "true" gia level for that sample. \code{varAAgia} is the sample level
#' variance.
#'
#' @format A data frame with six variables: \code{parasite}, \code{assay},
#'    \code{elisa}, \code{gia}, \code{sample}, \code{meanAAgia},
#'   and \code{varAAgia}
#' @aliases gia

"gia"
