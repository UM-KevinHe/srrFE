#=== HYPOTHESIS TESTS FUNCTION =================================
#' Conduct hypothesis testing for identifying outlying providers
#'
#' @param fit an object as output of \code{logis_fe} function.
#'
#' @param parm specify a subset of which providers are to be given confidence intervals. All providers are included by default.
#'
#' @param level confidence level used for constructing confidence intervals. Defaulting to 0.95.
#'
#' @param test a character string specifying the type of testing method to be conducted. Defaulting to "exact.poisbinom".
#'   \itemize{
#'   \item"exact.poisbinom": two-sided exact test based on Poisson-binomial distribution of \eqn{O_i|Z_i}
#'   \item"exact.bootstrap": two-sided exact test based on bootstrap procedure
#'   \item"wald": wald test
#'   \item"score": score test
#'   }
#'
#' @param null a character string or real number specifying null hypotheses of fixed provider effects.
#'
#' @param n resample size for bootstrapping. Defaulting to 10,000.
#'
#' @param ...
#'
#'
#'
#' @details
#' By default, the function uses the `"exact.poisbinom"` method.
#' The wald test is invalid for extreme providers (i.e. when provider effect goes to infinity).
#'
#'
#'
#' @return a dataframe containing:
#'
#' \item{flag}{a vector of flagging indicator. "1" means statistically higher than expected,
#' and "-1" means statistically lower than expected}
#'
#' \item{p}{p-value}
#'
#' \item{stat}{z-score}
#'
#'
#'
#' @examples
#' data(data_FE)
#' data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID, message = FALSE)
#' fit_fe <- logis_fe(data.prep)
#' test_fe(fit_fe)
#'
#' @importFrom stats plogis qnorm pnorm rbinom
#' @importFrom poibin ppoibin
#'
#'
#' @references
#' \itemize{
#' \item Wu, W, Yang, Y, Kang, J, He, K. (2022) Improving large-scale estimation and inference for profiling health care providers.
#' \emph{Statistics in Medicine}, \strong{41(15)}: 2840-2853.
#' }
#'
#' @export

test_fe <- function(fit, parm, level = 0.95, test = "exact.poisbinom", null = "median", n = 10000) {
  if (missing(fit)) stop ("Argument 'fit is required!", call.=F)
  if (!class(fit) %in% c("logis_fe")) stop("Object fit is not of the classes 'logis_fe'!", call.=F)
  if (!(test %in% c("exact.binom", "exact.poisbinom", "exact.bootstrap", "score", "wald")))
    stop("Argument 'test' NOT as required!",call.=F)
  alpha <- 1 - level

  Y.char <- fit$char_list$Y.char
  Z.char <- fit$char_list$Z.char
  prov.char <- fit$char_list$prov.char
  data <- fit$data_include[, c(Y.char, Z.char, prov.char)] #already sorted by provider ID

  gamma <- fit$df.prov$gamma_est #not use the potential Inf of gamma here
  beta <- fit$beta
  gamma.null <- ifelse(null=="median", median(gamma),
                       ifelse(is.numeric(null), null[1],
                              stop("Argument 'null' NOT as required!",call.=F)))
  if (missing(parm)) {
    # pass
  } else if (class(parm)==class(data[,prov.char]) & test!="wald") {
    data <- data[data[,prov.char] %in% parm,]
  } else if (class(parm)==class(data[,prov.char]) & test=="wald") {
    indices <- which(unique(data[,prov.char]) %in% parm)
  } else {
    stop("Argument 'parm' includes invalid elements!")
  }

  if (test=="exact.bootstrap") { #should consistent to exact.poisbinom method
    if (n<=0 | as.integer(n)!=n) stop("Argument 'n' NOT a positive integer!",call.=F)
    exact.bootstrap <- function(df, n) {
      probs <- plogis(gamma.null + unname(as.matrix(df[, Z.char])) %*% beta)
      obs <- sum(df[,Y.char])
      sums <- colSums(matrix(rbinom(n=length(probs)*n, size=1, prob=rep(probs,times=n)), ncol=n))
      p <- (sum(sums>obs)+ 0.5*sum(sums==obs))/n
      z.score <- qnorm(p, lower=F)
      flag <- ifelse(p<alpha/2, 1, ifelse(p<=1-alpha/2, 0, -1))
      p.val <- 2 * min(p, 1-p)
      return(c(flag, p.val, z.score))
    }
    results <- sapply(by(data, data[,prov.char],identity),
                      FUN=function(x) exact.bootstrap(x, n))
    return(data.frame(flag=factor(results[1,]),
                      p=results[2,],
                      stat=results[3,],
                      row.names=unique(data[, prov.char])))
  } else if (test=="score") {
    probs <- plogis(gamma.null+unname(as.matrix(data[,Z.char]))%*%beta)
    z.score <- sapply(split(data[,Y.char]-probs,data[,prov.char]),sum) /
      sqrt(sapply(split(probs*(1-probs),data[,prov.char]),sum))
    p <- pnorm(z.score, lower=F)
    flag <- ifelse(p<alpha/2, 1, ifelse(p<=1-alpha/2, 0, -1))
    p.val <- 2 * pmin(p, 1-p)
    return(data.frame(flag=factor(flag),
                      p=p.val,
                      stat=z.score,
                      row.names=unique(data[, prov.char])))
  } else if (test=="exact.poisbinom") {
    exact.poisbinom <- function(df) {
      probs <- plogis(gamma.null + unname(as.matrix(df[, Z.char])) %*% beta)
      obs <- sum(df[,Y.char])
      p <- 1 - poibin::ppoibin(obs, probs) + 0.5*poibin::dpoibin(obs, probs)  #"ppoibin": probability of "#event <= obs"
      z.score <- qnorm(p, lower=F)
      flag <- ifelse(p<alpha/2, 1, ifelse(p<=1-alpha/2, 0, -1))
      p.val <- 2 * min(p, 1-p)
      return(c(flag, p.val, z.score))
    }
    results <- sapply(by(data, data[,prov.char],identity),
                      FUN=function(x) exact.poisbinom(x))
    return(data.frame(flag=factor(results[1,]),
                      p=results[2,],
                      stat=results[3,],
                      row.names=unique(data[, prov.char])))
  } else if (test=="wald") { # invalid in presence of outlying providers
    if(!missing(parm)){
      if (sum(!is.finite(gamma[indices])) != 0){
        stop("wald test cannot be performed on providers with zero or all events!!")
      }
    } else {
      if (sum(!is.finite(gamma)) != 0){
        stop("wald test cannot be performed on providers with zero or all events!!")
      }
    }

    gamma.obs <- rep(gamma, sapply(split(data[,Y.char],data[,prov.char]),length)) #find gamma-hat
    probs <- as.numeric(plogis(gamma.obs+as.matrix(data[,Z.char])%*%beta))
    info.gamma.inv <- 1/sapply(split(probs*(1-probs), data[,prov.char]),sum) #I_11^-1
    info.betagamma <- sapply(by(probs*(1-probs)*as.matrix(data[,Z.char]),data[,prov.char],identity),colSums) #I_21
    info.beta <- t(as.matrix(data[,Z.char]))%*%(probs*(1-probs)*as.matrix(data[,Z.char]))
    schur.inv <- solve(info.beta-info.betagamma%*%(info.gamma.inv*t(info.betagamma))) # inv of Schur complement; S^-1
    if (missing(parm)) {
      mat.tmp <- info.gamma.inv*t(info.betagamma) #J_1^T
      names <- unique(data[, prov.char])
    } else {
      mat.tmp <- info.gamma.inv[indices]*t(info.betagamma[,indices])
      info.gamma.inv <- info.gamma.inv[indices]
      gamma <- gamma[indices]
      names <- unique(data[, prov.char])[indices]
    }
    se.gamma <- sqrt(info.gamma.inv+apply(mat.tmp, 1, FUN=function(x) t(matrix(x))%*%schur.inv%*%matrix(x))) #only diagonal elements of (1,1) block of I^-1(\theta)
    stat <- (gamma-gamma.null)/se.gamma
    p <- pnorm(stat, lower=F)
    flag <- ifelse(p<alpha/2, 1, ifelse(p<=1-alpha/2, 0, -1))
    p.val <- 2 * pmin(p, 1-p)
    return(data.frame(flag=factor(flag),
                      p=p.val,
                      stat=stat,
                      row.names=names))
  }
}









