#=== SUMMARY STATISTICS FOR COVARIATE ESTIMATES =================================
#' Provide the summary statistics for covariate estimates
#'
#' @param fit an object as output of \code{logis_fe} function.
#'
#' @param parm a character vector specifies a subset of covariates. All covariates are included by default.
#'
#' @param level confidence level used for constructing confidence intervals. Defaulting to 0.95.
#'
#' @param test a character string specifying the type of testing method. Defaulting to "wald.cpp".
#'   \itemize{
#'     \item "wald": wald test
#'     \item "wald.cpp": wald test using cpp function
#'     \item "lr": likelihood ratio test
#'     \item "score": score test
#'   }
#'
#' @param null the null value of the covariate estimate that requires testing. (e.g. test \eqn{H_0: \beta = 0})
#'
#' @param ...
#'
#'
#' @return a dataframe containing summary statistics for covariate estimates
#'
#' @examples
#' data(data_FE)
#' data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID, message = FALSE)
#' fit_fe <- logis_fe(data.prep)
#' summary.wald.cpp <- summary_fe_covar(fit_fe, level = 0.05, test = "wald")
#' summary.wald.cpp
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats plogis pnorm pchisq
#'
#'
#' @export

summary_fe_covar <- function(fit, parm, level = 0.95, test = "wald.cpp", null = 0) {
  Y.char <- fit$char_list$Y.char
  Z.char <- fit$char_list$Z.char
  prov.char <- fit$char_list$prov.char
  alpha <- 1 - level

  if (missing(parm)) {
    ind <- 1:length(Z.char)
  } else if (is.character(parm)) {
    ind <- which(Z.char %in% parm)
  } else if (is.numeric(parm) & max(abs(as.integer(parm)-parm))==0 & !(0 %in% parm)) {
    ind <- parm
  } else {
    stop("Argument 'parm' includes invalid elements!")
  }
  if (test=="wald") { # Wald test and test-based CIs
    data <- fit$data_include[, c(Y.char,Z.char,prov.char)]
    gamma.obs <- rep(fit$df.prov$gamma_est, sapply(split(data[,Y.char],data[,prov.char]),length))
    probs <- as.numeric(plogis(gamma.obs+as.matrix(data[,Z.char]) %*% fit$beta))
    probs <- pmin(pmax(probs,1e-10),1-1e-10)
    info.gamma.inv <- 1/sapply(split(probs*(1-probs), data[,prov.char]),sum)
    info.betagamma <- sapply(by(probs*(1-probs)*as.matrix(data[,Z.char]),data[,prov.char],identity),colSums)
    info.beta <- t(as.matrix(data[,Z.char]))%*%(probs*(1-probs)*as.matrix(data[,Z.char]))
    se.beta <- sqrt(diag(solve(info.beta-info.betagamma%*%(info.gamma.inv*t(info.betagamma)))))[ind] #S^-1
    p <- pnorm((fit$beta[ind]-null)/se.beta, lower=F)
    df <- data.frame(beta = fit$beta[ind],
                     se.beta = se.beta,
                     p=2*pmin(p,1-p),
                     CI.lower=fit$beta[ind]-qnorm(1-alpha/2)*se.beta,
                     CI.upper=fit$beta[ind]+qnorm(1-alpha/2)*se.beta)
    row.names(df) <- Z.char[ind]
    return(df)
  } else if (test=="wald.cpp") {  #use cpp function to calculate Cov(\beta)
    data <- fit$data_include[, c(Y.char,Z.char,prov.char)]
    n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length)
    ls <- wald_covar(data[,Y.char], as.matrix(data[,Z.char]), n.prov, fit$df.prov$gamma_est, as.numeric(fit$beta), ind, null, alpha)
    #ls$p <- c(ls$p); ls$stat <- c(ls$stat); ls$se.beta <- c(ls$se.beta)
    #ls$beta.lower <- c(ls$beta.lower); ls$beta.upper <- c(ls$beta.upper)
    df <- data.frame(beta = fit$beta[ind],
                     se.beta=c(ls$se.beta),
                     p=c(ls$p),
                     #stat=ls$stat,
                     CI.lower=c(ls$beta.lower),
                     CI.upper=c(ls$beta.upper))
    row.names(df) <- Z.char[ind]
    return(df)
  } else if (test=="lr") { # Note: only null=0 allowed for now and no CIs
    if (null!=0) stop("Argument 'null' is invalid!")
    data <- fit$data_include
    gamma.obs <- rep(pmax(pmin(fit$df.prov$gamma_est,median(fit$df.prov$gamma_est)+10),median(fit$df.prov$gamma_est)-10), sapply(split(data[,Y.char],data[,prov.char]),length))
    neg2Loglkd <- -2*sum((gamma.obs+as.matrix(data[,Z.char])%*%fit$beta)*data[,Y.char]-log(1+exp(gamma.obs+as.matrix(data[,Z.char])%*%fit$beta)))
    lr <- function(index) {
      data.null <- as.data.frame(cbind(data[,Y.char], data[,prov.char], data[,Z.char[-index]], data[, (ncol(data) - 2):ncol(data)]))
      char_list.null <- fit$char_list
      char_list.null$Z.char <- char_list.null$Z.char[-index]
      colnames(data.null)[1:2] <- c(char_list.null$Y.char, char_list.null$prov.char)
      data.prep.null <- list(data = data.null,
                             char_list = char_list.null)
      fe.null <- logis_fe(data.prep.null)
      Z.null <- as.matrix(data[,Z.char[-index]])
      gamma.obs.null <- rep(pmax(pmin(fe.null$df.prov$gamma,median(fe.null$df.prov$gamma)+10),median(fe.null$df.prov$gamma)-10), sapply(split(data[,Y.char],data[,prov.char]),length))
      neg2Loglkd.null <- -2*sum((gamma.obs.null+Z.null%*%fe.null$beta)*data[,Y.char]-log(1+exp(gamma.obs.null+Z.null%*%fe.null$beta)))
      p <- pchisq(neg2Loglkd.null-neg2Loglkd, 1, lower=F)
      return(p)
    }
    df <- data.frame(beta = fit$beta[ind],
                     p = sapply(ind, lr))
    rownames(df) <- Z.char[ind]
    return(df)
  } else if (test=="score") { # Note: only null=0 allowed for now and no CIs
    if (null!=0) stop("Argument 'null' is invalid!")
    data <- fit$data_include
    score <- function(index) {
      data.null <- as.data.frame(cbind(data[,Y.char], data[,prov.char], data[,Z.char[-index]], data[, (ncol(data) - 2):ncol(data)]))
      char_list.null <- fit$char_list
      char_list.null$Z.char <- char_list.null$Z.char[-index]
      colnames(data.null)[1:2] <- c(char_list.null$Y.char, char_list.null$prov.char)
      data.prep.null <- list(data = data.null,
                             char_list = char_list.null)
      fe.null <- logis_fe(data.prep.null)
      Z.null <- as.matrix(data[,Z.char[-index]])
      gamma.obs.null <- rep(fe.null$df.prov$gamma, sapply(split(data[,Y.char],data[,prov.char]),length))
      probs.null <- as.numeric(plogis(gamma.obs.null+as.matrix(data[,Z.char[-index]])%*%fe.null$beta))
      probs.null <- pmin(pmax(probs.null,1e-10),1-1e-10)
      info.gamma.inv.null <- 1/sapply(split(probs.null*(1-probs.null), data[,prov.char]),sum)
      info.betagamma.null <- sapply(by(probs.null*(1-probs.null)*Z.null,data[,prov.char],identity),colSums)
      info.beta.null <- t(Z.null)%*%(probs.null*(1-probs.null)*Z.null)
      mat.tmp <- info.gamma.inv.null*t(info.betagamma.null)
      schur.inv.null <- solve(info.beta.null-info.betagamma.null%*%mat.tmp)
      info.inv.11 <- mat.tmp%*%schur.inv.null%*%t(mat.tmp)
      diag(info.inv.11) <- diag(info.inv.11) + info.gamma.inv.null
      info.inv.21 <- -schur.inv.null%*%t(mat.tmp)
      info.beta.add.gamma <- sapply(by(probs.null*(1-probs.null)*data[,Z.char[index]],data[,prov.char],identity),sum)
      info.beta.add.beta <- t(as.matrix(data[,Z.char[index]]))%*%(probs.null*(1-probs.null)*Z.null)
      info.beta.add <- t(data[,Z.char[index]])%*%(probs.null*(1-probs.null)*data[,Z.char[index]]) -
        t(info.beta.add.gamma)%*%info.inv.11%*%info.beta.add.gamma -
        2*info.beta.add.beta%*%info.inv.21%*%info.beta.add.gamma -
        info.beta.add.beta%*%schur.inv.null%*%t(info.beta.add.beta)
      score.beta.add <- t(as.matrix(data[,Z.char[index]]))%*%(data[,Y.char]-probs.null)
      score.stat <- score.beta.add^2/info.beta.add
      p <- pchisq(score.stat, 1, lower=F)
      return(p)
    }
    df <- data.frame(beta = fit$beta[ind],
                     p = sapply(ind, score))
    rownames(df) <- Z.char[ind]
    return(df)
  }
}
