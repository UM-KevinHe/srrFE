#=== CONFIDENCE INTERVAL FUNCTION =================================
#' Provide confidence interval for provider effects or standardization ratios/rates
#'
#' @param fit an object as output of \code{logis_fe} function.
#'
#' @param parm specify a subset of which providers are to be given confidence intervals. All providers are included by default.
#'
#' @param level confidence level used for constructing confidence intervals. Defaulting to 0.95.
#'
#' @param option the confidence interval for the function's output, whether it is for gamma or standardization ratios/rates.
#'   \itemize{
#'   \item "gamma": provider effect
#'   \item "SR": standardization ratios/rates
#'   }
#'
#' @param test a character string specifying the type of testing method. Defaulting to "exact".
#'   \itemize{
#'   \item "exact": exact test
#'   \item "wald": wald test
#'   \item "score": score test
#'   }
#'
#' @param stdz if option = 'SR', a character string specifying the standardization method. Defaulting to "indirect".
#'   \itemize{
#'   \item "indirect": using indirect standardized method
#'   \item "direct": using direct standardized method
#'   }
#'
#' @param measure if option = 'SR', a boolean indicating whether the output measure is "ratio" or "rate". Both "rate" and "ratio" will be provided by default.
#'   \itemize{
#'   \item "rate": output the standardized rate. The "rate" has been restricted to 0% - 100%.
#'   \item "ratio":  output the standardized ratio
#'   }
#'
#' @param ...
#'
#' @return A dataframe containing the point estimate, and lower and upper bounds of the estimate.
#'
#' @examples
#' data(data_FE)
#' data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID, message = FALSE)
#' fit_fe <- logis_fe(data.prep)
#' confint(fit_fe, option = "gamma")
#' confint(fit_fe, option = "SR", stdz = "direct", measure = "rate")
#'
#' @importFrom stats plogis
#'
#' @exportS3Method confint logis_fe

confint.logis_fe <- function(fit, parm, level = 0.95, test = "exact",
                             option = c("gamma", "SR"),
                             stdz = "indirect", measure = c("rate", "ratio")) {
  if (missing(fit)) stop ("Argument 'fit' is required!",call.=F)
  if (!class(fit) %in% c("logis_fe")) stop("Object fit is not of the classes 'logis_fe'!",call.=F)
  if (! "gamma" %in% option & !"SR" %in% option) stop("Object fit is not of the classes 'logis_fe'!", call.=F)
  alpha <- 1 - level

  Y.char <- fit$char_list$Y.char
  Z.char <- fit$char_list$Z.char
  prov.char <- fit$char_list$prov.char
  gamma <- fit$df.prov$gamma_est
  beta <- fit$beta
  df.prov <- fit$df.prov
  prov.order <- rownames(fit$gamma)

  #confidence of gamma
  confint_fe_gamma <- function(fit, test, parm, alpha) {
    data <- fit$data_include
    if (missing(parm)) {
      # pass
    } else if (class(parm)==class(data[,prov.char]) & test!="wald") {
      data <- data[data[,prov.char] %in% parm,]
    } else if (class(parm)==class(data[,prov.char]) & test=="wald") {
      indices <- which(unique(data[,prov.char]) %in% parm)
    } else {
      stop("Argument 'parm' includes invalid elements!")
    }

    if (test %in% c("score", "exact")) {
      if (test=="score") {
        qnorm.halfalpha <- qnorm(alpha/2, lower=F)
        qnorm.alpha <- qnorm(alpha, lower=F)
        CL.finite <- function(df) {
          prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                         stop("Number of providers involved NOT equal to one!"))
          UL.gamma <- function(Gamma) { #increasing function w.r.t gamma_null
            p <- plogis(Gamma+Z.beta)
            return((Obs - sum(p)) / sqrt(sum(p*(1-p))) + qnorm.halfalpha)
          }
          LL.gamma <- function(Gamma) {
            p <- plogis(Gamma+Z.beta)
            return((Obs-sum(p)) / sqrt(sum(p*(1-p))) - qnorm.halfalpha)
          }
          Obs <- df.prov[prov, "Obs_facility"]
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          gamma.lower <- uniroot(LL.gamma, gamma[prov]+c(-5,0))$root
          gamma.upper <- uniroot(UL.gamma, gamma[prov]+c(0,5))$root
          return_mat <- c(gamma[prov], gamma.lower, gamma.upper)
          return(return_mat)
        }
        CL.no.readm <- function(df) { #only upper bound
          prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          max.Z.beta <- norm(Z.beta, "I")
          UL.gamma <- function(Gamma) {
            p <- plogis(Gamma+Z.beta)
            return(qnorm.alpha-sum(p)/sqrt(sum(p*(1-p))))
          }
          gamma.upper <- uniroot(UL.gamma,(10+max.Z.beta)*c(-1,1))$root
          return_mat <- c(gamma[prov], -Inf, gamma.upper)
          return(return_mat)
        }
        CL.all.readm <- function(df) { #only lower bound
          prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          max.Z.beta <- norm(Z.beta, "I")
          LL.gamma <- function(Gamma) {
            p <- plogis(Gamma+Z.beta)
            return(sum(1-p)/sqrt(sum(p*(1-p)))-qnorm.alpha)
          }
          gamma.lower <- uniroot(LL.gamma,(10+max.Z.beta)*c(-1,1))$root
          return_mat <- c(gamma[prov], gamma.lower, Inf)
          return(return_mat)
        }
      } else {
        CL.finite <- function(df) {
          UL.gamma <- function(Gamma)
            poibin::ppoibin(Obs-1,plogis(Gamma+Z.beta))+0.5*poibin::dpoibin(Obs,plogis(Gamma+Z.beta))-alpha/2
          LL.gamma <- function(Gamma)
            1-poibin::ppoibin(Obs,plogis(Gamma+Z.beta))+0.5*poibin::dpoibin(Obs,plogis(Gamma+Z.beta))-alpha/2
          prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Obs <- df.prov[prov, "Obs_facility"]
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          gamma.lower <- uniroot(LL.gamma, gamma[prov]+c(-5,0))$root
          gamma.upper <- uniroot(UL.gamma, gamma[prov]+c(0,5))$root
          return_mat <- c(gamma[prov], gamma.lower, gamma.upper)
          return(return_mat)
        }
        CL.no.readm <- function(df) {
          prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          max.Z.beta <- norm(Z.beta, "I")
          gamma.upper <- uniroot(function(x) prod(plogis(-x-Z.beta))/2-alpha,
                                 (10+max.Z.beta)*c(-1,1))$root
          return_mat <- c(gamma[prov], -Inf, gamma.upper)
          return(return_mat)
        }
        CL.all.readm <- function(df) {
          prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                         stop("Number of providers involved NOT equal to one!"))
          Z.beta <- as.matrix(df[,Z.char])%*%beta
          max.Z.beta <- norm(Z.beta, "I")
          gamma.lower <- uniroot(function(x) prod(plogis(x+Z.beta))/2-alpha,
                                 (10+max.Z.beta)*c(-1,1))$root
          return_mat <- c(gamma[prov], gamma.lower, Inf)
          return(return_mat)
        }
      }
      confint.finite <- sapply(by(data[(data$no.readm==0) & (data$all.readm==0),],
                                  data[(data$no.readm==0) & (data$all.readm==0),prov.char],identity),
                               FUN=function(df) CL.finite(df))
      confint.no.readm <- sapply(by(data[data$no.readm==1,], data[data$no.readm==1,prov.char],identity),
                                 FUN=function(df) CL.no.readm(df))
      confint.all.readm <- sapply(by(data[data$all.readm==1,], data[data$all.readm==1,prov.char],identity),
                                  FUN=function(df) CL.all.readm(df))
      confint_df <- as.numeric(cbind(confint.finite, confint.no.readm, confint.all.readm))
      confint_df <- as.data.frame(matrix(confint_df, ncol = 3, byrow = T))
      colnames(confint_df) <- c("gamma", "gamma.lower", "gamma.upper")
      return(confint_df[order(match(rownames(confint_df), prov.order)),])
    } else if (test=="wald") {
      if(!missing(parm)){
        if (sum(!is.finite(gamma[indices])) != 0){
          stop("wald test cannot be performed on facilities with zero or all readmissions!!")
        }
      } else {
        if (sum(!is.finite(gamma)) != 0){
          stop("wald test cannot be performed on facilities with zero or all readmissions!!")
        }
      }
      n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length)
      gamma.obs <- rep(gamma, n.prov)
      probs <- as.numeric(plogis(gamma.obs+as.matrix(data[,Z.char])%*%beta))
      info.gamma.inv <- 1/sapply(split(probs*(1-probs), data[,prov.char]),sum)
      info.betagamma <- sapply(by(probs*(1-probs)*as.matrix(data[,Z.char]),data[,prov.char],identity),colSums)
      info.beta <- t(as.matrix(data[,Z.char]))%*%(probs*(1-probs)*as.matrix(data[,Z.char]))
      schur.inv <- solve(info.beta-info.betagamma%*%(info.gamma.inv*t(info.betagamma))) # inv of Schur complement
      if (missing(parm)) {
        mat.tmp <- info.gamma.inv*t(info.betagamma)
        names <- unique(data[, prov.char])
        Z.beta <- as.matrix(data[data[,prov.char], Z.char]) %*% beta
      } else {
        mat.tmp <- info.gamma.inv[indices]*t(info.betagamma[,indices])
        info.gamma.inv <- info.gamma.inv[indices]
        gamma <- gamma[indices]
        Z.beta <- as.matrix(data[data[,prov.char] %in% parm, Z.char]) %*% beta
        n.prov <- n.prov[indices]
        names <- unique(data[, prov.char])[indices]
      }
      se.gamma <- sqrt(info.gamma.inv+apply(mat.tmp, 1, FUN=function(x) t(matrix(x))%*%schur.inv%*%matrix(x)))
      gamma.lower <- gamma - qnorm(alpha/2, lower=F)*se.gamma
      gamma.upper <- gamma + qnorm(alpha/2, lower=F)*se.gamma
      return_mat <- data.frame(gamma, gamma.lower, gamma.upper,
                               row.names=names)
      return(return_mat)
    }
  }
  if (option == "gamma"){
    return_mat <- confint_fe_gamma(fit, test, parm, alpha)
    return(return_mat)
  } else if (option == "SR"){
    data.ori <- fit$data_include
    population_rate <- sum(data.ori[,Y.char])/nrow(data.ori) * 100  #sum(O_i)/N *100%
    return_ls <- list()
    if ("indirect" %in% stdz) {
      SR.indirect <- SR_output(fit, stdz = c("indirect"), measure = c("ratio", "rate"))
      if (missing(parm)) {
        OE_df.indirect <- SR.indirect$OE$OE_indirect
        indirect.ratio_df <- SR.indirect$indirect.ratio
        indirect.rate_df <- SR.indirect$indirect.rate
        data <- data.ori
      } else if (class(parm) == class(data.ori[,prov.char])) {
        OE_df.indirect <- SR.indirect$OE$OE_indirect[rownames(SR.indirect$OE$OE_indirect) %in% parm,]
        data <- data.ori[data.ori[,prov.char] %in% parm,]
        indirect.ratio_df <- SR.indirect$indirect.ratio[rownames(SR.indirect$indirect.ratio) %in% parm,]
        indirect.rate_df <- SR.indirect$indirect.rate[rownames(SR.indirect$indirect.rate) %in% parm,]
      } else {
        stop("Argument 'parm' includes invalid elements!")
      }
      #functions for calculate CI of SRs
      qnorm.halfalpha <- qnorm(alpha/2, lower=F)
      qnorm.alpha <- qnorm(alpha, lower=F)
      SR_indirect.finite <- function(df) {
        prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                       stop("Number of providers involved NOT equal to one!"))
        Z.beta <- as.matrix(df[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = unique(df$ID), alpha = alpha)
        gamma.lower <- confint_gamma$gamma.lower
        gamma.upper <- confint_gamma$gamma.upper
        EXP.i <- OE_df.indirect[rownames(OE_df.indirect) == unique(df[,prov.char]), "Exp.indirect_facility"]
        SR.lower <- sum(plogis(gamma.lower+Z.beta)) / EXP.i
        SR.upper <- sum(plogis(gamma.upper+Z.beta)) / EXP.i
        return(c(SR.lower, SR.upper))
      }
      SR_indirect.no.readm <- function(df) {
        prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                       stop("Number of providers involved NOT equal to one!"))
        Z.beta <- as.matrix(df[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = unique(df$ID), alpha = alpha)
        gamma.upper <- confint_gamma$gamma.upper
        EXP.i <- OE_df.indirect[rownames(OE_df.indirect) == unique(df[,prov.char]), "Exp.indirect_facility"]
        SR.upper <- sum(plogis(gamma.upper+Z.beta)) / EXP.i
        return(c(0, SR.upper))
      }
      SR_indirect.all.readm <- function(df) {
        prov <- ifelse(length(unique(df[,prov.char]))==1, unique(df[,prov.char]),
                       stop("Number of providers involved NOT equal to one!"))
        Z.beta <- as.matrix(df[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = unique(df$ID), alpha = alpha)
        gamma.lower <- confint_gamma$gamma.lower
        EXP.i <- OE_df.indirect[rownames(OE_df.indirect) == unique(df[,prov.char]), "Exp.indirect_facility"]
        SR.lower <- sum(plogis(gamma.lower+Z.beta)) / EXP.i
        SR.upper <- nrow(df) / EXP.i
        return(c(SR.lower, SR.upper))
      }

      confint.finite <- sapply(by(data[(data$no.readm==0) & (data$all.readm==0),],
                                  data[(data$no.readm==0) & (data$all.readm==0),prov.char],identity),
                               FUN=function(df) SR_indirect.finite(df))
      confint.no.readm <- sapply(by(data[data$no.readm==1,], data[data$no.readm==1,prov.char],identity),
                                 FUN=function(df) SR_indirect.no.readm(df))
      confint.all.readm <- sapply(by(data[data$all.readm==1,], data[data$all.readm==1,prov.char],identity),
                                  FUN=function(df) SR_indirect.all.readm(df))


      CI.indirect_ratio <- as.numeric(rbind(t(indirect.ratio_df),
                                            cbind(confint.finite,
                                                  confint.no.readm,
                                                  confint.all.readm)))
      CI.indirect_ratio <- as.data.frame(matrix(CI.indirect_ratio, ncol = 3, byrow = T))
      colnames(CI.indirect_ratio) <- c("indirect_ratio", "CI_ratio.lower", "CI_ratio.upper")
      rownames(CI.indirect_ratio) <- names(indirect.rate_df)

      if ("ratio" %in% measure){
        return_ls$CI.indirect_ratio <- CI.indirect_ratio

      }

      if ("rate" %in% measure){
        rate.lower <- pmax(pmin(CI.indirect_ratio$CI_ratio.lower * population_rate, 100), 0)
        rate.upper <- pmax(pmin(CI.indirect_ratio$CI_ratio.upper * population_rate, 100), 0)
        CI.indirect_rate <- as.data.frame(cbind(indirect.rate_df, rate.lower, rate.upper))
        colnames(CI.indirect_rate) <- c("indirect_rate", "CI_rate.lower", "CI_rate.upper")
        return_ls$CI.indirect_rate <- CI.indirect_rate
      }
    }


    if ("direct" %in% stdz) {
      SR.direct <- SR_output(fit, stdz = c("direct"), measure = c("ratio", "rate"))
      if (missing(parm)) {
        OE_df.direct <- SR.direct$OE$OE_direct[1,1]
        direct.ratio_df <- SR.direct$direct.ratio
        direct.rate_df <- SR.direct$direct.rate
        data <- data.ori
      } else if (class(parm) == class(data[,prov.char])) {
        OE_df.direct <- SR.direct$OE$OE_direct[1,1]
        direct.ratio_df <- SR.direct$direct.ratio[rownames(SR.direct$direct.ratio) %in% parm,]
        direct.rate_df <- SR.direct$direct.rate[rownames(SR.direct$direct.rate) %in% parm,]
        data <- data.ori[data.ori[,prov.char] %in% parm,]
      } else {
        stop("Argument 'parm' includes invalid elements!")
      }
      #funcitons for calculate CI of SRs
      qnorm.halfalpha <- qnorm(alpha/2, lower=F)
      qnorm.alpha <- qnorm(alpha, lower=F)

      SR_direct.finite <- function(ID) {
        Z.beta.all <- as.matrix(data.ori[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = ID, alpha = alpha)
        gamma.lower <- confint_gamma$gamma.lower
        gamma.upper <- confint_gamma$gamma.upper
        SR.lower <- sum(plogis(gamma.lower+Z.beta.all)) / OE_df.direct
        SR.upper <- sum(plogis(gamma.upper+Z.beta.all)) / OE_df.direct
        return(c(SR.lower, SR.upper))
      }
      SR_direct.no.readm <- function(ID) {
        Z.beta.all <- as.matrix(data.ori[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = ID, alpha = alpha)
        gamma.upper <- confint_gamma$gamma.upper
        SR.upper <- sum(plogis(gamma.upper+Z.beta.all)) / OE_df.direct
        return(c(0, SR.upper))
      }
      SR_direct.all.readm <- function(ID) {
        Z.beta.all <- as.matrix(data.ori[,Z.char])%*%beta
        confint_gamma <- confint_fe_gamma(fit, test = test, parm = ID, alpha = alpha)
        gamma.lower <- confint_gamma$gamma.lower
        SR.lower <- sum(plogis(gamma.lower+Z.beta.all)) / OE_df.direct
        SR.upper <- nrow(data.ori) / OE_df.direct
        return(c(SR.lower, SR.upper))
      }

      confint.finite <- sapply(unique(data[(data$no.readm==0) & (data$all.readm==0),]$ID),
                               FUN = function(ID) SR_direct.finite(ID))
      confint.no.readm <- sapply(unique(data[(data$no.readm==1) & (data$all.readm==0),]$ID),
                                 FUN = function(ID) SR_direct.no.readm(ID))
      confint.all.readm <- sapply(unique(data[(data$no.readm==0) & (data$all.readm==1),]$ID),
                                  FUN = function(ID) SR_direct.all.readm(ID))
      CI.direct_ratio <- as.numeric(rbind(t(direct.ratio_df),
                                          cbind(confint.finite,
                                                confint.no.readm,
                                                confint.all.readm)))

      CI.direct_ratio <- as.data.frame(matrix(CI.direct_ratio, ncol = 3, byrow = T))
      colnames(CI.direct_ratio) <- c("direct_ratio", "CI_ratio.lower", "CI_ratio.upper")
      rownames(CI.direct_ratio) <- names(direct.rate_df)

      if ("ratio" %in% measure){
        return_ls$CI.direct_ratio <- CI.direct_ratio

      }

      if ("rate" %in% measure){
        rate.lower <- pmax(pmin(CI.direct_ratio$CI_ratio.lower * population_rate, 100), 0)
        rate.upper <- pmax(pmin(CI.direct_ratio$CI_ratio.upper * population_rate, 100), 0)
        CI.direct_rate <- as.data.frame(cbind(direct.rate_df, rate.lower, rate.upper))
        colnames(CI.direct_rate) <- c("direct_rate", "CI_rate.lower", "CI_rate.upper")
        return_ls$CI.direct_rate <- CI.direct_rate
      }
    }
    return(return_ls)
  }


}
