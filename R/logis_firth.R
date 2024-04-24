#=== FUNCTIONS FOR FITTING THE FIXED EFFECTS MODEL WITH FIRTH'S CORRECTION =================================
#' Main function for fitting fixed effects model with firth's bias correction
#'
#' @param data.prep an object from the `fe_data_prep()` function.
#'
#' @param max.iter maximum number of iterations. Defaulting to 10,000.
#'
#' @param tol a small positive number specifying stopping criterion of Newton-Raphson algorithm. Defaulting to 1e-5.
#'
#' @param bound a positive number to avoid inflation of provider effect. Defaulting to 10.
#'
#' @param backtrack a boolean indicating whether backtracking line search is implemented. Defaulting to FALSE.
#'
#' @param Rcpp a Boolean indicating whether the Rcpp function is used. Defaulting to TRUE.
#'
#' @param AUC a Boolean indicating whether report AUC. Defaulting to FALSE.
#'
#' @param message a Boolean indicating whether track the fitting process. Defaulting to TRUE.
#'
#' @param ...
#'
#'
#' @return An object with S3 class \code{logis_fe}.
#'
#' \item{beta}{a vector of fixed effects estimates of covariates}
#'
#' \item{gamma}{a vector of estimates of provider effects}
#'
#' \item{obs}{a vector of patients-level outcome}
#'
#' \item{neg2Loglkd}{minus two times log likelihood}
#'
#' \item{AIC}{Akaike info criterion}
#'
#' \item{BIC}{Bayesian info criterion}
#'
#' \item{AUC}{area under the ROC curve}
#'
#'
#' @examples
#' data(data_FE)
#' data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID)
#' fit_firth <- logis_firth(data.prep)
#'
#' @importFrom Rcpp evalCpp
#' @importFrom pROC auc
#' @importFrom RcppParallel RcppParallelLibs
#'
#'
#' @keywords Fixed Provider Effects, Firth Correction
#'
#' @export
#'
#'
logis_firth <- function(data.prep, max.iter = 10000, tol = 1e-5, bound = 10,
                        backtrack = TRUE, Rcpp = TRUE, AUC = FALSE, message = FALSE){
  if (missing(data.prep)) stop ("Argument 'data.prep' is required!", call.=F)
  if (!class(data.prep) %in% c("data_prep")) stop("Object 'data.prep' should be generated from 'fe_data_prep' function!", call.=F)

  if (!is.logical(backtrack)) stop("Argument 'backtrack' NOT as required!", call.=F)

  data <- data.prep$data
  Y.char <- data.prep$char_list$Y.char
  prov.char <- data.prep$char_list$prov.char
  Z.char <- data.prep$char_list$Z.char

  data <- data[data$included==1,]
  n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length) # provider-specific number of discharges
  n.events.prov <- sapply(split(data[, Y.char], data[, prov.char]), sum) # provider-specific number of events
  Z <- as.matrix(data[,Z.char])
  gamma.prov <- rep(log(mean(data[,Y.char])/(1-mean(data[,Y.char]))), length(n.prov))
  beta <- rep(0, NCOL(Z))
  m <- length(n.prov)
  n.obs <- nrow(data)


  if (Rcpp) { #Rcpp always use "backtrack"
    ls <- logis_firth_prov(as.matrix(data[, Y.char]), Z, n.prov, gamma.prov, beta,
                           n.obs, m, 0, 1, tol, max.iter,
                           bound, message, backtrack)
    gamma.prov <- as.numeric(ls$gamma)
    beta <- as.numeric(ls$beta)
  } else {
    iter <- 0
    beta.crit <- 100 # initialize stop criterion
    if (message == TRUE){
      message("Implementing firth's bias-corrected fixed provider effects model ...")
    }
    split.Z <- by(Z, data[,prov.char],identity)

    if (backtrack){ # initialize parameters for backtracking line search
      s <- 0.01
      t <- 0.6
      Loglkd <- function(gamma.obs, beta) {
        sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
      }
    }

    while (iter<=max.iter & beta.crit>=tol) {
      iter <- iter + 1
      gamma.obs <- rep(gamma.prov, n.prov)
      p <- c(plogis(gamma.obs+Z%*%beta))
      q <- p*(1-p)
      info.gamma.inv <- 1/sapply(split(q, data[,prov.char]),sum) #I_11^(-1)
      info.betagamma <- sapply(by(q*Z,data[,prov.char],identity),colSums) #I_21
      info.beta <- t(Z)%*%(q*Z) #I_22
      mat.tmp1 <- info.gamma.inv*t(info.betagamma) #J_1^T
      schur.inv <- solve(info.beta-info.betagamma%*%mat.tmp1) #S^-1
      mat.tmp2 <- mat.tmp1%*%schur.inv #J_2^T

      h <- q * (rep(info.gamma.inv + diag(mat.tmp1 %*% schur.inv %*% t(mat.tmp1)), n.prov) + #A_1 B_11 A_1^T
                  2 * do.call(c, lapply(1:m, function(i) {  #A_2 B_21 A_1^T
                    as.matrix(split.Z[[i]]) %*% -t(mat.tmp2)[, i, drop = FALSE]
                  })) +
                  sapply(1:n.obs, function(i) {  #A_2 B_22 A_2^T
                    Z[i, , drop = F] %*% schur.inv %*% t(Z[i, , drop = F])
                  }))
      A <- h * (0.5 - p) #a vector
      score.beta <- t(Z)%*%(data[,Y.char]-p + A)
      score.gamma <- sapply(split(data[,Y.char]-p + A, data[,prov.char]), sum)

      d.gamma.prov <- info.gamma.inv*score.gamma +
        mat.tmp2%*%(t(mat.tmp1)%*%score.gamma-score.beta)
      d.beta <- -t(mat.tmp2)%*%score.gamma+schur.inv%*%score.beta

      v <- 1 # initialize step size
      if (backtrack) {
        loglkd <- Loglkd(rep(gamma.prov, n.prov), beta)
        d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta+v*d.beta) - loglkd
        lambda <- c(score.gamma,score.beta)%*%c(d.gamma.prov,d.beta)
        while (d.loglkd < s*v*lambda) {  #update step size
          v <- t * v
          d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta+v*d.beta) - loglkd
        }
      }
      gamma.prov <- gamma.prov + v * d.gamma.prov
      gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
      beta.new <- beta + v * d.beta
      beta.crit <- norm(matrix(beta-beta.new),"I") # stopping criterion
      beta <- beta.new

      if (message){
        cat(paste0("Iter ",iter,": Inf norm of running diff in est reg parm is ",
                   formatC(beta.crit,digits=3,format="e"),";\n"))
      }
    }
    if (message == TRUE){
      message("\n Algorithm converged after ",iter," iterations!")
    }
  }

  gamma.obs <- rep(gamma.prov, n.prov)
  neg2Loglkd <- -2*sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
  AIC <- neg2Loglkd + 2 * (length(gamma.prov)+length(beta))
  BIC <- neg2Loglkd + log(nrow(data)) * (length(gamma.prov)+length(beta))

  df.prov <- data.frame(Obs_provider = sapply(split(data[,Y.char],data[,prov.char]),sum),
                        gamma_est = gamma.prov) #original gamma-hat, for internal using
  pred <- as.numeric(plogis(gamma.obs+Z%*%beta))

  # modify "outlier provider" to Inf or -Inf
  # if (sum(n.events.prov==n.prov) != 0 | sum(n.events.prov==0) != 0) {
  #   gamma.prov[n.events.prov==n.prov] <- Inf
  #   gamma.prov[n.events.prov==0] <- -Inf
  # }

  #change output format
  beta <- matrix(beta)
  gamma.prov <- matrix(gamma.prov)
  dimnames(beta) <- list(Z.char, "beta")
  dimnames(gamma.prov) <- list(names(n.prov), "gamma")

  char_list <- data.prep$char_list

  return_ls <- structure(list(beta = beta,
                              gamma = gamma.prov, #provider effect
                              pred = pred, #predicted probability
                              obs = data[, Y.char], #patient-level obs
                              neg2Loglkd = neg2Loglkd,
                              AIC = AIC,
                              BIC = BIC),
                         class = "logis_fe")
  if (AUC) {
    AUC <- pROC::auc(data[,Y.char], pred)
    return_ls$AUC <- AUC[1]
  }
  return_ls$df.prov <- df.prov
  return_ls$char_list <- char_list
  return_ls$data_include <- data
  return(return_ls)
}
