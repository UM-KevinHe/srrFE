#=== FUNCTIONS FOR FITTING THE Fixed Effects MODEL =================================
#' Main function for fitting fixed effects model
#'
#' @param data.prep an object from the `fe_data_prep()` function.
#'
#' @param algorithm a string specifying the algorithm to be used. Defaulting to "SerBIN".
#'   \itemize{
#'   \item "SerBIN": using the Serial blockwise inversion Newton algorithm to fit the model (See reference: Wu et al.).
#'   \item "BAN": using the block ascent Newton algorithm to fit the model (See reference: He et al.).
#'   }
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
#' @details
#'
#' The default algorithm is based on Serial blockwise inversion Newton (SerBIN) proposed by Wu et al. (2022),
#' but users can also choose to use the block ascent Newton (BAN) algorithm proposed by He et al. (2013) to fit the model.
#' We suggest using the default "SerBIN" option as it typically converges much faster for most datasets.
#' However, in rare cases where the SerBIN algorithm encounters second-order derivative irreversibility leading to an error,
#' users can consider using the "BAN" option as an alternative.
#'
#'
#'
#' @return An object with S3 class \code{logis_fe}.
#'
#' \item{beta}{a vector of fixed effects estimates of covariates}
#'
#' \item{gamma}{a vector of estimates of facility effects}
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
#' data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID, message = FALSE)
#' fit_fe <- logis_fe(data.prep, message = TRUE)
#'
#' @importFrom Rcpp evalCpp
#' @importFrom pROC auc
#'
#' @references
#' \itemize{
#' \item Wu, W, Yang, Y, Kang, J, He, K. (2022) Improving large-scale estimation and inference for profiling health care providers.
#' \emph{Statistics in Medicine}, \strong{41(15)}: 2840-2853.
#'
#' \item He K, Kalbfleisch, J, Li, Y, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
#' \emph{Lifetime Data Analysis}, \strong{19}: 490-512.
#' }
#'
#' @keywords Cost-Efficient Newton-Raphson (CENR) Algorithm, Fixed Provider Effects
#'
#' @export
#'


logis_fe <- function(data.prep, algorithm = "SerBIN", max.iter = 10000, tol = 1e-5, bound = 10,
                     backtrack = FALSE, Rcpp = TRUE, AUC = FALSE, message = FALSE){ #BAN #SerBIN
  if (missing(data.prep)) stop ("Argument 'data.prep' is required!", call.=F)
  if (!class(data.prep) %in% c("data_prep")) stop("Object 'data.prep' should be generated from 'fe_data_prep' function!", call.=F)

  if (!is.logical(backtrack)) stop("Argument 'backtrack' NOT as required!", call.=F)

  data <- data.prep$data
  Y.char <- data.prep$char_list$Y.char
  prov.char <- data.prep$char_list$prov.char
  Z.char <- data.prep$char_list$Z.char

  data <- data[data$included==1,]
  n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length) # provider-specific number of discharges
  n.readm.prov <- sapply(split(data[, Y.char], data[, prov.char]), sum) # provider-specific number of readmissions
  Z <- as.matrix(data[,Z.char])
  gamma.prov <- rep(log(mean(data[,Y.char])/(1-mean(data[,Y.char]))), length(n.prov))
  beta <- rep(0, NCOL(Z))


  if (algorithm == "SerBIN") {
    if (Rcpp) { #Rcpp always use "backtrack"
      ls <- logis_BIN_fe_prov(as.matrix(data[,Y.char]),Z,n.prov,gamma.prov,beta,0,1,tol,max.iter, bound, message)
      gamma.prov <- as.numeric(ls$gamma)
      beta <- as.numeric(ls$beta)
    } else {
      iter <- 0
      beta.crit <- 100 # initialize stop criterion
      if (message == TRUE){
        message("Implementing SerBIN algorithm for fixed provider effects model ...")
      }

      if (backtrack) {
        s <- 0.01; t <- 0.6 # initialize parameters for backtracking line search
        Loglkd <- function(gamma.obs, beta) {
          sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
        }
        loglkd <- Loglkd(rep(gamma.prov, n.prov), beta)
        while (iter<=max.iter & beta.crit>=tol) {
          iter <- iter + 1
          gamma.obs <- rep(gamma.prov, n.prov)
          p <- c(plogis(gamma.obs+Z%*%beta))
          q <- p*(1-p)
          score.gamma <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum)
          score.beta <- t(Z)%*%(data[,Y.char]-p)
          info.gamma.inv <- 1/sapply(split(q, data[,prov.char]),sum) #I_11^(-1)
          info.betagamma <- sapply(by(q*Z,data[,prov.char],identity),colSums) #I_21
          info.beta <- t(Z)%*%(q*Z) #I_22
          mat.tmp1 <- info.gamma.inv*t(info.betagamma) #J_1^T
          schur.inv <- solve(info.beta-info.betagamma%*%mat.tmp1) #S^-1
          mat.tmp2 <- mat.tmp1%*%schur.inv #J_2^T
          d.gamma.prov <- info.gamma.inv*score.gamma +
            mat.tmp2%*%(t(mat.tmp1)%*%score.gamma-score.beta)
          d.beta <- -t(mat.tmp2)%*%score.gamma+schur.inv%*%score.beta
          v <- 1 # initialize step size
          d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta+v*d.beta) - loglkd
          lambda <- c(score.gamma,score.beta)%*%c(d.gamma.prov,d.beta)
          while (d.loglkd < s*v*lambda) {
            v <- t * v
            d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta+v*d.beta) - loglkd
          }
          gamma.prov <- gamma.prov + v * d.gamma.prov
          gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
          beta.new <- beta + v * d.beta
          beta.crit <- norm(matrix(beta-beta.new),"I") # stopping criterion
          beta <- beta.new
          loglkd <- loglkd + d.loglkd
          if (message == TRUE){
            cat(paste0("Iter ",iter,": Inf norm of running diff in est reg parm is ",
                       formatC(beta.crit,digits=3,format="e"),";\n"))
          }
        }
      } else {
        while (iter<=max.iter & beta.crit>=tol) {
          iter <- iter + 1
          if (message == TRUE){
            cat(paste0("\n Iter ",iter,":"))
          }
          gamma.obs <- rep(gamma.prov, n.prov)
          p <- c(plogis(gamma.obs+Z%*%beta))
          score.gamma <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum)
          score.beta <- t(Z)%*%(data[,Y.char]-p)
          info.gamma.inv <- 1/sapply(split(p*(1-p), data[,prov.char]),sum)
          info.betagamma <- sapply(by(p*(1-p)*Z,data[,prov.char],identity),colSums)
          info.beta <- t(Z)%*%(p*(1-p)*Z)
          schur.inv <- solve(info.beta-info.betagamma%*%(info.gamma.inv*t(info.betagamma)))
          mat.tmp1 <- info.gamma.inv*t(info.betagamma)
          mat.tmp2 <- mat.tmp1%*%schur.inv
          gamma.prov <- gamma.prov +
            info.gamma.inv*score.gamma+mat.tmp2%*%(t(mat.tmp1)%*%score.gamma-score.beta)
          gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
          beta.new <- beta - t(mat.tmp2)%*%score.gamma+schur.inv%*%score.beta
          beta.crit <- norm(matrix(beta-beta.new),"I") # stopping criterion
          beta <- beta.new
          if (message == TRUE){
            cat(paste0("Inf norm of running diff in est reg parm is ",formatC(beta.crit,digits=3,format="e"),";"))
          }
        }
      }
      if (message == TRUE){
        message("\n SerBIN algorithm converged after ",iter," iterations!")
      }
    }
  } else if (algorithm == "BAN"){
    if (Rcpp) {
      ls <- logis_fe_prov(as.matrix(data[,Y.char]),Z,n.prov,gamma.prov,beta,backtrack,max.iter,bound,tol)
      gamma.prov <- as.numeric(ls$gamma); beta <- as.numeric(ls$beta)
    } else {
      iter <- 0
      beta.crit <- 100 # initialize stop criterion
      if (message == TRUE){
        message("Implementing BAN algorithm for fixed provider effects model ...")
      }

      if (backtrack) {
        s <- 0.01; t <- 0.8 # initialize parameters for backtracking line search
        Loglkd <- function(gamma.obs, beta) {
          sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
        }
        while (iter<=max.iter & beta.crit>=tol) {
          iter <- iter + 1
          # provider effect update
          gamma.obs <- rep(gamma.prov, n.prov)
          Z.beta <- Z%*%beta
          p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
          score.gamma.prov <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum)
          d.gamma.prov <- score.gamma.prov / sapply(split(pq, data[,prov.char]), sum)
          v <- 1 # initialize step size
          loglkd <- Loglkd(rep(gamma.prov, n.prov), beta)
          d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta) - loglkd
          lambda <- score.gamma.prov%*%d.gamma.prov
          while (d.loglkd < s*v*lambda) {
            v <- t * v
            d.loglkd <- Loglkd(rep(gamma.prov+v*d.gamma.prov, n.prov), beta) - loglkd
          }
          gamma.prov <- gamma.prov + v * d.gamma.prov
          gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
          gamma.obs <- rep(gamma.prov, n.prov)

          # regression parameter update
          p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
          score.beta <- t(Z)%*%(data[,Y.char]-p)
          info.beta <- t(Z)%*%(c(pq)*Z)
          d.beta <- as.numeric(solve(info.beta)%*%score.beta)
          v <- 1 # initialize step size
          loglkd <- Loglkd(gamma.obs, beta)
          d.loglkd <- Loglkd(gamma.obs, beta+v*d.beta) - loglkd
          lambda <- t(score.beta)%*%d.beta
          while (d.loglkd < s*v*lambda) {
            v <- t * v
            d.loglkd <- Loglkd(gamma.obs, beta+v*d.beta) - loglkd
          }
          beta.new <- beta + v * d.beta
          beta.crit <- norm(matrix(beta-beta.new),"I") # stopping criterion
          beta <- beta.new
          if (message == TRUE){
            cat(paste0("Iter ",iter,": Inf norm of running diff in est reg parm is ",
                       formatC(beta.crit,digits=3,format="e"),";\n"))
          }
        }
      } else {
        while (iter<=max.iter & beta.crit>=tol) {
          iter <- iter + 1
          if (message == TRUE){
            cat(paste0("\n Iter ",iter,":"))
          }
          # provider effect update
          gamma.obs <- rep(gamma.prov, n.prov)
          Z.beta <- Z%*%beta
          p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
          gamma.prov <- sapply(split(data[,Y.char]-p, data[,prov.char]), sum) /
            sapply(split(pq, data[,prov.char]), sum) + gamma.prov
          gamma.prov <- pmin(pmax(gamma.prov, median(gamma.prov)-bound), median(gamma.prov)+bound)
          gamma.obs <- rep(gamma.prov, n.prov)

          # regression parameter update
          p <- c(plogis(gamma.obs+Z.beta)); pq <- p*(1-p)
          score.beta <- t(Z)%*%(data[,Y.char]-p)
          info.beta <- t(Z)%*%(c(pq)*Z)
          beta.new <- beta + as.numeric(solve(info.beta)%*%score.beta)
          beta.crit <- norm(matrix(beta-beta.new),"I") # stopping criterion
          beta <- beta.new
          if (message == TRUE){
            cat(paste0("Inf norm of running diff in est reg parm is ",formatC(beta.crit,digits=3,format="e"),";"))
          }
        }

      }
      if (message == TRUE){
        message("\n BAN algorithm converged after ",iter," iterations!")
      }
    }
  } else {
    stop("Argument 'algorithm' NOT as required!")
  }

  gamma.obs <- rep(gamma.prov, n.prov)
  neg2Loglkd <- -2*sum((gamma.obs+Z%*%beta)*data[,Y.char]-log(1+exp(gamma.obs+Z%*%beta)))
  AIC <- neg2Loglkd + 2 * (length(gamma.prov)+length(beta))
  BIC <- neg2Loglkd + log(nrow(data)) * (length(gamma.prov)+length(beta))

  df.prov <- data.frame(Obs_facility = sapply(split(data[,Y.char],data[,prov.char]),sum),
                        gamma_est = gamma.prov) #original gamma-hat, for internal using

  # modify "outlier facility"
  if (sum(n.readm.prov==n.prov) != 0 | sum(n.readm.prov==0) != 0) {
    gamma.prov[n.readm.prov==n.prov] <- Inf
    gamma.prov[n.readm.prov==0] <- -Inf
  }

  #change output format
  beta <- matrix(beta)
  gamma.prov <- matrix(gamma.prov)
  dimnames(beta) <- list(Z.char, "beta")
  dimnames(gamma.prov) <- list(names(n.prov), "gamma")

  char_list <- data.prep$char_list

  return_ls <- structure(list(beta = beta,
                              gamma = gamma.prov, #facility effect
                              obs = data[, Y.char], #patient-level obs
                              neg2Loglkd = neg2Loglkd,
                              AIC = AIC,
                              BIC = BIC),
                         class = "logis_fe")
  if (AUC) {
    Pred <- as.numeric(plogis(gamma.obs+Z%*%beta))
    AUC <- pROC::auc(data[,Y.char], Pred)
    return_ls$AUC <- AUC[1]
  }
  return_ls$df.prov <- df.prov
  return_ls$char_list <- char_list
  return_ls$data_include <- data
  return(return_ls)
}






