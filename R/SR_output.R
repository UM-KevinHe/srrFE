#=== MEASURES OUTPUT FUNCTION =================================
#' Provide direct/indirect standardization ratio or rate
#'
#' @param fit an object as output of \code{logis_fe} function.
#'
#' @param stdz a character string specifying the standardization method. Defaulting to "indirect".
#' \itemize{
#'   \item{"indirect":} using indirect standardized method
#'   \item{"direct":} using direct standardized method
#' }
#'
#' @param measure a boolean indicating whether the output measure is "ratio" or "rate". Both "rate" and "ratio" will be provided by default.
#' \itemize{
#'   \item{"rate":} output the standardized rate. The "rate" has been restricted to 0% - 100%.
#'   \item{"ratio":} output the standardized ratio.
#' }
#'
#' @param null if "stdz = indirect", a character string or real number specifying null hypotheses of fixed provider effects for calculating standardized rate/ratio. Defaulting to "median".
#'
#' @param Rcpp a boolean indicating whether to use Rcpp. Defaulting to TRUE.
#'
#' @param threads an integer specifying the number of threads to use if "Rcpp = T". Defaulting to 1.
#'
#' @param ...
#'
#'
#' @details
#'
#' The `"stdz"` and `"measure"` arguments must be explicitly provided.
#' Users are allowed to specify `"stdz = c("indirect", "direct")"` or `"measure = c("rate", "ratio")"` to get multiple measures.
#'
#'
#' @return The return values depend on the user's choice of standardization method and measure type
#'
#' \item{indirect.ratio}{a vector of standardization ratio using indirect method}
#'
#' \item{direct.ratio}{a vector of standardization ratio using direct method}
#'
#' \item{indirect.rate}{a vector of standardization rate using indirect method}
#'
#' \item{direct.rate}{a vector of standardization rate using direct method}
#'
#'
#'
#' @examples
#' data(data_FE)
#' data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID, message = FALSE)
#' fit_fe <- logis_fe(data.prep)
#' SR <- SR_output(fit_fe, stdz = "direct", measure = "rate")
#' SR$direct.rate
#'
#' @references
#' \itemize{
#' \item He K, Kalbfleisch, J, Li, Y, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
#' \emph{Lifetime Data Analysis}, \strong{19}: 490-512.
#' }
#'
#' @export


SR_output <- function(fit, stdz = "indirect", measure = c("rate", "ratio"), null = "median",
                      Rcpp = TRUE, threads = 1, ...){
  if (missing(fit)) stop ("Argument 'fit' is required!", call.=F)
  if (!class(fit) %in% c("logis_fe")) stop("Object fit is not of the classes 'logis_fe'!", call.=F)

  if (!"indirect" %in% stdz & !"direct" %in% stdz){
    stop("Argument 'stdz' NOT as required!",call.=F)
  }
  if (!"rate" %in% measure & !"ratio" %in% measure){
    stop("Argument 'measure' NOT as required!",call.=F)
  }

  gamma.prov <- fit$gamma
  Z_beta <- fit$linear_pred


  return_ls <- list()
  OE_list <- list()

  if ("indirect" %in% stdz) {
    gamma.null <- ifelse(null=="median", median(gamma.prov),
                         ifelse(class(null)=="numeric", null[1],
                                stop("Argument 'null' NOT as required!",call.=F)))
    Exp <- as.numeric(plogis(gamma.null + Z_beta)) # expected prob of events under null

    df.prov <- data.frame(Obs_provider = fit$df.prov$Obs_provider,
                          Exp.indirect_provider = sapply(split(Exp, fit$prov), sum))
    OE_list$OE_indirect <- df.prov
    df.prov$IS_Ratio <- df.prov$Obs_provider / df.prov$Exp.indirect_provider #indirect standardized ratio: O_i/E_i
    if ("ratio" %in% measure){
      indirect_stdz.ratio <- matrix(df.prov$IS_Ratio)
      dimnames(indirect_stdz.ratio) <- list(rownames(fit$gamma), "Indirect_standardized.ratio")
      return_ls$indirect.ratio <- indirect_stdz.ratio
    }
    if ("rate" %in% measure) {
      population_rate <- sum(fit$obs)/length(Z_beta) * 100  #sum(O_i)/N *100%
      df.prov$IS_Rate <- pmax(pmin(df.prov$IS_Ratio * population_rate, 100), 0)  #restricted to 0%-100%
      indirect_stdz.rate <- matrix(df.prov$IS_Rate)
      dimnames(indirect_stdz.rate) <- list(rownames(fit$gamma), "Indirect_standardized.rate")
      return_ls$indirect.rate <- indirect_stdz.rate
    }
  }

  if ("direct" %in% stdz) {
    if (Rcpp) {
      Exp <- computeDirectExp(gamma.prov, Z_beta, threads)
    } else {
      exp_ZB <- exp(Z_beta)
      Exp.direct <- function(gamma){
        numer <- exp(gamma) * exp_ZB
        sum(1/(1 + 1/numer))
      }
      Exp <- sapply(gamma.prov, Exp.direct)
    }

    df.prov <- data.frame(Obs_all = rep(sum(fit$obs), length(gamma.prov)), #denominator
                          Exp.direct_all = Exp)  #numerator
    OE_list$OE_direct <- df.prov
    df.prov$DS_ratio <- df.prov$Exp.direct_all / df.prov$Obs_all #direct standardized ratio: E/sum(O)
    if ("ratio" %in% measure){
      direct_stdz.ratio <- matrix(df.prov$DS_ratio)
      dimnames(direct_stdz.ratio) <- list(rownames(fit$gamma), "Direct_standardized.ratio")
      return_ls$direct.ratio <- direct_stdz.ratio
    }
    if ("rate" %in% measure) {
      population_rate <- sum(fit$obs)/length(Z_beta) * 100  #sum(O_i)/N *100%
      df.prov$DS_Rate <- pmax(pmin(df.prov$DS_ratio * population_rate, 100), 0)  #restricted to 0%-100%
      direct_stdz.rate <- matrix(df.prov$DS_Rate)
      dimnames(direct_stdz.rate) <- list(rownames(fit$gamma), "Direct_standardized.rate")
      return_ls$direct.rate <- direct_stdz.rate
    }
  }
  return_ls$OE <- OE_list
  return(return_ls)
}








