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
#' SR <- SR_output(fit_fe, stdz = "indirect", measure = "ratio")
#' SR$indirect.ratio
#'
#' @references
#' \itemize{
#' \item He K, Kalbfleisch, J, Li, Y, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
#' \emph{Lifetime Data Analysis}, \strong{19}: 490-512.
#' }
#'
#' @export


SR_output <- function(fit, stdz = "indirect", measure = c("rate", "ratio"), null = "median"){
  if (missing(fit)) stop ("Argument 'fit' is required!", call.=F)
  if (!class(fit) %in% c("logis_fe")) stop("Object fit is not of the classes 'logis_fe'!", call.=F)



  if (!"indirect" %in% stdz & !"direct" %in% stdz){
    stop("Argument 'stdz' NOT as required!",call.=F)
  }
  if (!"rate" %in% measure & !"ratio" %in% measure){
    stop("Argument 'measure' NOT as required!",call.=F)
  }
  Y.char <- fit$char_list$Y.char
  Z.char <- fit$char_list$Z.char
  prov.char <- fit$char_list$prov.char
  data <- fit$data_include

  n.prov <- sapply(split(data[, Y.char], data[, prov.char]), length) # provider-specific number of discharges
  n.events.prov <- sapply(split(data[, Y.char], data[, prov.char]), sum) # provider-specific number of events
  Z <- as.matrix(data[,Z.char])
  gamma.prov <- fit$gamma
  beta <- fit$beta
  gamma.obs <- rep(gamma.prov, n.prov)
  population_rate <- sum(data[,Y.char])/length(gamma.obs) * 100  #sum(O_i)/N *100%

  if ("indirect" %in% stdz) {
    gamma.null <- ifelse(null=="median", median(gamma.prov),
                         ifelse(class(null)=="numeric", null[1],
                                stop("Argument 'null' NOT as required!",call.=F)))
    Exp <- as.numeric(plogis(gamma.null+Z%*%beta)) # expected prob of events under null

    df.prov <- data.frame(Obs_provider = sapply(split(data[,Y.char],data[,prov.char]),sum),
                          Exp.indirect_provider = sapply(split(Exp,data[,prov.char]),sum))
    OE_indirect <- df.prov

    df.prov$IS_Ratio <- df.prov$Obs_provider / df.prov$Exp.indirect_provider #indirect standardized ratio: O_i/E_i
    df.prov$IS_Rate <- pmax(pmin(df.prov$IS_Ratio * population_rate, 100), 0)  #restricted to 0%-100%

    indirect_stdz.ratio <- matrix(df.prov$IS_Ratio) #output measure
    dimnames(indirect_stdz.ratio) <- list(names(n.prov), "Indirect_standardized.ratio")
    indirect_stdz.rate <- matrix(df.prov$IS_Rate)
    dimnames(indirect_stdz.rate) <- list(names(n.prov), "Indirect_standardized.rate")
  }

  if ("direct" %in% stdz) {
    Exp.direct <- function(gamma, Z_beta) {
      sum(as.numeric(plogis(gamma + Z_beta)))
    }
    Z_beta <- Z %*% beta #common across providers
    Exp <- sapply(gamma.prov, Exp.direct, Z_beta = Z_beta)  #numerator
    df.prov <- data.frame(Obs_all = rep(sum(data[,Y.char]), length(gamma.prov)), #denominator
                          Exp.direct_all = Exp)  #numerator
    OE_direct <- df.prov
    df.prov$DS_ratio <- df.prov$Exp.direct_all / df.prov$Obs_all #direct standardized ratio: E/sum(O)
    df.prov$DS_Rate <- pmax(pmin(df.prov$DS_ratio * population_rate, 100), 0)  #restricted to 0%-100%

    direct_stdz.ratio <- matrix(df.prov$DS_ratio)
    dimnames(direct_stdz.ratio) <- list(names(n.prov), "Direct_standardized.ratio")
    direct_stdz.rate <- matrix(df.prov$DS_Rate)
    dimnames(direct_stdz.rate) <- list(names(n.prov), "Direct_standardized.rate")
  }


  # # modify "outlier provider"
  # if (sum(n.events.prov==n.prov) != 0 | sum(n.events.prov==0) != 0) {
  #   if ("direct" %in% stdz) { #only "direct standardized rates/ratios" are affected
  #     if ("ratio" %in% measure){
  #       direct_stdz.ratio[n.events.prov==n.prov,] <- length(gamma.obs) / df.prov$Obs_all[1]
  #       direct_stdz.ratio[n.events.prov==0,] <- 0
  #     }
  #     if ("rate" %in% measure){
  #       direct_stdz.rate[n.events.prov==n.prov,] <- 1
  #       direct_stdz.rate[n.events.prov==0,] <- 0
  #     }
  #   }
  # }

  return_ls <- list()
  OE_list <- list()

  if ("indirect" %in% stdz){
    if ("ratio" %in% measure){
      return_ls$indirect.ratio <- indirect_stdz.ratio
    }
    if ("rate" %in% measure) {
      return_ls$indirect.rate <- indirect_stdz.rate
    }
    OE_list$OE_indirect <- OE_indirect
  }

  if ("direct" %in% stdz) {
    if ("ratio" %in% measure){
      return_ls$direct.ratio <- direct_stdz.ratio
    }
    if ("rate" %in% measure) {
      return_ls$direct.rate <- direct_stdz.rate
    }
    OE_list$OE_direct <- OE_direct
  }
  return_ls$OE <- OE_list
  return(return_ls)
}








