#=== FUNCTIONS FOR EXTRACTING COEFFICIENTS =================================
#' Return the model coefficients of a \code{logis_fe} object
#'
#' @param fit an object as output of \code{logis_fe} function
#'
#' @param ...
#'
#' @exportS3Method coef logis_fe
#'
#' @examples
#' data(data_FE)
#' data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID)
#' fit_fe <- logis_fe(data.prep)
#' coef(fit_fe)  #covariate coefficient

coef.logis_fe <- function(fit, ...) {
  if (missing(fit)) stop ("Argument 'fit' is required!", call.=F)
  if (!class(fit) %in% c("logis_fe")) stop("Object fit is not of the classes 'logis_fe'!", call.=F)

  coef <- list(gamma = fit$gamma,
               beta = fit$beta)
  return(coef)

}
