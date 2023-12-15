#' Example dataset for fixed effects model
#'
#' A simulated data set containing response variable, provider information and 4 continuous covariates.
#' @name data_FE
#' @docType data
#' @usage data(data_FE)
#' @keywords datasets
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{Y}{a vector represents the binary outcome}
#'   \item{ID}{a vector represents the facility indicator (10 facilities in total)}
#'   \item{Z}{a data frame contains 4 continuous variables}
#' }
#' 
#' @examples
#' data(data_FE)
#' head(data_FE$Y) 
#' head(data_FE$ID)
#' head(data_FE$Z)
#' 
"data_FE"