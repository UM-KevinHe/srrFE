#=== FUNCTIONS FOR DATA PREPARATION FOR FIXED EFFECTS MODEL =================================
#' Data preparation before modelling
#'
#' @param Y a numerical vector, with values of 0 or 1, indicating the outcome variable.
#'
#' @param Z a matrix or data frame containing covariates.
#'
#' @param ID a vector representing the provider id. Its elements can be either numeric values or characters.
#'
#' @param cutoff an integer as cutoff of provider size with 10 as default. Providers with observations fewer than the "cutoff" value will be labeled as "include = 0".
#'
#' @param check a Boolean indicating whether checking missingness, variation, correlation and VIF of variables in the data. Defaulting to "TRUE".
#'
#' @param message a Boolean indicating whether printing out the information of the data preparation process. Defaulting to "TRUE".
#'
#' @param ...
#'
#'
#' @details
#' Major steps in this stage include (in order):
#'   \itemize{
#'   \item checking missingness, variation, correlation and VIF of variables in the data,
#'   \item provider screening based on the provider size, 10 as default,
#'   \item reporting proportions of providers with no events (i.e. all "0" outcomes) and with all events (i.e. all "1" outcomes),
#'   \item sorting data by provider identifiers.
#'   }
#'
#' The `fe_data_prep()` function returns data sorted by the provider identifiers,
#' accompanied by additional provider-related information indicating whether the provider's size exceeds a specified "cutoff"
#' and whether the respective provider has experienced either zero or all events.
#' The reason behind introducing a "cutoff" lies in findings from both simulated and real data studies,
#' revealing the instability of coefficient estimates for providers with small sizes.
#' Consequently, we strongly recommend excluding small providers during the model fitting process.
#' It is important to note that the resultant data frame retains all providers, including small ones,
#' but utilizes the "included = 0" label to signify the small providers. Subsequently, during the model fitting stage,
#' the `logis_fe()` function disregards records marked with "included = 0".
#'
#'
#' @return
#'
#' \item{data}{a sorted data frame including response, provider identifiers, covariates, and additional provider information.}
#'
#' \item{char_list}{a list including variable names.}
#'
#'
#' @examples
#' data(data_FE)
#' data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID)
#' head(data.prep$data)
#' data.prep$char_list
#'
#' @importFrom caret nearZeroVar
#' @importFrom olsrr ols_vif_tol
#'
#' @keywords Data Preparation, Fixed Provider Effects
#'
#' @export

fe_data_prep <- function(Y, Z, ID, cutoff = 10, check = TRUE, message = TRUE) {
  data <- as.data.frame(cbind(Y, ID, Z))
  Y.char <- colnames(data)[1]
  prov.char <- colnames(data)[2]
  Z.char <- colnames(Z)

  #check dimensions of the input data
  if (length(Y) != length(ID) | length(ID) != nrow(Z)){
    stop("Dimensions of the input data do not match!!", call.=F)
  }

  if (check) {
    ## check missingness of variables
    if (message == TRUE) message("Checking missingness of variables ... ")
    if (sum(complete.cases(data[,c(Y.char,Z.char,prov.char)]))==NROW(data)) {
      if (message == TRUE) message("Missing values NOT found. Checking missingness of variables completed!")
    } else {
      check.na <- function(name) {
        if (sum(is.na(data[,name])) > 0) {
          warning(sum(is.na(data[,name]))," out of ",NROW(data[,name])," in '",name,"' missing!",immediate.=T,call.=F)
        }
      }
      invisible(sapply(c(Y.char,Z.char,prov.char), check.na))
      missingness <- (1 - sum(complete.cases(data[,c(Y.char,Z.char,prov.char)])) / NROW(data)) * 100
      stop(paste(round(missingness,2), "% of all observations are missing!",sep=""),call.=F)
    }
    ## check variation in covariates
    if (message == TRUE) message("Checking variation in covariates ... ")
    nzv <- caret::nearZeroVar(data[,Z.char], saveMetrics=T)
    if (sum(nzv$zeroVar==T) > 0) {
      stop("Covariate(s) '", paste(row.names(nzv[nzv$zeroVar==T,]), collapse="', '"),
           "' with zero variance(s)!", call.=F)
    } else if (sum(nzv$nzv==T) > 0) {
      warning("Covariate(s) '",paste(row.names(nzv[nzv$nzv==T,]), collapse="', '"),
              "' with near zero variance(s)!",immediate.=T,call.=F)
    }
    if (message == TRUE) message("Checking variation in covariates completed!")
    ## check correlation
    if (message == TRUE) message("Checking pairwise correlation among covariates ... ")
    cor <- cor(data[,Z.char])
    threshold.cor <- 0.9
    if (sum(abs(cor[upper.tri(cor)])>threshold.cor) > 0) {
      cor[lower.tri(cor,diag=T)] <- 0
      ind <- which(abs(cor)>threshold.cor)
      pairs <- sapply(ind, function(ind) c(rownames(cor)[ind%%NROW(cor)],
                                           colnames(cor)[ind%/%NROW(cor)+1]))
      warning("The following ", NCOL(pairs),
              " pair(s) of covariates are highly correlated (correlation > ",
              threshold.cor,"): ", immediate.=T, call.=F)
      invisible(apply(pairs,2,function(col) message('("',paste(col, collapse='", "'),'")')))
    }
    if (message == TRUE) message("Checking pairwise correlation among covariates completed!")
    ## check VIF
    if (message == TRUE) message("Checking VIF of covariates ... ")
    m.lm <- lm(as.formula(paste(Y.char,"~",paste(Z.char, collapse="+"))), data=data)
    vif <- olsrr::ols_vif_tol(m.lm)
    if(sum(vif$VIF >= 10) > 0){
      warning("Covariate(s) '",
              paste(as.data.frame(vif)[vif$VIF>=10,"Variables"], collapse="', '"),
              "' with serious multicollinearity!",immediate.=T,call.=F)
    }
    if (message == TRUE) message("Checking VIF of covariates completed!")
  }

  data <- data[order(factor(data[,prov.char])),] # sort data by provider ID
  prov.size <- as.integer(table(data[,prov.char])) # provider sizes
  prov.size.long <- rep(prov.size,prov.size) # provider sizes assigned to patients
  data$included <- 1 * (prov.size.long > cutoff) # create variable 'included' as an indicator
  if (message == TRUE) warning(sum(prov.size<=cutoff)," out of ",length(prov.size),
          " providers considered small and filtered out!",immediate.=T,call.=F)


  prov.list <- unique(data[data$included==1,prov.char])   # a reduced list of provider IDs
  prov.no.events <-      # providers with no events
    prov.list[sapply(split(data[data$included==1,Y.char], factor(data[data$included==1,prov.char])),sum)==0]
  data$no.events <- 0
  data$no.events[data[,prov.char]%in%c(prov.no.events)] <- 1
  if (message == TRUE) message(paste(length(prov.no.events),"out of",length(prov.list),
                "remaining providers with no events."))
  prov.all.events <-     # providers with all events
    prov.list[sapply(split(1-data[data$included==1,Y.char],factor(data[data$included==1,prov.char])),sum)==0]
  data$all.events <- 0
  data$all.events[data[,prov.char]%in%c(prov.all.events)] <- 1
  if (message == TRUE) message(paste(length(prov.all.events),"out of",length(prov.list),
                "remaining providers with all events."))
  if (message == TRUE) message(paste0("After screening, ", round(sum(data[data$included==1,Y.char])/length(data[data$included==1,Y.char])*100,2),
                 "% of all records exhibit occurrences of events (Y = 1)"))


  char_list <- list(Y.char = Y.char,
                    prov.char = prov.char,
                    Z.char = Z.char)

  return_ls <- list(data = data,
                    char_list = char_list)


  return_ls <- structure(list(data = data,
                              char_list = char_list),
                         class = "data_prep")
  return(return_ls)
}









