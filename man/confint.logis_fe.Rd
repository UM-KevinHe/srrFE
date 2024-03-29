% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint.logis_fe.R
\name{confint.logis_fe}
\alias{confint.logis_fe}
\title{Provide confidence interval for provider effects or standardization ratios/rates}
\usage{
\method{confint}{logis_fe}(
  fit,
  parm,
  level = 0.95,
  test = "exact",
  option = c("gamma", "SR"),
  stdz = "indirect",
  measure = c("rate", "ratio")
)
}
\arguments{
\item{fit}{an object as output of \code{logis_fe} function.}

\item{parm}{specify a subset of which providers are to be given confidence intervals. All providers are included by default.}

\item{level}{confidence level used for constructing confidence intervals. Defaulting to 0.95.}

\item{test}{a character string specifying the type of testing method. Defaulting to "exact".
\itemize{
\item "exact": exact test
\item "wald": wald test
\item "score": score test
}}

\item{option}{the confidence interval for the function's output, whether it is for gamma or standardization ratios/rates.
\itemize{
\item "gamma": provider effect
\item "SR": standardization ratios/rates
}}

\item{stdz}{if option = 'SR', a character string specifying the standardization method. Defaulting to "indirect".
\itemize{
\item "indirect": using indirect standardized method
\item "direct": using direct standardized method
}}

\item{measure}{if option = 'SR', a boolean indicating whether the output measure is "ratio" or "rate". Both "rate" and "ratio" will be provided by default.
\itemize{
\item "rate": output the standardized rate. The "rate" has been restricted to 0\% - 100\%.
\item "ratio":  output the standardized ratio
}}

\item{...}{}
}
\value{
A dataframe containing the point estimate, and lower and upper bounds of the estimate.
}
\description{
Provide confidence interval for provider effects or standardization ratios/rates
}
\examples{
data(data_FE)
data.prep <- fe_data_prep(data_FE$Y, data_FE$Z, data_FE$ID, message = FALSE)
fit_fe <- logis_fe(data.prep)
confint(fit_fe, option = "gamma")
confint(fit_fe, option = "SR")

}
