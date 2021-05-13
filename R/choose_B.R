#' Choose B for permutation tests
#' 
#' Computes the value of B for a permutation test required to obtain a specified accuracy when
#' approximating the permutation p-values using B random permutations. 
#' 
#' @param p0 A guess for the p-value. Can be based e.g. on a small number of simulations. The default is 0.05.
#' @param width The desired width of the Clopper-Pearson interval. The default is 0.01.
#' @param conf.level The confidence level of the Clopper-Pearson interval. The default is 0.95.
#' 
#' @return B
#' 
#' @details Computations are based on the Clopper-Pearson interval, using a formula from
#' Thulin (2014). The procedure is described in Section 3.3 in Persson et al. (2019).
#' 
#' @examples 
#' 
#' 
#' @export
choose_B <- function(p0 = 0.05, width = 0.01, conf.level = 0.95)
{
  alpha <- 1-conf.level
  ceiling((2*qnorm(1-alpha/2)^2*p0*(1-p0)+2*qnorm(1-alpha/2)*sqrt(qnorm(1-alpha/2)^2*p0^2*(1-p0)^2+width*p0*(1-p0))+width)/width^2)
}