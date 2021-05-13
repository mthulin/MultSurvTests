### Gehan
#' Multivariate permutation Gehan test
#' 
#' Computes the p-value of the multivariate permutation Gehan test described in Persson et al. (2019).
#' 
#' @param B An integer specifying the number of permutations to perform. The default is 999. It is recommended to use \code{\link{choose_B}} for choosing \code{B}.
#' @param z A matrix containing the observed (possibly censored) survival times for the two groups. The observations for the first group should be one the first \code{n1} rows.
#' @param delta.z A matrix containing the censoring status of each observation in \code{z}.
#' @param n1 An integer specifying the sample size of the first group.
#' 
#' @return A p-value.
#' 
#' @details Multivariate version of the logrank and Gehan tests were described by
#' Wei & Lachin (1984). Persson et al. (2019) described permutation versions of 
#' these tests, with improved performance.
#' 
#' @examples
#' # Diabetes data:
#' ?diabetes
#' # Survival times for the two groups:
#' x <- as.matrix(subset(diabetes, LASER==1)[c(6,8)])
#' y <- as.matrix(subset(diabetes, LASER==2)[c(6,8)])
#' # Censoring status for the two groups:
#' delta.x <- as.matrix(subset(diabetes, LASER==1)[c(7,9)])
#' delta.y <- as.matrix(subset(diabetes, LASER==2)[c(7,9)])
#' # Create the input for the test:
#' z <- rbind(x, y)
#' delta.z <- rbind(delta.x, delta.y)
#' # Run the test with 99 permutations:
#' perm_gehan(B = 99, z, delta.z, n1 = nrow(x))
#' 
#' @export
perm_gehan <- function(B = 999, z, delta.z, n1, g.test)
{
  p <- ncol(z)
  n2 <- nrow(z)-n1
  n <- nrow(z)
  
  # Compute the test statistic for the original data:
  x.rows <- 1:n1
  x.new <- z[x.rows,]
  y.new <- z[-x.rows,]
  delta.x.new <- delta.z[x.rows,]
  delta.y.new <- delta.z[-x.rows,]
  g.test <- gehan(x.new, y.new, delta.x.new, delta.y.new, n1, n2, p)
  
  # Progress Bar
  pb = txtProgressBar(min = 0, max = B, initial = 0, 
                      char = "=", width = getOption("width"))
  T.perm <- rep(NA, B)

  for(i in 1:B)
  {
    # Update progress bar at the beginning of each loop:
    setTxtProgressBar(pb, i)
    
    # Random indices:
    x.rows <- sample(n,n1)
    # Subset the full matrix:
    x.new <- z[x.rows,]
    y.new <- z[-x.rows,]
    delta.x.new <- delta.z[x.rows,]
    delta.y.new <- delta.z[-x.rows,]
    # Compute the test statistic for the permuted sample:
    T.perm[i]<- gehan(x.new, y.new, delta.x.new, delta.y.new, n1, n2, p)
  }
  cat("\n")
  return(sum(T.perm >= as.numeric(g.test))/B)
}

### Mvlogrank
#' Multivariate permutation logrank test
#' 
#' Computes the p-value of the multivariate permutation logrank test described in Persson et al. (2019).
#' 
#' @param B An integer specifying the number of permutations to perform. The default is 999. It is recommended to use \code{\link{choose_B}} for choosing \code{B}.
#' @param z A matrix containing the observed (possibly censored) survival times for the two groups. The observations for the first group should be one the first \code{n1} rows.
#' @param delta.z A matrix containing the censoring status of each observation in \code{z}.
#' @param n1 An integer specifying the sample size of the first group.
#' 
#' @return A p-value.
#' 
#' @details Multivariate version of the logrank and Gehan tests were described by
#' Wei & Lachin (1984). Persson et al. (2019) described permutation versions of 
#' these tests, with improved performance.
#' 
#' @examples
#' # Diabetes data:
#' ?diabetes
#' # Survival times for the two groups:
#' x <- as.matrix(subset(diabetes, LASER==1)[c(6,8)])
#' y <- as.matrix(subset(diabetes, LASER==2)[c(6,8)])
#' # Censoring status for the two groups:
#' delta.x <- as.matrix(subset(diabetes, LASER==1)[c(7,9)])
#' delta.y <- as.matrix(subset(diabetes, LASER==2)[c(7,9)])
#' # Create the input for the test:
#' z <- rbind(x, y)
#' delta.z <- rbind(delta.x, delta.y)
#' # Run the test with 99 permutations:
#' perm_mvlogrank(B = 99, z, delta.z, n1 = nrow(x))
#' 
#' @export
perm_mvlogrank <- function(B, z, delta.z, n1)
{
  p <- ncol(z)
  n2 <- nrow(z)-n1
  n <- nrow(z)
  
  # Compute the test statistic for the original data:
  x.rows <- 1:n1
  x.new <- z[x.rows,]
  y.new <- z[-x.rows,]
  delta.x.new <- delta.z[x.rows,]
  delta.y.new <- delta.z[-x.rows,]
  mv.test <- mvlogrank(x.new, y.new, delta.x.new, delta.y.new, n1, n2, p)
  
  # Progress Bar
  pb = txtProgressBar(min = 0, max = B, initial = 0, 
                      char = "=", width = getOption("width"))
  T.perm <- rep(NA, B)
  
  for(i in 1:B)
  {
    # Update progress bar at the beginning of each loop:
    setTxtProgressBar(pb, i)
    
    # Random indices:
    x.rows <- sample(n,n1)
    # Subset the full matrix:
    x.new <- z[x.rows,]
    y.new <- z[-x.rows,]
    delta.x.new <- delta.z[x.rows,]
    delta.y.new <- delta.z[-x.rows,]
    # Compute the test statistic for the permuted sample:
    T.perm[i]<- mvlogrank(x.new,y.new,delta.x.new,delta.y.new,n1,n2,p)
  }
  cat("\n")
  return(sum(T.perm >= as.numeric(mv.test))/B)
}
