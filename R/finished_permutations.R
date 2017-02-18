### Gehan
#' Permutations of gehan test
#' 
#' Performs permutations of the gehan test through random sub sampling 
#' 
#' @param B Integer specifying the number of permutations to perform
#' @param z A matrix which combines the x and y matricies
#' @param delta.z A matrix which combines the delta.x and delta.y matricies
#' @param n1 Integer. Should be specified as number of rows in the x matrix
#' @param n2 Integer. Shouled be specified as number of rows in the y matrix
#' @param p Integer. No clue what it does
#' @param g.test Numeric from gehan
#' 
#' @return Numeric
#' 
#' @examples
#' See ?MultSurvTests
#' 
#' @export
perm_gehan <- function(B, z, delta.z, n1, n2, p, g.test)
{
  # Progress Bar
  pb = txtProgressBar(min = 0, max = B, initial = 0, 
                      char = "=", width = getOption("width"))
  
  T.perm <- rep(NA, B)
  n <- n1+n2
  
  for(i in 1:B)
  {
    # Update progress bar at the beginning of each loop
    setTxtProgressBar(pb, i)
    
    # Random indicies
    x.rows<-sample(n,n1)
    # Subset the full matrix
    x.new<-z[x.rows,]
    y.new<-z[-x.rows,]
    delta.x.new<-delta.z[x.rows,]
    delta.y.new<-delta.z[-x.rows,]
    # Add each gehan result to vector T.perm
    T.perm[i]<- gehan(x.new,y.new,delta.x.new,delta.y.new,n1,n2,p)
  }
  cat("\n")
  return(sum(T.perm >= as.numeric(g.test))/B)
}

### Mvlogrank
#' Permutations of mvlogrank test
#' 
#' Performs permutations of the mvlogrank test through random sub sampling 
#' 
#' @param B Integer specifying the number of permutations to perform
#' @param z A matrix which combines the x and y matricies
#' @param delta.z A matrix which combines the delta.x and delta.y matricies
#' @param n1 Integer. Should be specified as number of rows in the x matrix
#' @param n2 Integer. Shouled be specified as number of rows in the y matrix
#' @param p Integer. No clue what it does
#' @param g.test Numeric from gehan
#' 
#' @return Numeric
#' 
#' @example 
#' See ?MultSurvTests
#' 
#' @export
perm_mvlogrank <- function(B, z, delta.z, n1, n2, p, mv.test)
{
  # Progress Bar
  pb = txtProgressBar(min = 0, max = B, initial = 0, 
                      char = "=", width = getOption("width"))
  # Vector to hold results
  T.perm <- rep(NA, B)
  n <- n1+n2
  
  for(i in 1:B)
  {
    # Update progress bar at the beginning of each loop
    setTxtProgressBar(pb, i)
    
    # Random indicies
    x.rows<-sample(n,n1)
    # Subset the full matrix by random indicies
    x.new<-z[x.rows,]
    y.new<-z[-x.rows,]
    delta.x.new<-delta.z[x.rows,]
    delta.y.new<-delta.z[-x.rows,]
    # Add each mvlogrank result to vector T.perm[i]
    T.perm[i]<- mvlogrank(x.new,y.new,delta.x.new,delta.y.new,n1,n2,p)
  }
  cat("\n")
  return(sum(T.perm >= as.numeric(mv.test))/B)
}
