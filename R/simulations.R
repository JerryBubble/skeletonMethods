#' Generating Ring data
#'
#' @param n_c a number indicating the number of points in the middle circle
#' @param n_r a number indicating the number of points in the outer ring
#' @param sd_c a number indicating the standard deviation of shifts for points on the middle circle
#' @param sd_r a number indicating the standard deviation of shifts for points on the outer ring
#' @param d a number indicating the dimension of the generated data. d>2 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = ring_data()$data
#' plot(X[,1], X[,2])
ring_data = function(n_c = 200, n_r=1000, sd_c = 0.1, sd_r = 0.1, d=2, sd_high = 0.1){
  c_lab = c(rep(1,n_r), rep(2,n_c))

  # the ring
  th = stats::runif(n_r,0,2*pi)
  x  = cos(th)
  y  = sin(th)
  X  = cbind(x,y) + matrix(stats::rnorm(2*length(x), sd=sd_r), ncol=2)

  # the node
  u = matrix(stats::rnorm(2*n_c, sd=sd_c), ncol=2)
  X = rbind(X, u)

  if(d==2){
    X0 = X
  }
  if(d>2){
    n = nrow(X)
    d_add = d-2
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  Y0 = c(rep(1,n_r), rep(2,n_c))
  return(list(data = X0, clus = Y0))
}

#' Generating Mickey data
#'
#' @param n1 a number indicating the number of points in the larger circle
#' @param n2 a number indicating the number of points in each of the two smaller circles
#' @param s1 a number indicating the base radius for the larger circle
#' @param s2 a number indicating the base radius for each of the two smaller circles
#' @param d a number indicating the dimension of the generated data. d>2 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = Mickey_data()$data
#' plot(X[,1], X[,2])
Mickey_data = function(n1=1000,n2=100,s1=1,s2=0.2, d=2, sd_high=0.1){
  n  = n1
  th = stats::runif(n,0,2*pi)
  r =  s1*stats::rbeta(n,2,1)
  x1 = r*cos(th)
  x2 = r*sin(th)

  n  = n2
  th = stats::runif(n,0,2*pi)
  r =  s2*stats::rbeta(n,2,1)
  a1 = r*cos(th) - 1.2
  a2 = r*sin(th) + 1.2

  n  = n2
  th = stats::runif(n,0,2*pi)
  r =  s2*(stats::rbeta(n,2,1))
  b1 = r*cos(th) + 1.2
  b2 = r*sin(th) + 1.2

  x1 = c(x1,a1,b1)
  x2 = c(x2,a2,b2)
  X = cbind(x1,x2)
  if(d==2){
    X0 = X
  }
  if(d>2){
    n = nrow(X)
    d_add = d-2
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  Y0 = c(rep(1,n1), rep(2,n2), rep(3,n2))
  return(list(data = X0, clus = Y0))
}


#' Generating MixMickey data
#'
#' @param n1 a number indicating the number of points in the larger Gaussian cluster
#' @param n2 a number indicating the number of points in each of the two smaller Gaussian clusters
#' @param s1 a number indicating the base radius for the larger Gaussian cluster
#' @param s2 a number indicating the base radius for each of the two smaller Gaussian clusters
#' @param d a number indicating the dimension of the generated data. d>2 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = MixMickey_data()$data
#' plot(X[,1], X[,2])
MixMickey_data = function(n1=2000,n2=600,s1=2,s2=1, d=2, sd_high=0.1){
  X1 = mvtnorm::rmvnorm(n1, mean = c(0,0), sigma = s1*diag(2))
  X2 = mvtnorm::rmvnorm(n2, mean = c(3,3), sigma = s2*diag(2))
  X3 = mvtnorm::rmvnorm(n2, mean = c(-3,3), sigma = s2*diag(2))
  X = rbind(X1,X2,X3)
  if(d==2){
    X0 = X
  }
  if(d>2){
    n = nrow(X)
    d_add = d-2
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  Y0 = c(rep(1,n1), rep(2,n2), rep(3,n2))
  return(list(data = X0, clus = Y0))
}


#' Generating MixStar data
#'
#' @param n a number indicating the number of points in each of the Gaussian clusters
#' @param d a number indicating the dimension of the generated data. d>2 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = MixStar_data()$data
#' plot(X[,1], X[,2])
MixStar_data = function(n=1000, d=2, sd_high=0.1){
  X1 = mvtnorm::rmvnorm(n, mean = c(4,0), sigma = matrix(c(5,0,0,0.3), ncol = 2, byrow = TRUE))
  rotate1 = matrix(c( sqrt(3)/2, 1/2, -1/2, sqrt(3)/2), ncol = 2, byrow = TRUE)
  X1 = X1 %*% rotate1
  X2 = mvtnorm::rmvnorm(n, mean = c(-4,0),  sigma = matrix(c(5,-0,-0,0.3), ncol = 2, byrow = TRUE))
  rotate2 = matrix(c( sqrt(3)/2, -1/2, 1/2, sqrt(3)/2), ncol = 2, byrow = TRUE)
  X2 = X2 %*% rotate2
  X3 = mvtnorm::rmvnorm(n, mean = c(0,-4),  sigma = matrix(c(0.3,-0,-0,5), ncol = 2, byrow = TRUE))
  X = rbind(X1,X2,X3)
  if(d==2){
    X0 = X
  }
  if(d>2){
    n = nrow(X)
    d_add = d-2
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  Y0 = c(rep(1,n), rep(2,n), rep(3,n))
  return(list(data = X0, clus = Y0))
}


#' Generating Yinyang data
#'
#' @param n_m number of points in each of the two moon shapes. Positive integer value or vector with length=2
#' @param n_c number of points in each of the two small circles
#' @param n_r number of points in the outer ring
#' @param sd_c the standard deviation of shifts for points in the two small circles
#' @param sd_r the standard deviation of shifts for points in the outer ring
#' @param d a number indicating the dimension of the generated data. d>2 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = Yinyang_data()$data
#' plot(X[,1], X[,2])
Yinyang_data = function(n_m=400,n_c=200,n_r=2000,sd_c=0.1,sd_r=0.1, d=2, sd_high=0.1){
  X_m = clusterSim::shapes.two.moon(n_m)

  x1 = cbind(stats::rnorm(n_c, mean= 0.5, sd=sd_c), stats::rnorm(n_c, mean= -1.5, sd=sd_c))
  x2 = cbind(stats::rnorm(n_c, mean= -1, sd=sd_c), stats::rnorm(n_c, mean= 0.5, sd=sd_c))

  th = stats::runif(n_r,0,2*pi)
  x  = 2.5*cos(th) - 0.25
  y  = 2.5*sin(th) -0.5
  x3  = cbind(x,y) + matrix(stats::rnorm(2*length(x), sd=sd_r), ncol=2)

  X = rbind(X_m$data, x1, x2, x3)

  if(d==2){
    X0 = X
  }
  if(d>2){
    n = nrow(X)
    d_add = d-2
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  Y0 = c(X_m$clusters, rep(3,n_c), rep(4,n_c), rep(5,n_r))
  return(list(data = X0, clus = Y0))
}

#' Generating Manifold Mixture data
#'
#' @param n1 number of points in the ring structure
#' @param n2 number of points in the perturbed line segment structure
#' @param n3 number of points in the 2-dimensional plane structure
#' @param sd1 the standard deviation of shifts for points in the ring structure
#' @param sd2 the standard deviation of shifts for points in the perturbed line segment structure
#' @param sd3 the standard deviation of shifts for points in the 2-dimensional plane structure
#' @param d a number indicating the dimension of the generated data. d>3 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = MM_data()$data
MM_data = function(n1=800,n2=400,n3=2000,sd1=0.04,sd2=0.04,sd3=0.04, d=3, sd_high=0.1){
  th = stats::runif(n1,0,3*pi)
  x  = cos(th)
  y  = sin(th)
  z  = 1
  X1  = cbind(x,y,z) + matrix(stats::rnorm(3*n1, sd=sd1), ncol=3)

  x = 0
  y = 0
  z = rep(0.5, n2)
  X2  = cbind(x,y,z) + matrix(stats::rnorm(3*n2, sd=sd2), ncol=3)


  x = stats::runif(n = n3, min= -1, max=1)
  y = stats::runif(n = n3, min= -1, max=1)
  z = 0
  X3 = cbind(x,y,z) + matrix(stats::rnorm(3*n3, sd=sd3), ncol=3)

  X = rbind(X1,X2,X3)


  if(d==3){
    X0 = X
  }
  if(d>3){
    n = nrow(X)
    d_add = d-3
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  Y0 = c(rep(1,n1), rep(2,n2), rep(3,n3))
  return(list(data = X0, clus = Y0))
}


#' Generating SwissRoll Regression data
#'
#' @param N the sample size
#' @param Height the maximum radius
#' @param start the starting angle for the Swiss Roll
#' @param roll the number of rolls
#' @param d a number indicating the dimension of the generated data. d>3 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @param sd_e a number indicating the standard deviation for the additive Gaussian noises in the response
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = MM_data()$data
SwissRollReg <- function(N=2000, Height=4*pi,start = pi, roll=1, d = 3, sd_high=0.1, sd_e = 0.3){

  ## build manifold
  c = log((start + 2*pi*roll)/start) #need start > 0
  p = start * exp(c *runif(N, 0, 1))
  y = Height * runif(N, 0 , 1)
  X = cbind(p * cos(p), y, p * sin(p))
  X = scale(X)

  if(d==3){
    X0 = X
  }
  if(d>3){
    n = nrow(X)
    d_add = d-3
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  p0 = p - start -  pi * roll
  Y0 = 0.1*p0^3*( (y<pi) + (y > 2*pi & y<3*pi) ) + rnorm(N, mean = 0, sd = sd_e)

  return(list(data = X0, Y = Y0, angles = p, heights = y))

}



#' Generating Yinyang regression data
#'
#' @param n_m number of points in each of the two moon shapes. Positive integer value or vector with length=2
#' @param n_c number of points in each of the two small circles
#' @param n_r number of points in the outer ring
#' @param sd_c the standard deviation of shifts for points in the two small circles
#' @param sd_r the standard deviation of shifts for points in the outer ring
#' @param d a number indicating the dimension of the generated data. d>2 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = Yinyang_data()$data
#' plot(X[,1], X[,2])
Yinyang_reg_data = function(n_m=400,n_c=200,n_r=2000,sd_c=0.1,sd_r=0.1, d=2, sd_high=0.1, sd_e = 0.1){

  X_m = clusterSim::shapes.two.moon(n_m)

  x1 = cbind(stats::rnorm(n_c, mean= 0.5, sd=sd_c), stats::rnorm(n_c, mean= -1.5, sd=sd_c))
  x2 = cbind(stats::rnorm(n_c, mean= -1, sd=sd_c), stats::rnorm(n_c, mean= 0.5, sd=sd_c))

  th = stats::runif(n_r,0,2*pi)
  x  = 2.5*cos(th) - 0.25
  y  = 2.5*sin(th) -0.5
  x3  = cbind(x,y) + matrix(stats::rnorm(2*length(x), sd=sd_r), ncol=2)
  Y_r = sin(th*4) + 1.5

  X = rbind(X_m$data, x1, x2, x3)

  if(d==2){
    X0 = X
  }
  if(d>2){
    n = nrow(X)
    d_add = d-2
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  Y0 = c(X_m$clusters, rep(0,n_c), rep(3,n_c), Y_r) + rnorm(nrow(X_m$data)+2*n_c+n_r, mean = 0, sd = sd_e)
  return(list(data = X0, Y = Y0))

}

#' Generating Yinyang regression data
#'
#' @param n_m number of points in each of the two moon shapes. Positive integer value or vector with length=2
#' @param n_c number of points in each of the two small circles
#' @param n_r number of points in the outer ring
#' @param sd_c the standard deviation of shifts for points in the two small circles
#' @param sd_r the standard deviation of shifts for points in the outer ring
#' @param d a number indicating the dimension of the generated data. d>2 dimensions are filled with independent Gaussian noises
#' @param sd_high a number indicating the standard deviation for Gaussian noises in the higher dimensions
#' @return A list of generated data points and the true labels
#' \itemize{
#'   \item data - The generated data matrix
#'   \item clus - The vector of true cluster labels for each data point
#' }
#' @export
#' @examples
#' X = Yinyang_data()$data
#' plot(X[,1], X[,2])
NoiseYinyang_reg_data = function(nnoise = 800, n_m=400,n_c=200,n_r=2000,sd_c=0.1,sd_r=0.1, d=2, sd_high=0.1, sd_e = 0.1){

  X_m = clusterSim::shapes.two.moon(n_m)

  x1 = cbind(stats::rnorm(n_c, mean= 0.5, sd=sd_c), stats::rnorm(n_c, mean= -1.5, sd=sd_c))
  x2 = cbind(stats::rnorm(n_c, mean= -1, sd=sd_c), stats::rnorm(n_c, mean= 0.5, sd=sd_c))

  th = stats::runif(n_r,0,2*pi)
  x  = 2.5*cos(th) - 0.25
  y  = 2.5*sin(th) -0.5
  x3  = cbind(x,y) + matrix(stats::rnorm(2*length(x), sd=sd_r), ncol=2)
  Y_r = sin(th*4) + 1.5

  x4 = cbind( stats::runif(nnoise, -3.5,3.5), stats::runif(nnoise, -3.5,3.5))

  X = rbind(X_m$data, x1, x2, x3, x4)

  if(d==2){
    X0 = X
  }
  if(d>2){
    n = nrow(X)
    d_add = d-2
    X0 = cbind(X, matrix(stats::rnorm(d_add*(n),  sd = sd_high), nrow = n))
  }
  Y0 = c(X_m$clusters, rep(0,n_c), rep(3,n_c), Y_r, rep(1.5,nnoise)) + rnorm(nrow(X_m$data)+2*n_c+n_r+nnoise, mean = 0, sd = sd_e)
  return(list(data = X0, Y = Y0))

}
