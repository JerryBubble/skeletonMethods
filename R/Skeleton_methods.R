#' Vectorized pairwise Euclidean distance calculation
#'
#' @param A,B Matrices of data points with same number of columns
#' @return A matrix of pairwise Euclidean distances
#' @export
#' @examples
#' A = matrix(c(1,2,3,4,5,6,7,8), ncol = 2, byrow = TRUE)
#' vectorized_pdist(A,A)
vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}


#' Construct the skeleton representation of a dataset
#'
#' @param data a matrix/dataframe of data points.
#' @param centers an optional matrix/dataframe of the knots in the skeleton representation.
#' @param labels an optional numeric vector of skeleton membership labels for each data points.
#' @param k a number indication the numer of knots to use for skeleton representation.
#' @param rep a number indicating the repetition when overfitting kmeans to choose knots.
#' @param wedge a character or a vector of characters indicating the types of edge measure weights to include.
#' Can take 'all', 'none', 'voronoi', 'face', 'frustum', 'avedist'.
#' @param h a number for the bandwidth when using KDE to calculate Face or Frustum density.
#' @param kernel a character string giving the smoothing kernel to be used in KDE. Same as in the density function. This must partially match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "gaussian", and may be abbreviated to a unique prefix.
#' @param R0 a number indicating the disk radius for Frustum density
#' @param idx_frustum  logical; if TRUE, use same disk radius for different pairs of knots. If FALSE, use different within cluster variance as radius
#' @return A list with following components:
#' \describe{
#'   \item{centers}{The matrix of the skeleton knots.}
#'   \item{cluster}{The numeric vector indicating the neareast knot for each point.}
#'   \item{nknots}{The number indicating the number of knots.}
#'   \item{voron_weights}{The matrix of Voronoi density weights for edges in the skeleton.}
#'   \item{face_weights}{The matrix of Face density weights for edges in the skeleton.}
#'   \item{frustum_weights}{The matrix of Tube/Frustum density weights for edges in the skeleton.}
#'   \item{avedist_weights}{The matrix of Average Distance weights for edges in the skeleton.}
#'   \item{bw}{The bandwidth used for Face and Frustum density calculation. }
#'   \item{R}{The disk radius used for Frustum density calculation. }
#' }
#' @export
#' @examples
#' ###create ring data
#' n_c = 200
#' n_r=1000
#' sd_c = 0.1
#' sd_r = 0.1
#' d=2
#' sd_high = 0.1
#' c_lab = c(rep(1,n_r), rep(2,n_c))
#' th = runif(n_r,0,2*pi)
#' x  = cos(th)
#' y  = sin(th)
#' X  = cbind(x,y) + matrix(rnorm(2*length(x), sd=sd_r), ncol=2)
#' u = matrix(rnorm(2*n_c, sd=sd_c), ncol=2)
#' X0 = rbind(X, u) #created Ring data
#' ###Construct skeleton representation
#' skeleton = skeletonCons(X0, rep = 1000)
#' ### Voronoi density
#' F_NNden_tmp = max(skeleton$voron_weights) - skeleton$voron_weights
#' diag(F_NNden_tmp) = 0
#' F_NNden_tmp = as.dist(F_NNden_tmp)
#' F_NNden_hclust = hclust(F_NNden_tmp, method="single") #used average distance instead
#' plot(F_NNden_hclust, main = 'Cluster Dendrogram with Single Linkage')
#' F_NNden_lab = cutree(F_NNden_hclust, k=2)
#' X_lab_F_NNden = F_NNden_lab[skeleton$cluster]
#' plot(X[,1], X[,2], col = X_lab_F_NNden+1, pch = 19)
skeletonCons = function(data, centers=NULL, labels=NULL,
                        k=NA, #the number of knots for skeleton
                        rep=1000, #number of repetition for overfitting kmeans choosing knots
                        wedge='all', #types of edge weights to include
                        #wedge can take 'all', 'none', 'voronoi', 'face', 'frustum', 'avedist'
                        h=NA, #bandwidth for KDE for Face and Frustum density
                        hadj = 1/5, #adjust the rate of bandwidth with sample size
                        kernel = "gaussian", #kernel for KDE for Face and Frustum density
                        R0=NA, #disk radius for Frustum density
                        idx_frustum = T #use same disk radius for different pairs of knots
                        #otherwise use different within cluster variance as radius
){
  X_km = NULL
  n = nrow(data)
  d = ncol(data)
  
  X0 = as.matrix(data)
  
  #construct knots
  if(is.null(centers)&is.null(labels)){
    
    # Overfitting k-means
    #setting the number of knots k
    if(is.na(k)){
      k = round(sqrt(n))
      # sqrt-n rule
    }
    X_km = stats::kmeans(as.matrix(X0), center=k, nstart = rep)
    centers = X_km$centers
    labels = X_km$cluster
  }else if(is.null(labels)){#centers provided but not labels
    centers = as.matrix(centers)
    labels = RANN::nn2(centers, X0, k=1)$nn.idx[,1]
    k = nrow(centers)
    withinss = numeric(k)
    for (i in 1:k) {
      withinss[i] = sum(apply(t(X0[labels == i,]) - centers[i,], 1, FUN=function(x) sum(x^2)))
    }
    X_km$k = k
    X_km$withinss = withinss
    X_km$size = as.numeric(table(labels))
  }else if(is.null(centers)){#labels provided but not centers
    k = sum(table(labels) >3)
    centers = matrix(nrow = k, ncol = d)
    withinss = numeric(k)
    validlabels = as.numeric(levels(as.factor(labels))[table(labels) >3])
    rowid = 1
    for (i in validlabels) {
      centers[rowid,] = apply(X0[labels == i,], 2, mean)
      withinss[rowid] = sum(apply(t(X0[labels == i,]) - centers[rowid,], 1, FUN=function(x) sum(x^2)))
      rowid = rowid+1
    }
    X_km$k = k
    X_km$withinss = withinss
    X_km$size = as.numeric(table(labels))
  }else{#both centers and labels provided
    centers = as.matrix(centers)
    labels = as.numeric(labels)
    k = length(levels(as.factor(labels)))
    withinss = numeric(k)
    for (i in 1:k) {
      withinss[i] = sum(apply(t(X0[labels == i,]) - centers[i,], 1, FUN=function(x) sum(x^2)))
    }
    X_km$k = k
    X_km$withinss = withinss
    X_km$size = as.numeric(table(labels))
  }
  
  output = list()
  output$centers = centers
  output$cluster = labels
  output$nknots = k
  
  #identify what edge weights to include
  edge_include = c(F,F,F,F)
  if(length(wedge) == 1){
    if(wedge == 'all'){
      edge_include = c(T,T,T,T)
    }else{
      edge_include = c('voronoi', 'face', 'frustum', 'avedist') %in% wedge
    }
  }else if(length(wedge) > 1){
    edge_include = c('voronoi', 'face', 'frustum', 'avedist') %in% wedge
  }
  
  if(sum(edge_include) == 0){#stop if no edge weight method specified
    return(output)
  }
  
  
  m = max(labels)
  
  #edge weight matrices
  if(edge_include[1]){voron_weights = matrix(0,m,m)}
  if(edge_include[2]){face_weights = matrix(0,m,m)}
  if(edge_include[3]){
    frustum_weights = matrix(0,m,m)
    #calculate disk radius
    if(is.na(R0)){ #fill in average within cluster variance if not specified
      R_cluster = sqrt(X_km$withinss/(X_km$size-1))
      R0 = mean(R_cluster)
      if(idx_frustum){#same radius for all pairs
        R0_lv = rep(R0, k)
      }else{
        R0_lv = R_cluster
      }
    }else if(length(R0)==1){#one specified R0
      R0_lv = rep(R0, k)
      idx_frustum = T
    }else if(length(R0) == k){
      R0_lv = R0
      idx_frustum = F
    }else{#catch bad cases
      print('wrong R0 specified, use default')
      R_cluster = sqrt(X_km$withinss/(x_km$size-1))
      R0 = mean(R_cluster)
      R0_lv = rep(R0, k)
      idx_frustum = T
    }
  }
  if(edge_include[4]){avedist_weights = matrix(0,m,m)}
  
  
  X_nn = RANN::nn2(centers, X0, k=2) #2-nearest neighbor calculation
  
  for(i in 1:(m-1)){
    center1 = centers[i,]
    wi1 = which(X_nn$nn.idx[,1]==i)
    wi2 = which(X_nn$nn.idx[,2]==i)
    for(j in (i+1):m){
      center2 = centers[j,]
      wj1 = which(X_nn$nn.idx[,1]==j)
      wj2 = which(X_nn$nn.idx[,2]==j)
      nn2ij = union(intersect(wi1, wj2), intersect(wi2, wj1))#2nn neighborhood
      
      
      if(length(nn2ij) <= 1 ){#not in Delaunnay Triangulation
        if(edge_include[1]){voron_weights[i,j] = 0}
        if(edge_include[2]){
          face_weights[i,j] = 0}
        if(edge_include[3]){
          frustum_weights[i,j] = 0}
        if(edge_include[4]){avedist_weights[i,j] = 0}
        next
      }
      
      d12 = sqrt(sum((center1-center2)^2))
      
      if(edge_include[1]){#Voronoi density
        voron_weights[i,j] = length(nn2ij)/d12
      }
      
      
      if(edge_include[2]){#face density
        v0 = (center2-center1)/d12 #direction vector
        p0_length = colSums((t(X0)-(center1+center2)/2)*v0) #projected distance to middle point
        p_dot = p0_length/d12 #standardize the projected distances
        
        if(!is.numeric(h)){ #bandwidth selection
          if(is.na(h)){
            h1 = ks::hns(p_dot[nn2ij])
          }else if(h=="hns"){
            h1 = ks::hns(p_dot[nn2ij])
          }else if(h=="hlscv"){
            h1 = ks::hlscv(p_dot[nn2ij])
          }else{
            h1 = ks::hpi(p_dot[nn2ij])
          }
          h1 = h1*(length(nn2ij)^{1/5-hadj}) #adjusting the rate of h with sample size
        }else{h1 = h}
        
        
        den_proj = stats::density(p_dot[nn2ij], kernel = kernel, bw=abs(h1), from = -1, to = 1) #KDE with projected points
        face_weights[i,j] = (den_proj$y[256]+den_proj$y[257])/2 #interpolated density at middle point
      }#end face density calculation
      
      if(edge_include[3]){#frustum density
        v0 = (center2-center1)/d12 #direction vector
        p0_length = colSums((t(X0)-center1)*v0) #projected distance to center1
        p_dot = p0_length/d12 #standardize the projected distances
        c1_length = sqrt(colSums((t(X0)-center1)^2)) #total distance to center1
        perp_length = sqrt(c1_length^2-p0_length^2) #orthogonal distance to the center-passing line
        
        ### setting the threshold for each disk
        R0_threshold = R0_lv[i] + (R0_lv[j]-R0_lv[i])*p_dot
        R0_threshold[which(p_dot<= 0)] = R0_lv[i]
        R0_threshold[which(p_dot> 1)] = R0_lv[j]
        
        w_edge = which(p_dot>= -1.5 &p_dot<= 1.5 &perp_length<R0_threshold)#points that can be used in KDE
        
        if(length(w_edge)>1){#KDE works with more than 1 data points
          if(!is.numeric(h)){ #bandwidth selection
            if(is.na(h)){
              h1 = ks::hns(p_dot[nn2ij])
            }else if(h=="hns"){
              h1 = ks::hns(p_dot[nn2ij])
            }else if(h=="hlscv"){
              h1 = ks::hlscv(p_dot[nn2ij])
            }else{
              h1 = ks::hpi(p_dot[nn2ij])
            }
            h1 = h1*(length(nn2ij)^{1/5-hadj}) #adjusting the rate of h with sample size
          }else{h1 = h}
          
          
          den_proj = stats::density(p_dot[w_edge], kernel = kernel, bw= abs(h1), from = -0.5, to = 1.5) #KDE with projected points
          # finding the minimal density
          min_den_proj = min(den_proj$y[which(den_proj$x>0&den_proj$x< 1)])
          frustum_weights[i,j] = min_den_proj
        }else{frustum_weights[i,j] = 0}
      }#end frustum density calculation
      
      if(edge_include[4]){#avedist density
        dists = vectorized_pdist(X0[wi1,], X0[wj1,])
        avedist_weights[i,j] = mean(dists)
      }#end avedist density
    } #end for j
  } #end for i
  
  if(edge_include[1]){output$voron_weights = voron_weights + t(voron_weights)}
  if(edge_include[2]){
    output$face_weights = face_weights + t(face_weights)
    output$bw = h1
  }
  if(edge_include[3]){
    output$frustum_weights = frustum_weights + t(frustum_weights)
    output$R = R0_lv
    output$bw = h1
  }
  if(edge_include[4]){output$avedist_weights = avedist_weights+t(avedist_weights)}
  
  return(output)
}




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
    ntol = nrow(X)
    d_add = d-2
    X0 = cbind(X, matrix(stats::rnorm(d_add*(ntol),  sd = sd_high), nrow = ntol))
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





