#' skeleton distance calculation
#' , , , , ,
#' @param nnc1 a vector of the indexes of the two nearest neighbors of one point.
#' @param nnc2 a vector of the indexes of the two nearest neighbors of another point.
#' @param px1 the projection proportion of point 1 onto the edge by its two nearest knots, with the lower index knot as the origin.
#' @param px2 the projection proportion of point 2 onto the edge by its two nearest knots, with the lower index knot as the origin.
#' @param skeleton a list returned by the voronSkeleton method constructing the skeleton graph.
#' @return skeleton-based distance between two points on the skeleton
#' @export
dskeleton = function(nnc1, nnc2, px1, px2, skeleton) {
  g = skeleton$g
  edgedists = skeleton$edgedists

  #simplify calculation by only focusing on pair of points sharing at least one knot
  if(length(intersect(nnc1,nnc2))<1){ return(Inf)}

  if(is.na(px1) & is.na(px2)){
    #case when both data points cannot be projected onto an edge
    return(igraph::distances(g, v=nnc1[1], to =nnc2[1]))

  }else if(!is.na(px1) & is.na(px2) ){
    #case when point id1 can be projected onto an edge
    btnodes = c(igraph::distances(g, v=nnc1[1], to =nnc2[1]),
                igraph::distances(g, v=nnc1[2], to =nnc2[1]))
    if(!is.finite(min(btnodes))){
      return(Inf)
    }
    if(which.min(btnodes) == 1){#shortest path from nnc1[1]
      if(nnc1[1] < nnc1[2]){ #nnc1[1] is the reference point for projection
        return(min(btnodes) + edgedists[nnc1[1], nnc1[2]]*px1 )
      }else{
        return(min(btnodes) + edgedists[nnc1[1], nnc1[2]]*(1-px1) )
      }
    }else{
      if(nnc1[1] < nnc1[2]){ #nnc1[1] is the reference point for projection
        return(min(btnodes) + edgedists[nnc1[1], nnc1[2]]*(1-px1) )
      }else{
        return(min(btnodes) + edgedists[nnc1[1], nnc1[2]]*px1 )
      }
    }

  }else if(is.na(px1) & !is.na(px2)){
    #case when point id2 can be projected onto an edge
    btnodes = c(igraph::distances(g, v=nnc1[1], to =nnc2[1]),
                igraph::distances(g, v=nnc1[1], to =nnc2[2]))
    if(!is.finite(min(btnodes))){
      return(Inf)
    }
    if(which.min(btnodes) == 1){#shortest path from nnc2[1]
      if(nnc2[1] < nnc2[2]){#nnc2[1] is the reference point for projection
        return(min(btnodes) + edgedists[nnc2[1], nnc2[2]]*px2 )
      }else{
        return(min(btnodes) + edgedists[nnc2[1], nnc2[2]]*(1-px2) )
      }
    }else{
      if(nnc2[1] < nnc2[2]){#nnc2[1] is the reference point for projection
        return(min(btnodes) + edgedists[nnc2[1], nnc2[2]]*(1-px2) )
      }else{
        return(min(btnodes) + edgedists[nnc2[1], nnc2[2]]*px2 )
      }
    }

  }else{ #case when both points are projected onto edges

    btnodes = c(igraph::distances(g, v=nnc1[1], to =nnc2[1]),
                igraph::distances(g, v=nnc1[1], to =nnc2[2]),
                igraph::distances(g, v=nnc1[2], to =nnc2[1]),
                igraph::distances(g, v=nnc1[2], to =nnc2[2]))
    if(!is.finite(min(btnodes))){
      return(Inf)
    }

    temp = min(btnodes)

    if(which.min(btnodes) == 1){ #nnc1[1] and nnc2[1] are closest
      if(nnc1[1] < nnc1[2]){#nnc1[1] is the reference point for id1 projection
        temp = temp + edgedists[nnc1[1], nnc1[2]]*px1
      }else{
        temp = temp + edgedists[nnc1[1], nnc1[2]]*(1-px1)
      }
      if(nnc2[1] < nnc2[2]){#nnc2[1] is the reference point for id2 projection
        temp = temp + edgedists[nnc2[1], nnc2[2]]*px2
      }else{
        temp = temp + edgedists[nnc2[1], nnc2[2]]*(1-px2)
      }
      return(temp)
    }else if(which.min(btnodes) == 2){ #nnc1[1] and nnc2[2] are closest
      if(nnc1[1] < nnc1[2]){#nnc1[1] is the reference point for id1 projection
        temp = temp + edgedists[nnc1[1], nnc1[2]]*px1
      }else{
        temp = temp + edgedists[nnc1[1], nnc1[2]]*(1-px1)
      }
      if(nnc2[1] < nnc2[2]){#nnc2[1] is the reference point for id2 projection
        temp = temp + edgedists[nnc2[1], nnc2[2]]*(1-px2)
      }else{
        temp = temp + edgedists[nnc2[1], nnc2[2]]*px2
      }
      return(temp)
    }else if(which.min(btnodes) == 3){ #nnc1[2] and nnc2[1] are closest
      if(nnc1[1] < nnc1[2]){#nnc1[1] is the reference point for id1 projection
        temp = temp + edgedists[nnc1[1], nnc1[2]]*(1-px1)
      }else{
        temp = temp + edgedists[nnc1[1], nnc1[2]]*px1
      }
      if(nnc2[1] < nnc2[2]){#nnc2[1] is the reference point for id2 projection
        temp = temp + edgedists[nnc2[1], nnc2[2]]*px2
      }else{
        temp = temp + edgedists[nnc2[1], nnc2[2]]*(1-px2)
      }
      return(temp)
    }else{ #nnc1[2] and nnc2[2] are closest
      if(nnc1[1] < nnc1[2]){#nnc1[1] is the reference point for id1 projection
        temp = temp + edgedists[nnc1[1], nnc1[2]]*(1-px1)
      }else{
        temp = temp + edgedists[nnc1[1], nnc1[2]]*px1
      }
      if(nnc2[1] < nnc2[2]){#nnc2[1] is the reference point for id2 projection
        temp = temp + edgedists[nnc2[1], nnc2[2]]*(1-px2)
      }else{
        temp = temp + edgedists[nnc2[1], nnc2[2]]*px2
      }
      return(temp)
    }
  }

}

#' Construct the skeleton graph of a dataset based on Voronoi denisty weights
#'
#' @param trainX a matrix/dataframe of training data points.
#' @param centers an optional matrix/dataframe of the knots in the skeleton representation.
#' @param labels an optional numeric vector of skeleton membership labels for each data points.
#' @param numknots a number indicating the number of knots to use for skeleton representation.
#' @param k_cut a number indicating the number of resulting disconnected components in the skeleton representation.
#' @param rep a number indicating the repetition when overfitting kmeans to choose knots.
#' @param seed the random seed for the construction
#' @return A list with following components:
#' \describe{
#'   \item{centers}{The matrix of the skeleton knots.}
#'   \item{labels}{The numeric vector indicating the cluser membership label of each center.}
#'   \item{kdists}{The matrix of the pairwise Euclidean distances between the centers.}
#'   \item{edgedists}{The matrix of edges in the skeleton with lengths. Entries in kdists not an edge set to 0.}
#'   \item{nn}{The matrix where each row records the indices of the two closest centers of the corresponding data point.}
#'   \item{g}{The igraph object for the constructed skeleton graph.}
#' }
#' @export
voronSkeleton = function(trainX,centers=NULL, labels=NULL,
                         numknots, k_cut=1, rep=1000, seed=NULL){
  tol = 1e-10
  if(!is.null(seed)){set.seed(seed)}
  skeleton = skeletonCons(data = trainX, centers=centers, labels=labels,
                          k= numknots , rep = rep,wedge = "voronoi")

  centers = as.matrix(skeleton$centers)
  labels = as.numeric(skeleton$cluster)

  #Voronoi density similarities
  voronweights = skeleton$voron_weights

  #choose correct cut
  F_NNden_tmp = max(skeleton$voron_weights) - skeleton$voron_weights
  diag(F_NNden_tmp) = 0
  F_NNden_tmp = stats::as.dist(F_NNden_tmp)
  F_NNden_hclust = stats::hclust(F_NNden_tmp, method="single")
  cutheight = F_NNden_hclust$height[length(F_NNden_hclust$height)-k_cut+2]
  #cut some edges
  voronweights[voronweights<(max(voronweights)-cutheight+tol)] = 0

  #Euclidean distance between centers
  kdists = vectorized_pdist(centers,centers)

  edgedists = kdists #use euclidean distance as edge weights
  edgedists[voronweights == 0] = 0

  # ###############################################
  # # get the graph based on Euclidean similarities
  # ###############################################
  g = igraph::graph.adjacency(edgedists, mode='undirected', weighted=TRUE)
  # coords_fr = layout.fruchterman.reingold(g, weights=E(g)$weight)
  # # igraph plot options
  # igraph.options(vertex.size=8, edge.width=0.75)
  # # plot network
  # plot(g, layout=coords_fr)
  return(list(centers = centers,
              labels = labels,
              kdists = kdists,
              edgedists = edgedists,
              g = g,
              nn = skeleton$nn))
}



#' Project data points onto the skeleton
#'
#' @param nn the matrix where each row records the indices of the two closest centers of the corresponding data point.
#' @param X the matrix of the data points to be projected.
#' @param skeleton a list returned by the voronSkeleton method constructing the skeleton graph.
#' @return A vector where each entry records the projection proportion of the corresponding data point on the edge by the two closest knots, from the smaller index knot. NA if the two closest knots not connected.
#' @export
skelProject = function(nn = NULL, X, skeleton){
  if(is.null(nn)){
    nn = RANN::nn2(centers, X, k=2)$nn.idx}
  centers = skeleton$centers
  kdists = skeleton$kdists
  g = skeleton$g
  #####
  #calculate the projection for each data point
  ######
  #a proportion between [0,1] is recorded
  px = numeric(nrow(X))
  for (i in 1:nrow(X)) {
    #For case x cannot be projected to an edge, project to nearest center
    if(igraph::get.edge.ids(g, nn[i,])==0){
      px[i] = NA
    }else{
      k0 = min(nn[i,])
      k1 = max(nn[i,])
      #directional vector
      #always use the smaller index as the starting point
      v = (centers[k1,]-centers[k0,])/kdists[k0,k1]
      prop =(X[i,]-centers[k0,])%*%v/kdists[k0,k1] #projected proportion
      px[i] = max(min(prop, 1),0) #confine to [0,1]
    }
  }
  #####
  #End calculate the projection for each data point
  ######
  return(px)
}



#' Calculating the base bandwidth for the S-kernel regressor
#'
#' @param centernn the matrix where each row records the indices of the two closest knots of each knot.
#' @param px the projection proportion of the training set returned by the skelProject function.
#' @param skeleton a list returned by the voronSkeleton method constructing the skeleton graph.
#' @return A numeric value for the global bandwidth to be used.
#' @export
skelBand = function(centernn=NULL, px, skeleton){
  centers = skeleton$centers
  kdists = skeleton$kdists
  g = skeleton$g
  nn = skeleton$nn
  if(is.null(centernn)){
    centernn = RANN::nn2(centers, centers, k=2)$nn.idx}
  centerpx=skelProject(nn = centernn, centers, skeleton = skeleton)

  centerSkeldists = matrix(Inf,nrow = nrow(centers), ncol = nrow(nn))
  for (i in 1:nrow(centers)) {
    #calculate graph distances between new obs and sample points
    for (j in 1:nrow(nn)) {
      centerSkeldists[i,j] = dskeleton(nn[j,], centernn[i,], px[j], centerpx[i], skeleton = skeleton)
    }
  }

  centerbands = numeric(nrow(centers))
  for (i in 1:nrow(centers)) {
    centerbands[i] = ks::hns(centerSkeldists[i,is.finite(centerSkeldists[i,])])
  }
  centerband = mean(centerbands, na.rm = T)
  return(centerband)
}



#' Fitting the S-kernel regressor with fixed bandwidth
#'
#' @param h1 the fixed bandwidth for kernel regression
#' @param testX the matrix of the test data points.
#' @param testSkeldists the matrix of skeleton-based distances between test and training points.
#' @param trainY the responses of the training points.
#' @return A vector of predicted responses for the test points.
#' @export
skelKernel = function(h1, testX,testSkeldists, trainY){
  tol = 1e-10
  testfit_tmp2 = numeric(nrow(testX))
  #predictions on test data
  for (i in 1:nrow(testX)) {
    x = testSkeldists[i,]
    useid = is.finite(x)
    x = x[useid]
    if(length(x)==1){ #happens when there is only 1 finite distance training point
      testfit_tmp2[i] = trainY[useid]
    }else if(stats::sd(x) < tol){ # case when all x are 0, which happens on singly separated knot
      testfit_tmp2[i] = mean(trainY[useid])
    }else{
      testfit_tmp2[i] = stats::ksmooth(x,trainY[useid],"normal", bandwidth = h1, x.points = c(0))$y
      if(is.na(testfit_tmp2[i])){
        testfit_tmp2[i] = stats::ksmooth(x,trainY[useid],"normal", bandwidth = stats::quantile(x,0.25), x.points = c(0))$y
      }
    }
  }
  return(testfit_tmp2)
}



#' Fitting the S-kernel regressor with adaptive bandwidth
#'
#' @param hrate the scaling constant to apply to the adaptive calculated bandwidth
#' @param testX the matrix of the test data points.
#' @param testSkeldists the matrix of skeleton-based distances between test and training points.
#' @param trainY the responses of the training points.
#' @return A vector of predicted responses for the test points.
#' @export
skelKerneladopt = function(hrate, testX,testSkeldists, trainY){
  testfit_tmp = numeric(nrow(testX))
  #predictions on test data
  for (i in 1:nrow(testX)) {
    x = testSkeldists[i,]
    useid = is.finite(x)
    x = x[useid]
    h1 = (ks::hns(x))*hrate
    if(is.na(h1)){ # happens when the only 1 finite distance training point
      testfit_tmp[i] = trainY[useid]
    }
    if(h1 == 0){ # case when all x are 0, which happens on singly separated knot
      testfit_tmp[i] = mean(trainY[useid])
    }else{
      testfit_tmp[i] = stats::ksmooth(x,trainY[useid],"normal", bandwidth = h1, x.points = c(0))$y
      if(is.na(testfit_tmp[i])){
        testfit_tmp[i] = stats::ksmooth(x,trainY[useid],"normal", bandwidth = stats::quantile(x,0.25), x.points = c(0))$y
      }
    }
  }
  return(testfit_tmp)
}


#' Fitting the S-kNN regressor
#'
#' @param k the number of neighbors used in k-nearest-neighbor regression.
#' @param testX the matrix of the test data points.
#' @param testSkeldists the matrix of skeleton-based distances between test and training points.
#' @param trainY the responses of the training points.
#' @return A vector of predicted responses for the test points.
#' @export
skelknn = function(k=3, testX,testSkeldists, trainY){
  gnntestfit_tmp = numeric(nrow(testX))
  for (i in 1:nrow(testX)) {
    x = testSkeldists[i,]
    nnid = order(x)[1:k]
    gnntestfit_tmp[i] = mean(trainY[nnid])
  }
  return(gnntestfit_tmp)
}



#' Fitting the S-Lspline regressor
#'
#' @param newnn the matrix where each row records the indices of the two closest knots of the new test points.
#' @param newpx the projection proportion for the new test points.
#' @param px the projection proportion for the training points.
#' @param trainY the responses of the training points.
#' @param skeleton a list returned by the voronSkeleton method constructing the skeleton graph.
#' @return A list with following components:
#' \describe{
#'   \item{trainZ}{The transformed training data matrix.}
#'   \item{testZ}{The transformed test data matrix.}
#'   \item{model}{The lm object for the resulting regression model.}
#'   \item{pred}{The vector of predicted responses for the test points.}
#' }
#' @export
skelLinear = function(newnn, newpx, px, trainY, skeleton){
  numknots = nrow(skeleton$centers)
  nn = skeleton$nn
  #modified data matrix for training data based projected values
  trainZ = matrix(0, nrow = length(trainY), ncol = numknots)
  for (i in 1:nrow(nn)) {
    if(is.na(px[i])){#where the data point is projected to nearest knot
      trainZ[i,nn[i,1]] = 1
    }else{
      k0 = min(nn[i,])
      k1 = max(nn[i,])
      #Y0 + px(Y1-Y0), hence px(Y1), (1-px)Y0
      trainZ[i, k0] = 1-px[i] #*edgedists[k0,k1]
      trainZ[i, k1] = px[i]
    }
  }

  #modified data matrix for test data based projected values
  testZ = matrix(0, nrow = length(newpx), ncol = numknots)
  for (i in 1:nrow(newnn)) {
    if(is.na(newpx[i])){#where the data point is projected to nearest knot
      testZ[i,newnn[i,1]] = 1
    }else{
      k0 = min(newnn[i,])
      k1 = max(newnn[i,])
      #Y0 + px(Y1-Y0), hence px(Y1), (1-px)Y0
      testZ[i, k0] = 1-newpx[i] #*edgedists[k0,k1]
      testZ[i, k1] = newpx[i]
    }
  }

  colnames(trainZ) = paste0("k", 1:numknots)
  colnames(testZ) = paste0("k", 1:numknots)

  ZtrainDat = as.data.frame(cbind(trainY, trainZ))
  testZ = as.data.frame(testZ)
  lmod = stats::lm(trainY ~.-1, data = ZtrainDat) #no intercept

  lmpred = stats::predict(lmod, newdata = testZ)

  return(list(trainZ = trainZ, testZ = testZ, model = lmod, pred = lmpred))
}

#' Fitting the S-Qspline regressor
#'
#' @param newnn the matrix where each row records the indices of the two closest knots of the new test points.
#' @param newpx the projection proportion for the new test points.
#' @param px the projection proportion for the training points.
#' @param trainY the responses of the training points.
#' @param skeleton a list returned by the voronSkeleton method constructing the skeleton graph.
#' @return A list with following components:
#' \describe{
#'   \item{trainZ}{The transformed training data matrix.}
#'   \item{testZ}{The transformed test data matrix.}
#'   \item{model}{The lm object for the resulting regression model.}
#'   \item{pred}{The vector of predicted responses for the test points.}
#' }
#' @export
skelQuadratic = function(newnn, newpx, px, trainY,skeleton){
  numknots = nrow(skeleton$centers)
  nn = skeleton$nn
  #modified data matrix for training data based projected values
  trainZ = matrix(0, nrow = length(trainY), ncol = (2*numknots))
  for (i in 1:nrow(nn)) {
    if(is.na(px[i])){#where the data point is projected to nearest knot
      trainZ[i,nn[i,1]] = 1
    }else{
      k0 = min(nn[i,])
      k1 = max(nn[i,])
      trainZ[i, k0] = 2*px[i]^3 - 3*px[i]^2 + 1
      trainZ[i, k1] = -2*px[i]^3 + 3*px[i]^2
      trainZ[i, (k0+numknots)] = px[i]^3 - 2*px[i]^2 + px[i]
      trainZ[i, (k1+numknots)] = px[i]^3 - px[i]^2
    }
  }

  #modified data matrix for test data based projected values
  testZ = matrix(0, nrow = length(newpx), ncol = (2*numknots))
  for (i in 1:nrow(newnn)) {
    if(is.na(newpx[i])){#where the data point is projected to nearest knot
      testZ[i,newnn[i,1]] = 1
    }else{
      k0 = min(newnn[i,])
      k1 = max(newnn[i,])
      testZ[i, k0] = 2*newpx[i]^3 - 3*newpx[i]^2 + 1
      testZ[i, k1] = -2*newpx[i]^3 + 3*newpx[i]^2
      testZ[i, (k0+numknots)] = newpx[i]^3 - 2*newpx[i]^2 + newpx[i]
      testZ[i, (k1+numknots)] = newpx[i]^3 - newpx[i]^2
    }
  }

  colnames(trainZ) = c(paste0("k", 1:numknots),paste0("g", 1:numknots))
  colnames(testZ) = c(paste0("k", 1:numknots),paste0("g", 1:numknots))

  ZtrainDat = as.data.frame(cbind(trainY, trainZ))
  testZ = as.data.frame(testZ)
  lmod2 = stats::lm(trainY ~.-1, data = ZtrainDat) #no intercept

  lmpred2 = stats::predict(lmod2, newdata = testZ)

  return(list(trainZ=trainZ, testZ=testZ, model = lmod2, pred = lmpred2))
}


#' Fitting the S-Cspline regressor
#'
#' @param newnn the matrix where each row records the indices of the two closest knots of the new test points.
#' @param newpx the projection proportion for the new test points.
#' @param px the projection proportion for the training points.
#' @param trainY the responses of the training points.
#' @param skeleton a list returned by the voronSkeleton method constructing the skeleton graph.
#' @return A list with following components:
#' \describe{
#'   \item{trainZ}{The transformed training data matrix.}
#'   \item{testZ}{The transformed test data matrix.}
#'   \item{model}{The lm object for the resulting regression model.}
#'   \item{pred}{The vector of predicted responses for the test points.}
#' }
#' @export
skelCubic = function(newnn, newpx, px,trainY, skeleton){
  numknots = nrow(skeleton$centers)
  nn = skeleton$nn
  #modified data matrix for training data based projected values
  trainZ = matrix(0, nrow = length(trainY), ncol = (3*numknots))
  for (i in 1:nrow(nn)) {
    if(is.na(px[i])){#where the data point is projected to nearest knot
      trainZ[i,nn[i,1]] = 1
    }else{
      k0 = min(nn[i,]) #index j as in formula
      k1 = max(nn[i,]) #index l as in formula
      trainZ[i, k0] = 1 - 10*px[i]^3 + 15*px[i]^4 -6*px[i]^5
      trainZ[i, k1] = 10*px[i]^3 -15*px[i]^4 + 6 *px[i]^5
      trainZ[i, (k0+numknots)] = px[i] - 6*px[i]^3 + 8*px[i]^4 - 3*px[i]^5
      trainZ[i, (k1+numknots)] = -4*px[i]^3 + 7*px[i]^4 - 3*px[i]^5
      trainZ[i,(k0+ 2*numknots)] = 0.5*px[i]^2 - 1.5*px[i]^3 + 1.5*px[i]^4 - 0.5*px[i]^5
      trainZ[i,(k1+ 2*numknots)] = 0.5*px[i]^3 - px[i]^4 + 0.5*px[i]^5
    }
  }

  #modified data matrix for test data based projected values
  testZ = matrix(0, nrow = length(newpx), ncol = (3*numknots))
  for (i in 1:nrow(newnn)) {
    if(is.na(newpx[i])){#where the data point is projected to nearest knot
      testZ[i,newnn[i,1]] = 1
    }else{
      k0 = min(newnn[i,])
      k1 = max(newnn[i,])

      testZ[i, k0] = 1 - 10*newpx[i]^3 + 15*newpx[i]^4 -6*newpx[i]^5
      testZ[i, k1] = 10*newpx[i]^3 -15*newpx[i]^4 + 6 *newpx[i]^5
      testZ[i, (k0+numknots)] = newpx[i] - 6*newpx[i]^3 + 8*newpx[i]^4 - 3*newpx[i]^5
      testZ[i, (k1+numknots)] = -4*newpx[i]^3 + 7*newpx[i]^4 - 3*newpx[i]^5
      testZ[i,(k0+ 2*numknots)] = 0.5*newpx[i]^2 - 1.5*newpx[i]^3 + 1.5*newpx[i]^4 - 0.5*newpx[i]^5
      testZ[i,(k1+ 2*numknots)] = 0.5*newpx[i]^3 - newpx[i]^4 + 0.5*newpx[i]^5
    }
  }

  colnames(trainZ) = c(paste0("k", 1:numknots),paste0("g", 1:numknots), paste0("h", 1:numknots))
  colnames(testZ) = c(paste0("k", 1:numknots),paste0("g", 1:numknots), paste0("h", 1:numknots))

  ZtrainDat = as.data.frame(cbind(trainY, trainZ))
  testZ = as.data.frame(testZ)
  lmod3 = stats::lm(trainY ~.-1, data = ZtrainDat) #no intercept

  lmpred3 = stats::predict(lmod3, newdata = testZ)

  return(list(trainZ = trainZ, testZ = testZ, model = lmod3, pred = lmpred3))
}


#' Run all the skeleton-based methods with cross-validation on a dataset
#'
#' @param X0 a matrix/dataframe of the covariates.
#' @param Y0 a vector of the responses.
#' @param seed a random seed for skeleton construction.
#' @param nfold the number of folds for cross validation
#' @param hrate_seq a vector of scaling factors for bandwidth used for kernel regression.
#' @param k_seq a vector of the number of neighbors used for S-knn regressor.
#' @param numknots the number of knots to use for skeleton representation.
#' @param k_cut the number of disconnected components of the resulting skeleton graph.
#' @param rep number of repetitions for k-Means during skeleton construction process.
#' @return A list with following components:
#' \describe{
#'   \item{SEE}{A list recording the Sum of Squared Errors of the tested methods.}
#'   \item{fits}{A list of vectors recording the predicted values of the data points by different methods.}
#' }
#' @export
#' @examples
#' \dontrun{
#' d = 100 #dimension of covariates
#' set.seed(1234)
#' traindata = Yinyang_reg_data(d = d) #simulate the Yinyang Regression data example
#' X0  = as.matrix(traindata$data)
#' Y0  = traindata$Y
#' results = skelReg_runAllCV(X0,Y0, seed = NULL, nfold = 5,
#'                            hrate_seq = c(1,2),
#'                            k_seq = c(3,6,9,12),
#'                            numknots=NULL,
#'                            k_cut=5,
#'                            rep = 1000)
#' }
skelReg_runAllCV = function(X0,Y0, seed = NULL, nfold = 5,
                            hrate_seq = c(1/8, 1/4,1/2,1,2,4,6, 8,10,12, 14,16),
                            k_seq = seq(3,36, by=3),
                            numknots=NULL,
                            k_cut=1,
                            rep = 1000
){
  #split the data for cross validation
  if(!is.null(seed)){set.seed(seed)}
  n = length(Y0)
  fsplit = c(rep( floor(n/nfold), nfold-1), n - floor(n/nfold)*(nfold-1))
  folds = split(1:nrow(X0), sample(rep(1:nfold, fsplit)))
  #choices for number of knots
  ntrain = nrow(X0)*(nfold-1)/nfold
  if(is.null(numknots)){ numknots= round(sqrt(ntrain))}
  ####################
  #Allocate space for saving results
  ####################
  skelResults = matrix(nrow = (nfold*length(hrate_seq)), ncol = 3)#skeleton kernel regression with adaptive bandwidth
  skelResults2 = matrix(nrow = (nfold*length(hrate_seq)), ncol = 3)#skeleton kernel regression with global bandwidth
  skelnnResults = matrix(nrow = (nfold*length(k_seq)), ncol = 3)#skeleton kNN regressor
  lsplineResults = matrix(nrow = nfold, ncol = 2) #skeleton linear spline
  qsplineResults = matrix(nrow = nfold, ncol = 2) #skeleton quadratic spline
  csplineResults = matrix(nrow = nfold, ncol = 2) #skeleton cubic spline
  skelpred = list()
  for (j in 1:length(hrate_seq)) {
    skelpred[[j]] = rep(NA,n)
  }
  skelpred2 = list()
  for (j in 1:length(hrate_seq)) {
    skelpred2[[j]] = rep(NA,n)
  }
  skelnnpred = list()
  for (j in 1:length(k_seq)) {
    skelnnpred[[j]] = rep(NA,n)
  }
  lsplinefit = rep(NA,n)
  qsplinefit = rep(NA,n)
  csplinefit = rep(NA,n)
  ################################################################################
  #Perform cross validation
  ################################################################################
  for(ifold in 1:nfold){
    testIndexes = folds[[ifold]]
    testX = X0[testIndexes, ,drop = FALSE]
    trainX = X0[-testIndexes, ]
    testY = Y0[testIndexes,drop = FALSE]
    trainY = Y0[-testIndexes]
    ######
    #construct skeleton
    #####
    voronSkel = voronSkeleton(trainX,numknots = numknots, k_cut = k_cut, rep=rep, seed=seed)
    centers = voronSkel$centers
    labels = voronSkel$labels
    kdists = voronSkel$kdists
    edgedists = voronSkel$edgedists
    g = voronSkel$g
    nn = voronSkel$nn
    #####
    #calculate the projection for training data point
    ######
    px = skelProject(nn = nn, X = trainX, skeleton = voronSkel)
    #####
    #situate test observations on the graph
    #####
    newnn = RANN::nn2(centers, testX, k=2)$nn.idx
    #calculate the projection for new data points
    newpx = skelProject(nn = newnn, X = testX, skeleton = voronSkel)
    #graph distances between test observations and training points
    testSkeldists = matrix(Inf,nrow = nrow(testX), ncol = nrow(trainX))
    for (i in 1:nrow(testX)) {
      #calculate graph distances between new obs and sample points
      for (j in 1:nrow(trainX)) {
        if(length(intersect(newnn[i,], nn[j,]))>0){ #only calculate distance between nearby data points
          testSkeldists[i,j] = dskeleton(nn[j,], newnn[i,], px[j], newpx[i], voronSkel)
        }
      }
    }
    #######
    #S-kernel regression with global bandwidth
    #######
    #calculate the based bandwidth for cross validation
    centerband = skelBand(centernn = NULL, px = px, skeleton = voronSkel)
    testfit_tmp2 = numeric(nrow(testX))
    for (ihrate in 1:length(hrate_seq)) {
      hrate_tmp = hrate_seq[ihrate]
      h1 = centerband*hrate_tmp
      #regression on test data
      testfit_tmp2 = skelKernel(h1, testX,testSkeldists, trainY)
      skelpred2[[ihrate]][testIndexes] = testfit_tmp2
      skelResults2[(ihrate-1)*nfold+ifold,1] = sum((testfit_tmp2 - testY)^2)
      skelResults2[(ihrate-1)*nfold+ifold,2] = hrate_tmp
      skelResults2[(ihrate-1)*nfold+ifold,3] = ifold

    }
    #######
    #S-kernel regression with varying bandwidth
    #######
    testfit_tmp = numeric(nrow(testX))
    for (ihrate in 1:length(hrate_seq)) {
      hrate_tmp = hrate_seq[ihrate]
      #regression on test data
      testfit_tmp = skelKerneladopt(hrate_tmp, testX,testSkeldists, trainY)
      skelpred[[ihrate]][testIndexes] = testfit_tmp
      skelResults[(ihrate-1)*nfold+ifold,1] = sum((testfit_tmp - testY)^2)
      skelResults[(ihrate-1)*nfold+ifold,2] = hrate_tmp
      skelResults[(ihrate-1)*nfold+ifold,3] = ifold
    }
    ######
    #S-kNN regressor, using graph distance for kNN
    ######
    gnntestfit_tmp = numeric(nrow(testX))
    for (ik in 1:length(k_seq)) {
      k_tmp = k_seq[ik]
      #fit on test data
      gnntestfit_tmp = skelknn(k_tmp, testX,testSkeldists, trainY)
      skelnnpred[[ik]][testIndexes] = gnntestfit_tmp
      skelnnResults[(ik-1)*nfold+ifold,1] = sum((gnntestfit_tmp - testY)^2)
      skelnnResults[(ik-1)*nfold+ifold,2] = k_tmp
      skelnnResults[(ik-1)*nfold+ifold,3] = ifold
    }
    ######
    #linear spline on graph with closed-form solution
    ######
    slinearmodel = skelLinear(newnn, newpx,  px, trainY, voronSkel)
    trainZ = slinearmodel$trainZ  #modified data matrix for training data from projection
    testZ = slinearmodel$testZ   #modified data matrix for test data from projection
    lmod = slinearmodel$model #the model
    lmpred = slinearmodel$pred #the predicted values
    lsplinefit[testIndexes] =  lmpred
    lsplineResults[ifold,1] = sum((testY -lmpred)^2)
    lsplineResults[ifold,2] = ifold
    ######
    #quadratic spline on graph with closed-form solution
    ######
    squadraticmodel = skelQuadratic(newnn, newpx,  px, trainY, voronSkel)
    trainZ = squadraticmodel$trainZ
    testZ = squadraticmodel$testZ
    lmod2 = squadraticmodel$model
    lmpred2 = squadraticmodel$pred
    qsplinefit[testIndexes] =  lmpred2
    qsplineResults[ifold,1] = sum((testY -lmpred2)^2)
    ######
    #cubic spline on general graph with closed-form solution (2p+1 polynomial)
    ######
    scubicmodel = skelCubic(newnn, newpx,  px, trainY, voronSkel)
    trainZ = scubicmodel$trainZ
    testZ = scubicmodel$testZ
    lmod3 = scubicmodel$model
    lmpred3 = scubicmodel$pred
    csplinefit[testIndexes] =  lmpred3
    csplineResults[ifold,1] = sum((testY -lmpred3)^2)
    csplineResults[ifold,2] = ifold
  }### end cross validation
  #gather up results over folds
  `%>%` <- magrittr::`%>%`
  skelSSE = as.data.frame(skelResults) %>% dplyr::group_by(V2) %>% dplyr::summarise(skelSSE = sum(V1))
  names(skelSSE)[1] = "bandrate"
  skelSSE2 = as.data.frame(skelResults2) %>% dplyr::group_by(V2) %>% dplyr::summarise(skelSSE2 = sum(V1))
  names(skelSSE2)[1] = "bandrate"
  skelknnSSE = as.data.frame(skelnnResults) %>% dplyr::group_by(V2) %>% dplyr::summarise(skelknnSSE = sum(V1))
  names(skelknnSSE)[1] = "nneighbor"
  SSEresults = list(skelSSE= as.data.frame(skelSSE),
                    skelSSE2= as.data.frame(skelSSE2),
                    skelknnSSE=as.data.frame(skelknnSSE),
                    lsplineSSE = sum(lsplineResults[,1]),
                    qsplineSSE = sum(qsplineResults[,1]),
                    csplineSSE = sum(csplineResults[,1])
  )
  fits = list(skelpred= skelpred,
              skelpred2= skelpred2,
              skelnnpred= skelnnpred,
              lsplinefit= lsplinefit,
              qsplinefit= qsplinefit,
              csplinefit= csplinefit
  )
  return(list(SSE = SSEresults, fits = fits))
}


#' Run all the skeleton-based methods with training and test datasets
#'
#' @param trainX a matrix/dataframe of the training covariates.
#' @param trainY a vector of the responses on training data.
#' @param testX a matrix/dataframe of the test covariates.
#' @param seed a random seed for skeleton construction.
#' @param hrate_seq a vector of scaling factors for bandwidth used for kernel regression.
#' @param k_seq a vector of the number of neighbors used for S-knn regressor.
#' @param numknots the number of knots to use for skeleton representation.
#' @param k_cut the number of disconnected components of the resulting skeleton graph.
#' @param rep number of repetitions for k-Means during skeleton construction process.
#' @return A list with following components:
#' \describe{
#'   \item{fits}{A list of vectors recording the predicted values of the data points by different methods.}
#' }
#' @export
#' @examples
#' \dontrun{
#' fits = skelReg_runAll(trainX,trainY,testX, seed = NULL,
#'                       hrate_seq = c(1,2),
#'                       k_seq = seq(3,6, by=3),
#'                       numknots=NULL,
#'                       k_cut=5,
#'                       rep = 1000)
#' }
skelReg_runAll = function(trainX,trainY,testX, seed = NULL,
                          hrate_seq = c(1/8, 1/4,1/2,1,2,4,6, 8,10,12, 14,16),
                          k_seq = seq(3,36, by=3),
                          numknots=NULL,
                          k_cut=1,
                          rep = 1000
){
  #split the data for cross validation
  if(!is.null(seed)){set.seed(seed)}
  n = nrow(testX)
  #choices for number of knots
  ntrain = nrow(trainX)
  if(is.null(numknots)){ numknots= round(sqrt(ntrain))}
  ####################
  #Allocate space for saving results
  ####################
  skelpred = list()
  for (j in 1:length(hrate_seq)) {
    skelpred[[j]] = rep(NA,n)
  }
  skelpred2 = list()
  for (j in 1:length(hrate_seq)) {
    skelpred2[[j]] = rep(NA,n)
  }
  skelnnpred = list()
  for (j in 1:length(k_seq)) {
    skelnnpred[[j]] = rep(NA,n)
  }
  lsplinefit = rep(NA,n)
  qsplinefit = rep(NA,n)
  csplinefit = rep(NA,n)
  ################################################################################
  #Perform fitting
  ################################################################################
  ######
  #construct skeleton
  #####
  voronSkel = voronSkeleton(trainX,numknots = numknots, k_cut = k_cut, rep=rep, seed=seed)
  centers = voronSkel$centers
  labels = voronSkel$labels
  kdists = voronSkel$kdists
  edgedists = voronSkel$edgedists
  g = voronSkel$g
  nn = voronSkel$nn
  #####
  #calculate the projection for training data point
  ######
  px = skelProject(nn = nn, X = trainX, skeleton = voronSkel)
  #####
  #situate test observations on the graph
  #####
  newnn = RANN::nn2(centers, testX, k=2)$nn.idx
  #calculate the projection for new data points
  newpx = skelProject(nn = newnn, X = testX, skeleton = voronSkel)
  #graph distances between test observations and training points
  testSkeldists = matrix(Inf,nrow = nrow(testX), ncol = nrow(trainX))
  for (i in 1:nrow(testX)) {
    #calculate graph distances between new obs and sample points
    for (j in 1:nrow(trainX)) {
      if(length(intersect(newnn[i,], nn[j,]))>0){ #only calculate distance between nearby data points
        testSkeldists[i,j] = dskeleton(nn[j,], newnn[i,], px[j], newpx[i], voronSkel)
      }
    }
  }
  #######
  #S-kernel regression with global bandwidth
  #######
  #calculate the based bandwidth for cross validation
  centerband = skelBand(centernn = NULL, px = px, skeleton = voronSkel)
  testfit_tmp2 = numeric(nrow(testX))
  for (ihrate in 1:length(hrate_seq)) {
    hrate_tmp = hrate_seq[ihrate]
    h1 = centerband*hrate_tmp
    #regression on test data
    testfit_tmp2 = skelKernel(h1, testX,testSkeldists, trainY)
    skelpred2[[ihrate]]= testfit_tmp2
  }
  #######
  #S-kernel regression with varying bandwidth
  #######
  testfit_tmp = numeric(nrow(testX))
  for (ihrate in 1:length(hrate_seq)) {
    hrate_tmp = hrate_seq[ihrate]
    #regression on test data
    testfit_tmp = skelKerneladopt(hrate_tmp, testX,testSkeldists, trainY)
    skelpred[[ihrate]] = testfit_tmp
  }
  ######
  #S-kNN regressor, using graph distance for kNN
  ######
  gnntestfit_tmp = numeric(nrow(testX))
  for (ik in 1:length(k_seq)) {
    k_tmp = k_seq[ik]
    #fit on test data
    gnntestfit_tmp = skelknn(k_tmp, testX,testSkeldists, trainY)
    skelnnpred[[ik]] = gnntestfit_tmp
  }
  ######
  #linear spline on graph with closed-form solution
  ######
  slinearmodel = skelLinear(newnn, newpx,  px, trainY, voronSkel)
  trainZ = slinearmodel$trainZ  #modified data matrix for training data from projection
  testZ = slinearmodel$testZ   #modified data matrix for test data from projection
  lmod = slinearmodel$model #the model
  lmpred = slinearmodel$pred #the predicted values

  ######
  #quadratic spline on graph with closed-form solution
  ######
  squadraticmodel = skelQuadratic(newnn, newpx,  px, trainY, voronSkel)
  trainZ = squadraticmodel$trainZ
  testZ = squadraticmodel$testZ
  lmod2 = squadraticmodel$model
  lmpred2 = squadraticmodel$pred

  ######
  #cubic spline on general graph with closed-form solution (2p+1 polynomial)
  ######
  scubicmodel = skelCubic(newnn, newpx,  px, trainY, voronSkel)
  trainZ = scubicmodel$trainZ
  testZ = scubicmodel$testZ
  lmod3 = scubicmodel$model
  lmpred3 = scubicmodel$pred

  fits = list(skelpred= skelpred,
              skelpred2= skelpred2,
              skelnnpred= skelnnpred,
              lsplinefit= lmpred,
              qsplinefit= lmpred2,
              csplinefit= lmpred3
  )
  return(fits)
}
