#' skeleton distance calculation
#' , , , , ,
#' @param g an igraph for skeleton knots and edges.
#' @param nnc1 a vector of the indeces of the two nearest neighbors of one point.
#' @param nnc2 a vector of the indeces of the two nearest neighbors of another point.
#' @param px1 the projection proportion of point 1 onto the edge by its two nearest knots, with the lower index knot as the origin.
#' @param px2 the projection proportion of point 2 onto the edge by its two nearest knots, with the lower index knot as the origin.
#' @param edgedists a matrix of Euclidean distances between knots, which are used as edge distances.
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
dskeleton = function(g, nnc1, nnc2, px1, px2, edgedists) {

  #simplify calculation by only focusing on pair of points sharing at least one knot
  if(length(intersect(nnc1,nnc2))<1){ return(Inf)}

  if(is.na(px1) & is.na(px2)){
    #case when both data points cannot be projected onto an edge
    return(distances(g, v=nnc1[1], to =nnc2[1]))

  }else if(!is.na(px1) & is.na(px2) ){
    #case when point id1 can be projected onto an edge
    btnodes = c(distances(g, v=nnc1[1], to =nnc2[1]),
                distances(g, v=nnc1[2], to =nnc2[1]))
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
    btnodes = c(distances(g, v=nnc1[1], to =nnc2[1]),
                distances(g, v=nnc1[1], to =nnc2[2]))
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

    btnodes = c(distances(g, v=nnc1[1], to =nnc2[1]),
                distances(g, v=nnc1[1], to =nnc2[2]),
                distances(g, v=nnc1[2], to =nnc2[1]),
                distances(g, v=nnc1[2], to =nnc2[2]))
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
