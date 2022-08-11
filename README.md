# skeletonMethods

R package for Skeleton-based Analysis Methods

[![LifeCycle Status](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://deploy-preview-81--tidyverse-org.netlify.com/lifecycle/#maturing)


`skeletonMethods` is [R package](https://www.r-project.org/package) implementing the **Skeleton Clustering** and **Skeleton Clustering**  methods in : 
> Zeyu Wei, and Yen-Chi Chen. "Skeleton Clustering: Dimension-Free Density-based Clustering." 2021.
> Zeyu Wei, and Yen-Chi Chen. "Skeleton Regression." 2021.

The manuscript of the paper can be found at [here](https://arxiv.org/abs/2104.10770).


## Installation 
The package currently only exists on github. Installed directly from this repo with help from the [devtools](https://github.com/hadley/devtools) package. Run the following in R:

	```R
	install.packages("devtools")
	devtools::install_github("JerryBubble/skeletonClus") 
	```
	
#### Development note 

The package is actively developing.


## Description of Clustering Framework
The follwoing figure illustrates the overall procedure of the skeleton clustering method.
<div style = "text-align:center" align="center"> <img src="https://jerrybubble.github.io/images/skeletonWorkFlow.jpg" width = "660"/> </div>

Starting with a collection of observations (panel (a)),
we first find knots, the representative points of the entire data (panel (b)). By default we use the k-means algorithm to choose knots. 
Then we compute the corresponding Voronoi cells induced by the knots (panel (c))
and the edges associating the Voronoi cells (panel (d), this is the Delaunay triangulation). 
For each edge in the graph, we compute a density-based similarity measure that quantifies the closeness of each pair of knots.
For the next step we segment knots into groups based on a linkage criterion (single linkage in this example), leading to the dendrogram in panel (e). 
Finally, we choose a threshold that cuts the dendrogram into $S = 2$ clusters (panel (f))
and assign cluster label to each observation according to the knot-cluster that it belongs to (panel (g)).



## Code Example for Clustering

Load R-package and simulate Yinyang data of dimension 100
```R
library("skeletonClus")

dat = Yinyang_data(d = 100)
X0 = dat$data
y0 = dat$clus
```
Other simulation data used in the skeleton clustering paper can be genereted similarly: `Mickey_data()` for Mickey data, `MM_data()` for Manifold Mixture data, `MixMickey_data()` for Mix Mickey data, `MixStar_data()` for mixture of three Gaussian cluster giving a star shape, and `ring_data()` for Ring data

Also, with `set.seed(1234)` the 1000-dimensional simulated datasets are also built into the skeletonClus package, with names `Yinyang1000, Ring1000, Mickey1000, MixMickey1000, MixStar1000, ManifoldMix1000`. Those built in datasets can be load into R working environment by
```R
data(Yinyang1000)
```

We proceed with constructing the weighted skeleton with overfitting k-means, wihch can be donw with the `skeletonCons` function:
```R
skeleton = skeletonCons(X0, rep = 1000)
```

For more flexible usage of skeletonCons see the help page at
```R
?skeletonCons
```
The knots and the Voronoi density edge weights can be accessed as
```R
skeleton$centers
skeleton$voron_weights
```
The Face density, Tube density, and Average Distance density weights can be accessed in similar fashion. We then segment the knots based on the edge weights with hierarchical clustering (taking Voronoi density for example):

```R
##turn similarity into distance
VD_tmp = max(skeleton$voron_weights) - skeleton$voron_weights 
diag(VD_tmp) = 0
VD_tmp = as.dist(VD_tmp)

##perform hierarchical clustering on the distance matrix
VD_hclust = hclust(VD_tmp, method="single")
```
If the final number of clusters is unknown, we can plot the dendrogram of hierarchical clustering to see the structure of the data and pick a cutting threshold. Here we know the data have 5 components and hence we cut the tree into 5 groups.
```R
##plot dendrogram
plot(VD_hclust)
##cut tree 
VD_lab = cutree(VD_hclust, k=5)

```

Then we assign the label of each observation to have the same group membership as the closest knot

```R
X_lab_VD = VD_lab[skeleton$cluster]
```
With the __mclust__ R package we can calculate the adjusted Rand index:
```R
library("mclust")
voron_rand = adjustedRandIndex(X_lab_VD, y0)
```
FInally we can get a plot to visualize the clustering result:

```R
 plot(X0[,1], X0[,2], col = cols[X_lab_VD], main = "Skeleton Clustering with Voronoi Density", xlab = "X", ylab = "Y", cex.main=1.5)
```


## Description of Regression Framework
Like in the skeleton clustering framework, for regression puporse we also first construct a skeleton graph to represent the covariates. Particularly, for regression we use the density-aided similarity measure called Voronoi density to help cut skeleton if needed, and such construction can be done with the function `voronSkeleton`. Then we project the data covariates onto the skeleton graph with function `skelProject` and compute the skeleton-based distances between them with function `dskeleton`. Then we implement different non-parametric regression techiniques on the skeleton graph:
### Skeleton Kernel Regression
Let $K_h(.) = K(./h)$ be a non-negative kernel function with bandwidth $h > 0$, let $d_{\mathcal S} and let $s_1, \dots, s_n$ denote the skeleton-projected covariates. the corresponding skeleton-based kernel (S-kernel) regressor for a point $s$ is
$$
    \hat{m}(s) = \frac{\sum_{j=1}^N K_h(d_{\mathcal S}(s_j, s)) Y_j}{\sum_{j=1}^N K_h(d_{\mathcal S}(s_j, s))} 
\end{align}
An example kernel function  is the Gaussian kernel that
\begin{align}
    K_h(d_{\mathcal S}(s_j, s_\ell)) = \exp\left(- \frac{d_{\mathcal S}(s_j, s_\ell)^2}{h^2}\right)
$$
and the function `skelKernel` fits this regressor.