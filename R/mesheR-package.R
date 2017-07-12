#' Meshing tools
#' 
#' A toolbox for meshing operations: elastic mesh mapping, Iteratively closest point matching, vertex selection, ...,
#' 
#' \tabular{ll}{
#' Package: \tab mesheR\cr
#' Type: \tab Package\cr
#' Version: \tab 0.4.160301\cr
#' Date: \tab 2016-03-01\cr
#' License: \tab GPL\cr
#' LazyLoad: \tab yes\cr }
#' 
#' @name mesheR-package
#' @aliases mesheR-package mesheR
#' @docType package
#' @author Stefan Schlager
#' 
#' Maintainer: Stefan Schlager <zarquon42@@gmail.com>
#' @references To be announced
#' @keywords package
#' @useDynLib mesheR, .registration=TRUE
#' @importFrom Morpho angle.calc closemeshKD invertFaces facenormals file2mesh mcNNindex meshcube meshDist meshres projRead rmUnrefVertex rmVertex rotmesh.onto rotonmat rotonto tangentPlane tps3d unrefVertex vert2points tps3d
#' @importFrom Rvcg vcgBorder vcgClost vcgGetEdge vcgRaySearch vcgNonBorderEdge vcgSmooth vcgClostKD vcgCreateKDtreeFromBarycenters vcgClostOnKDtreeFromBarycenters
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom graphics plot
#' @importFrom grDevices col2rgb colorRampPalette extendrange rgb
#' @importFrom methods as
#' @importFrom stats optim pnorm prcomp qchisq quantile

NULL
