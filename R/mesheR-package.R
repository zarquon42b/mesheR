#' Meshing tools
#' 
#' A toolbox for meshing operations: elastic mesh mapping, Iteratively closest point matching, vertex selection, ...,
#' 
#' \tabular{ll}{
#' Package: \tab mesheR\cr
#' Type: \tab Package\cr
#' Version: \tab 0.4-00\cr
#' Date: \tab 2013-07-19\cr
#' License: \tab GPL\cr
#' LazyLoad: \tab yes\cr }
#' 
#' @name mesheR-package
#' @aliases mesheR-package mesheR
#' @docType package
#' @author Stefan Schlager
#' 
#' Maintainer: Stefan Schlager <zarquon42@gmail.com>
#' @references To be announced
#' @keywords package
#' @useDynLib mesheR
#' @importClassesFrom spam spam spam.chol.NgPeyton
#' @importFrom Morpho vert2points closemeshKD tanplan tps3d conv2backf facenormals adnormals meshDist rmUnrefVertex rotmesh.onto file2mesh meshres warp.mesh unrefVertex rotonto angle.calc rotonmat meshcube rmVertex projRead mcNNindex
#' @importFrom Rvcg vcgNonBorderEdge vcgClost vcgBorder vcgSmooth vcgIntersect vcgGetEdge
#' @importFrom parallel detectCores mclapply
#' @importFrom rgl translate3d

NULL
