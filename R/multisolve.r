#' Fortran routine for fast inversion of multiple 3x3 matrices
#' 
#' Fortran routine for fast inversion of multiple 3x3 matrices
#' 
#' 
#' Be m=k*3. The matrix will then be interpreted as k 3x3 matrices, which are
#' subsequently inverted
#' 
#' @param A m x 3 matrix, where m is a multiple of 3 %% ~~Describe \code{A}
#' here~~
#' @param trans calculate inverse of transposed subblocks %% ~~Describe
#' \code{trans} here~~
#' @return Returns matrix containing inverted submatrices.
#' @author Stefan Schlager
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' A <- matrix(rnorm(18),6,3)
#' B <- multisolve3(A)
#' @export multisolve3
multisolve3 <- function(A,trans=FALSE)
{
    A <- as.matrix(A)
    n <- ncol(A)
    m <- nrow(A)
    l <- (m/n)-1
    if (m %% n !=0)
        stop("number of rows must be multiple of columns")
    storage.mode(A) <- "double"
    storage.mode(m) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(l) <- "integer"
    trans <- as.integer(trans)
    out <- .Fortran("multisolve3",A,m,l,trans)
    return(out[[1]])
}

testsolve3 <- function(A)
    {
        A <- as.matrix(A)
         n <- ncol(A)
         storage.mode(n) <- "integer"
         storage.mode(A) <- "double"
         B <- A
         out <- .Fortran("M33INV",A,A,TRUE)
         return(out)
     }
multisolve <- function(A,trans=FALSE)
{
    A <- as.matrix(A)
    n <- ncol(A)
    m <- nrow(A)
    l <- (m/n)-1
    if (m %% n !=0)
        stop("number of rows must be multiple of columns")
    storage.mode(A) <- "double"
    storage.mode(m) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(l) <- "integer"
    trans <- as.integer(trans)
    out <- .Fortran("multisolve",A,m,n,l,trans)
    return(out[[1]])
}
