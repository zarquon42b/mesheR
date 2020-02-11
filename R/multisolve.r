#' routine for fast inversion of multiple 3x3 matrices
#' 
#' routine for fast inversion of multiple 3x3 matrices
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
#' 
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
    trans <- as.integer(trans)
    out <- .Call("multisolve3Cpp",A,trans)
    return(out)
}

