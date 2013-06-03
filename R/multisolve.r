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
