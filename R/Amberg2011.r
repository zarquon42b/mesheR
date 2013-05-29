edgeSelectAmberg <- function(mesh)
    {
        verts <- vert2points(mesh)
        norms <- t(facenormals(mesh)$normals[1:3,])##calculate face normals
        allcoo <- rbind(verts,norms)##combine vertices and face normals
        rowind <- (1:ncol(mesh$it))
        rowind3 <- rowind*3
        rowind2 <- rowind3-1
        rowind1 <- rowind3-2
        whichsel <- cbind(rowind1,rowind2,rowind3,t(mesh$it),nrow(verts)+rowind)
        sel <- Matrix(0,3*ncol(mesh$it),nrow(allcoo))
        sel[whichsel[,c(1,4)]] <- -1
        sel[whichsel[,c(1,5)]] <- 1
        sel[whichsel[,c(2,4)]] <- -1
        sel[whichsel[,c(2,6)]] <- 1
        sel[whichsel[,c(3,7)]] <- 1
        return(list(sel=sel,allcoo=allcoo))
    }

mesh2sparse2 <- function (mat) 
{
    nvb <- nrow(mat)
    j <- matrix(1:(nvb),3,nvb/3)
    j1 <- cbind(j,j,j)*0
    j1[,c(1:ncol(j))*3] <- j
    j1[,(c(1:ncol(j))*3)-1] <- j
    j1[,(c(1:ncol(j))*3)-2] <- j
    j1 <- as.vector(j1)
    i1 <- as.vector(rbind(1:(nvb),1:(nvb),1:(nvb))[,1:nvb])
    out <- sparseMatrix(j = j1, i = i1, x = as.vector(t(mat)))
    invisible(t(out))
}
invertBD <- function(bd,cores=detectCores())
{
    matfun <- function(x,bd)
        {
            out <- t(solve(as.matrix(bd[((x*3)+1):((x*3)+3),((x*3)+1):((x*3)+3)])))
            return(out)
        }
    ## for(i in 0:(ncol(bd)/3-1)){matlist[[i+1]] <- as.matrix(bd[((i*3)+1):((i*3)+3),((i*3)+1):((i*3)+3)])}
    la <- as.list(0:(ncol(bd)/3-1))
    solvit <- mclapply(la,matfun,bd,mc.cores=cores)
    ##solvit <- mclapply(matlist,solve,mc.cores = cores)
    ##solvit <- mclapply(solvit,t,mc.cores = cores)
    solvemat <- matrix(unlist(solvit),ncol(bd),3,byrow = T)
    solvemat <- mesh2sparse2(solvemat)
    return(solvemat)
}

