buildAffineAmberg <- function(mesh)
    {
        verts <- vert2points(mesh)
        norms <- t(facenormals(mesh)$normals[1:3,])##calculate face normals
        allcoo <- rbind(verts,norms)##combine vertices and face normals
        rowind <- (1:ncol(mesh$it))
        rowind3 <- rowind*3
        rowind2 <- rowind3-1
        rowind1 <- rowind3-2
        whichsel <- cbind(rowind1,rowind2,rowind3,t(mesh$it),nrow(verts)+rowind)
        sel <- Matrix::Matrix(0,3*ncol(mesh$it),nrow(allcoo))
        sel[whichsel[,c(1,4)]] <- -1
        sel[whichsel[,c(1,5)]] <- 1
        sel[whichsel[,c(2,4)]] <- -1
        sel[whichsel[,c(2,6)]] <- 1
        sel[whichsel[,c(3,7)]] <- 1
        return(list(sel=sel,allcoo=allcoo))
    }

mesh2sparse2 <- function (mat) 
{
    mat <- as.matrix(mat)
    nvb <- nrow(mat)
    data <- t(mat)
    data <- as.vector(data)
    j <- matrix(1:(nvb),3,nvb/3)
    j1 <- cbind(j,j,j)*0
    j1[,c(1:ncol(j))*3] <- j
    j1[,(c(1:ncol(j))*3)-1] <- j
    j1[,(c(1:ncol(j))*3)-2] <- j
    j1 <- as.vector(j1)
    i1 <- as.vector(rbind(1:(nvb),1:(nvb),1:(nvb))[,1:nvb])
    out <- Matrix::sparseMatrix(j = j1, i = i1, x = data)
    invisible((out))
}

createArcNode <- function(mesh)
    {
        edges <- vcgNonBorderEdge(mesh)
        n <- dim(edges)[1]
        dimit <- ncol(mesh$it)
        arcnode <- Matrix(0,n,dimit)
        edges$ind <- 1:n
        arcnode[cbind(edges$ind,edges$face1)] <- -1
        arcnode[cbind(edges$ind,edges$face2)] <- 1
        invisible(arcnode)
    }

createS <- function(mesh)
    {
        out <- list()
        sel <- buildAffineAmberg(mesh)
        transel <- sel$sel %*% sel$allcoo
        invtransel <- multisolve3(transel,trans=FALSE)
        out$sparseAijk <- mesh2sparse2(invtransel)
        out$sel <- sel
        dia3 <- sparseMatrix(i=1:3,j=1:3,x=rep(1,3))
        arcnodes <- createArcNode(mesh) 
        arcnodes <- kronecker(arcnodes,dia3)
        out$process <- arcnodes%*%out$sparseAijk
        return(out)
    }


