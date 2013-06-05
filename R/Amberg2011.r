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
        return(list(sel=sel,allcoo=allcoo,verts=verts,normals=norms))
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
        out <- list()
        out$arcnode <- arcnode
        out$edges <- edges
        invisible(out)
    }
AijklWeights <- function(mesh,arcnode)
    {
        verts <- vert2points(mesh)
        bary <-t( facenormals(mesh)$vb[1:3,])
        xijk <- bary[arcnode$edges$face1,]-bary[arcnode$edges$face2,]
        xijk <- apply(xijk,1,function(x){x <- sqrt(sum(x^2))})
        nedge <- length(xijk)
        arraw <- cbind(verts[arcnode$edges$vert1,],verts[arcnode$edges$vert2,],bary[arcnode$edges$face1,],bary[arcnode$edges$face2,])

        aijkl <- sqrt(areafun(arraw))
        diagmat <- aijkl/xijk
        diatmp <- sparseMatrix(i=1:nedge,j=1:nedge,x=diagmat)
        dia3 <- sparseMatrix(i=1:3,j=1:3,x=rep(1,3))
        diatmp <- kronecker(diatmp,dia3)
        return(diatmp)
        
    }
areafun <- function(x)
    {

        m <- nrow(x)
        storage.mode(x) <- "double"
        out <- as.double(rep(0,m))
        out <- .Fortran("areafun",x,m,out)
        return(out[[3]])
    }



createS <- function(mesh)
    {
        out <- list()
        sel <- buildAffineAmberg(mesh)
        transel <- sel$sel %*% sel$allcoo
        invtransel <- multisolve3(transel,trans=TRUE)
        arcnodes <- createArcNode(mesh) 
        weight <- AijklWeights(mesh,arcnodes)
        out$sparseAijk <- Matrix::t(mesh2sparse2(invtransel))
        out$sel <- sel
        dia3 <- sparseMatrix(i=1:3,j=1:3,x=rep(1,3))
        arcnodes3 <- kronecker(arcnodes$arcnode,dia3)
        out$process <- weight%*%arcnodes3%*%out$sparseAijk%*%sel$sel
        
        return(out)
    }


createC <- function(lm1,mesh)
    {
        proj <- vcgClost(lm1,mesh,barycentric=TRUE)
        C <- Matrix(0,nrow(lm1),ncol(mesh$vb))
        vertptr <- t(mesh$it[,proj$faceptr])
        faceptr <- cbind(1:nrow(lm1),vertptr)
        for(i in 2:4)
            C[faceptr[,c(1,i)]] <- proj$barycoords[i-1,]
        return(C)
    }

AmbergDeformSpam <- function(mesh,lm1,lm2,k0=1,lambda=1)
    {
        t0 <- Sys.time()
        out <- list()
        S <- createS(mesh)
        ##convert S to class spam
        spamS <- as.spam.dgCMatrix(S$process)
        fn <- S$sel$normals
        
        ##Jc <- Matrix(0,nrow(lm1),dim(S$process)[2])
        ##Jc[1:nrow(lm1),1:ncol(C)] <- C
        ## setup Jc Jacobian submatrix
        
        Jc <- createJc(lm1,ncol(S$process),mesh)
        Jc <- lambda*Jc
        ## setup Jacobian
        J <- rbind(spamS,Jc)

        ## setup Jn Jacobian submatrix
        Jn0 <- sparseMatrix(i=1:ncol(mesh$it),j=1:ncol(mesh$it),x=rep(1,ncol(mesh$it)))
        Jn <- Matrix(0,ncol(mesh$it),ncol(J)-ncol(Jn0))
                                        #Jn[,-c(1:(ncol(J)-ncol(Jn0)))] <- Jn0
        Jn <- cBind(Jn,Jn0)
        Jn <- k0*Jn
        Jn <- as.spam.dgCMatrix(Jn)
        ##
        t1 <-  Sys.time();print(t1-t0)
        ## convert the matrices to class spam and let spam handle it (much faster)
        #J <- as.spam.dgCMatrix(J)
        #Jc <- as.spam.dgCMatrix(Jc)
       
        ## append Jn to Jacobian
        J <- rbind(J,Jn)
        H <- spam::t(J)%*%J
        Jtc <- spam::t(Jc)%*%lm2
        Hchol <- chol(H)
        print(system.time(k <- solve.spam(Hchol,lambda*Jtc)))
        v <- S$sel$allcoo
        v[-c(1:dim(vert2points(mesh))[1]),] <- 0
        v <- Matrix(v)
        k_v <- k-v
        Jtnv <- spam::t(Jn)%*%(fn)
        print(system.time(deltav <- solve(H,k0*Jtnv)+k_v))
        vert_new <- as.matrix(spam::t((v+deltav)[1:ncol(mesh$vb),]))
        mesh_new <- mesh
        mesh_new$vb[1:3,] <- (vert_new)
        meshout <- adnormals(mesh_new)
        return(list(mesh=meshout,Jn=Jn,Jc=Jc,J=J,H=H,Hchol=Hchol))
    }
createJc <- function(lm1,ncol,mesh)
    {##function to setup Jc Jacobian submatrix
        C <- createC(lm1,mesh)
        Jc <- Matrix(0,nrow(lm1),ncol)
        Jc[1:nrow(lm1),1:ncol(C)] <- C
        Jc <- as.spam.dgCMatrix(Jc)
        return(Jc)
    }

AmbergDeform <- function(mesh,lm1,lm2,k0=1,lambda=1)
    {
        t0 <- Sys.time()
        ##
        S <- createS(mesh)
        fn <- S$sel$normals
        C <- createC(lm1,mesh)
        Jc <- Matrix(0,nrow(lm1),dim(S$process)[2])
        Jc[1:nrow(lm1),1:ncol(C)] <- C
        Jc <- lambda*Jc
        J <- rBind(S$process,Jc)
        Jn0 <- sparseMatrix(i=1:ncol(mesh$it),j=1:ncol(mesh$it),x=rep(1,ncol(mesh$it)))
        Jn <- Matrix(0,ncol(mesh$it),ncol(J))
        Jn[,-c(1:(ncol(J)-ncol(Jn0)))] <- Jn0
        Jn <- k0*Jn
        ##
        t1 <-  Sys.time();print(t1-t0)
        ##
        J <- rBind(J,Jn)

        
        H <- Matrix::t(J)%*%J
        Jtc <- Matrix::t(Jc)%*%lm2
        print(system.time(k <- solve(H,lambda*Jtc)))
        v <- S$sel$allcoo
        v[-c(1:dim(vert2points(mesh))[1]),] <- 0
        v <- Matrix(v)
        k_v <- k-v
        Jtnv <- t(Jn)%*%(fn)
        print(system.time(deltav <- solve(H,k0*Jtnv)+k_v))
        vert_new <- as.matrix(t((v+deltav)[1:ncol(mesh$vb),]))
        mesh_new <- mesh
        mesh_new$vb[1:3,] <- (vert_new)
        meshout <- adnormals(mesh_new)
        return(meshout)
    }
