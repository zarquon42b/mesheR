colorTrans <- function(mesh1,mesh2,tol=1)
    {
        if (is.null(mesh2$material$color))
            stop("no color to transfer")
        vn <- ncol(mesh2$vb)
        clost <- closemeshKD(mesh1,mesh2)
        col <-  rep("#FFFFFF", vn)
        tmp1 <- data.frame(it = as.vector(mesh2$it))
        tmp1$rgb <- as.vector(mesh2$material$color)
        tmp1 <- unique(tmp1)
        col[tmp1$it] <- tmp1$rgb
        colout <- hex2RGB(col,gamma=T)
        hits <- which(clost$quality <= tol)
        if (length(hits) == 0)
            stop("no hits within the given tolerance - increase tol!")
        colmesh1 <- matrix(0,ncol(mesh1$vb),4)
        colmesh1[hits,1] <- clost$ptr[hits]
        colmesh1[hits,2:4] <- t(mesh2$it[,clost$ptr[hits]])
        vert1col <- as(colout[colmesh1[hits,2],],"LAB")
        vert2col <- as(colout[colmesh1[hits,3],],"LAB")
        vert3col <- as(colout[colmesh1[hits,4],],"LAB")

        #vert1col <- as(vert1col,"LAB")
        #vert2col <- as(vert2col,"LAB")
        #vert3col <- as(vert3col,"LAB")
        outcol <- vert1col
        outcol@coords <- (vert1col@coords+vert2col@coords+vert3col@coords)/3

        goodcoll <- as(outcol,"RGB")
        blackit <- rep(max(outcol@coords),ncol(mesh1$vb))
        colmat <- RGB(blackit,blackit,blackit)
        #colmat <- matrix(1,ncol(mesh1$vb),3)
        colmat@coords[hits,] <- goodcoll@coords
        #colmat <- rgb(colmat[, 1], colmat[, 2], colmat[,3], maxColorValue = max(colmat))
        colmat <- hex(colmat)
        colfun <- function(x) {
            x <- colmat[x]
            return(x)
        }
        material <- list()
        material$color <- matrix(colfun(mesh1$it), dim(mesh1$it))
        mesh1$material$color <- material$color
        return(list(mesh=mesh1,clost=clost,colmat=colmat,outcol=outcol))
    }
