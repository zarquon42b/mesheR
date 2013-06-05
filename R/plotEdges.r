plotEdges <- function(mesh,col=2,lwd=2)
    {
        edges <- vcgGetEdge(mesh)
        bord <- which(edges$border == 1)
        bordmesh <- mesh
        bordmesh$it <- rbind(edges$vert1[bord],edges$vert2[bord],edges$vert1[bord])
        wire3d(bordmesh,col=col,lwd=lwd,lit=F)
    }
