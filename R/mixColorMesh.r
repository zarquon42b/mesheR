mixColorMesh <- function(mesh1,mesh2,alpha=0.5)
    {
       col1 <- mesh1$material$color
       col2 <- mesh2$material$color

       nas1 <- which(is.na(col1))
       nas2 <- which(is.na(col2))
       if (length(nas1) > 0)
            col1[nas1] <- "#FFFFFF"
        if (length(nas2) > 0)
            col1[nas2] <- "#FFFFFF"

       mx <- hex(mixcolor(alpha=alpha,hex2RGB(col1),hex2RGB(col2)))
       mesh1$material$color <- matrix(mx,dim(mesh1$material$color))
       return(mesh1)
   }
       
