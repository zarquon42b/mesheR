selectVertex <- function(mesh,col=3,visible=TRUE)
  {
    open3d()
    wire3d(mesh, col = col, specular="black")
    selcheck <- 0
    run <- 0
    cat("select a region using the right mouse button\n")
    while (run == 0) { #initial selection
      rgl.bringtotop(stay = FALSE)
      if (interactive()) {
        f <- select3d("right")
        subset <- t(mesh$vb[1:3, ])
        tmpsel <- which(f(subset))
        if(visible)
          {
            visi <- which(glVisible(mesh))
            tmpsel <- visi[which(visi%in%tmpsel)]
          }
        selected <- tmpsel
        iter=0   
        while (selcheck == 0) {
          if (iter == 0)
            view <- points3d(subset[selected,], col = 2, cex = 2)
          answer <- readline("do more? (q/a/r/i/s)(yes|add|remove|invert|switch view)\n")
          if (answer == "q") {
            selcheck <- 1
            run <- 1
          }
          if (answer == "s") {
            visible <- !visible
          }
          if (substr(answer,1L,1L) == "a") {
            rgl.bringtotop(stay = FALSE)
            f <- select3d("right")
            tmpsel <- which(f(subset))
            if(visible)
              {
                visi <- which(glVisible(mesh))
                tmpsel <- visi[which(visi%in%tmpsel)]
              }
            selected <- unique(c(selected,tmpsel))
            rgl.pop("shapes", id = view)
            view <- points3d(subset[selected,], col = 2, cex = 2)
          }
          if (substr(answer,1L,1L) == "i") {
            rgl.bringtotop(stay = FALSE)
            #f <- select3d("right")
            
            selected <- (1:ncol(mesh$vb))[-selected]
            rgl.pop("shapes", id = view)
            view <- points3d(subset[selected,], col = 2, cex = 2)
          }
          if (substr(answer,1L,1L) == "r") {
            rgl.bringtotop(stay = FALSE)
            f <- select3d("right")
            tmpsel <- which(f(subset))
            if(visible)
              {
                visi <- which(glVisible(mesh))
                tmpsel <- visi[which(visi%in%tmpsel)]
              }
            remov <- which(selected%in%tmpsel)
            rgl.pop("shapes", id = view)
            if (length(remov) > 0 )
              selected <- selected[-remov]
            view <- points3d(subset[selected,], col = 2, cex = 2)
          }
          iter <- iter+1
        }
      }
    }
    return(selected)
  }

cutMesh <- function(mesh,visible=TRUE,invert=TRUE,col=3)
  {
    removal <- selectVertex(mesh,col=col,visible=visible)
    vb <- 1:ncol(mesh$vb)

    if (invert)
      vb <- vb[-removal]
    else
      vb <- removal
    outmesh <- rmVertex(mesh,vb)
    rgl.clear()
    wire3d(outmesh,col=col,specular="black")
    cat(paste(length(vb)," vertices removed\n"))
    invisible(outmesh)
    
  }
    

  
