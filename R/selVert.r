selectVertex <- function(mesh,col=3,visible=TRUE)
  {
    open3d()
    wire3d(mesh, col = col, specular="black")
    selcheck <- 0
    run <- 0
     cat("select a region using the right mouse button\n")
            while (run == 0) {
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
                      
                  #selcheck <- length(selected)
                  while (selcheck == 0) {
                    view <- points3d(subset[selected,], col = 2, cex = 2)
                    answer <- readline("do you like the view? (y/add/remove)\n")
                    if (answer == "y") {
                      selcheck <- 1
                      run <- 1
                      #rgl.pop("shapes", id = view)
                    }
                    if (substr(answer,1L,1L) == "a") {
                      #selcheck <- 1
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
                     if (substr(answer,1L,1L) == "r") {
                     # selcheck <- 1
                      f <- select3d("right")
                      tmpsel <- which(f(subset))
                       if(visible)
                        {
                          visi <- which(glVisible(mesh))
                          visi <- visi[which(visi%in%selected)]
                        }
                      remov <- which(selected%in%which(f(subset)))
                     
                      if (length(remov) >0 )

                        selected <- selected[-remov]
                      rgl.pop("shapes", id = view)
                      if (length(selected) > 0)
                        view <- points3d(subset[selected,], col = 2, cex = 2)
                    }
                  }
                  
                }
            }
    invisible(selected)
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
    

  
