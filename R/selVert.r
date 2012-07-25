selectVertex <- function(mesh,col=3)
  {
    wire3d(mesh, col = col, specular="black")
    selcheck <- 0
    run <- 0
     cat("select a region using the right mouse button\n")
            while (run == 0) {
                rgl.bringtotop(stay = FALSE)
                if (interactive()) {
                  f <- select3d("right")
                  subset <- t(mesh$vb[1:3, ])
                  selected <- which(f(subset))
                  #selcheck <- length(selected)
                  while (selcheck == 0) {
                    view <- points3d(subset[selected,], col = 2, cex = 2)
                    answer <- readline("do you like the view? (y/add/remove)\n")
                    if (answer == "y") {
                      selcheck <- 1
                      run <- 1
                      rgl.pop("shapes", id = view)
                    }
                    if (substr(answer,1L,1L) == "a") {
                      #selcheck <- 1
                      f <- select3d("right")
                      selected <- unique(c(selected,which(f(subset))))
                      rgl.pop("shapes", id = view)
                      view <- points3d(subset[selected,], col = 2, cex = 2)
                    }
                     if (substr(answer,1L,1L) == "r") {
                     # selcheck <- 1
                      f <- select3d("right")
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
