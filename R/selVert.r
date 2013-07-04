#' select and crop triangular surface meshes
#' 
#' select a subset of a surface mesh interactively.
#' 
#' select vertices or trim a triangular mesh by selecting a subset - either
#' restricted to visible from the current point of view or by rectangular
#' region. Visibility is determined by the fact that the vector from the
#' viewpoint to the selected vertex does not intersect the mesh.
#' 
#' @param mesh triangular mesh of class "mesh3d"
#' @param col color to render the mesh. 
#' @param visible select only vertices visible (from the present point of view)
#' @param add logical: add the surface to an existing window.
#' @param render character: how to render the surface. Possible values are
#' "shade" or "wire".
#' @param offset initial offset to move vertex slightly away from the surface.
#' @param \dots additional arguments passed to the rendering functions
#' shade3d and wire3d from package "shapes". 
#' @return selectVertex returns the indices of the selected vertices.
#' 
#' cutMesh returns the trimmed mesh.
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead}}, \code{\link{glVisible}}, \code{\link{cutMesh}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' \dontrun{
#' data(nose)
#' selection <- selectVertex(shortnose.mesh)
#' }
#' 
#' @export selectVertex 
selectVertex <- function(mesh,col=3,visible=TRUE,add=FALSE,render=c("shade","wire"), offset=1e-3, ...)
  {
    visifun <- function()
      {
        visi <- which(glVisible(mesh, offset=offset))
        tmpsel <- visi[which(visi%in%tmpsel)]
        return(tmpsel)
      }
        
    render <- substr(render[1],1L,1L)
    do3d <- wire3d
    if (render == "s")
      do3d <- shade3d
    if (!add)
      open3d()
    do3d(mesh, col = col, specular="black",...)
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
          tmpsel <- visifun()
        selected <- tmpsel
        iter=0   
        while (selcheck == 0) {
          if (iter == 0)
            view <- points3d(subset[selected,1],subset[selected,2],subset[selected,3], col = 2, cex = 2)
          
          if (visible)
            visiquestion <- ("do more? (q/a/r/i/s)(quit|add|remove|invert|switch selection mode)\ncurrent selection mode: visible only\n")
          else
            visiquestion <-("do more? (q/a/r/i/s)(quit|add|remove|invert|switch selection mode)\ncurrent selection mode: all vertices\n")
          answer <- readline(visiquestion)
          
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
               tmpsel <- visifun()
            selected <- unique(c(selected,tmpsel))
            rgl.pop("shapes", id = view)
            view <- points3d(subset[selected,1],subset[selected,2],subset[selected,3], col = 2, cex = 2)
          }
          if (substr(answer,1L,1L) == "i") {
            rgl.bringtotop(stay = FALSE)
            #f <- select3d("right")
            
            selected <- (1:ncol(mesh$vb))[-selected]
            rgl.pop("shapes", id = view)
            view <- points3d(subset[selected,1],subset[selected,2],subset[selected,3], col = 2, cex = 2)
          }
          if (substr(answer,1L,1L) == "r") {
            rgl.bringtotop(stay = FALSE)
            f <- select3d("right")
            tmpsel <- which(f(subset))
            if(visible)
             tmpsel <- visifun()
            remov <- which(selected%in%tmpsel)
            rgl.pop("shapes", id = view)
            if (length(remov) > 0 )
              selected <- selected[-remov]
            view <- points3d(subset[selected,1],subset[selected,2],subset[selected,3],col = 2, cex = 2)
          }
          iter <- iter+1
        }
      }
    }
    return(sort(selected))
  }

#' crop triangular surface meshes
#' 
#' crop a surface mesh interactively.
#' 
#' select vertices or trim a triangular mesh by selecting a subset - either
#' restricted to visible from the current point of view or by rectangular
#' region. Visibility is determined by the fact that the vector from the
#' viewpoint to the selected vertex does not intersect the mesh.
#' 
#' @param mesh triangular mesh of class "mesh3d"
#' @param visible select only vertices visible (from the present point of view)
#' @param keep.selected logical: determines if the selected vertices or their
#' complement are deleted.
#' @param col color to render the mesh.
#' @param add logical: add the surface to an existing window.
#' @param render character: how to render the surface. Possible values are
#' "shade" or "wire".
#' @param offset initial offset to move vertex slightly away from the surface. 1e-3 seems to be a good threshold for objects from the macroscopical world.
#' @param \dots additional arguments passed to the rendering functions
#' shade3d and wire3d from package "rgl".
#' @return selectVertex returns the indices of the selected vertices.
#' 
#' cutMesh returns the trimmed mesh.
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead}},  \code{\link{glVisible}}, \code{\link{selectVertex}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' \dontrun{
#' data(nose)
#' selection <- selectVertex(shortnose.mesh)
#' }
#' @export cutMesh
cutMesh <- function(mesh,visible=TRUE,keep.selected=TRUE,col=3,add=FALSE,render=c("shade","wire"),offset=1e-3,...)
  {
      mesh <- adnormals(mesh)
    render <- substr(render[1],1L,1L)
    do3d <- wire3d
    if (render == "s")
      do3d <- shade3d
    removal <- selectVertex(mesh,col=col,visible=visible,add=add,render=render,offset=offset,...)
    vb <- 1:ncol(mesh$vb)
    
    if (keep.selected)
      vb <- vb[-removal]
    else
      vb <- removal
    outmesh <- rmVertex(mesh,vb)
    rgl.clear()
    do3d(outmesh,col=col,specular="black")
    cat(paste(length(vb)," vertices removed\n"))
    invisible(outmesh)
    
  }
    
