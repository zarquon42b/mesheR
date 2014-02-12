#' Place Dowels representing tissue depth on a skull
#'
#' Place Dowels representing tissue depth on a skull by using predefined tissue-values
#'
#' @param lm k x 3landmarks representing dowel locations
#' @param mesh triangular mesh (object of class 'mesh3d')
#' @param ldowel vector of length k giving a length for each dowel
#' @param render logical: if TRUE, the result will be rendered in an rgl window
#' @param col color of the dowels
#' @param radius numeric: diameter of dowels
#' @param meshcol mesh color
#' @param fine integer: amount of vertices to generate the cylinders representing tissue thickness.
#' @param smooth logical: use smoothed normals for dowel orientation
#' @param dowelcol specify dowelcolor.
#'
#' @details For orientation of the dowels, the (angle weighted) normal vectors of the surface is used. 
#' @return an invisible list with
#' \item{dowels }{a list containing the meshes for all dowels}
#' \item{endpoints }{a matrix containing the coordinates of the dowels endpoints}
#' 
#' @examples
#' require(Rvcg)
#' data(humface)
#' lms <- matrix(c(17.6061 , 9.2072 , -6.9917 , 44.7959 , 37.1135 , 76.0469 , 3.7734 , 12.8234 , 81.3138 , 63.4207 , 47.5009 , 44.6468), 4, 3)
#' #place markers of 5mm length
#' placeDowels(lms, humface, rep(5, 4))
#' @seealso \code{\link{cylinder}}
#' @importFrom rgl shade3d
#' @export placeDowels
placeDowels <- function(lm, mesh, ldowel, render=TRUE,col=1,radius=1,meshcol=3, fine=50, smooth=TRUE, dowelcol=NULL)
    {

        ldo <- FALSE
        colvec <- FALSE
        projLM <- vcgClost(lm,mesh,smoothNormals=smooth)
                                        ##projLM <- closemeshKD(lm,mesh)
        #projLM <- projRead(lm,mesh,smooth=TRUE)
        if (length(ldowel) > 1)
            ldo <- TRUE
                  
        if (length(col) > 1)
            colvec <- TRUE

        dowels <- list()
        endpoints <- matrix(NA, dim(lm)[1] ,3)
        endpoints <- vert2points(projLM)+t(projLM$normals[1:3,])*ldowel
### create dowels and render them if required
        for (i in 1:dim(lm)[1])
            {
                if (ldo)
                    ltmp <-ldowel[i]
                else
                    ltmp=ldowel

                dowels[[i]] <- cylinder(projLM$vb[1:3,i],projLM$normals[1:3,i],length=ltmp,radius=radius,fine=fine,adNormals=FALSE)
                
                if (!is.null(dowelcol))
                    dowels[[i]]$material$color <- matrix(dowelcol,dim(dowels[[i]]$it)[1],dim(dowels[[i]]$it)[2])
                    
                if (render)
                    {
                        if (colvec)
                            coltmp <- col[i]
                        else
                            coltmp <- col

                        shade3d(dowels[[i]],col=coltmp)
                    }
            }
        if (render) ## render mesh ##
            shade3d(mesh,col=meshcol)

        invisible(list(dowels=dowels, endpoints=endpoints))
    }
