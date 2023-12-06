#' write landmarks to json file readable by NORA software
#'
#' write landmarks to json file readable by NORA software
#' @param x  \code{k x 3} matrix containing landmarks
#' @param filename character: output file name
#' @param labels character vector of length \code{k} containing landmark names
#' @param size integer: voxel size around landmark
#' @export
#' @importFrom jsonlite write_json
write.ano.json <- function(x,filename=dataname,labels=dataname,size=1) {
    dataname <- deparse(substitute(x))
    if (!grepl("*.json$", filename)) 
        filename <- paste0(filename,".ano.json")
    nrx <- nrow(x)

    if (labels[[1]] == dataname)
        mylabels <- paste0(dataname,"-",1:nrx)
    else if (length(labels) == nrx)
        mylabels <- labels
    
       
   
    ## setup markups
  
    position <- lapply(1:nrx, function(y) y <- x[y,] )
    poschk <- which(as.logical(sapply(lapply(position,is.na),sum)))
    if (length(poschk))
        position[poschk] <- ""
    cp <- data.frame(name=mylabels,color=rainbow(nrx),size=as.integer(size))
    cp$coords <- position
    

    markups <- data.frame(type="pointset", name="detections")
    markups$state <- list(NULL)
    markups$points=list(cp)

    
    out <- list("annotations"=markups)
    
    #out$points <- markups
    write_json(out,pretty=T,auto_unbox = T,filename,digits=NA,always_decimal=TRUE,null="list")
    
}
