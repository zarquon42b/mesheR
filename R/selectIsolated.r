#' interactively select isolated pieces from a mesh
#'
#' interactively select isolated pieces from a mesh
#' @param mesh triangular mesh of class mesh3d
#' @param maxpiece integer: the n-largest (number of vertices) pieces to consider
#' @importFrom Rvcg vcgIsolated
#' @return returns a mesh with all selected parts merged and colored
#' @export
selectIsolated <- function(mesh,maxpiece=10) {
    pieces <- vcgIsolated(mesh,split=TRUE)
    ll <- length(pieces)
    nverts <- sapply(pieces,function(x) x <- ncol(x$vb))
    pieces <- pieces[order(nverts,decreasing = T)]
    pieces <- pieces[1:min(maxpiece,ll)]
    ll <- length(pieces)
    cols <- colorRampPalette(c("red","green","blue"))
    cols <- cols(ll)
    answerargs <- c("yes","no","back","stop")
    cat(paste0("there are ",ll," pieces\n"))
    open3d()
    out <- list()
    answerList <- list()
    if (interactive()) {
        i <- 1
        while (i < ll) {
            print(i)
            wire3d(pieces[[i]],col=cols[i])
            answer <- list();class(answer) <- "try-error"
            while(inherits(answer,"try-error"))
                answer <- try(match.arg(tolower(readline("select (yes/no/back/stop - abbreviations allowed)\n")),answerargs),silent = TRUE)
            answerList <- append(answerList, answer)
            if (answer == "yes") {
                out <- append(out,list(colorMesh(pieces[[i]],cols[i])))
            }
            else if (answer == "no")
                rgl.pop()
            else if (answer == "stop") {
                rgl.pop()
                break
            }
            else if ((answer == "back") && i > 1) {

                if (answerList[[i-1]] == "yes"){
                    out <- out[1:length(out)-1]
                    rgl.pop()
                    rgl.pop()
                }
                else if (answerList[[i-1]] %in%("no")) {
                    rgl.pop()
                }
                answerList <- answerList[1:(length(answerList)-2)]
                i <- i - 2
            }
            else {
                out <- out[1:length(out)-1]
                if ((i > 1) && answerList[[i - 1]] == "yes")
                   rgl.pop()
                answerList <- answerList[1:(length(answerList)-1)]
                i <- i - 1
            }
            i <- i + 1
        }
    if (length(out) > 1)
        out <- mergeMeshes(out)
    }
    return(out)
}
