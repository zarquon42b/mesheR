#' create an object of class BayesDeform to be used in gaussMatch and AmbergRegister
#'
#' create an object of class BayesDeform to be used in gaussMatch and AmbergRegister
#'
#' @param model a model of class pPCA created with package RvtkStatismo
#' @param sdmax numeric or vector of numerics (length must be of value =< \code{iterations} as to be specified in gaussMatch or AmbergRegister. If \code{length(sd) =n < iterations}, then the first n iterations will be restricted to the according standard deviations
#' @param mahalanobis parameter according to mahaprob in function PredictSample from RvtkStatismo
#' @param ptValueNoise parameter according to ptValueNoise in function statismoConstrainModel from RvtkStatismo
#' @param wt numeric: a weight put on the closest model configuration. E.g. if wt = 2 the deformation will be 0.33*closestPoints + 0.66*most probable model shape. This can either be a constant scalar (that might be modified in each iteration by a shrinkage function - see below) or a vector of length sdmax.
#' @param align logical: if the target mesh is already aligned to the model, align can be set to FALSE.
#' @param shrinkfun a two-valued function for the generation of the weight in the i-th iteration. See examples below.
#' @return returns a list of class "BayesDeform"
#' @examples
#' \dontrun{
#' shrinkfun <- function(x,i){x <- x*0.9^i; return(x)}
#' # for wt=0.5 this will produce 
#' # 
#' }
#' @export 
createBayes <- function(model, sdmax=numeric(0),mahalanobis="chisq",ptValueNoise=2,wt=NULL,align=TRUE,shrinkfun=NULL) {
    Bayes <- list(); class(Bayes) <- "BayesDeform"
    Bayes$model <- model
    Bayes$sdmax <- sdmax
    Bayes$mahaprob <- mahalanobis
    Bayes$ptValueNoise <- ptValueNoise
    if (length(wt) == 1) {
        Bayes$wt <- rep(wt,length(sdmax))
        if (!is.null(shrinkfun)) {
            for (i in 2:length(sdmax))
                Bayes$wt[i] <- shrinkfun(Bayes$wt[i],i)
        }
    } else if ( length(wt) != length(sdmax)) {
        stop("wt must be of same length as sdmax, or single value")
    } else
        Bayes$wt <- wt
            
    Bayes$align <- align
    return(Bayes)
}
