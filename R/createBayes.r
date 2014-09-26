#' create an object of class BayesDeform to be used in gaussMatch and AmbergRegister
#'
#' create an object of class BayesDeform to be used in gaussMatch and AmbergRegister
#'
#' @param model a model of class pPCA created with package RvtkStatismo
#' @param sdmax numeric or vector of numerics (length must be of value =< \code{iterations} as to be specified in gaussMatch or AmbergRegister. If \code{length(sd) =n < iterations}, then the first n iterations will be restricted to the according standard deviations
#' @param mahalanobis parameter according to mahaprob in function PredictSample from RvtkStatismo
#' @param ptValueNoise parameter according to ptValueNoise in function statismoConstrainModel from RvtkStatismo
#' @param wt weight put on the closest model configuration. E.g. if wt = 2 the deformation will be 0.33*closestPoints + 0.66*most probable model shape.
#' @param align logical: if the target mesh is already aligned to the model, align can be set to FALSE.
#' @return returns a list of class "BayesDeform"
#' @export 
createBayes <- function(model, sdmax=numeric(0),mahalanobis="chisq",ptValueNoise=2,wt=NULL,align=TRUE) {
    Bayes <- list(); class(Bayes) <- "BayesDeform"
    Bayes$model <- model
    Bayes$sdmax <- sdmax
    Bayes$mahaprob <- mahalanobis
    Bayes$ptValueNoise <- ptValueNoise
    Bayes$wt <- wt
    Bayes$align <- align
    return(Bayes)
}
