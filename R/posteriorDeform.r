## find correspondences
getCorrespondences <- function(mesh,targetmesh,distance,silent=TRUE,slide=ifelse(bending,3,10),bending=TRUE,partsample=partsample,ray=TRUE,tol=pi/5,k=200,meanmod,modlm=NULL, tarlm=NULL,forceLM=FALSE) {
    myslide <- NULL
    customLM <- FALSE
    parttofixed <- vcgClostKD(transferPoints(partsample,meanmod,mesh,tolwarn = 5),mesh,k=50)
    if (ray)
        part2raw    <- vcgRaySearch(parttofixed,targetmesh,mindist=T)
    else {
        
        part2raw <- vcgClostKD(parttofixed,targetmesh,angdev = tol,k=k,sign=FALSE)
        part2raw$distance <- part2raw$quality
        part2raw$quality <- rep(1,length(part2raw$quality))
    }
    goodclost   <- which(as.logical((normcheck(parttofixed,part2raw) < tol) * (abs(part2raw$distance) < distance)*(part2raw$quality==1)))
    
    referencepoints <- vert2points(parttofixed)[goodclost,]
    targetpoints <- vert2points(part2raw)[goodclost,]
    nref <- nrow(referencepoints)
    stepsize=1
    #iterations=3
    if (!bending) {
        stepsize=0.1
        #iterations=10
    }
    if (!is.null(modlm) && !is.null(tarlm)) {
        customLM <- TRUE
        ## print("using landmarks")
        modlm <- transferPoints(modlm,meanmod,mesh,tolwarn = 5)
        nlm <- nrow(modlm)
        npos <- (1:nlm)+nref
        referencepoints <- rbind(referencepoints,modlm)
        targetpoints <- rbind(targetpoints,tarlm)
    }
    if (slide > 0) {
        if (forceLM) {
            SMvector <- surp <- 1:nref
        } else {
            SMvector <- surp <- 1:nrow(referencepoints)
        }
                
                myslide    <- relaxLM(referencepoints,targetpoints,mesh=mesh,iterations = slide,SMvector = SMvector,surp=surp,silent=silent,bending=bending,stepsize =stepsize ,tol=0)
    }
        
    return(list(goodclost=goodclost,myslide=myslide,targetpoints=targetpoints))
}


#' Fit an SSM to a target based on subsampling corresponding points
#'
#'  Fit an SSM to a target based on subsampling corresponding points and compute the posterior mean
#'
#' @param model statismo shape model
#' @param target target mesh
#' @param reference model instance other than the mean (class mesh3d)
#' @param partsample predetermined corresponding points on the sample mean
#' @param samplenum integer: if partsample=NULL, this specifies the number of coordinates sampled on the model mean
#' @param distance numeric: constrain maximum distance to mark target point as appropriate
#' @param slide integer: if > 0 the valid correspondences on the model instance will be relaxed minimizing bending energy/procrustes distance.
#' @param bending logical: if TRUE, the coordinates on the model instance are relaxed using bending energy, Procrustes distance otherwise
#' @param ray logical: if TRUE, the closest point search will be performed along the normals only
#' @param deform logical if TRUE, the posterior mean will also be deformed to the target using an elastic deformation
#' @param Amberg if TRUE the deformation will use the function \code{\link{AmbergDeformSpam}} and \code{\link{tps3d}} otherwise
#' @param rhotol maximal tolerated angle between normals to be considered a valid match
#' @param modlm matrix containing 3D landmarks on the model mean (not for alignment)
#' @param tarlm  matrix containing 3D landmarks on the target surface
#' @param align2mod logical: if TRUE, the prediction step will perform an alignment to the model using the valid correspondences.
#' @param alignbymesh logical: if TRUE, the alignment to the SSM will be computed by the entire mesh, if FALSE only the valid correspondences are used.
#' @param forceLM if TRUE, predfined landmarks \code{modlm} are not allowed to slide. For cases with high uncertainty, this can lead to unwanted mesh distortions.
#' @param silent logical: supress debug output
#' @param threads integer: number of threads to use for tps interpolation (set to 1 if using openblas, or otherwise it can become instable)
#' @param ... additional parameters passed to \code{\link{AmbergDeformSpam}}.
#' @return returns a deformed version of a model instance fitted to the target
#' @note Please note that it is required to align the target mesh to the model mean beforehand. This can be performed using the function \code{\link{icp}}, for example.
#' @examples
#' \dontrun{
#' require(RvtkStatismo)
#'download.file(url="https://github.com/marcelluethi/statismo-shaperegistration/raw/master/data/VSD001_femur.vtk","./VSD001_femur.vtk",method = "w")
#' download.file(url="https://github.com/marcelluethi/statismo-shaperegistration/raw/master/data/VSD002_femur.vtk","./VSD002_femur.vtk",method = "w")
#' download.file(url="https://github.com/marcelluethi/statismo-shaperegistration/raw/master/data/VSD001-lm.csv","./VSD001-lm.csv",method = "w")
#' download.file(url="https://github.com/marcelluethi/statismo-shaperegistration/raw/master/data/VSD002-lm.csv","./VSD002-lm.csv",method = "w")
#' ref <- read.vtk("VSD001_femur.vtk")
#' tar <- read.vtk("VSD002_femur.vtk")
#' ref.lm <- as.matrix(read.csv("VSD001-lm.csv",row.names=1,header = FALSE))
#' tar.lm <- as.matrix(read.csv("VSD002-lm.csv",row.names=1,header = FALSE))
#' Kernels <- SumKernels(GaussianKernel(50,50),IsoKernel(0.1,ref))
#' mymod <- statismoModelFromRepresenter(ref,kernel=Kernels,ncomp=100)
#' postDef <- posteriorDeform(mymod,tar,modlm=ref.lm,tarlm = tar.lm,samplenum = 1000)
#' ## run a loop redoing that step using the result of the previous step as input
#' for (i in 1:5)
#'    postDef <- posteriorDeform(mymod,tar,modlm=ref.lm,tarlm = tar.lm,samplenum = 1000,reference=postDef)
#'
#' ## now we leave the model space for a final deform involving a TPS deform
#' postDefFinal <- postDef
#' for (i in 1:3)
#'     postDefFinal <- posteriorDeform(mymod,tar,modlm=ref.lm,tarlm = tar.lm,samplenum = 3000,reference=postDefFinal,deform=T,distance=3)
#' 
#' Morpho::meshDist(postDefFinal,tar,from=-2,to=2,tol=.5)
#' rgl::wire3d(tar,col="white")
#' }
#' @importFrom Rvcg vcgSample
#' @importFrom Morpho relaxLM
#' @export
posteriorDeform <- function(model,target,reference=NULL,partsample=NULL,samplenum=1000,distance=1e10,slide=3,bending=TRUE,ray=FALSE,deform=FALSE, Amberg=FALSE,rhotol=pi/2,modlm=NULL,tarlm=NULL,align2mod=TRUE,alignbymesh=FALSE,forceLM=FALSE,silent=FALSE,threads=1,...) {
      if (!requireNamespace("RvtkStatismo"))
          stop("for using the option Bayes, please install RvtkStatismo from https://github.com/zarquon42b/RvtkStatismo")
    meanmod <- RvtkStatismo::DrawMean(model)
    if (is.null(reference))
        reference <- meanmod
## cat("Deformation step without curvature\n")
    if (Amberg)
        deformfun <- function(mesh,lm1,lm2,...) { return( AmbergDeformSpam(mesh,lm1,lm2,...)$mesh)}
    else
        deformfun <- function(mesh,lm1,lm2,...)  {return(tps3d(mesh,lm1,lm2,threads=threads))}
    if (is.null(partsample))
        partsample <- vcgSample(meanmod,samplenum)
    corrs   <- getCorrespondences(reference,target,distance,silent=silent,bending=bending,partsample = partsample,ray=ray,meanmod = meanmod,tol = rhotol,modlm = modlm,tarlm = tarlm,forceLM=forceLM)
    myslide <- corrs$myslide
    targetpoints <- corrs$targetpoints
    back2mod <- transferPoints(myslide,reference,meanmod,tolwarn = 5)
    if (!alignbymesh)
        deformed <- RvtkStatismo::PredictSample(model,lmDataset = targetpoints,lmModel=back2mod,sdmax=7,mahaprob="dist",align=align2mod)
    else {
        if (align2mod)
            trafo <- computeTransform(meanmod,reference)
        else
            trafo <- diag(4)
        targetpoints2mod <- applyTransform(targetpoints,trafo)
        deformed <- applyTransform(RvtkStatismo::PredictSample(model,lmDataset = targetpoints2mod,lmModel=back2mod,sdmax=7,mahaprob="dist",align=FALSE),trafo,inverse = TRUE)        
        
        
    }
      if (deform) {
        myslide <- transferPoints(myslide,reference,deformed,tolwarn = 5)
        deformed <- deformfun(deformed,myslide,targetpoints,...)
    }
    count <- 0
    return(deformed)
    
}

#' Deforms a reference to a target based on a TPS or AmbergDeform
#'
#' Deforms a (pre-aligned) reference to a target based on a TPS/AmbergDeform and automatically sampled sliding semi-landmarks
#'
#' @param reference reference mesh
#' @param target target mesh
#' @param partsample predetermined sample points on the reference mesh
#' @param samplenum integer: if partsample=NULL, this specifies the number of coordinates sampled on the model mean
#' @param distance numeric: constrain maximum distance to mark target point as appropriate
#' @param slide integer: if > 0 the valid correspondences on the model instance will be relaxed minimizing bending energy/procrustes distance.
#' @param bending logical: if TRUE, the coordinates on the model instance are relaxed using bending energy, Procrustes distance otherwise
#' @param ray logical: if TRUE, the closest point search will be performed along the normals only
#' @param Amberg if TRUE the deformation will use the function \code{\link{AmbergDeformSpam}} and \code{\link{tps3d}} otherwise
#' @param rhotol maximal tolerated angle between normals to be considered a valid match
#' @param reflm matrix containing 3D landmarks on the reference
#' @param tarlm  matrix containing 3D landmarks on the target surface
#' @param forceLM if TRUE, predfined landmarks \code{modlm} are not allowed to slide. For cases with high uncertainty, this can lead to unwanted mesh distortions.
#' @param silent logical: supress debug output
#' @param threads integer: number of threads to use for tps interpolation (set to 1 if using openblas, or otherwise it can become instable)
#' @param ... additional parameters passed to \code{\link{AmbergDeformSpam}}.
#' @return returns a deformed version of a model instance fitted to the target
#' @note Please note that it is required to align the target mesh to the reference mesh beforehand. This can be performed using the function \code{\link{icp}}, for example.
#' @export 
subsampleDeform <- function(reference,target,partsample=NULL,samplenum=1000,distance=1e10,slide=3,bending=TRUE,ray=FALSE, Amberg=FALSE,rhotol=pi/2,reflm=NULL,tarlm=NULL,forceLM=FALSE,silent=FALSE,threads=1,...) {
    if (Amberg)
        deformfun <- function(mesh,lm1,lm2,...) { return( AmbergDeformSpam(mesh,lm1,lm2,...)$mesh)}
    else
        deformfun <- function(mesh,lm1,lm2,...)  {return(tps3d(mesh,lm1,lm2,threads=threads))}
    if (is.null(partsample))
        partsample <- vcgSample(reference,samplenum)
    corrs   <- getCorrespondences(reference,target,distance,silent=silent,bending=bending,partsample = partsample,ray=ray,meanmod = reference,tol = rhotol,modlm = reflm,tarlm = tarlm,forceLM=forceLM)
    myslide <- corrs$myslide
    targetpoints <- corrs$targetpoints  
    message("Relaxing landmarks\n")
    deformed <- deformfun(reference,myslide,targetpoints,...)
    count <- 0
    return(deformed)
    

    
}
