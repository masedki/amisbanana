
## Classe S4 pour une population pondérées amis à revoir si package R !!!!!!!!!!!!!!
amispop <- setClass(
  # Set the name for the class
  Class = "amispop",
  # Define the slots
  representation = representation(x = "matrix",w = "numeric", denom="numeric", numera = "numeric", d = "numeric",n = "numeric"),
  # Set the default values for the slots. (optional)
  prototype=prototype(x = matrix(0,0,0),w = numeric(), denom=numeric() , numera=numeric(), d = numeric(), n = numeric()),
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity=function(object)
  {
    if(sum(object@w < 0)>0){
      stop("A negative weight")
    }
    return(TRUE)
  }
)

amispop <- function(x)
{
return(new("amispop", x=x, w=rep(1, nrow(x)), d=ncol(x), n=nrow(x)))
}






setGeneric(name= "samplefrom",  def = function(obj){ standardGeneric("samplefrom")})
setMethod(f = "samplefrom",
           signature(obj="amispop"),
           definition = function(obj){
             #idx = sample.int(obj@n, size = obj@n, prob = as.vector(obj@w) , rep=TRUE)
             #obj1 <- new("amispop", x=matrix(0,0,0), w=rep(1,obj@n), denom=as.vector(obj@denom), numera=as.vector(obj@numera), d=ncol(obj@x), n=nrow(obj@x))
             idx = sample.int(obj@n, size = 50000, prob = as.vector(obj@w) , rep=TRUE)
             obj1 <- new("amispop", x=matrix(0,0,0), w=rep(1,50000), denom=as.vector(obj@denom), numera=as.vector(obj@numera), d=ncol(obj@x), n=50000)
             obj1@x <- obj@x[idx,]
             return(obj1)
           })

setGeneric(name= "initsample",  def = function(obj, n, s, b){ standardGeneric("initsample")})
setMethod(f = "initsample", signature(obj="amispop", n="numeric", s="numeric", b="numeric"),
            definition = function(obj, n, s, b){
            obj@x <- matrix(NA, n, length(s))
            obj@d <- ncol(obj@x)
            obj@n <- nrow(obj@x)
            for(j in 1:obj@d)
              obj@x[,j] <- rlogis(n, loc = 0, scale = s[j])

            obj@numera  <- as.vector(dbanana(obj@x, b))
            obj@denom <- as.vector(dlogistic(obj@x, s))
            obj@w <- obj@numera/obj@denom
            obj@w <- obj@w/sum(obj@w)
            return(obj)
          })


########################################################################################################################
## !!!!!!!!!!!!!!!! Classe S4 amisstrategy provisoire, à remvoir si transformation en package !!!!!!!!!!!!!!!!!!!!!!!!!!
########################################################################################################################
setClass(
  Class = "amisstrategy",
  representation = representation(nu="numeric", g="numeric", nbSmall="numeric", iterSmall="numeric", nbKeep="numeric", iterKeep="numeric", tolKeep="numeric"),
  prototype = prototype(nu=numeric(),g=numeric(), nbSmall=numeric(), iterSmall=numeric(), nbKeep=numeric(), iterKeep=numeric(), tolKeep=numeric())
)

## Constructeur de la classe S4 amisstrategy
amisstrategy <- function(nu=3,g=6, nbSmall=150, iterSmall=50, nbKeep=25, iterKeep=10**3, tolKeep=10**(-3)){
  if( nbKeep > nbSmall)
    nbKeep <- nbSmall
  new("amisstrategy",nu=nu,g=g,nbSmall=nbSmall, iterSmall=iterSmall, nbKeep=nbKeep,iterKeep=iterKeep, tolKeep=tolKeep)
}

########################################################################################################################
## !!!!!!!!!!!!!!!! Classe S4 amisparameter provisoire, à remvoir si transformation en package !!!!!!!!!!!!!!!!!!!!!!!!!
########################################################################################################################
amisparameter <- setClass(
  # Set the name for the class
  Class = "amisparameter",
  # Define the slots
  representation = representation(mu = "matrix", pi = "numeric", sigma = "array", nu = "integer"),
  # Set the default values for the slots. (optional)
  prototype=prototype(mu = matrix(0,0,0), pi = numeric(), sigma = array(data = NA, dim = length(data), dimnames = NULL), nu=integer())
)

########################################################################################################################
## !!!!!!!!!!!!!!!! Classe S4 amiscriteria provisoire, à remvoir si transformation en package !!!!!!!!!!!!!!!!!!!!!!!!!!
########################################################################################################################
amiscriteria <- setClass(
  # Set the name for the class
  Class = "amiscriteria",
  # Define the slots
  representation = representation(loglik = "numeric", degeneracyrate = "numeric"),
  # Set the default values for the slots. (optional)
  prototype=prototype(loglik = numeric(), degeneracyrate= numeric())
)

########################################################################################################################
## Classe S4 amisXEM qui contient tout le reste
########################################################################################################################
setClass(
  Class = "amisXEM",
  representation = representation(pop="amispop", criteria="amiscriteria", strategy="amisstrategy", param="amisparameter"),
  prototype = prototype(pop=new("amispop"), criteria=new("amiscriteria"), strategy=new("amisstrategy"), param=new("amisparameter"))
)

########################################################################################################################
## Classe S4 amisSystem qui contient tout le reste
########################################################################################################################
setClass(
  Class = "amissystem",
  representation = representation(pop="amispop", param="amisparameter"),
  prototype = prototype(pop=new("amispop"), param=new("amisparameter"))
)

generatesystem <- function(obj,n,b)
{
  x<-matrix(NA, n, nrow(obj@param@mu))
  idx <- sample(obj@strategy@g, size=n,prob = obj@param@pi, rep=TRUE)
  for(k in 1:obj@strategy@g)
     x[which(idx==k),]  <- rmt(n=length(which(idx==k)), mean=obj@param@mu[,k], S=obj@param@sigma[,,k], df=obj@param@nu)
  pop <- new("amispop", x=x, w=numeric(), denom=numeric(), numera=numeric(), d=ncol(x), n=n)
  mysystem <- new("amissystem", pop=pop, param=obj@param)
  computedensities(mysystem, b)
  return(mysystem)
}
