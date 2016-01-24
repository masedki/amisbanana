amisbanana <- function(nsec=seq(4000, 13000, len=10), update=FALSE,b=0.03, g=6, nu=3, d=2, n0=3000,s=c(1.73*14.14/pi, 17.3/pi, rep(1.73/pi, (d-2))))
{
  if(update)
  {
    s4vec <- lapply( rep("amisparameter", length(nsec)), new)
    strategy <- amisstrategy(nu=nu,g=g, nbSmall=500)
    #strategy <- amisstrategy(nu=nu,g=g)
    if(d>2)
    {
      pop <- initsample(new("amispop"), n=100000, s=s,b=b)
      #strategy <- amisstrategy(nu=nu,g=g, nbSmall=500)
    }
    else
    {
      pop <- initsample(new("amispop"), n=n0, s=s,b=b)
      #strategy <- amisstrategy(nu=nu,g=g)
    }
    poptmp <- samplefrom(pop)
    param <- new("amisparameter")
    xem <- new("amisXEM", pop=poptmp, criteria=new("amiscriteria"), strategy=strategy, param=param)
    runXEM(xem)
    system <- generatesystem(xem, nsec[1],b)
    s4vec[1] <-xem@param
    pop <-system@pop
    computedmixture(pop, s4vec[[1]], nsec[1])
    pop@w <- as.vector(pop@numera/pop@denom)
    for(t in 2:length(nsec))
    {
      poptmp <- samplefrom(pop)
      xem <- new("amisXEM", pop=poptmp, criteria=new("amiscriteria"), strategy=strategy, param=new("amisparameter"))
      runXEM(xem)
      system <- generatesystem(xem, nsec[t],b)
      s4vec[t] <-xem@param
      pop@x<- rbind(pop@x, system@pop@x)
      pop@numera <- c(pop@numera, system@pop@numera)
      pop@denom <- rep(0, sum(nsec[1:t]))
      pop@n <- sum(nsec[1:t])
      for(h in 1:t)
        computedmixture(pop, s4vec[[h]], nsec[h])
      pop@w <- as.vector(pop@numera/pop@denom)
    }
    return(list(pop=pop, param=s4vec))
  }
  else
  {
    s4vec <- lapply( rep("amissystem", length(nsec)), new)
    pop <- initsample(new("amispop"), n=n0, s=s,b=b)
    strategy <- amisstrategy(nu=nu,g=g)
    xem <- new("amisXEM", pop=pop, criteria=new("amiscriteria"), strategy=strategy, param=new("amisparameter"))
    runXEM(xem)

    for(t in 1:length(nsec))
    {
      system <- generatesystem(xem, nsec[t],b)
      computedensities(system, b)
      s4vec[t] <- system
      if(t < length(nsec))
      {
        pop <- samplefrom(system@pop)
        xem <- new("amisXEM", pop=pop, criteria=new("amiscriteria"), strategy=strategy, param=new("amisparameter"))
        runXEM(xem)
      }
    }

    return(s4vec)
  }
}
