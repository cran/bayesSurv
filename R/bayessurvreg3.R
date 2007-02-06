################################################
#### AUTHOR:     Arnost Komarek             ####
####             (2005)                     ####
####                                        ####
#### FILE:       bayessurvreg3.R            ####
####                                        ####
#### FUNCTIONS:  bayessurvreg3              ####
################################################

### ======================================
### bayessurvreg3
### ======================================
## Survival regression with G-spline random intercept
## and G-spline error
##
## 01/02/2005: start working on it
## 02/02/2005: finished
## 07/12/2006: extension allowing inclulsion of the estimated correlation between the onset and time-to-event random intercepts
##
bayessurvreg3 <- function
(  formula,
   random,
   formula2,
   random2,
   data = parent.frame(),
   na.action = na.fail,
   onlyX = FALSE,
   nsimul = list(niter = 10, nthin = 1, nburn = 0, nwrite = 10),   
   prior,
   prior.beta,
   prior.b,
   init = list(iter = 0),
   mcmc.par = list(type.update.a = "slice", k.overrelax.a = 1, k.overrelax.sigma = 1, k.overrelax.scale = 1,
                   type.update.a.b = "slice", k.overrelax.a.b = 1, k.overrelax.sigma.b = 1, k.overrelax.scale.b = 1),
   prior2,
   prior.beta2,
   prior.b2,
   init2,
   mcmc.par2 = list(type.update.a = "slice", k.overrelax.a = 1, k.overrelax.sigma = 1, k.overrelax.scale = 1,
                    type.update.a.b = "slice", k.overrelax.a.b = 1, k.overrelax.sigma.b = 1, k.overrelax.scale.b = 1),
   priorinit.Nb,
   rho = list(type.update = "fixed.zero", init=0, sigmaL=0.1),
   store = list(a = FALSE, a2 = FALSE, y = FALSE, y2 = FALSE, r = FALSE, r2 = FALSE, b = FALSE, b2 = FALSE,
                a.b = FALSE, a.b2 = FALSE, r.b = FALSE, r.b2 = FALSE), 
   dir = getwd())
{
  thispackage = "bayesSurv"
  #thispackage = NULL
  store <- bayessurvreg3.checkStore(store)
  nsimul <- bayessurvreg.checknsimul(nsimul)
  
  transform = function(t){log(t)}
  dtransform = function(t){1/t}
  transform2 = function(t){t}
  dtransform2 = function(t){1}
  
  ## Give a function call to be recorded in a resulting object.
  call <- match.call(expand.dots = TRUE)

  ## Extract all the design information from the function call
  ## Do not transform the response at this moment
  doubly <- ifelse(missing(formula2), FALSE, TRUE)
  m <- match.call(expand.dots = FALSE)
  des <- bayessurvreg.design(m=m, formula=formula, random=random, data=data, transform=transform2, dtransform=transform2)
  if (doubly) des2 <- bayessurvreg.design(m=m, formula=formula2, random=random2, data=data, transform=transform2, dtransform=transform2)
  else        des2 <- list(nX = 0, n = des$n, nrandom = 0, randomInt = FALSE)
  
  if (onlyX){
    if (doubly) return(list(X=des$X, X2=des2$X))
    else        return(des$X)
  } 
  
  if (!des$nrandom){
    store$b <- FALSE
    store$a.b <- FALSE
    store$r.b <- FALSE
  }
  if (!des2$nrandom){
    store$b2 <- FALSE
    store$a.b2 <- FALSE
    store$r.b2 <- FALSE
  }  

  ## Perform some checks
  if (des$nrandom > 1) stop("Only intercept may appear in 'random'")
  if (des$nrandom & !des$randomInt) stop("Only intercept may appear in 'random'")
  
  if (des2$nrandom > 1) stop("Only intercept may appear in 'random2'")
  if (des2$nrandom & !des2$randomInt) stop("Only intercept may appear in 'random2'")

  nobs <- des$n
  if (doubly){
    nobs2 <- des2$n
    if (nobs != nobs2) stop("Inconsistent formula and formula2 (different number of observations indicated)")
  }

  ## Priors and inits for beta parameters
  if (missing(init))        init <- list()
  if (missing(init2))       init2 <- list()
  if (missing(mcmc.par))    mcmc.par <- list()
  if (missing(mcmc.par2))   mcmc.par2 <- list()
  if (!doubly)              mcmc.par2 <- list()

  if (missing(prior.beta)) prior.beta <- list()
  betadi <- bayessurvreg3.priorBeta(prior.beta, init, des)
  init$beta <- attr(betadi, "init")
  prior.beta <- attr(betadi, "prior.beta")
  
  if (missing(prior.beta2)) prior.beta2 <- list()
  betadi2 <- bayessurvreg3.priorBeta(prior.beta2, init2, des2)
  init2$beta <- attr(betadi2, "init")
  prior.beta2 <- attr(betadi2, "prior.beta")

  ## Priors and inits for rho (correlation coefficient between two random intercepts)
  if (missing(priorinit.Nb)){
    rho <- bayessurvreg3.checkrho(rho=rho, doubly=doubly)
    if (rho$type.update == "fixed.zero") version <- 3
    else                                 version <- 31
  }
  else{
    rho <- list(type.update = "fixed.zero", init=0, sigmaL=0.1)
    rho <- bayessurvreg3.checkrho(rho=rho, doubly=doubly)    
    version <- 32
  }  
  
  ## Priors and inits for random intercept
  if (version %in% c(3, 31)){
    if (missing(prior.b)) prior.b <- list()  
    reffdi <- bayessurvreg3.priorb(prior.b=prior.b, init=init, design=des, mcmc.par=mcmc.par)
    prior.b  <- attr(reffdi, "prior.b")
    init     <- attr(reffdi, "init")
    mcmc.par <- attr(reffdi, "mcmc.par")
  
    ## Priors and inits for random intercept
    if (missing(prior.b2)) prior.b2 <- list()  
    reffdi2 <- bayessurvreg3.priorb(prior.b=prior.b2, init=init2, design=des2, mcmc.par=mcmc.par2)
    prior.b2  <- attr(reffdi2, "prior.b")
    init2     <- attr(reffdi2, "init")
    mcmc.par2 <- attr(reffdi2, "mcmc.par")
  }
  else{
    if (version == 32){
      Reff <- bayessurvreg3.priorinitNb(priorinit.Nb=priorinit.Nb, init=init, init2=init2, design=des, design2=des2, doubly=doubly)
      reffdi  <- Reff$reffdi
      reffdi2 <- Reff$reffdi2      
      priorinit.Nb <- attr(Reff, "priorinit.Nb")
      init         <- attr(Reff, "init")
      init2        <- attr(Reff, "init2")
      rm(list="Reff")
    }
    else{
      stop("It is strange but I cannot determine the version to be used")
    }  
  }  
  
  ## Priors and inits for error G-spline, (censored) observations and observational allocations
  if (!doubly) prior2 <- list()
  prinit <- bayessurvreg3.priorInit(prior, init, des, mcmc.par, prior2, init2, des2, mcmc.par2, doubly)
  init      <- attr(prinit, "init")
  prior     <- attr(prinit, "prior")  
  mcmc.par  <- attr(prinit, "mcmc.par")
  init2     <- attr(prinit, "init2")
  prior2    <- attr(prinit, "prior2")  
  mcmc.par2 <- attr(prinit, "mcmc.par2")
    
  ## Compute quantities to determine the space needed to be allocated
  ##   and numbers of iterations in different phases
  if (nsimul$nburn >= nsimul$niter) nsimul$nburn <- nsimul$niter - 1
  if (nsimul$nburn < 0) nsimul$nburn <- 0
 
  if (nsimul$nburn == 0) nruns <- 1
  else                   nruns <- 2

  nrun <- numeric(2)
  nrun[2] <- nsimul$niter - nsimul$nburn
  nrun[1] <- nsimul$nburn

  nwrite.run <- nrun
  nwrite.run[nsimul$nwrite <= nrun] <- nsimul$nwrite   
  max.nwrite <- max(nwrite.run)

  ## Write headers to files with stored values  
  bayessurvreg3.writeHeaders(dir=dir, doubly=doubly, prior.init=prinit,
                             priorb.di=reffdi, priorb2.di=reffdi2, store=store, design=des, design2=des2, version=version)  

  ## Combine similar parameters into one vector  
  dims <- c(nobs, as.numeric(doubly))
  storeV <- c(store$a, store$y, store$r, store$b, store$a2, store$y2, store$r2, store$b2, store$a.b, store$r.b, store$a.b2, store$r.b2)
  nsimul.run1 <- c(nrun[1], nsimul$nthin, nwrite.run[1])
  nsimul.run2 <- c(nrun[2], nsimul$nthin, nwrite.run[2])
  names(nsimul.run1) <- names(nsimul.run2) <- c("niter", "nthin", "nwrite")
  nsample1 <- nsimul.run1["niter"] %/% nsimul.run1["nthin"]
  nsample2 <- nsimul.run2["niter"] %/% nsimul.run2["nthin"]  

  cat("Simulation started on                       ", date(), "\n", sep = "")
  fit <- .C("bayessurvreg2", as.character(dir),
                             dims = as.integer(dims),
                             X1 = as.double(if(des$nX) t(des$X) else 0),
                             X2 = as.double(if(des2$nX) t(des2$X) else 0),
                             y1.left = as.double(prinit$y.left),
                             y1.right = as.double(prinit$y.right),
                             status1 = as.integer(prinit$status),
                             t2.left = as.double(prinit$t2.left),
                             t2.right = as.double(prinit$t2.right),
                             status2 = as.integer(prinit$status2),
                             Ys1 = as.double(prinit$y),
                             Ys2 = as.double(prinit$y2),
                             r1 = as.integer(prinit$r),
                             r2 = as.integer(prinit$r2),
                             specif = as.integer(prinit$specification),
                             r1.b = as.integer(reffdi$r),
                             r2.b = as.integer(reffdi2$r),
                             specif.b = as.integer(c(reffdi$specification, reffdi2$specification)),
                             GsplI1 = as.integer(prinit$Gparmi),
                             GsplD1 = as.double(prinit$Gparmd),
                             GsplI2 = as.integer(prinit$Gparmi2),
                             GsplD2 = as.double(prinit$Gparmd2),
                             priorBetaI1 = as.integer(betadi$parmI),
                             priorBetaD1 = as.double(betadi$parmD),
                             priorBetaI2 = as.integer(betadi2$parmI),
                             priorBetaD2 = as.double(betadi2$parmD),
                             priorbI1 = as.integer(reffdi$bparmI),
                             priorbD1 = as.double(reffdi$bparmD),
                             priorbI2 = as.integer(reffdi2$bparmI),
                             priorbD2 = as.double(reffdi2$bparmD),
                             priorCovMatI1 = as.integer(reffdi$GsplI),
                             priorCovMatD1 = as.double(reffdi$GsplD),
                             priorCovMatI2 = as.integer(reffdi2$GsplI),
                             priorCovMatD2 = as.double(reffdi2$GsplD),
                             rho             = as.double(rho$init),
                             rho.accept      = integer(nsample1),
                             rho.type.update = as.integer(rho$typeI),
                             rho.sigmaL      = as.double(rho$sigmaL),           
                             iter = as.integer(prinit$iter),
                             nsimul = as.integer(nsimul.run1),
                             store = as.integer(storeV),
                             version = as.integer(version),
                             mainSimul = as.integer(0),            
                             err = integer(1),
           PACKAGE = thispackage)
  
  if (fit$err != 0) stop ("Something went wrong during the simulation.")
  cat("Burn-up finished on                         ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")
  
  ## Rewrite sampled values by new files
  bayessurvreg3.writeHeaders(dir=dir, doubly=doubly, prior.init=prinit,
                             priorb.di=reffdi, priorb2.di=reffdi2, store=store, design=des, design2=des2, version=version)  
  
  ## Main simulation
  fit <- .C("bayessurvreg2", as.character(dir),
                             dims = as.integer(dims),
                             X1 = as.double(fit$X1),
                             X2 = as.double(fit$X2),
                             y1.left = as.double(fit$y1.left),
                             y1.right = as.double(fit$y1.right),
                             status1 = as.integer(fit$status1),
                             t2.left = as.double(fit$t2.left),
                             t2.right = as.double(fit$t2.right),
                             status2 = as.integer(fit$status2),
                             Ys1 = as.double(fit$Ys1),
                             Ys2 = as.double(fit$Ys2),
                             r1 = as.integer(fit$r1),
                             r2 = as.integer(fit$r2),
                             specif = as.integer(fit$specif),
                             r1.b = as.integer(fit$r1.b),
                             r2.b = as.integer(fit$r2.b),
                             specif.b = as.integer(fit$specif.b),            
                             GsplI1 = as.integer(fit$GsplI1),
                             GsplD1 = as.double(fit$GsplD1),
                             GsplI2 = as.integer(fit$GsplI2),
                             GsplD2 = as.double(fit$GsplD2),
                             priorBetaI1 = as.integer(fit$priorBetaI1),
                             priorBetaD1 = as.double(fit$priorBetaD1),
                             priorBetaI2 = as.integer(fit$priorBetaI2),
                             priorBetaD2 = as.double(fit$priorBetaD2),
                             priorbI1 = as.integer(fit$priorbI1),
                             priorbD1 = as.double(fit$priorbD1),
                             priorbI2 = as.integer(fit$priorbI2),
                             priorbD2 = as.double(fit$priorbD2),
                             priorCovMatI1 = as.integer(fit$priorCovMatI1),
                             priorCovMatD1 = as.double(fit$priorCovMatD1),
                             priorCovMatI2 = as.integer(fit$priorCovMatI2),
                             priorCovMatD2 = as.double(fit$priorCovMatD2),
                             rho             = as.double(fit$rho),
                             rho.accept      = integer(nsample2),
                             rho.type.update = as.integer(fit$rho.type.update),
                             rho.sigmaL      = as.double(fit$rho.sigmaL),                       
                             iter = as.integer(fit$iter),
                             nsimul = as.integer(nsimul.run2),
                             store = as.integer(fit$store),
                             version = as.integer(fit$version),
                             mainSimul = as.integer(1),
                             err = integer(1),
           PACKAGE = thispackage)
  
  if (fit$err != 0) stop ("Something went wrong during the simulation.")  
  cat("Simulation finished on                      ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")     
  
  toreturn <- fit$iter
  attr(toreturn, "call") <- call

  attr(toreturn, "init") <- init
  attr(toreturn, "prior") <- prior
  attr(toreturn, "prior.beta") <- prior.beta
  attr(toreturn, "mcmc.par") <- mcmc.par

  if (doubly){
    attr(toreturn, "init2") <- init2
    attr(toreturn, "prior2") <- prior2
    attr(toreturn, "prior.beta2") <- prior.beta2
    attr(toreturn, "mcmc.par2") <- mcmc.par2
  }

  if (version %in% c(3, 31)){
    attr(toreturn, "prior.b") <- prior.b
    if (doubly) attr(toreturn, "prior.b2") <- prior.b2
  }  
  if (version == 31){
    attr(toreturn, "rho.accept") <- fit$rho.accept
  }  

  if (version == 32){
    attr(toreturn, "priorinit.Nb") <- fit$priorinit.Nb
  }  
    
  class(toreturn) <- "bayessurvreg3"
  return(toreturn)    
}
