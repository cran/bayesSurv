## Subfunction for bayessurvreg1.R
##  -> just to make it more readable
##
## Manipulation with initial values and prior specifications

bayessurvreg1.priorInit <- function(prior, init, Yinit, Xinit, n, nX, nrandom, ncluster, indb, randomInt, toler.chol){

   ## ============================================================================================
   ## Prior parameters (calculate these that were not given by the user and change the notation)
   ## part 1
   ## ============================================================================================
   prior.pari <- numeric(3)
   names(prior.pari) <- c("kmax", "k.prior", "Eb0.depend.mix")
   if (is.null(prior$kmax)) prior$kmax <- 5
   if (is.null(prior$k.prior)) prior$k.prior <- "poisson"
   if (is.null(prior$Eb0.depend.mix)) prior$Eb0.depend.mix <- FALSE
   if (is.null(prior$poisson.k)) prior$poisson.k <- 3
   prior.pari["k.prior"] <- pmatch(prior$k.prior, c("poisson", "uniform", "fixed"), nomatch = -1) - 1
               ## 0 = Poisson, 1 = Uniform, 2 = Fixed
   if (prior.pari["k.prior"] < 0) stop("Prior for k (number of mixture components) must be either poisson, uniform or fixed.")
   prior.pari["kmax"] <- prior$kmax
   prior.pari["Eb0.depend.mix"] <- 1*(prior$Eb0.depend.mix)

   prior.pard <- numeric(2*prior.pari["kmax"] + 7)
   names(prior.pard) <- c(paste("pi.split", 1:prior.pari["kmax"], sep = ""), paste("pi.birth", 1:prior.pari["kmax"], sep = ""),
                          "lambda", "delta", "xi", "kappa", "zeta", "g", "h")
   prior.pard["lambda"] <- prior$poisson.k


   ## =========================================================
   ## Get initial estimates
   ## =========================================================
   init.error <- "Something is wrong with your initials."
   fit.init <- survreg(Yinit ~ Xinit - 1, dist = "lognormal")          ## intercept is already included in Xinit

     ## index of the first iteration           
   if (is.null(init$iter) | is.na(init$iter)) init$iter <- 0                                 
   else                                       init$iter <- init$iter[1]
   
     ## initial mixture   
   if (is.null(init$mixture)){
     if (prior.pari["k.prior"] == 2) stop("init$mixture must be given when prior$k.prior is 'fixed'.")
     init$mixture <- numeric(1 + 3*prior$kmax)
     init$mixture[1] <- 1                                        ## initial k
     init$mixture[2] <- 1.0                                      ## initial weight of the first mixture component   
     init$mixture[2 + prior$kmax] <- fit.init$coefficients[1]    ## initial mean of the first mixture component
     init$mixture[2 + 2*prior$kmax] <- fit.init$scale^2          ## initial variance of the first mixture component
   }
   else{
     if (length(init$mixture) != 1 + 3*prior$kmax) stop("Incorrect init$mixture parameter supplied.")
     if (is.na(init$mixture[1])) stop("Incorrect init$mixture parameter supplied.")
     wi <- init$mixture[2:(1 + init$mixture[1])]
     mui <-init$mixture[(2 + prior$kmax):(1 + prior$kmax + init$mixture[1])]       
     sig2i <- init$mixture[(2 + 2*prior$kmax):(1 + 2*prior$kmax + init$mixture[1])]   
     if (sum(is.na(wi)) | sum(is.na(mui)) | sum(is.na(sig2i))) stop("Incorrect init$mixture parameter supplied.")
     if (sum(wi < 0)) stop("Incorrect init$mixture parameter supplied.")
     if (sum(sig2i <= 0)) stop("Incorrect init$mixture parameter supplied.")
     wi <- wi/sum(wi)                              ## to make sure that the sum is 1
     ordermu <- order(mui)
     k.temp <- init$mixture[1]
     init$mixture <- numeric(1 + 3*prior$kmax)     
     wi <- wi[ordermu]
     mui <- mui[ordermu]
     sig2i <- sig2i[ordermu]
     init$mixture[1] <- k.temp
     init$mixture[2:(1 + init$mixture[1])] <- wi
     init$mixture[(2 + prior$kmax):(1 + prior$kmax + init$mixture[1])] <- mui
     init$mixture[(2 + 2*prior$kmax):(1 + 2*prior$kmax + init$mixture[1])] <- sig2i     
   }

     ## initial beta parameter
   if (!nX) init$beta <- 0
   else{
     if (is.null(init$beta)){
       init$beta <- fit.init$coefficients[-1]                  ## remove the intercept
     }
     else{
       if (length(init$beta) < nX) stop("Incorrect init$beta parameter supplied.")
       init$beta <- init$beta[1:nX]
     }
     if (sum(is.na(init$beta))) stop("Incorrect init$beta parameter supplied.")
   }     
   
     ## initial values of the random effects
   if (!nrandom) init$b <- 0
   else{
     if (is.null(init$b)){
       bb <- fit.init$coefficients[-1][indb > 0]
       if (randomInt) bb <- c(0, bb)
       init$b <- rep(bb, ncluster)
     } 
     else{
       if (length(init$b) < nrandom*ncluster) stop("Incorrect init$b parameter supplied.")
       init$b <- init$b[1:(nrandom*ncluster)]
     }
     if (sum(is.na(init$b))) stop("Incorrect init$b parameter supplied.")
   }     
   
     ## initial values of the (transformed) latent response
   if (is.null(init$y)){
     init$y <- as.numeric(log(Yinit[,1]))          
   } 
   else{
     if (length(init$y) < n) stop("Incorrect init$y parameter supplied.")    
     init$y <- init$y[1:n]
   }
   if (sum(is.na(init$y))) stop("Incorrect init$y parameter supplied.")
   
     ## initial values of the component pertinences
   if (is.null(init$r)){
     init$r <- numeric(n) + 1         ## initially, everyone belongs to the first mixture component    
   }
   else{
     if (length(init$r) < n) stop("Incorrect init$r parameter supplied.")
     init$r <- init$r[1:n]
   } 
   if (sum(is.na(init$r)) | sum(init$r <= 0) | sum(init$r > init$mixture[1]))
     stop("Incorrect init$r parameter supplied.")
   
     ## initial values of the matrix D (covariance matrix of the random effects)
     ## (lower triangle of the matrix D in column major order)
   if (!nrandom) init$D <- 0   
   else{
     if (is.null(init$D)){
       init$D <- diag(nrandom)[lower.tri(diag(nrandom), diag = TRUE)]   ## identity matrix as initial D
     }
     else{
       if (length(init$D) < 0.5*nrandom*(1+nrandom)) stop("Incorrect init$D parameter supplied.")
       init$D <- init$D[1:(0.5*nrandom*(1+nrandom))]
     }
     if (sum(is.na(init$D))) stop("Incorrect init$D parameter supplied.")
   }     
   
     ## initial values of the remaining parameters (only eta at this moment)
   if (is.null(init$otherp)){
     init$otherp <- 1              ## initial of eta
   }
   else{
     init$otherp <- init$otherp[1]
   }
   if (sum(is.na(init$otherp))) stop("Incorrect init$otherp parameter supplied.")
   
     ## initial proposal vector for a split-combine move
   if (is.null(init$u)){
     init$u <- c(runif(1), 0, 0, runif(3*(prior$kmax - 1)))
   }
   else{
     if (length(init$u) < 3*prior$kmax) stop("Incorrect init$u parameter supplied.")  
     init$u <- init$u[1:(3*prior$kmax)]
   }     
   if (sum(is.na(init$u))) stop("Incorrect init$u parameter supplied.")  
   if (sum(init$u < 0 | init$u > 1)) stop("Incorrect init$u parameter supplied.")  

   
   ## =========================================================================================
   ## Prior parameters (calculate these that were not given by the user and change the notation)
   ## part 2
   ## =========================================================================================
   if (is.null(prior$dirichlet.w)) prior$dirichlet.w <- 1
   if (prior$dirichlet.w < 1) stop ("prior$dirichlet.w must be at least 1.")
   prior.pard["delta"] <- prior$dirichlet.w

   if (is.null(prior$mean.mu)) prior$mean.mu <- init$mixture[2 + prior.pari["kmax"]]    ## mean of the first comp.
   if (is.null(prior$var.mu)) prior$var.mu <- 2*init$mixture[2 + 2*prior.pari["kmax"]]  ## 2*variance of the first comp.
   if (prior$var.mu <= 0) stop("prior$var.mu must be positive.")
   prior.pard["xi"] <- prior$mean.mu
   prior.pard["kappa"] <- prior$var.mu

   if (is.null(prior$shape.invsig2)) prior$shape.invsig2 <- 1.5
   if (is.null(prior$shape.hyper.invsig2)) prior$shape.hyper.invsig2 <- 0.8
   if (is.null(prior$rate.hyper.invsig2)) prior$rate.hyper.invsig2 <- prior$var.mu
   if (prior$shape.invsig2 <= 0) stop("prior$shape.invsig2 must be positive.")
   if (prior$shape.hyper.invsig2 <= 0) stop("prior$shape.hyper.invsig2 must be positive.")
   if (prior$rate.hyper.invsig2 <= 0) stop("prior$rate.hyper.invsig2 must be positive.")      
   prior.pard["zeta"] <- prior$shape.invsig2
   prior.pard["g"] <- prior$shape.hyper.invsig2
   prior.pard["h"] <- prior$rate.hyper.invsig2

   if (is.null(prior$pi.split)) prior$pi.split <- c(1, rep(0.5, prior.pari["kmax"]-2), 0)
   if (is.null(prior$pi.birth)) prior$pi.birth <- c(1, rep(0.5, prior.pari["kmax"]-2), 0)
   if (prior$pi.split[1] != 1) stop("prior$pi.split[1] must be equal to 1.")
   if (prior$pi.birth[1] != 1) stop("prior$pi.birth[1] must be equal to 1.")
   if (prior$pi.split[prior.pari["kmax"]] != 0) stop("prior$pi.split[kmax] must be equal to 0.")
   if (prior$pi.birth[prior.pari["kmax"]] != 0) stop("prior$pi.birth[kmax] must be equal to 0.")      
   if (length(prior$pi.split) != prior.pari["kmax"]) stop("Incorrect length of a vector prior$pi.split.")
   if (length(prior$pi.birth) != prior.pari["kmax"]) stop("Incorrect length of a vector prior$pi.birth.")      
   prior.pard[paste("pi.split", 1:prior.pari["kmax"], sep = "")] <- prior$pi.split
   prior.pard[paste("pi.birth", 1:prior.pari["kmax"], sep = "")] <- prior$pi.birth
   
   priordi <- list(integer = prior.pari, double = prior.pard)
   attr(priordi, "init") <- init
   attr(priordi, "prior") <- prior
   
   return(priordi)
   
 }
