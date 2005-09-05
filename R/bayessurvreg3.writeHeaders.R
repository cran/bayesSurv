####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2005)                         ####
####                                            ####
#### FILE:       bayessurvreg3.writeHeaders.R   ####
####                                            ####
#### FUNCTIONS:  bayessurvreg3.writeHeaders     ####
####################################################

### ======================================
### bayessurvreg3.writeHeaders
### ======================================
## Subfunction for bayessurvreg3.R
##  -> just to make it more readable
##
## Write headers to files where simulated values will be stored
##
bayessurvreg3.writeHeaders <- function(dir, doubly, prior.init, priorb.di, priorb2.di, store, design, design2)
{
  bayesBisurvreg.writeHeaders(dir=dir, dim=1, nP=design$n, doubly=doubly,
                              prior.init=prior.init, store=store, design=design, design2=design2)

  ## Files with sampled values of random effects
  if (store$b){
    sink(paste(dir, "/b.sim", sep = ""), append = FALSE)
    cat(paste(rep(design$names.random, design$ncluster), ".",
              rep(unique(design$cluster), rep(design$nrandom, design$ncluster)), sep = ""), "\n", sep = "    ")
    sink()
  }
  else{
    file.remove(paste(dir, "/b.sim", sep = ""))
  }

  if (doubly){
    if (store$b2){
      sink(paste(dir, "/b_2.sim", sep = ""), append = FALSE)
      cat(paste(rep(design2$names.random, design2$ncluster), ".",
                rep(unique(design2$cluster), rep(design2$nrandom, design2$ncluster)), sep = ""), "\n", sep = "    ")
      sink()
    }
    else{
      file.remove(paste(dir, "/b_2.sim", sep = ""))
    }    
  }
  else{
    file.remove(paste(dir, "/b_2.sim", sep = ""))
  }

  ## Files with sampled G-splines - distribution of the random intercept
  if (design$nrandom){
    write.headers.Gspline(dir=dir, dim=1, nP=design$ncluster, label="_b", gparmi=priorb.di$GsplI,
                          store.a=store$a.b, store.y=FALSE, store.r=store$r.b, care.of.y=FALSE)
  }
  else{
    clean.Gspline(dir=dir, label="_b", care.of.y=FALSE)
  }    

  if (doubly){
    if (design$nrandom){
      write.headers.Gspline(dir=dir, dim=1, nP=design2$ncluster, label="_b2", gparmi=priorb2.di$GsplI,
                            store.a=store$a.b2, store.y=FALSE, store.r=store$r.b2, care.of.y=FALSE)
    }
    else{
      clean.Gspline(dir=dir, label="_b2", care.of.y=FALSE)
    }     
  }
  else{
    clean.Gspline(dir=dir, label="_b2", care.of.y=FALSE)
  }    
}
