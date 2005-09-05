###############################################
#### AUTHOR:    Arnost Komarek             ####
####            (2004)                     ####
####                                       ####
#### FILE:      zzz.R                      ####
####                                       ####
#### FUNCTIONS: .First.lib                 ####
###############################################

### =============================================
### .First.lib
### =============================================
.First.lib <- function(lib, pkg)
{
   require(survival)
   require(coda)
   require(smoothSurv)
   library.dynam("bayesSurv", pkg, lib)

   invisible()
}

