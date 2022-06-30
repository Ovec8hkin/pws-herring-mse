# fun_obsm.r
# Created by John Trochta
# Date created:  07/28/2020
# Summary:
# Generate survey observations from pop dynamics of fun_operm
library(gtools)

fun_obsm <- function(pop_dyn, waa, perc.female, y){
  
  waa <- apply(waa, 2, mean)+c(0, 0, 0, rnorm(7))
  female <- rnorm(1, mean=mean(perc.female), sd=sd(perc.female))
  
  # Read in predicted survey/fishery values
  mdm.predicted             <- pop_dyn$mdm
  PWSSC.hydro.predicted     <- pop_dyn$PWSSC.hydro
  ADFG.hydro.predicted      <- pop_dyn$ADFG.hydro
  egg.dep.predictied        <- pop_dyn$egg
  spawn.age.comp.predicted  <- pop_dyn$spawn.age.comp
  seine.age.comp.predicted  <- pop_dyn$seine.age.comp

  # Sampling variances from Muradian et al. 2019
  mdm.sd          <- 0.32
  pwssc.hydro.sd  <- 0.35
  adfg.hydro.sd   <- 0.29
  egg.sd          <- 0.35

  # Generate observations with error
  mdm <- mdm.predicted[y]*exp(rnorm(1, 0, mdm.sd)-mdm.sd^2/2)
  pwssc.hydro <- PWSSC.hydro.predicted[y]*exp(rnorm(1, 0, pwssc.hydro.sd)-pwssc.hydro.sd^2/2)
  
  #spawn.age.comp <- gtools::rdirichlet(1, 39*spawn.age.comp.predicted[y, ])
  spawn.age.comp <- rmultinom(1, 1500, spawn.age.comp.predicted[y, ]/sum(spawn.age.comp.predicted[y, ]))
  spawn.age.comp <- t(spawn.age.comp/sum(spawn.age.comp))

  if(any(seine.age.comp.predicted==-9)){
    seine.age.comp <- rep(-9, ncol(seine.age.comp.predicted))
  }else{
    seine.age.comp <- gtools::rdirichlet(1, 119*seine.age.comp.predicted[y, ])
  }

  if(any(ADFG.hydro.predicted==-9)){
    adfg.hydro <- -9
    adfg.hydro.sd <- -9
  }else{
    adfg.hydro <- ADFG.hydro.predicted[y]*exp(rnorm(1, 0, adfg.hydro.sd)-adfg.hydro.sd^2/2)
  }

  if(any(egg.dep.predictied==-9)){
    egg <- -9
    egg.sd <- -9
  }else{
    egg <- egg.dep.predictied[y]*exp(rnorm(1, 0, egg.sd)-egg.sd^2/2)
  }

  aerial.juv <- -9
  vhsv.antibody <- rep(-9, 20)
  ich.antibody  <- rep(-9, 20)

  return(listN(waa, female, mdm, mdm.sd, pwssc.hydro, pwssc.hydro.sd, adfg.hydro, adfg.hydro.sd, egg, egg.sd, spawn.age.comp, seine.age.comp, aerial.juv, vhsv.antibody, ich.antibody))
}

# Function for simultaneously creating names of variables/elements within a list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}