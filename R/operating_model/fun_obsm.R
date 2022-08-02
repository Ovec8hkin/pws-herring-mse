# fun_obsm.r
# Created by John Trochta
# Date created:  07/28/2020
# Summary:
# Generate survey observations from pop dynamics of fun_operm
library(gtools)

fun_obsm <- function(survey.indices, waa, fec, perc.female, juv.overdisp, sample.sizes, y, sim.seed,
                     survey.controls = list(mdm=TRUE, egg=FALSE, pwssc.hydro=TRUE, adfg.hydro=FALSE, spac=TRUE, seac=TRUE, juv=FALSE, vhsv=FALSE, ich=FALSE)){
  
  set.seed(sim.seed)

  waa <- apply(waa, 2, mean)+c(0, 0, 0, rnorm(7))
  fec <- c(0, 0, 0, apply(fec[fec[,5]!=-9,4:10], 2, mean))+c(0, 0, 0, rnorm(7))
  female <- rnorm(1, mean=mean(perc.female), sd=sd(perc.female))
  
  # Read in predicted survey/fishery values
  mdm.predicted             <- survey.indices$mdm
  PWSSC.hydro.predicted     <- survey.indices$PWSSC.hydro
  ADFG.hydro.predicted      <- survey.indices$ADFG.hydro
  egg.dep.predicted         <- survey.indices$egg
  spawn.age.comp.predicted  <- survey.indices$spawn.age.comp
  seine.age.comp.predicted  <- survey.indices$seine.age.comp
  vhsv.antibody.predicted   <- survey.indices$vhsv.antibody
  ich.antibody.predicted    <- survey.indices$ich.antibody

  # Sampling variances from Muradian et al. 2019
  mdm.sd          <- 0.32
  pwssc.hydro.sd  <- 0.35
  adfg.hydro.sd   <- 0.29
  egg.sd          <- 0.35

  # Generate observations with error
  mdm <- mdm.predicted[y]*exp(rnorm(1, 0, mdm.sd)-mdm.sd^2/2)
  pwssc.hydro <- PWSSC.hydro.predicted[y]*exp(rnorm(1, 0, pwssc.hydro.sd)-pwssc.hydro.sd^2/2)
  
  #spawn.age.comp <- gtools::rdirichlet(1, 39*spawn.age.comp.predicted[y, ])
  spawn.age.comp <- rmultinom(1, sample.sizes$spac, spawn.age.comp.predicted[y, ]/sum(spawn.age.comp.predicted[y, ]))
  spawn.age.comp <- t(spawn.age.comp/sum(spawn.age.comp))

  if(any(seine.age.comp.predicted==-9)){
    seine.age.comp <- rep(-9, ncol(seine.age.comp.predicted))
  }else{
    #seine.age.comp <- gtools::rdirichlet(1, 119*seine.age.comp.predicted[y, ])
    seine.age.comp <- rmultinom(1, sample.sizes$seac, seine.age.comp.predicted[y, ]/sum(seine.age.comp.predicted[y, ]))
    seine.age.comp <- t(seine.age.comp/sum(seine.age.comp))
  }

  adfg.hydro <- ifelse(survey.controls$adfg.hydro, ADFG.hydro.predicted[y]*exp(rnorm(1, 0, adfg.hydro.sd)-adfg.hydro.sd^2/2), -9)
  adfg.hydro.sd <- ifelse(survey.controls$adfg.hydro, adfg.hydro.sd, -9)

  egg <- ifelse(survey.controls$egg, egg.dep.predicted[y]*exp(rnorm(1, 0, egg.sd)-egg.sd^2/2), -9)
  egg.sd <- ifelse(survey.controls$egg, egg.sd, -9)

  #aerial.juv <- -9 # negative binomial error distribution (check formula)
  aerial.juv <- ifelse(survey.controls$juv, rnbinom(1, size=1/juv.overdisp, mu=survey.indices$juv.schools), -9)
  vhsv.antibody <- ifelse(survey.controls$vhsv, vhsv.antibody.predicted[y, ], rep(-9, 20)) # ask John for simulation code for these
  ich.antibody  <- ifelse(survey.controls$ich,  ich.antibody.predicted[y, ], rep(-9, 20))

  return(listN(waa, fec, female, 
               mdm, mdm.sd, 
               pwssc.hydro, pwssc.hydro.sd, 
               adfg.hydro, adfg.hydro.sd, 
               egg, egg.sd, 
               spawn.age.comp, seine.age.comp, 
               aerial.juv, 
               vhsv.antibody, ich.antibody
              )
        )
}

# Function for simultaneously creating names of variables/elements within a list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}