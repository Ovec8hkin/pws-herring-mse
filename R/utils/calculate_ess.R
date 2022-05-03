calculate.ess <- function(age.comps, start.ess, samp.size, ess.its, nyr.fit){
    its=1
    for(i in 1:2){
        # Create "PWS_ASA(ESS_estimate).ctl" with sample sizes (the original sample sizes on the first iteration)
        write.table(
            rbind("# PWS age comp effective sample sizes", 
                "# Seine ESS",      start.ess$seine,  " ", 
                "# Spawn ESS",      start.ess$spawn,  " ", 
                "# VHS sero ESS",   start.ess$vhsv,   " ", 
                "# Ich prev ESS",   start.ess$ich
            ),
            file = "PWS_ASA(ESS_estimate).ctl", 
            append = F, 
            sep = " ",
            row.names=FALSE, 
            col.names=FALSE, 
            quote=F
        )
        
        # Compile and Run PWS_ASA
        if(OS=="MAC"){
            # system("./PWS_ASA -pinwrite -nohess")
            system("./PWS_ASA -pinwrite", intern=TRUE)
        }else if(OS=="PC"){
            shell("PWS_ASA  -pinwrite -nohess")
        }
        
        
        # Read in the estimated seine and spawner age comps
        seine.age.comp.est <- read.table("rep_out/SeAC_pd.rep",     header = FALSE) 
        spawn.age.comp.est <- read.table("rep_out/SpAC_pd.rep",     header = FALSE)
        vhsv.age.comp.est  <- read.table("rep_out/vhscomp_pd.rep",  header = FALSE)
        ich.age.comp.est   <- read.table("rep_out/ichcomp_pd.rep",  header = FALSE)
        
        # Calculate the ESS
        seine.ess <- rowSums(seine.age.comp.est*(1-seine.age.comp.est))/rowSums((age.comps$seine[1:nyr.fit, ]-seine.age.comp.est)^2)
        spawn.ess <- rowSums(spawn.age.comp.est*(1-spawn.age.comp.est))/rowSums((age.comps$spawn[1:nyr.fit, ]-spawn.age.comp.est)^2)
        vhsv.ess  <- rowSums(vhsv.age.comp.est*(1-vhsv.age.comp.est))/rowSums((age.comps$vhsv[1:nyr.fit, ]-vhsv.age.comp.est)^2)
        ich.ess   <- rowSums(ich.age.comp.est*(1-ich.age.comp.est))/rowSums((age.comps$ich[1:nyr.fit, ]-ich.age.comp.est)^2)
        
        # Remove the missing years of age comps
        seine.ess.rem <- seine.ess[!(age.comps$seine[1:nyr.fit, 1]==-9)]
        spawn.ess.rem <- spawn.ess[!(age.comps$spawn[1:nyr.fit, 1]==-9)]
        vhsv.ess.rem <- vhsv.ess[!(age.comps$vhsv[1:nyr.fit, 1]==-9)]
        ich.ess.rem <- ich.ess[!(age.comps$ich[1:nyr.fit, 1]==-9)]
        
        # Calculate the ratio of ESS to original sample sizes
        seine.ratio <- seine.ess.rem/samp.size$seine[1:nyr.fit, 1][age.comps$seine[1:nyr.fit, 1]!=-9]
        spawn.ratio <- spawn.ess.rem/samp.size$spawn[1:nyr.fit, 1][age.comps$spawn[1:nyr.fit, 1]!=-9]
        vhsv.ratio  <- vhsv.ess.rem/samp.size$vhsv[1:nyr.fit, 1][age.comps$vhsv[1:nyr.fit, 1]!=-9]
        ich.ratio   <- ich.ess.rem/samp.size$ich[1:nyr.fit, 1][age.comps$ich[1:nyr.fit, 1]!=-9]
        
        # Calculate the harmonic means (see Muradian et al. 2017 and Stewart & Hamel 2014)
        seine.hm <- 1/mean(1/seine.ratio)
        spawn.hm <- 1/mean(1/spawn.ratio)
        vhsv.hm <- 1/mean(1/vhsv.ratio)
        ich.hm <- 1/mean(1/ich.ratio)
        
        # Compare this harmonic mean to the previous using a convergence criteria (WHAT AM I CONVERGING!!!!)
        if(its==1) {
            convergence <- 0
            seine.hmS <- seine.hm
            spawn.hmS <- spawn.hm
            vhsv.hmS  <- vhsv.hm
            ich.hmS   <- ich.hm
        }else{
            seine.test <- abs(seine.hm - seine.hmS[its-1])/seine.hmS[its-1]*100
            spawn.test <- abs(spawn.hm - spawn.hmS[its-1])/spawn.hmS[its-1]*100
            vhsv.test  <- abs(vhsv.hm - vhsv.hmS[its-1])/vhsv.hmS[its-1]*100
            ich.test   <- abs(ich.hm - ich.hmS[its-1])/ich.hmS[its-1]*100

            # This criteria was arbitrarily chosen (0.1% change)
            convergence <- (seine.test<0.1 & spawn.test<0.1 & vhsv.test<0.1) 
            
            seine.hmS <- rbind(seine.hmS, seine.hm)
            spawn.hmS <- rbind(spawn.hmS, spawn.hm) 
            vhsv.hmS  <- rbind(vhsv.hmS, vhsv.hm) 
            ich.hmS   <- rbind(ich.hmS, ich.hm) 
        }
        
        # Now multiply the harmonic mean by the sample size to get the new ESS 
        seine.ess <- round(seine.hm*samp.size$seine, 0)
        spawn.ess <- round(spawn.hm*samp.size$spawn, 0)
        vhsv.ess  <- round(vhsv.hm*samp.size$vhsv,   0)
        ich.ess   <- round(ich.hm*samp.size$ich,     0)
        
        # Use the average ESS for all years (each years obs weighted equally)
        # seine.ess[seine.ess>0] <- round(mean(seine.ess[seine.ess>0]), digits=0)
        # spawn.ess[spawn.ess>0] <- round(mean(spawn.ess[spawn.ess>0]), digits=0)
        
        # Denote the missing values
        seine.ess[(age.comps$seine[, 1] == -9), 1] <- -9
        spawn.ess[(age.comps$spawn[, 1] == -9), 1] <- -9
        vhsv.ess[(age.comps$vhsv[, 1]   == -9), 1] <- -9
        ich.ess[(age.comps$ich[, 1]     == -9), 1] <- -9
        
        vhsv.ess <- samp.size$vhsv
        ich.ess <- samp.size$ich
        
        # Fill in this iteration"s ESS
        ess.its$seine.ess.its <- cbind(ess.its$seine.ess.its, round(seine.ess, 0))
        ess.its$spawn.ess.its <- cbind(ess.its$spawn.ess.its, round(spawn.ess, 0))
        ess.its$vhsv.ess.its <- cbind(ess.its$vhsv.ess.its, round(vhsv.ess, 0))
        ess.its$ich.ess.its <- cbind(ess.its$ich.ess.its, round(ich.ess, 0))
        
        # Cease iterations if convergence hasn"t happened after so many...
        if(its==10) {
            break
        }
        its <- its+1
    }
    return(list(seine.ess=seine.ess, spawn.ess=spawn.ess, vhsv.ess=vhsv.ess, ich.ess=ich.ess))
}