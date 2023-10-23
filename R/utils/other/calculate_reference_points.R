library(tidyverse)
source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

total.mortality <- function(F, selectivity, mortality){
    return(-(F*selectivity + mortality))
}

calc.naa.pr <- function(mu, v, s){
    age.classes <- length(v)
    naa <- rep(0, age.classes)
    naa[1] <- 1
    for(a in 2:(age.classes-1)){
        naa[a] <- naa[a-1]*exp(total.mortality(mu, v[a-1], s[a-1]))
    }
    a <- a+1
    naa[a] <-  naa[a-1]*exp(total.mortality(mu, v[a-1], s[a-1]))/(1-exp(total.mortality(mu, v[a-1], s[a-1])))
    return(naa)
}

calc.ypr <- function(mu, v, s, waa){
    naa <- calc.naa.pr(mu, v, s)
    return(sum(naa * waa * (1-exp(total.mortality(mu, v, s))) * ((mu*v)/(mu*v + s))))
}

calc.ypr0 <- function(v, s, waa){
    return(calc.ypr(0, v, s, waa))
}

calc.sbpr <- function(mu, v, s, fec){
    naa <- calc.naa.pr(mu, v, s)
    return(sum(naa * fec))
}

calc.sbpr0 <- function(v, s, fec){
    return(calc.sbpr(0, v, s, fec))
}

calc.ssb0 <- function(v, s, f, R0){
    return(calc.sbpr0(v, s, f)*R0)
}

calc.msy <- function(){

}

compute.mortality <- function(model.dir, start.year=1980, end.year=NA, nyr=NA){

    if(is.na(end.year) & is.na(nyr)){
        warning("Neither end.year nor nyr specified. Please specify one.")
        stop()
    }

    if(!is.na(nyr)){
         end.year <- start.year+nyr
    }

    if(!is.na(end.year)){
        nyr <- end.year-start.year
    }
   
    survival.sum <- readr::read_csv(paste0(model.dir, "/mcmc_out/adult_survival_effects_summer.csv"), col_names = FALSE) %>% select_at(((start.year-1980)*10+1):((end.year-1980+1)*10))
    survival.win <- readr::read_csv(paste0(model.dir, "/mcmc_out/adult_survival_effects_winter.csv"), col_names = FALSE) %>% select_at(((start.year-1980)*10+1):((end.year-1980+1)*10))
    survival <- survival.sum*survival.win

    mortality <- tibble(survival) %>%
                    `colnames<-`(paste0(rep(paste0("age",0:9),times=nyr), "_", rep(start.year:end.year,each=10))) %>%
                    mutate(mcmc=1:nrow(survival)) %>%
                    pivot_longer(                               # Convert to long format
                        !mcmc,
                        names_to = c("age","year"),
                        names_sep="_",
                        values_to="survival"
                    ) %>%
                    mutate(
                        mortality=-log(survival),
                        year = as.integer(year)
                    ) %>%        # Compute mortality from survival
                    group_by(mcmc, age, year) %>%               # Group by age class and year
                    summarise(mortality=mean(mortality)) %>%    # mean mortality by age class and mcmc save
                    ungroup() %>%
                    pivot_wider(                                # Convert back to  
                        names_from = age,
                        values_from = mortality
                    ) %>%
                    select(!mcmc) %>%
                    print(n=20)

    return(mortality)

}

compute.waa <- function(model.dir, start.year=1980, end.year=NA, nyr=NA, wide=FALSE){

    if(is.na(end.year) & is.na(nyr)){
        warning("Neither end.year nor nyr specified. Please specify one.")
        stop()
    }

    if(!is.na(nyr)){
         end.year <- start.year+nyr
    }

    if(!is.na(end.year)){
        nyr <- end.year-start.year
    }

    waa <- as_tibble(read.data.files(paste0(model.dir, "/"))$PWS_ASA.dat$waa) %>%
                `colnames<-`(paste0("age",0:9)) %>%
                mutate(year=1980:(1980+nrow(.)-1)) %>%
                filter(year >= start.year & year <= end.year) %>%
                pivot_longer(!year,names_to = "age",values_to="waa") %>%
                mutate(year=as.integer(year))

    if(wide){
        waa <- waa %>% 
                pivot_wider(
                    names_from=age,
                    values_from = waa
                ) %>%
                mutate(count=5200) %>%
                uncount(count)
    }

    return(waa)
}

compute.fecundity <- function(model.dir, start.year=1980, end.year=NA, nyr=NA){

    if(is.na(end.year) & is.na(nyr)){
        warning("Neither end.year nor nyr specified. Please specify one.")
        stop()
    }

    if(!is.na(nyr)){
         end.year <- start.year+nyr
    }

    if(!is.na(end.year)){
        nyr <- end.year-start.year+1
        print(nyr)
    }

    return(
        as_tibble(read.data.files(paste0(model.dir, "/"))$PWS_ASA.dat$fecundity) %>%
            `colnames<-`(paste0("age",0:9)) %>%
            mutate(year=1980:(1980+nrow(.)-1)) %>%
            filter(year >= start.year & year <= end.year) %>%
            pivot_longer(!year,names_to = "age",values_to="fec") %>%
            mutate(
                year=as.integer(year),
                fec=ifelse(fec==-9, NA, fec)
            ) %>%
            group_by(age) %>%
            filter(!is.na(fec)) %>%
            summarise(fec=mean(fec, na.rm=TRUE)) %>%
            pivot_wider(
                names_from=age,
                values_from=fec
            ) %>%
            mutate(count=5200*nyr) %>%
            uncount(count) %>%
            mutate(year = rep(start.year:end.year, 5200))
    )
}

compute.naa <- function(model.dir, start.year=1980, end.year=NA, nyr=NA){

    if(is.na(end.year) & is.na(nyr)){
        warning("Neither end.year nor nyr specified. Please specify one.")
        stop()
    }

    if(!is.na(nyr)){
         end.year <- start.year+nyr
    }

    if(!is.na(end.year)){
        nyr <- end.year-start.year
    }

    return(
        readr::read_csv(paste0(model.dir, "/mcmc_out/Num_at_age.csv"),col_names = FALSE) %>%
            select_at(((start.year-1980)*10+1):((end.year-1980)*10)) %>%
            `colnames<-`(paste0(rep(paste0("age",0:9),times=nyr), "_", rep(start.year:(end.year-1),each=10))) %>%
            mutate(mcmc=1:nrow(.)) %>% 
            pivot_longer(
                !mcmc,
                names_to = c("age","year"),
                names_sep="_",
                values_to="naa"
            ) %>% 
            mutate(year = as.integer(year))
    )
}

compute.selectivity <- function(model.dir, start.year=1980, end.year=NA, nyr=NA){

    catches <- compute.catches(model.dir, start.year, end.year, nyr)

    catches <- catches %>%
                mutate(
                    exploit_tot = ifelse(is.nan(catch/biomass) | is.infinite(catch/biomass), 0, catch/biomass),
                    exploit_std = exploit_tot/sum(exploit_tot),
                    select_tot = exploit_std/max(exploit_std),
                    exploit_se = ifelse(is.nan(se/biomass) | is.infinite(se/biomass), 0, se/biomass),
                    exploit_se_std = exploit_se/sum(exploit_se),
                    select_se = exploit_se_std/max(exploit_se_std),
                    exploit_fb = ifelse(is.nan(fb/biomass) | is.infinite(fb/biomass), 0, fb/biomass),
                    exploit_fb_std = exploit_fb/sum(exploit_fb),
                    select_fb = exploit_fb_std/max(exploit_fb_std),
                    exploit_pk = ifelse(is.nan(pk/biomass) | is.infinite(pk/biomass), 0, pk/biomass),
                    exploit_pk_std = exploit_pk/sum(exploit_pk),
                    select_pk = exploit_pk_std/max(exploit_pk_std),
                    exploit_gn = ifelse(is.nan(gn/biomass) | is.infinite(gn/biomass), 0, gn/biomass),
                    exploit_gn_std = exploit_gn/sum(exploit_gn),
                    select_gn = exploit_gn_std/max(exploit_gn_std)
                ) %>%
                select(c(mcmc, age, year, waa, se, fb, pk, gn, naa, catch, biomass, select_tot, select_se, select_fb, select_pk, select_gn)) %>%
                print(n=10)
    
    selectivity <- catches %>% select(mcmc,year,age, select_tot, select_se, select_fb, select_pk, select_gn) %>%
                    group_by(mcmc,year) %>%
                    ungroup() %>%
                    group_by(mcmc,age,year) %>%
                    summarise(
                        sel_tot=mean(select_tot),
                        sel_se=mean(select_se),
                        sel_fb=mean(select_fb),
                        sel_pk=mean(select_pk),
                        sel_gn=mean(select_gn)
                    ) %>%
                    ungroup() %>%
                    pivot_wider(
                        names_from = age,
                        values_from = starts_with("sel")
                    ) %>%
                    select(!mcmc) %>%
                    replace(is.na(.), 0) %>%
                    print(n=30)

    return(
        list(
            full = selectivity %>% select(c(year, starts_with("sel_tot"))) %>% `colnames<-`(c("year", paste0("age",0:9))),
            se = selectivity %>% select(c(year, starts_with("sel_se"))) %>% `colnames<-`(c("year", paste0("age",0:9))),
            fb = selectivity %>% select(c(year, starts_with("sel_fb"))) %>% `colnames<-`(c("year", paste0("age",0:9))),
            pk = selectivity %>% select(c(year, starts_with("sel_pk"))) %>% `colnames<-`(c("year", paste0("age",0:9))),
            gn = selectivity %>% select(c(year, starts_with("sel_gn"))) %>% `colnames<-`(c("year", paste0("age",0:9)))
        )
    )
}

compute.catches <- function(model.dir, start.year=1980, end.year=NA, nyr=NA){

    if(is.na(end.year) & is.na(nyr)){
        warning("Neither end.year nor nyr specified. Please specify one.")
        stop()
    }

    if(!is.na(nyr)){
         end.year <- start.year+nyr
    }

    if(!is.na(end.year)){
        nyr <- end.year-start.year
    }

    raw.data <- read.data.files(paste0(model.dir, "/"))$PWS_ASA.dat

    fb <- as_tibble(raw.data$foodbait_catch)
    pk <- as_tibble(raw.data$pound_catch)
    gn <- as_tibble(raw.data$gillnet_catch)
    se <- as_tibble(raw.data$seine_yield)
    
    names(fb) <- names(pk) <- names(gn) <- paste0("age",0:9)
    names(se) <- "se"

    fb <- fb %>% mutate(year=1980:(1980+nrow(.)-1)) %>% filter(year >= start.year & year <= end.year) %>% pivot_longer(!year,names_to = "age",values_to="fb")
    pk <- pk %>% mutate(year=1980:(1980+nrow(.)-1)) %>% filter(year >= start.year & year <= end.year) %>% pivot_longer(!year,names_to = "age",values_to="pk")
    gn <- gn %>% mutate(year=1980:(1980+nrow(.)-1)) %>% filter(year >= start.year & year <= end.year) %>% pivot_longer(!year,names_to = "age",values_to="gn")
    se <- se %>% mutate(year=1980:(1980+nrow(.)-1)) %>% filter(year >= start.year & year <= end.year)

    waa <- compute.waa(model.dir, start.year, end.year, nyr)

    se <- readr::read_csv(paste0(model.dir, "/mcmc_out/SeAC.csv"), col_names = FALSE) %>%
                    select_at(((start.year-1980)*10+1):((end.year-1980+1)*10)) %>%
                    `colnames<-`(paste0(rep(paste0("age",0:9),times=nyr), "_", rep(start.year:end.year,each=10))) %>%
                    mutate(mcmc=1:nrow(.)) %>%
                    pivot_longer(
                        !mcmc,
                        names_to = c("age","year"),
                        names_sep="_",
                        values_to="ac"
                    ) %>% 
                    mutate(year = as.integer(year)) %>%
                    left_join(waa) %>% left_join(se) %>% 
                    group_by(mcmc,year) %>%
                    mutate(se=se/sum(waa*ac)*ac) %>%
                    ungroup() %>% select(!ac)

    naa <- compute.naa(model.dir, start.year, end.year, nyr)

    catches <- se %>% left_join(fb) %>% left_join(pk) %>% left_join(gn) %>% left_join(naa) %>%
                group_by(mcmc, year) %>%
                mutate(
                    catch = (se+fb+pk+gn),
                    biomass = naa*waa,
                ) %>%
                print(n=10)

    return(catches)
}

compute.exploit.rate <- function(model.dir, start.year=1980, end.year=NA, nyr=NA){

    if(is.na(end.year) & is.na(nyr)){
        warning("Neither end.year nor nyr specified. Please specify one.")
        stop()
    }

    if(!is.na(nyr)){
         end.year <- start.year+nyr
    }

    if(!is.na(end.year)){
        nyr <- end.year-start.year
    }

    exploit.rate <- read.exploit.rates(model.dir)

    exploit.rate <- as_tibble(exploit.rate) %>%
                        pivot_longer(
                            everything(),
                            names_to="year",
                            values_to="exploit"
                        ) %>%
                        mutate(year=as.integer(year)) %>%
                        filter(year >= start.year & year <= end.year)

    return(exploit.rate)

}

compute.recruitment <- function(model.dir, start.year=1980, end.year=NA, nyr=NA){
    
    if(is.na(end.year) & is.na(nyr)){
        warning("Neither end.year nor nyr specified. Please specify one.")
        stop()
    }

    if(!is.na(nyr)){
         end.year <- start.year+nyr
    }

    if(!is.na(end.year)){
        nyr <- end.year-start.year
    }

    return(
        readr::read_csv(paste0(model.dir, "/mcmc_out/Age3.csv"), col_names = FALSE) %>% 
            select_at(((start.year-1980)+1):((end.year-1980+1))) %>%
            `colnames<-`(as.integer(start.year:end.year)) %>%
            pivot_longer(everything(), names_to="year", values_to="age3")
    )

}

ref.points.ts <- function(data){

    mu   <- data[, grepl("exp", colnames(data))]
    sel  <- data[, grepl("sel", colnames(data))]
    mort <- data[, grepl("mor", colnames(data))]
    waa  <- data[, grepl("waa", colnames(data))]
    fec  <- data[, grepl("fec", colnames(data))]
    R0   <- data[, grepl("R0",  colnames(data))]
    
    ypr   <- rep(NA, length(mu))
    sbpr  <- rep(NA, length(mu))
    sbpr0 <- rep(NA, length(mu))
    ssb0  <- rep(NA, length(mu))
    for(i in 1:length(mu)){
        ypr[i] <- calc.ypr(mu[i], sel[i,], mort[i,], waa[i,])
        sbpr[i] <- calc.sbpr(mu[i], sel[i, ], mort[i, ], fec[i, ])
        sbpr0[i] <- calc.sbpr0(sel[i, ], mort[i, ], fec[i, ])
        ssb0[i] <- calc.ssb0(sel[i, ], mort[i, ], fec[i, ], R0[i])
    }

    years <- data[, "year"]
    
    ypr.data   <- tibble(year=years, ypr=ypr)     %>% group_by(year) %>% median_qi(ypr,   .width=c(0.5, 0.95))
    sbpr.data  <- tibble(year=years, sbpr=sbpr)   %>% group_by(year) %>% median_qi(sbpr,  .width=c(0.5, 0.95))
    sbpr0.data <- tibble(year=years, sbpr0=sbpr0) %>% group_by(year) %>% median_qi(sbpr0, .width=c(0.5, 0.95))
    ssb0.data  <- tibble(year=years, ssb0=ssb0)   %>% group_by(year) %>% median_qi(ssb0,  .width=c(0.5, 0.95))
    spr.data   <- tibble(year=years, spr=sbpr/sbpr0) %>% group_by(year) %>% median_qi(spr,  .width=c(0.5, 0.95))

    return(list(
        ypr = ypr.data,
        sbpr = sbpr.data,
        sbpr0 = sbpr0.data,
        ssb0 = ssb0.data,
        spr = spr.data
    ))

}

start=1980
end=2033

exploit.rates <- compute.exploit.rate(model.dir, start.year=start, end.year=end)
selectivity <- compute.selectivity(model.dir, start.year=start, end.year=end)
mortality <- compute.mortality(model.dir, start.year=start, end.year=end)
waa <- compute.waa(model.dir, start.year=start, end.year=end, wide=TRUE)
fec <- compute.fecundity(model.dir, start.year=start, end.year=end)
R0 <- compute.recruitment(model.dir, start.year=start, end.year=end)

big <- mortality %>%
            bind_cols(
                selectivity$full %>% select(starts_with("age"))
            ) %>%
            bind_cols(
                waa %>% select(starts_with("age"))
            ) %>%
            bind_cols(
                fec %>% select(starts_with("age"))
            ) %>%
            bind_cols(exploit.rates$exploit) %>%
            bind_cols(R0$age3) %>%
            `colnames<-`(c("year", paste0("mort.",0:9), paste0("sel.",0:9), paste0("waa.",0:9), paste0("fec.",0:9), "exploit", "R0")) %>%
            as.matrix

ref.points <- ref.points.ts(big)

ggplot(ref.points$ypr, aes(x=year, y=ypr, ymin=.lower, ymax=.upper))+
  geom_lineribbon(size=0.1) +
  scale_fill_brewer(palette="Blues")+
  labs(x="Year", y="Yield-per-Recruit", title="Yield-per-Recruit")

ggplot(ref.points$sbpr, aes(x=year, y=sbpr, ymin=.lower, ymax=.upper))+
  geom_lineribbon(size=0.1) +
  scale_fill_brewer(palette="Blues")+
  labs(x="Year", y="Spawning-Biomass-per-Recruit", title="Spawning-Biomass-per-Recruit")

ggplot(ref.points$sbpr0, aes(x=year, y=sbpr0, ymin=.lower, ymax=.upper))+
  geom_lineribbon(size=0.1) +
  scale_fill_brewer(palette="Blues")+
  labs(x="Year", y="Unfished SBPR", title="Unfished SBPR")

ggplot(ref.points$ssb0, aes(x=year, y=ssb0, ymin=.lower, ymax=.upper))+
  geom_lineribbon(size=0.1) +
  scale_fill_brewer(palette="Blues")+
  labs(x="Year", y="Unfished Spawning Biomass", title="Unfished Spawning Biomass")

ggplot(ref.points$spr, aes(x=year, y=spr, ymin=.lower, ymax=.upper))+
  geom_lineribbon(size=0.1) +
  scale_fill_brewer(palette="Blues")+
  coord_cartesian(ylim=c(0, 1))+
  labs(x="Year", y="Spawner Potential Ratio", title="Spawner Potential Ratio")


selectivity.long <- function(data, name){
    return(
        data %>% 
            filter(year > 1983) %>%
            select(!year) %>%
            filter_all(any_vars(. != 0)) %>%
            pivot_longer(everything(), names_to="age", values_to="sel") %>%
            mutate(
                age=rep(0:9, length.out=nrow(.)),
                fishery=rep(name, length.out=nrow(.))
            )
    )
}

selectivity.conf <- selectivity.long(selectivity$full, "full") %>%
                        bind_rows(selectivity.long(selectivity$se, "se")) %>%
                        bind_rows(selectivity.long(selectivity$fb, "fb")) %>%
                        bind_rows(selectivity.long(selectivity$pk, "pk")) %>%
                        bind_rows(selectivity.long(selectivity$gn, "gn")) %>%
                        group_by(fishery, age) %>%
                        median_qi(sel, .width=c(0.5, 0.95)) %>%
                        print(n=10)
sel.full <- selectivity.long(selectivity$full, "full") %>% group_by(age) %>%
                        median_qi(sel, .width=c(0.5, 0.95))

ggplot(sel.full, aes(x=age, y=sel, ymin=.lower, ymax=.upper)) +
    geom_lineribbon(size=0.9) +
    geom_point() +
    geom_line() +
    scale_fill_brewer(palette="Blues")+
    scale_x_continuous(breaks=0:9)+
    labs(x="Age", y="Selectivity", title="Selectivity by Age")

ggplot(selectivity.conf, aes(x=age, y=sel, color=fishery, ymin=.lower, ymax=.upper)) +
    geom_lineribbon(size=0.9) +
    geom_point() +
    geom_line() +
    scale_fill_brewer(palette="Blues")+
    scale_x_continuous(breaks=0:9)+
    labs(x="Age", y="Selectivity", title="Selectivity by Age")+
    facet_wrap(~fishery)


mean.mortality <- mortality %>% 
                        select(!year) %>%
                        filter_all(any_vars(. != 0)) %>%
                        pivot_longer(everything(), names_to="age", values_to="mort") %>%
                        mutate(age=rep(0:9, length.out=nrow(.))) %>%
                        group_by(age) %>%
                        summarise(mort=mean(mort))

mean.selectivity <- selectivity %>% 
                        select(!year) %>%
                        filter_all(any_vars(. != 0)) %>%
                        pivot_longer(everything(), names_to="age", values_to="sel") %>%
                        mutate(age=rep(0:9, length.out=nrow(.))) %>%
                        group_by(age) %>%
                        summarise(sel=mean(sel))

mean.waa <- waa %>% 
                        select(!year) %>%
                        filter_all(any_vars(. != 0)) %>%
                        pivot_longer(everything(), names_to="age", values_to="waa") %>%
                        mutate(age=rep(0:9, length.out=nrow(.))) %>%
                        group_by(age) %>%
                        summarise(waa=mean(waa))

mean.fecundity <- fec %>% 
                        select(!year) %>%
                        filter_all(any_vars(. != 0)) %>%
                        pivot_longer(everything(), names_to="age", values_to="fec") %>%
                        mutate(age=rep(0:9, length.out=nrow(.))) %>%
                        group_by(age) %>%
                        summarise(fec=mean(fec))

## Yield-per-Recruit and Spawning-Bimoass-per-Recruit Curve
## using mean mortality, selectivity, waa
exploitation.rates <- seq(0, 1, 0.01)
yprs <- rep(NA, length(exploitation.rates))
sbprs <- rep(NA, length(exploitation.rates))
for(i in 1:length(exploitation.rates)){
    e <- exploitation.rates[i]
    yprs[i] <- calc.ypr(e, mean.selectivity$sel, mean.mortality$mort, mean.waa$waa)
    sbprs[i] <- calc.sbpr(e, mean.selectivity$sel, mean.mortality$mort, mean.fecundity$fec)
}

ref.points.df <- tibble(exploit=exploitation.rates, ypr=yprs, sbpr=sbprs)
sbpr0 <- calc.sbpr0(mean.selectivity$sel, mean.mortality$mort, mean.fecundity$fec)
sbpr.40 <- sbpr0*0.4
f.40 <- optim(c(0.6), min.f40)$par

ggplot(ref.points.df, aes(x=exploit, y=ypr))+
    geom_point()+
    geom_line()+
    scale_y_continuous(breaks=seq(0, 30, 5))+
    labs(x="Exploitation Rate", y="Yield-per-Recruit", title="Yield-per-Recruit Curve")

ggplot(ref.points.df, aes(x=exploit, y=sbpr))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=sbpr.40)+
    geom_vline(xintercept=f.40)+
    scale_y_continuous(breaks=seq(0, 30000, 5000))+
    coord_cartesian(ylim=c(0, 30000))+
    labs(x="Exploitation Rate", y="Spawning-Biomass-per-Recruit", title="Spawning-Biomass-per-Recruit Curve")

min.f40 <- function(e){
    sbpr <- calc.sbpr(e, mean.selectivity$sel, mean.mortality$mort, mean.fecundity$fec)
    return((sbpr-12003)^2)
}

## Unfished spawning biomass
## Using entire recruitment timeseries as recruitment is assumed 
## independent of spawning biomass.

r.eq <- 492.7

ssb0 <- calc.ssb0(mean.selectivity$sel, mean.mortality$mort, mean.fecundity$fec, r.eq)


