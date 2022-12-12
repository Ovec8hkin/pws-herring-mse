library(tidyverse)
source(file=paste0(here::here("R/utils/"), "fun_read_dat.R"))

model.dir <- paste0(here::here("results/base/sim_197/year_0/model/"))

years <- 1980:2021
nyr <- length(years)

biomass.est <- read.biomass.estimates(model.dir, nyr=42)

annual.biomass <- as_tibble(biomass.est) %>%
                    pivot_longer(everything(), names_to="year", values_to="biomass") %>%
                    group_by(year) %>%
                    summarise(biomass=median(biomass)) %>%
                    print(n=10)
biomass <- annual.biomass$biomass

three.year.rel.change <- n.year.rel.biomass(biomass, n=3)
evenness <- c(NA, NA, as.vector(even.50)) # run age_struct_evenness.R for the even.50 object

n.year.rel.biomass <- function(biomass, n=3){
    rel.biomass.change <- rep(NA, length(biomass)-n+1)
    for(i in n:length(biomass)){
        rel.biomass.change[i-n+1] = biomass[i]/biomass[i-n+1]
    }
    return(rel.biomass.change)
}

compute.f.thresh <- function(biomass, b.lim=19956, b.tar=38556, f.max=0.2, scaling="linear"){
    f <- 0
    if(biomass <= b.lim){
        f <- 0
    }else if(biomass > b.lim && biomass <= b.tar){
        if(scaling == "linear"){
            f <- (biomass - b.lim) * f.max / (b.tar - b.lim)
        }else if(scaling == "logistic"){
            f <- f.max/(1+exp(-0.0007*(biomass - (b.lim+((b.tar-b.lim)/2)))))
        }else if(scaling == "exponential"){
            f <- (f.max*2)/(1+exp(-0.0004*(biomass - b.tar)))
        }else if(scaling == "logarithmic"){
            f <- (f.max*2)/(1+exp(-0.0004*(biomass - b.lim)))-f.max
        }else{
            f <- NA
        }
    }else{
        f <- f.max
    }

    return(f)
}

compute.f.gradient <- function(biomass, rel.ssb.change, b.lim=19958){
    f.def <- compute.f.thresh(biomass)
    f.grad <- f.def*rel.ssb.change
    p=0.5
    return(p*f.def + (1-p)*f.grad)
}

compute.f.age.struct <- function(biomass, evenness, b.lim=19956){
    f <- compute.f.thresh(biomass)
    if(is.na(evenness)) 
        return(NA)
    e <- compute.f.thresh(evenness, b.lim=0.4, b.tar=0.8, f.max=0.5)

    return(f*(e + 0.5))
}

fs.default <- apply(as.matrix(biomass), 1, compute.f.thresh)
fs.logistic     <- apply(matrix(biomass), 1, compute.f.thresh, scaling="logistic")
fs.exponential  <- apply(matrix(biomass), 1, compute.f.thresh, scaling="exponential")
fs.logarithmic  <- apply(matrix(biomass), 1, compute.f.thresh, scaling="logarithmic")
fs.low.thresh   <- apply(matrix(biomass), 1, compute.f.thresh, b.lim=10000)
fs.high.thresh  <- apply(matrix(biomass), 1, compute.f.thresh, b.lim=30000)
fs.high.f       <- apply(matrix(biomass), 1, compute.f.thresh, f.max=0.4)
fs.low.f        <- apply(matrix(biomass), 1, compute.f.thresh, f.max=0.1)
fs.gradient     <- apply(matrix(3:nyr), 1, function(i) compute.f.gradient(biomass=biomass[i], rel.ssb.change=three.year.rel.change[i-2]))
fs.evenness     <- apply(matrix(1:nyr), 1, function(i) compute.f.age.struct(biomass=biomass[i], evenness=evenness[i]))

fs.df <- data.frame(year=years, 
                    default=fs.default, 
                    logistic=fs.logistic,
                    exponential=fs.exponential,
                    logarithmic=fs.logarithmic,
                    low.thresh=fs.low.thresh,
                    high.thresh=fs.high.thresh,
                    low.f=fs.low.f,
                    high.f=fs.high.f,
                    gradient=c(NA, NA, fs.gradient),
                    age.struct=fs.evenness
                    )
round(fs.df, 2)

fs.df.relative <- fs.df/fs.df$default
fs.df.relative[is.na(fs.df.relative) | (fs.df.relative) > 1e5] <- 1.0
fs.df.relative$year <- years


fs.df.long <- as_tibble(fs.df) %>%
                pivot_longer(-year, names_to="control.rule", values_to="harvest.rate")
fs.df.long$control.rule <- factor(fs.df.long$control.rule, 
                                  levels=c("default", "logistic", "exponential", "logarithmic", "low.thresh", "high.thresh", "low.f", "high.f", "age.struct", "gradient"),
                                  labels=c("Default", "Logistic", "Exponential", "Logarithmic", "Low Threshold", "High Threshold", "Low F", "High F", "Age Structured", "Gradient")
                                 )

fs.df.relative.long <- as_tibble(fs.df.relative) %>%
                pivot_longer(-year, names_to="control.rule", values_to="harvest.rate")
fs.df.relative.long$control.rule <- factor(fs.df.relative.long$control.rule, 
                                  levels=c("default", "logistic", "exponential", "logarithmic", "low.thresh", "high.thresh", "low.f", "high.f", "age.struct", "gradient"),
                                  labels=c("Default", "Logistic", "Exponential", "Logarithmic", "Low Threshold", "High Threshold", "Low F", "High F", "Age Structured", "Gradient")
                                 )


p1 <- ggplot(fs.df.long, aes(x=year, y=harvest.rate, color=control.rule)) +
    geom_line(size=1)+
    scale_y_continuous(limits=c(0, 1.0), expand=c(0, 0))+
    scale_x_continuous(expand=c(0, 0))+
    labs(title="Annual Harvest Rates", y="Harvest Rate", x="Year", color="Control Rule")

p2 <- ggplot(fs.df.relative.long, aes(x=year, y=harvest.rate, color=control.rule)) +
    geom_line(size=1)+
    scale_y_continuous(limits=c(0, 10), expand=c(0, 0))+
    scale_x_continuous(expand=c(0, 0))+
    labs(title="Annual Harvest Rates Relative to Default Control Rule", y="Relative Harvest Rate", x="Year", color="Control Rule")

library(ggpubr)

ggarrange(p1, p2, common.legend = TRUE, legend="bottom")