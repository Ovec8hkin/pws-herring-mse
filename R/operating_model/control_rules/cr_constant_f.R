hcr.constant.f <- function(ssb.projected, naa.projected, 
                           options=list(
                               f.rate = 0.0
                           )){

    return(ssb.projected*options$f.rate)                

}