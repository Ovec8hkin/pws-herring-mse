hcr.constant.escapement <- function(ssb.projected, naa.projected, 
                                    options=list(
                                        threshold = 20000,
                                        proportion = 1.0
                                    )){

    f = 0
    if(ssb.projected > options$threshold){
        f <- ((ssb.projected-options$threshold)*options$proportion)/ssb.projected
    }

    return(f)

}