#' Function that orders all interaction pairs as L-R. Leaves unchanged the
#' R-R and L-L
#'
#' @param input.data uploaded data
#'
#' @return ordered input data
#' @importFrom dplyr arrange filter %>%
#'

updateInputLR <- function(input.data){
    RLint <- input.data %>%
        filter(typeA == "R" & typeB == "L")
    
    if(nrow(RLint)>0){
        otherINT <- input.data %>%
            filter(!(typeA == "R" & typeB == "L"))
        RLint.swap <- swap.RLint(RLint)
        updated.input <- rbind(otherINT, RLint.swap)
    } else{
        updated.input <- input.data
    }

    updated.input <- updated.input %>%
        arrange(int_pair)
    return(updated.input)
}


#' Swaps interaction pairs that are R-L to L-R
#'
#' @param RLint subset of R-L interactions
#'
#' @return input data with ordered L-R pairs and L-L/R-R
#'

swap.RLint <- function(RLint){
    RLint.swap <- RLint
    RLint.swap$int_pair <- paste0(sapply(strsplit(RLint$int_pair, " & "), 
                                         function(x) x[2]), " & ", 
                                  sapply(strsplit(RLint$int_pair, " & "), 
                                         function(x) x[1]))
    RLint.swap$geneA <- RLint$geneB
    RLint.swap$geneB <- RLint$geneA
    RLint.swap$typeA <- RLint$typeB
    RLint.swap$typeB <- RLint$typeA
    RLint.swap$clustA <- RLint$clustB
    RLint.swap$clustB <- RLint$clustA
    
    return(RLint.swap)                      
}
