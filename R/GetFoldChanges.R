#' GetFoldChanges
#'
#' Calculate the fold change equivalent for the final state of two simulated datasets e.g. with vs without influx.
#' 
#' @param ReferenceDataObject The reference data object (result from DataSimulateR function)
#' @param AlternativeDataObject The alternative situation data object (result from DataSimulateR function)
#' @param plot_out Whether to plot the Fold change distribution.
#' @param bw The bw for the density plot.
#' @param plot_title Optional plot title.
#'
#' @return A list with: the time vector and a matrix with the simulated data. (1 row per node)
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' Nmetabos <- 20L
#' Nrates <- 10L
#' 
#' Network <- NetworkCreateR(N = Nmetabos, BA_power = 0.5, BA_mValue = 4)
#' 
#' Rate_function <- RateFunctionBuildR(type = "sigmoid")
#' 
#' rate_vector <- round(5*runif(Nrates))
#' 
#' rate_mapping <- Network
#' active_rates <- which(Network == 1, arr.ind = TRUE)
#' for(rr in 1:nrow(active_rates)){
#'     rate_mapping[active_rates[rr,1], active_rates[rr,2]] <- sample(seq_along(rate_vector), size = 1)
#'     }
#'     
#'  No_influx <- DataSimulateR(NetworkMatrix = Network, dT = 0.01, Tstart = 0, Tstop = 3, 
#'                             T0_nodes = 100, rate_vector = rate_vector, rate_mapping = rate_mapping, 
#'                             RateFunctionObject = Rate_function, plot_out = FALSE)
#'                             
#'  influx_vector <- c(rep(1,10),rep(0,Nmetabos-10))                           
#'  With_influx <- DataSimulateR(NetworkMatrix = Network, dT = 0.01, Tstart = 0, Tstop = 3,
#'                               T0_nodes = 100, influx_vector = influx_vector, influx_Tframe = 0.5,
#'                               rate_vector = rate_vector, rate_mapping = rate_mapping, 
#'                               RateFunctionObject = Rate_function, plot_out = FALSE)
#' 
#' GetFoldChanges(ReferenceDataObject = No_influx, AlternativeDataObject = With_influx)
#' 
#' @export
#' 
#' @importFrom graphics plot 
#' @importFrom stats density
#' 
GetFoldChanges <- function(ReferenceDataObject, AlternativeDataObject, plot_out = TRUE, bw = 0.05, plot_title = NULL){
    
    if(length(AlternativeDataObject$time) != length(ReferenceDataObject$time)){
        warning("The data objects differ in simulation time length")
    }
    
    if(nrow(AlternativeDataObject$data) != nrow(ReferenceDataObject$data)){
        stop("Unequal number of nodes in the data objects. ")
    }
    
    FoldChanges <- AlternativeDataObject$data[,length(AlternativeDataObject$time)]/ReferenceDataObject$data[,length(ReferenceDataObject$time)]
    
    if(plot_out){
        plot(density(FoldChanges, bw = bw), main = plot_title)
    }
    
    return(FoldChanges)
}