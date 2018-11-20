#' NetworkCreateR 
#'
#' This function creates a fictional metabolomics network according to the Barabasi-Albert model. 
#'
#' @param N the number of nodes in the network (metabolites)
#' @param BA_power (igrap parameter) The power of the prefferential attachment. 
#' @param BA_mValue (igrap parameter) The number of edges to add in each time step
#' @param plothist Whether to plot a histogram of the network's connectivity distribution
#' @param histbreaks The number of breaks in the histogram 
#' @param ... Optional additional arguments to be passed along to the network generator funcion \code{\link[igraph]{sample_pa}}. 
#'
#' @return A matrix with the network structure
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' Network <- NetworkCreateR(N = 50, BA_power = 0.5, BA_mValue = 4)
#' image(Network)
#' 
#' @export
#' 
#' @importFrom igraph sample_pa
#' @importFrom graphics hist
#' 
NetworkCreateR <- function(N, BA_power = 0.5, BA_mValue = 4, plothist = TRUE, histbreaks = 50, ...) {
    
    BA <- sample_pa(n = N, power = BA_power, m = BA_mValue, directed = TRUE, ...)
    NetworkMatrix_original <- as.matrix(BA[])
    conns <- which(NetworkMatrix_original == 1, arr.ind = T)
    shuffle <- sample(1:nrow(conns),  round(nrow(conns)/2))
    conns_shuffled <- rbind(conns[shuffle,c(2,1)],conns[-shuffle,])
    
    NetworkMatrix_shuffled <- matrix(0, nrow = nrow(NetworkMatrix_original), ncol = ncol(NetworkMatrix_original))
    NetworkMatrix_shuffled[conns_shuffled] <- 1
    
    graph_connections_shuffled <- rep(0,N)
    for(m in 1:length(graph_connections_shuffled)){
        graph_connections_shuffled[m] <- sum(NetworkMatrix_shuffled[m,]) + sum(NetworkMatrix_shuffled[,m])
    }
    
    if(plothist){
        hist(graph_connections_shuffled, breaks = histbreaks)
    }
    
    NetworkMatrix <- NetworkMatrix_shuffled
    
    return(NetworkMatrix)
    
}