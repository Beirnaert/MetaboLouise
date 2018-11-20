#' DataSimulateR
#'
#' Simulated dynamic/longitudinal data based on an underlying network. The network is initialized with values 
#' for every node (e.g. concentrations in the case of metabolites). These values evolve over time caused by the
#' (variable) rates. 
#' 
#' @param NetworkMatrix The underlying network in matrix form (e.g. form the NetworkCreateR function).
#' @param dT Simulation time step (must be small enough to avoid approximation errors).
#' @param Tstart Starting time point.
#' @param Tstop Ending time point.
#' @param T0_nodes Vector (or single value) with the initial node value(s) (metabolite concentrations).
#' @param influx_vector Vector with the influx (per time unit) received by the corresponding metabolites.
#' @param influx_Tframe Vector of two values indicating the start and ending time of the influx (if only single ending time value supplied, assumption of influx start = Tstart is made). This can also be a data frame/matrix with 2 columns and 1 row per metabolite, or just a single ending value.       
#' @param rate_vector Vector with the initial rates. Number of rates can be between 1 and the number of edges in the network.
#' @param rate_mapping A matrix (same size as NetworkMatrix) with for every link/edge in the network (1 in NetworkMatrix) an index of rate_vector to be matched.
#' @param RateFunctionObject An object from RateFunctionBuildR describing the evolution of rates. If not provided, the rates are constant.
#' @param plot_out Whether to plot the simulated data.   
#' @param plot_title  Optional plot title.
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
#'                             RateFunctionObject = Rate_function, plot_out = TRUE)
#' 
#' @export
#' 
#' @importFrom graphics plot lines
#' 
DataSimulateR <- function(NetworkMatrix, dT, Tstart, Tstop, T0_nodes = 100, influx_vector = NULL,influx_Tframe = NULL, 
                          rate_vector = NULL, rate_mapping = NULL, RateFunctionObject = NULL, plot_out = TRUE, plot_title = NULL){
    
    if(dim(NetworkMatrix)[1] != dim(NetworkMatrix)[2]){
        stop("NetworkMatrix is not a square matrix.")
    }
    N <- dim(NetworkMatrix)[1]
    
    if(max(NetworkMatrix) > 1){
        warning("NetworkMatrix has elements larger than 1. This has no meaning and these are converted to 1.")
        NetworkMatrix[NetworkMatrix > 1] = 1
    }
    
    if(length(T0_nodes) != N){
        T0_nodes <- rep(T0_nodes[1], N)
    }
    
    if(length(T0_nodes) != 1 & length(T0_nodes) != N){
        stop("Initial node values not set correctly.")
    }
    if(length(T0_nodes) == 1){
        T0_nodes <- rep(T0_nodes, N)
    }
    
    if(is.null(rate_vector)){
        message("rate_vector not supplied. Assuming constant rates of 1.")
        rate_vector <- 1
    }
    if(is.null(rate_mapping)){
        message(paste("No rate_mapping supplied. Assuming all rates map to the first element of rate_vector, with a rate of ",rate_vector[1], sep = ""))
        rate_mapping <- NetworkMatrix
    }else{
        if(any(dim(rate_mapping) != dim(NetworkMatrix))){
            stop("rate_mapping is incorrect.")
        }
        if(max(rate_mapping) > length(rate_vector)){
            stop("rate_mapping is incorrect. There are larger indices in rate_mapping than there are elements in rate_vector")
        }
        
    } 
    if(is.null(influx_vector) | is.null(influx_Tframe)){
        Influx <- FALSE
    }else{
        if(is.null(nrow(influx_Tframe)) & length(influx_Tframe) == 1){
            influx_Tframe <- c(Tstart,influx_Tframe)
        }
        if(is.null(nrow(influx_Tframe))){
            influx_Tframe <- matrix(rep(influx_Tframe, N), ncol = 2, nrow = N, byrow = T)
        }
        if(nrow(influx_Tframe) != N){
            stop("The influx_Tframe is not properly defined. The amount of rows is not equal to the number of nodes")
        }
        
        Influx <- TRUE
        if(length(influx_vector) != N){
            stop("The length of the influx vector does not match with the number of nodes")
        }
        influx_nodes <- which(influx_vector > 0)
        Influx_df <- data.frame(node = influx_nodes, 
                                influx = influx_vector[influx_nodes],
                                influx_start = influx_Tframe[influx_nodes,1],
                                influx_end = influx_Tframe[influx_nodes,2])
    }
    
    # rate equation
    node_evolution <-  function(C, r) {
        # basic euler rate equation
        return(r*C )
    }

    # Constructing the Rate evolution function
    if(is.null(RateFunctionObject)){
        Rate_evolution <-  function(C) {
            return(1)
        }
    }else{
        if("linear" %in% RateFunctionObject$type){
            Rate_evolution <-  function(C) {
                Rfactor <- RateFunctionObject$lin_xy_start[2] + RateFunctionObject$lin_slope*(C-RateFunctionObject$lin_xy_start[2])
                return(Rfactor)
            }
        }else if("sigmoid" %in% RateFunctionObject$type){
            Rate_evolution <-  function(C) {
                Rfactor <- RateFunctionObject$sig_max / (1 + exp(-RateFunctionObject$sig_k*(C-RateFunctionObject$sig_C0)))
                return(Rfactor)
            }
        }else if("step" %in% RateFunctionObject$type){
            Rate_evolution <-  function(C) {
                Clarger <- which(RateFunctionObject$step_switchpoints > C)[1]
                if(!is.na(Clarger)){
                    Rfactor <- RateFunctionObject$step_levels[Clarger]
                }else{
                    Rfactor <- rev(RateFunctionObject$step_levels)[1]
                }
                return(Rfactor)
            }
        }else{
            warning("No correct type found in RateFunctionObject$type. Assuming no variable rates.")
            Rate_evolution <-  function(C) {
                return(1)
            }
        }
    }

    # simulation settings
    
    t_vector <- seq(Tstart, Tstop, by = dT)
    Tpoints <- length(t_vector)
    
    SimulatedData <- matrix(0, ncol = Tpoints, nrow = N) # population size
    if(length(t_vector) > ncol(SimulatedData) ){
        t_vector <- t_vector[1:ncol(SimulatedData)]
    }
    
    # initial conditions
    SimulatedData[,1] <- T0_nodes
    
    NetworkConnections <- which(NetworkMatrix == 1, arr.ind = TRUE)

    for (i in 2:Tpoints) {
        SimulatedData[,i] <- SimulatedData[,(i-1)]
        for (p in 1:nrow(NetworkConnections)){
            # update rates
            curr_node_from <- NetworkConnections[p,1] # from
            curr_node_to <- NetworkConnections[p,2] # to
            
            rate_vector[rate_mapping[curr_node_from,curr_node_to]] <- rate_vector[rate_mapping[curr_node_from,curr_node_to]] * ( 1 + dT*( Rate_evolution(C = SimulatedData[curr_node_from,i]) -1)  )
            
            dN <- node_evolution(C = SimulatedData[curr_node_from,(i)],
                                 r = rate_vector[rate_mapping[curr_node_from,curr_node_to]]) * dT
            
            if(dN > SimulatedData[curr_node_from,(i)]){
                warning("rate or time step too large")
                dN <- SimulatedData[curr_node_from,(i)]
                SimulatedData[curr_node_from,i] <- 0
                SimulatedData[curr_node_to,i] <- dN + SimulatedData[curr_node_to,(i)]
            } else{
                SimulatedData[curr_node_from,i] <- - dN + SimulatedData[curr_node_from,(i)]
                SimulatedData[curr_node_to,i] <- dN + SimulatedData[curr_node_to,(i)]
            }
            
        }
        
        if(Influx){
            for(fl in 1:nrow(Influx_df)){
                if(t_vector[i] >= Influx_df$influx_start[fl] & t_vector[i] <= Influx_df$influx_end[fl]){
                    SimulatedData[Influx_df$node[fl],i] <- SimulatedData[Influx_df$node[fl],i] * (1 + dT * Influx_df$influx[fl])
                }
            }
        }
    }

    if(plot_out){
        if(is.null(influx_vector)){
            plot( t_vector, SimulatedData[1,], 
                  type = "n", 
                  ylim = c(min(SimulatedData), max(SimulatedData)), 
                  main = plot_title )
            for(pl in 1:N){
                lines(t_vector, SimulatedData[pl,], lty = pl)
            }
        }else{
            plot( t_vector, SimulatedData[1,], 
                  type = "n", 
                  ylim = c(min(SimulatedData), max(SimulatedData)), 
                  main = plot_title )
            for(pl in 1:N){
                if(influx_vector[pl] >0){
                    lines(t_vector, SimulatedData[pl,], lty = pl, col = "red")
                } else{
                    lines(t_vector, SimulatedData[pl,], lty = pl, col = "black")
                }
            }
        }
    }
    
    return(list(time = t_vector, data = SimulatedData))
}

