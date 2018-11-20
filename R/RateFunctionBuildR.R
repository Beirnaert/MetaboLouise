#' RateFunctionBuildR
#'
#' Build and visualize an appropriate rate multiplier function. The flow between nodes in the network is governed by certain rates.
#' These rates can be made dependent on the values of the source node (flow from source node to receiver node). 
#' With this function you can visualize different types and receive the parameters to be submitted to the simulation function.
#' Note that this function returns an object containing the plot parameters. The simulation function can take only function type.
#' 
#' @param type Type of the function, any or multiple of "linear", "sigmoid" or "step" can be used.
#' @param C_range The range of node concentrations over which to plot the rate function. 
#' @param lin_xy_start Linear parameter: the starting point for the linear curve e.g. c(0,0). 
#' @param lin_slope Linear parameter: In case only a single point is provided, the slope is also necessary.
#' @param sig_C0 Sigmoid parameter: The midpoint concentration.
#' @param sig_k Sigmoid parameter: The curve steepness.
#' @param sig_max Sigmoid parameter: The maximal height of the function.
#' @param step_levels Step function parameter: the distinct levels of the function.
#' @param step_switchpoints Step function parameter: The points (nr. of levels minus 1) at which a switch is made.
#' @param plot.out Whether to plot the resulting functions.
#' @param Nplotpoints The number of plotpoints for the optional plot
#'
#' @return A plot with the visualize rate function and a list with the necessary parameters
#' 
#' @author Charlie Beirnaert, \email{charlie.beirnaert@@uantwerpen.be}
#'
#' @examples
#' RateFunctionBuildR()
#' 
#' @export
#' 
#' @importFrom graphics plot lines legend
#' 
RateFunctionBuildR <- function(type = c("linear", "sigmoid", "step"), C_range = c(0,200), 
                               lin_xy_start = c(0,0), lin_slope = 0.01, 
                               sig_C0 = 100, sig_k = 0.05, sig_max = 2,
                               step_levels = c(0,1,2), step_switchpoints = c(50, 150),
                               plot.out = TRUE, Nplotpoints = 100) {
    
    if(min(C_range) < 0){
        warning("In the simulation the concentration will never be negative.")
    }
    if(any(!type %in% c("linear", "sigmoid", "step"))){
        stop("Wrong function type selected.")
    }
    if("step" %in% type){
        if(length(step_switchpoints) != length(step_levels)-1){
            stop("Do not know how to generate the step function. The number of switchpoints must be 1 less than the number of levels.")
        }
        if(any(diff(step_switchpoints)<0)){
            stop("step_switchpoints not monotone increasing.")
            
        }
    }
    
    C_plotpoints <- seq(from = C_range[1], to = C_range[2], length.out = Nplotpoints)
    
    Rate_plotlist <- list()
    
    for(tt in 1:length(type)){
        
        if(type[tt] == "linear"){
            Rate_plot <- lin_xy_start[2] + lin_slope*(C_plotpoints-lin_xy_start[2])
        }
        if(type[tt] == "sigmoid"){
            Rate_plot <- sig_max / (1 + exp(-sig_k*(C_plotpoints-sig_C0)))
        }
        if(type[tt] == "step"){
            Rate_plot <- rep(step_levels[1], Nplotpoints)
            for(k in seq_along(step_switchpoints)){
                Rate_plot[C_plotpoints >= step_switchpoints[k]] <- step_levels[k+1]
            }
        }
        Rate_plotlist[[tt]] <- Rate_plot
        
    }
    
    if(plot.out){
        
        plot(x = C_plotpoints, 
             y = Rate_plot, 
             type = "n",
             ylim = c(min(unlist(Rate_plotlist)), max(unlist(Rate_plotlist))),
             xlab = "concentration",
             ylab = "rate multiplication factor",
             main = "Rate function")
        
        for(tt in 1:length(type)){
            lines(C_plotpoints, Rate_plotlist[[tt]], lty = tt)
        }
        
        legend(x = "topleft",
               legend = type,
               lty=1:length(type))
        
    }
    

    return(list(type = type, lin_xy_start = lin_xy_start, lin_slope = lin_slope, 
                sig_C0 = sig_C0, sig_k = sig_k, sig_max = sig_max,
                step_levels = step_levels, step_switchpoints = step_switchpoints ))
}