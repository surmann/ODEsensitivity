# Helper functions for plot.ODEmorris():

##### Auxiliary function: Plotting mu.star and sigma separately ###########
plotSep <- function(res, pars, state_name = NULL, colors_pars = NULL, 
                    common_title = NULL, legendPos, type = "l", ...) {
  t.vec <- res[1, ]
  k <- length(pars)
  if(is.null(colors_pars)){
    my.cols <- grDevices::rainbow(k)
  } else{
    my.cols <- colors_pars
  }
  if(is.null(common_title) && !is.null(state_name)){
    common_title <- paste0("Morris Screening for State Variable \"", state_name, "\"")
  } else if(is.null(common_title) && is.null(state_name)){
    common_title <- "Morris Screening"
  }
  if(legendPos == "outside"){
    if(k > 6){
      # Multirow legend:
      legend_ncol <- ceiling(k / 2)
      oldpar <- graphics::par(mfrow = c(1, 2),
                    oma = c(2.5, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
    } else{
      # One-row legend:
      legend_ncol <- k
      oldpar <- graphics::par(mfrow = c(1, 2),
                    oma = c(1.7, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
    }
  } else{
    oldpar <- graphics::par(mfrow = c(1, 2),
                  oma = c(0, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
  }
  
  # Plot mu.star:
  graphics::plot(t.vec, y = res[paste0("mu.star_", pars[1]), ], col = my.cols[1], 
       lwd = 1, type = type, xlab = "Time", ylab = expression(mu*"*"),
       ylim = range(res[paste0("mu.star_", pars), ], na.rm = TRUE), ...)
  if(k > 1) {
    for(i in 2:k) {
      graphics::lines(t.vec, y = res[paste0("mu.star_", pars[i]), ], col = my.cols[i], 
            lwd = 1, type = type, ...)
    }
  }
  if(legendPos != "outside"){
    graphics::legend(legendPos, legend = pars, lty = 1, col = my.cols)
  }
  
  # Plot sigma:
  graphics::plot(t.vec, y = res[paste0("sigma_", pars[1]), ], col = my.cols[1], 
       lwd = 1, type = type, xlab = "Time", ylab = expression(sigma), 
       ylim = range(res[paste0("sigma_", pars), ], na.rm = TRUE), ...)
  if(k > 1) {
    for(i in 2:k) {
      graphics::lines(t.vec, y = res[paste0("sigma_", pars[i]), ], col = my.cols[i], 
            lwd = 1, type = type, ...)
    }
  }
  if(legendPos != "outside"){
    graphics::legend(legendPos, legend = pars, lty = 1, col = my.cols)
  }
  
  # Create the big common title for the two plots:
  graphics::mtext(common_title, side = 3, line = 0, outer = TRUE, cex = 1.2, font = 2)
  
  # Legend outside of the plotting region:
  if(legendPos == "outside"){
    # Dummy plot:
    oldpar2 <- graphics::par(mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4), 
                   new = TRUE)
    graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    graphics::legend("bottom", legend = pars, lty = 1, col = my.cols, 
           bty = "n", xpd = TRUE, ncol = legend_ncol, inset = c(0, 0))
    graphics::par(oldpar2)
  }
  
  graphics::par(oldpar)
}

##### Auxiliary function: Plotting trajectories ######################
plotTrajectories <- function(res, pars, state_name = NULL, colors_pars = NULL, 
                             main_title = NULL, legendPos, type = "l", ...) {
  t.vec <- res[1, ]
  k <- length(pars)
  if(is.null(colors_pars)){
    my.cols <- grDevices::rainbow(k)
  } else{
    my.cols <- colors_pars
  }
  if(is.null(main_title) && !is.null(state_name)){
    main_title <- paste0("Morris Screening for State Variable \"", state_name, 
                         "\": Trajectories")
  } else if(is.null(main_title) && is.null(state_name)){
    main_title <- "Morris Screening: Trajectories"
  }
  if(legendPos == "outside"){
    if(k > 6){
      # Multirow legend:
      legend_ncol <- ceiling(k / 2)
      oldpar <- graphics::par(oma = c(2.5, 0, 0, 0))
    } else{
      # One-row legend:
      legend_ncol <- k
      oldpar <- graphics::par(oma = c(1.7, 0, 0, 0))
    }
  }
  
  # Plot the trajectories:
  graphics::plot(x = res[paste0("mu.star_", pars[1]), ], 
       y = res[paste0("sigma_", pars[1]), ], 
       col = my.cols[1], lwd = 1, type = type, main = main_title,
       xlim = range(res[paste0("mu.star_", pars), ], na.rm = TRUE),
       ylim = range(res[paste0("sigma_", pars), ], na.rm = TRUE),
       xlab = expression(mu*"*"), ylab = expression(sigma), ...)
  if(k > 1) {
    for(i in 2:k) {
      graphics::lines(x = res[paste0("mu.star_", pars[i]), ], 
            y = res[paste0("sigma_", pars[i]), ], 
            col = my.cols[i], lwd = 1, type = type, ...)
    }
  }
  if(legendPos == "outside"){
    # Dummy Plot:
    oldpar2 <- graphics::par(oma = rep(0, 4), mar = rep(0, 4), new = TRUE)
    graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    graphics::legend("bottom", legend = pars, lty = 1, col = my.cols, 
           bty = "n", xpd = TRUE, ncol = legend_ncol, inset = c(0, 0))
    
    graphics::par(oldpar2)
    graphics::par(oldpar)
  } else{
    graphics::legend(legendPos, legend = pars, lty = 1, col = my.cols)
  }
}
