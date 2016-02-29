# Subroutine for processing the results of soboljansen_list() or 
# sobolmartinez_list():

sobol_process <- function(x, pars, times){
  stopifnot(class(x) %in% c("soboljansen", "sobolmartinez"))
  
  k <- length(pars)
  
  ST_by_state <- lapply(1:dim(x$y)[3], function(i){
    S <- rbind(times, x$S[, , i])
    T <- rbind(times, x$T[, , i])
    rownames(S) <- rownames(T) <- c("time", pars)
    return(list(S = S, T = T))
  })
  names(ST_by_state) <- dimnames(x$y)[[3]]
  
  # Handling of invalid SA indices:
  # "minor" problem with first order indices: -0.05 <= first order index < 0,
  # "major" problem with first order indices: first order index < -0.05,
  # "minor" problem with total indices: 1 < total index <= 1.05,
  # "major" problem with total indices: 1.05 < total index
  # First order indices > 1 or total indices < 0 shouldn't occur. If they do,
  # further investigation will be needed.
  
  check_indices <- sapply(ST_by_state, function(ST){
    S <- ST$S[-1, ]
    T <- ST$T[-1, ]
    negative_S <- c(minor = any(-0.05 <= S & S < 0), 
                    major = any(S < -0.05))
    greater1_S <- any(1 < S)
    negative_T <- any(T < 0)
    greater1_T <- c(minor = any(1 < T & T <= 1.05), 
                    major = any(1.05 < T))
    return(c(neg_S = negative_S, 
             big_S = greater1_S,
             neg_T = negative_T,
             big_T = greater1_T))
  })
  
  if(any(check_indices["neg_S.minor", ])){
    # "Repair" minor negative first order indices by setting them to 0:
    for(i in seq_along(ST_by_state)[check_indices["neg_S.minor", ]]){
      check_neg_S.minor <- -0.05 <= ST_by_state[[i]]$S & ST_by_state[[i]]$S < 0
      ST_by_state[[i]]$S[check_neg_S.minor] <- 0
    }
  }
  
  if(any(check_indices["big_T.minor", ])){
    # "Repair" total indices being slightly bigger than 1 by setting them to 1:
    for(i in seq_along(ST_by_state)[check_indices["big_T.minor", ]]){
      check_big_T.minor <- 1 < ST_by_state[[i]]$T & ST_by_state[[i]]$T <= 1.05
      ST_by_state[[i]]$T[check_big_T.minor] <- 1
    }
  }
  
  if(any(check_indices[c("neg_S.major", "big_T.major"), ])){
    warning("Negative first order indices (< -0.05) and/or total indices ",
            "> 1.05 detected. Argument \"n\" might be too low. If using a ",
            "higher value for \"n\" does not help, please check if the ",
            "parameter distributions (\"rfuncs\") and their arguments ", 
            "(\"rargs\") generate valid parameter values.",
            call. = FALSE)
  }
  
  if(any(check_indices[c("big_S", "neg_T"), ])){
    warning("First order indices > 1 and/or negative total indices detected. ",
            "This shouldn't happen. Please contact the package author.",
            call. = FALSE)
  }
  
  return(ST_by_state)
}
