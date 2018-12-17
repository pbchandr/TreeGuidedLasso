setParams <- function(opts) {

  if (is.list(opts) && length(opts)==0) {
    print('empty data frame');
    opts <- list()
    
    # Starting point
    opts[["init"]] <- 2 
    
    # Termination point
    opts[["tFlag"]] <- 6
    opts[["maxIter"]] <- 100
    opts[["tol"]]=1e-3;
    
    # Normalization
    opts[["nFlag"]] <- 0
    
    # Regularization
    opts[["rFlag"]] <- 1 # the input parameter 'rho' is a ratio in (0, 1)
    
    # Method (Line Search)
    opts[["mFlag"]] <- 0
    opts[["lFlag"]] <- 0
    
  }
  else {
  
    ## Starting point
      if ("init" %in% names(opts)) {
        if ((opts[["init"]]!=0) && (opts[["init"]]!=1) && (opts[["init"]]!=2)){
          opts[["init"]] <- 0 # if .init is not 0, 1, or 2, then use the default 0
        }
        if ((!("x0" %in% opts)) && (opts[["init"]]==1)){
          opts[["init"]] <- 0 # if .x0 is not defined and .init=1, set .init=0
        }
      } else {
        opts[["init"]] <- 0 # if .init is not specified, use "0"
      }
      # print(opts[["init"]])
  
    
    ## Termination
      # Maximum Iteration
      if ("maxIter" %in% names(opts)) {
        if (opts[["maxIter"]]<1) {
          opts[["maxIter"]] <- 10000
        }
      } else {
        opts[["maxIter"]] <- 10000
      }
      # print(opts[["maxIter"]])
      
      # Tolerance
      if (!("tol" %in% names(opts))) {
        opts[["tol"]] <- 1e-3
      }

      # Flag for termination
      if ("tFlag" %in% names(opts)) {
        if (opts[["tFlag"]]< 1){
          opts[["tFlag"]] <- 1
        } else if (opts[["tFlag"]]> 6) {
          opts[["tFlag"]] <- 6
        } else {
          opts[["tFlag"]] <- floor(opts[["tFlag"]]);
        }
      }
      else {
        opts[["tFlag"]] <- 1
      }
  
  
    ## Normalization
    if ("nFlag" %in% names(opts)) {
      if ((opts[["nFlag"]]!=1) && (opts[["nFlag"]]!=2)) {
        opts[["nFlag"]] <- 0
      } 
    } else {
      opts[["nFlag"]] <- 0
    }
    # print(opts[["nFlag"]])
    
    ## Regularization
    if ("rFlag" %in% names(opts)) {
      if (opts[["rFlag"]]!=1) {
        opts[["rFlag"]] <- 0
      }
    } else {
      opts[["rFlag"]] <- 0
    }
    # print(opts[["rFlag"]])
    
    ## Method (Line Search)
    if ("lFlag" %in% names(opts)) {
      if (opts[["lFlag"]]!=1) {
        opts[["lFlag"]] <- 0
      }
    } else {
      opts[["lFlag"]] <- 0
    }
  
    if ("mFlag" %in% names(opts)) {
      if (opts[["mFlag"]]!=1){
        opts[["mFlag"]] <- 0
      }
    } else {
      opts[["mFlag"]] <- 0
    }
  }
  
  # print(opts)
  return(opts)
}
