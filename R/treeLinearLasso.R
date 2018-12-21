treeLinearLasso <- function(A, y, z, opts) {
  
  valueL <- c()
  funVal <- c()
  
# Verify and initialize the parameters ---------------------------------------------------------------------------
  
  # Verify the number of input parameters
    # if (missing(A) | missing(y) | missing(z)){
    if (nargs() < 3) { 
      stop("Inputs: A, y and z should be specified!")
    } else if (nargs() == 3) {
      opts <- list()
    }
  
  # Get the dimension of the matrix A
    m <- nrow(A)
    n <- ncol(A)
  
  # Verify the length of y
    if (length(y) !=m) {
      stop('\n Check the length of y!\n')
    }
  
  # Verify the value of z
    if (z <= 0) {
      stop('\n z should be positive!\n')
    }
  
  # Run sll_opts to set default values (flags)
    opts <- setParams(opts) # Run setParams to set default values (flags)
    
  
  # Restart the program for better efficiency. This is a newly added function
    if (!("rStartNum" %in% names(opts))){
      opts[["rStartNum"]] <- opts$maxIter
    } else if (opts[["rStartNum"]] <= 0){
      opts[["rStartNum"]] <- opts$maxIter
    }
  

# Normalization---------------------------------------------------------------------------------------------------
    if (opts[["nFlag"]] !=0) {
      if ("mu" %in% names(opts)) {
        mu <- opts[["mu"]]
        if(ncol(mu) !=n) {
          stop('\n Check the input .mu')
        }
      } else {
        mu <- colMeans(A);
      }
      
      if (opts[["nFlag"]]==1) {
        if ("nu" %in% names(opts)) {
          nu <- opts[["nu"]]
          if(ncol(nu) !=n) {
            stop('\n Check the input .nu!')
          }
        } else {
          nu <- (colSums(A^2)/m)^(0.5)
          nu <- t(nu)
        }
      } else {
        if ("nu" %in% names(opts)) {
          nu <- opts[["nu"]]
          if(ncol(nu) !=m) {
            stop('\n Check the input .nu!')
          } 
        } else {
          nu <- (colSums(A^2)/n)^(0.5)
        }
      }
      
      ind_zero <- find(abs(nu)<= 1e-10)    
      nu[ind_zero] <- 1
    }
  
  
# Group and Others -----------------------------------------------------------------------------------------------
  # Initialize ind
    if(!("ind" %in% names(opts))){
      stop('\n In tree_LeastR, the field .ind should be specified')
    } else {
      ind <- opts[["ind"]]
      if(nrow(ind) != 3){
        stop("Check opts.ind")
      }
    }
  
  # Setting GFlag
    GFlag <- 0
    # If GFlag = 1, We will use general_altra c function
    if("G" %in% names(opts)){
      GFlag <- 1
      
      G <- opts[["G"]]
      if(max(G) > n || max(G) < 1){
        stop(paste0('The input G is incorrect. It should be within ',1,' and ',n))
      }
    }
    
  
# Starting Point Initialization ----------------------------------------------------------------------------------
  # compute AT y
    if (opts[["nFlag"]] ==0) {
      ATy <- t(A)%*%y
    } else if(opts[["nFlag"]] ==1){
      ATy <- t(A)%*%y - sum(y)*t(mu)
      ATy <- ATy./nu
    } else {
      invNu <- y./nu
      ATy <- t(A)%*%invNu - sum(invNu)*t(mu);
    }
    
  # Process the regularization parameter
    if(opts[["rFlag"]] == 0) {
      lambda <- z
    } else {
      computedFlag <- 0
      if("lambdaMax" %in% names(opts) && opts[["lambdaMax"]] != -1) {
        lambda <- z * opts$lambdaMax
        computedFlag<- 1
      }
      
      if(computedFlag == 0) {
        if(GFlag == 0){
          lambda_max <- findLambdaMax(ATy, n, ind, ncol(ind))
        }
        else {
          
          lambda_max <- general_findLambdaMax(ATy, n, ind, ncol(ind))
        }
        lambda <- z*lambda_max
      }
    }
    
    
    # The following is for computing lambdaMax. we use opts.lambdaMax=-1 to show that we need the computation. 
      # One can use this for setting up opts.lambdaMax
    if("lambdaMax" %in% names(opts) && opts[["lambdaMax"]] != -1) {
      if(GFlag == 0){
        lambda_max <- findLambdaMax(ATy, n, ind, ncol(ind))
      }
      else {
        
        lambda_max <- general_findLambdaMax(ATy, n, ind, ncol(ind))
      }
      
      x <-lambda_max
      funVal <- lambda_max
      valueL <- lambda_max
      
      TreeLeastR_out <- list("beta" = x, "valueL" = valueL, "funVal" <- funVal)
      return(TreeLeastR_out)
    }
    
    
    # Initialize a starting point
    if (opts[["init"]] == 2) {
      x <- matrix(0, n, 1)  
    } else {
      if ("x0" %in% names(opts)) {
        x <- opts[["x0"]]
        if(length(x) != n) {
          stop("Check input x0")
        }
      } else {
        x <- ATy
      }
    }
    
    # Compute Ax
    if(opts[["nFlag"]] == 0) {
      Ax <- A%*%x
    } else if(opts[["nFlag"]] == 1) {
      invNu <- x./nu
      mu_invNu <- mu * invNu
      
      Ax <- A*invNu - rep(mu_invNu, m)
    } else {
      Ax <- A%*%x - rep(mu*x, m)
      Ax <- Ax./nu
    }
    
    # if(opts[["init"]] == 0) {
    #   x <- matrix(0, n, 1)
    # }
    
  
# The Armijo Goldstein line search scheme + accelearted gradient descent -----------------------------------------
  bFlag <- 0 # this flag tests whether the gradient step only changes a little
  L <- 1  # We assume that the maximum eigenvalue of A'A is over 1
  
  # Assign xp with x and Axp with xp
  xp <-x
  Axp <- Ax
  xxp <- matrix(0, ncol(A),1)
  
  # alphap and alpha are used for computing the weight in forming search point
  alphap <- 0
  alpha <- 1
  
  
  ind_work <- matrix(0, nrow(ind), ncol(ind))
  
  for(iterStep in 1: opts$maxIter){
    # Step 1 - Compute search point based on xp and x (with beta)
    beta <- (alphap - 1)/alpha
    s <- x + beta * xxp
    
    # Step 2 - line search for : and compute the new approximate solution x
    # Compute gradient (g) at s
    As <- Ax + beta * (Ax - Axp)
    
    # compute ATAs
    ATAs <- crossprod(A, As)
    
    # Obtain the gradient g
    g <- (ATAs - ATy)
    
    # Assign x and Ax to xp and Axp
    xp <- x
    Axp <- Ax
    
    while(TRUE) {
      # let s walk in a step in the antigradient of s to get v and then do the l1-norm regularized projection
      v <- s - g/L
      
      # Tree overlapping group Lasso projection
      ind_work[1:2,] <- ind[1:2,]
      ind_work[3,] <- ind[3, ]* (lambda/L)
      
      if(GFlag == 0){
        x <- altra(v, n, ind_work, ncol(ind_work))
      }else {
        x <- general_altra(v, n, G, ind_work, ncol(ind_work))
      }
      
      # Difference between the new approximate solution x and the search point s
      v <- x - s
      
      # Compute Ax
      Ax <- tcrossprod(A, t(x))
      Av <- Ax - As
      rSum <- crossprod(v) #t(v)%*%v
      lSum <- crossprod(Av) #t(Av)%*%Av
      
      if(rSum <= 1e-20){
        bFlag <- 1 # Shows the gradient step makes little improvement
        break
      }
      
      # the condition is ||Av||_2^2 <= L * ||v||_2^2
      # cat(iterStep, lSum, rSum*L, '\n')
      # cat(typeof(lSum), typeof(rsum_L), '\n')
      if(((rSum*L) - lSum) >= 1e-5){
        break
      } 
      else {
        L <- max(2*L, lSum/rSum)
        
      }
    }
    valueL <- c(valueL, L)
    
    # Step 3 - Update alpha and alphap and check for convergence
    alphap <- alpha
    alpha <- (1+ sqrt(4*alpha*alpha +1))/2
    xxp <- x - xp
    Axy <- Ax - y
    
    # Compute the regularization part
    if (GFlag==0) {
      tree_norm <- treeNorm(x, n, ind, ncol(ind))
    } else {
      tree_norm <- general_treeNorm(x, n, G, ind, ncol(ind))
    }
    
    # function value = loss + regularization
    funVal <- c(funVal, crossprod(Axy)/2 + lambda * tree_norm)
    
    # If there is only little improvement, ten break
    if (bFlag == 1) {
      break
    }
    
    # Look for the value of termination flag
    switch (opts[["tFlag"]],
            {
              if (iterStep>=2){
                if (abs(funVal(iterStep) - funVal(iterStep-1) ) <= opts$tol) {
                  break
                }
              }
            },
            
            {
              if (iterStep>=2) {
                if (abs(funVal(iterStep) - funVal(iterStep-1) ) <= opts$tol * funVal(iterStep-1)){
                  break
                }
              }
            },
            
            {
              if (funVal[iterStep] <= opts$tol) {
                break
              }
            },
            
            {
              norm_xxp <- sqrt(crossprod(xxp))
              if (norm_xxp <= opts$tol) {
                break
              }
            },
            
            {
              norm_xp <- sqrt(crossprod(xp))
              norm_xxp <- sqrt(crossprod(xxp))
              if (norm_xxp <= opts$tol * max(norm_xp,1)) {
                break
              }
            },
            
            {
              if ((iterStep >= opts$maxIter) | (abs(funVal(iterStep) - funVal(iterStep-1) ) <= opts$tol)) {
                break
              }
            }
    )
    
    if(("rStartNum" %in% names(opts)) && (iterStep %% opts$rStartNum == 0)) {
      alphap <- 0
      alpha <- 1
      xp <- xAxp <- Ax
      xxp <- matrix(0, n, 1)
      L <- L/2;
    }
  }
  
  tll_out <- list("beta" = x, "valueL" = valueL, "funVal" <- funVal)
  return(tll_out)
  
}