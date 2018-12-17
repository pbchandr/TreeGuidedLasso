##################################################################################################################
#                      Function linearLasso : Regression with l1-ball constraint (aka LASSO)                          #
#                                   Problem : min  1/2 || A x - y||^2                                            #
#                                             s.t. ||z||_1 <= z                                                  #
##################################################################################################################
linearLasso <- function(A, y, z, opts) {
  
  # z <- rho
  valueL <- c()
  funVal <- c()
# Verify and initialize the parameters ---------------------------------------------------------------------------

  # Verify the number of input parameters
  # if (missing(A) | missing(y) | missing(z)){
  if (nargs() < 3) {
    stop("Inputs: A, y and z should be specified!")
  } else if (nargs() == 3) {
    opts <- list()
    print("Empty opts")
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
  opts <- setParams(opts)

# Normalization --------------------------------------------------------------------------------------------------
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

# Starting Point Initialization ----------------------------------------------------------------------------------
  # compute AT y
  if (opts[["nFlag"]] ==0) {
    # print("nFlag = 0")
    ATy <- t(A)%*%y
  } else if(opts[["nFlag"]] ==1){
    ATy <- t(A)%*%y - sum(y)*t(mu)
    ATy <- ATy./nu
  } else {
    invNu <- y./nu
    ATy <- t(A)%*%invNu - sum(invNu)*t(mu);
  }

  # Process the regularization parameter
  # L2 norm regularization
  if("rsL2" %in% names(opts)) {
    rsL2 <- opts[["rsL2"]]
    if(rsL2 < 0) {
      stop("rsL2 should not be non-negative!")
    }
  } else {
    rsL2 <- 0
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

  if(opts[["init"]] == 0) {
    x_norm <- sum(abs(x))
    x_2norm <-  as.numeric(t(x)%*%x)
    if (x_norm >= 1e-6) {
      ratio_max <- z/x_norm;
      ratio_optimal <-  as.numeric((t(Ax)%*%y)/ (t(Ax)%*%Ax + rsL2 * x_2norm))

      if (abs(ratio_optimal) <= ratio_max) {
        ratio <- ratio_optimal
      }
      else if (ratio_optimal < 0){
        ratio <- -ratio_max
      }
      else {
        ratio <- ratio_max
      }

      x <- ratio * x
      Ax <- ratio*Ax
    }
  }

# The Armijo Goldstein line search scheme + accelearted gradient descent -----------------------------------------
  bFlag <- 0 # this flag tests whether the gradient step only changes a little
  L <- 1 + rsL2 # % We assume that the maximum eigenvalue of A'A is over 1

  # Assign xp with x and Axp with xp
  xp <-x
  Axp <- Ax
  xxp <- matrix(0, ncol(A),1)

  # alphap and alpha are used for computing the weight in forming search point
  alphap <- 0
  alpha <- 1

  lambda0 <- 0

  for(iterStep in 1: opts$maxIter){
    # Step 1 - Compute search point based on xp and x (with beta)
    beta <- (alphap - 1)/alpha
    s <- x + beta * xxp
    
    # Step 2 - line search for : and compute the new approximate solution x
    
    # Compute gradient (g) at s
    As <- Ax + beta * (Ax - Axp)
    
    # compute ATAs
    ATAs <- t(A)%*%As
    
    # Obtain the gradient g
    g <- (ATAs - ATy) + rsL2*s
    
    # Assign x and Ax to xp and Axp
    xp <- x
    Axp <- Ax
    
    while(TRUE) {
      # let s walk in a step in the antigradient of s to get v and then do the l1-norm regularized projection
      v <- s - g/L
      
      # Projection - eplb
      
      out <- eplb(v, n, z, lambda0)
      x <- out[1:n]
      zfStep <- out[n + 1]
      lambda0 <- out[n + 2]
      
      # Difference between the new approximate solution x and the search point s
      v <- x - s
      
      # Compute Ax
      Ax <- tcrossprod(A, t(x))
      Av <- Ax - As
      rSum <- crossprod(v) #t(v)%*%v
      lSum <- crossprod(Av) #t(Av)%*%Av
      
      if(rSum <= 1e-200){
        bFlag <- 1 # Shows the gradient step makes little improvement
        break
      }
      
      # the condition is ||Av||_2^2 <= (L - rsL2) * ||v||_2^2
      if(lSum <= rSum*(L - rsL2)){
        break
      } else {
        L <- max(2*L, lSum/rSum + rsL2)
      }
    }
    
    valueL <- c(valueL, L)
    
    # Step 3 - Update alpha and alphap and check for convergence
    alphap <- alpha
    alpha <- (1+ sqrt(4*alpha*alpha +1))/2
    
    xxp <- x - xp
    xxp[1:5]
    Axy <- Ax - y
    Axy[1:5]
    
    funVal <- c(funVal, (t(Axy)%*%Axy)/2 + rsL2/2 * t(x)%*%x)
    
    
    # If there is only little improvement, ten break
    if (bFlag == 1) {
      break
    }
    
    # Look for the value of termination flag
    switch (opts[["tFlag"]],
            {
              if (iterStep>=2){
                if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts$tol) {
                  break
                }
              }
            },
            
            {
              if (iterStep>=2) {
                if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts$tol * funVal(iterStep-1)){
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
              norm_xxp <- sqrt(t(xxp)%*%xxp)
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
              if (iterStep >= opts$maxIter) {
                break
              }
            }
    )
  }
  # print(iterStep)

  ll_out <- list("beta" = x, "valueL" = valueL, "funVal" <- funVal)
  return(ll_out)
}
