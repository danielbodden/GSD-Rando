source("standards.R") 
library(randomizeR)
library(rpact)
library(gsDesign) 
library(mvtnorm)    
library(openxlsx)
library(parallel)
library(readxl)
library(ggplot2)
library(dplyr)
library(rpact)              

# Generates randomization sequences for power calculation
Sequence_generation <- function(n, K, RP, n_sim,  rb = 4, mti =3, p=2/3) {
  if (!((n/K) %% 1 == 0)) {                                                 
    stop("The amount of stages is not divisible by the sample size.")
  }
  
  # Generate randomization object
  randobj <- switch(as.character(RP),
                    "CR" = crPar(n, K = 2),
                    "RAR" = rarPar(n, K = 2, groups = c("0", "1")),
                    "BSD" = bsdPar(n, mti = mti, groups = c("0", "1")),
                    "EBC" = ebcPar(n, p, groups = c("0", "1")),
                    "CHEN" = chenPar(n, mti = mti, p = p, groups = c("0", "1")),
                    "PBR" = rpbrPar(n, rb = rb, groups = c("0", "1")),
                    "MP" = mpPar(n, mti = mti, ratio = c(1, 1)),
                    stop("Invalid RP parameter.")
  )
  # Generates n_sim randomization sequences from a randomization procedure with 
  # the actual probability of the randomization sequence happening for the randomization procedure
  seq = genSeq(randobj, n_sim, seed = 42) 
  return(seq)
}

###########################################################
# Calculates the power for
# - Pocock / OF
# - LanDeMets alpha spending Pocock / OF
# - Pocock / Of with correction for actual information
# based on MVN distribution conditioned on randomization sequence
# by Monte Carlo simulations over randomization sequences generated from a specific randomization procedure
# Note: alpha is always one-sided!
Power_condMVN <- function(n, n_sim, K, RP, sfu, sides = 1, alpha =0.025,  rb = 4, mti =3, p=2/3, delta, futility=FALSE, futility_binding=FALSE) {
  Power = rep(999999, n_sim)   
  seq = Sequence_generation(n=n, K=K, RP = RP, n_sim= n_sim, rb=rb,mti=mti,p=p) # Generate randomization sequences
  counter_zero_allocations = 0      # Counter for skipped sequences
  
#      seq@M = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1), nrow=1, byrow=TRUE)              #   Hardcoded sequence for testing

  
  # Set up parallel computing
  #  cl <- makeCluster(detectCores() - 1) # Create a cluster using one less than the number of cores available
  #cl <- makeCluster(4) # Create a cluster using one less than the number of cores available
  # Export needed variables and functions to the cluster
  #clusterExport(cl, varlist=c("seq", "n", "K", "RP", "sfu", "sides", "alpha", "rb", "mti", "p", "delta", "futility"), envir=environment())
  # Load required libraries on all cluster nodes. Also source the standards.R script.
#  clusterEvalQ(cl, {
 #   source("standards.R")                                                       # Make sure the path to standards.R is correct and accessible from the cluster nodes
 #   library(randomizeR)
 #   library(RPact)
 #   library(gsDesign)
#    library(mvtnorm)
#    library(openxlsx)
#  })
  # Execute the core of Power_condMVN inw parallel
  #Power <- unlist(parLapply(cl, 1:n_sim, function(j) {
   Power <- unlist(lapply(1:n_sim, function(j) {                               # Use this loop when no Cluster is in use
      
    # Initalize vectors for calculation
    current_seq = seq@M[j,]
    n_A = n_B =  numeric(K+1)
    I = Summe = resp = numeric(K)
    n_A[[1]] = n_B[[1]] = 0 
    k=n/K
    sigma = 1
    
    subseq = current_seq[1:k]                                   # Allocations in Stage i 
    n_A[2] = sum(subseq)                                        # total sample size in group A until stage 2
    n_B[2] = k-n_A[2]                                           # total sample size in group B until stage 2
    
    if (futility == FALSE && futility_binding== TRUE) {
      futility_binding = FALSE
      warning ("Futility_binding set to FALSE AS futility is FALSE in parameter settings for Power_condMVN.")
    }

    # When any group is empty in stage 1, skip the randomization sequence
    if(n_A[[2]] == 0 | n_B[2] == 0) {                
      Pow = -99
      counter_zero_allocations <<- counter_zero_allocations + 1
    }
    else {
        I[1] =  1 / ( sigma/ n_A[[2]]+ sigma /n_B[[2]] )
      for (i in (2:K)) { # For every stage
        subseq = current_seq[((i-1)*k+1):(i*k)]                                   # Allocations in Stage i 
        n_A[i+1] = n_A[i]+sum(subseq)                                             # total sample size in group A until stage i
        n_B[i+1] = (i*k)-n_A[i+1]
        I[i] =  1 / ( sigma/ n_A[[i+1]]+ sigma /n_B[[i+1]] )                    # Information for each stage 
      }
      if (!(K==1)) {
        if (futility_binding==FALSE) {
          if (sfu == "LDMOF") { 
            testdesign = gsDesign(k=K, test.type = sides , sfu = sfLDOF, alpha= alpha, n.I=I)  # alpha is always one-sided in gsDesign, even if 2-sided test design is selected.
          } else if (sfu == "LDMPocock") {
            testdesign = gsDesign(k=K, test.type = sides , sfu = sfLDPocock, alpha= alpha, n.I=I)  # alpha is always one-sided in gsDesign, even if 2-sided test design is selected.
          } else  {
          if (sfu=="coRPOC") {
            sfu = "Pocock"
            testdesign = gsDesign(k=K, test.type = sides , sfu = sfu, alpha= alpha, n.I=I)  
          } else if (sfu=="corOF") {
            sfu="OF"
            testdesign = gsDesign(k=K, test.type = sides , sfu = sfu, alpha= alpha, n.I=I)  
          } else {
          testdesign = gsDesign(k=K, test.type = sides , sfu = sfu, alpha= alpha)       # alpha is always one-sided in gsDesign, even if 2-sided test design is selected.
          }
        }   
        upper_bound = testdesign$upper$bound
        if (sides==2) {                       
          lower_bound =testdesign$lower$bound
        }
        } else if (futility_binding == TRUE) {   
          if (sfu=="coRPOC") {            # corPoc not supported with binding futility. Futility boundary set to non-binding.
            sfu = "Pocock"
            testdesign = gsDesign(k=K, test.type = sides , sfu = sfu, alpha= alpha, n.I=I)  
            upper_bound = testdesign$upper$bound
          } else if (sfu=="corOF") {
            sfu="OF"    # corOF not supported with binding futility. Futility boundary set to non-binding.
            testdesign = gsDesign(k=K, test.type = sides , sfu = sfu, alpha= alpha, n.I=I)   
            upper_bound = testdesign$upper$bound
          } else {
          I = I/I[[length(I)]]            # rpact uses information fraction instead of total information
          design_types_LDM <- c("LDMOF" = "asOF", "LDMPocock" = "asP")
          
          if (sfu %in% names(design_types_LDM)) {
            design <- getDesignGroupSequential(
              sided = 1, alpha = alpha,
              informationRates = I, typeOfDesign = design_types_LDM[[sfu]],
              futilityBounds = rep(0, K-1), bindingFutility = TRUE
            )
          }
          design_types_standard <- c("OF" = "OF", "Pocock" = "P")
          
          if (sfu %in% names(design_types_standard)) {
            design <- getDesignGroupSequential(
              sided = 1, alpha = alpha, typeOfDesign = design_types_standard[[sfu]],
              futilityBounds = rep(0, K-1), bindingFutility = TRUE
            )
          }
          
          
          upper_bound <- design$criticalValues
          }
        }
      } else {
        stop("Please input K>1.")
      }
      if (sides == 1) {                                                           # If one-sided design is selected, set the lower_bound to -infinity
        lower_bound=rep(-99, K)
      }

      if (futility == TRUE) {                                                    # (binding) futility boundary 
        lower_bound=rep(0, K)                                                    
      }
 
      r = 120  
      zbdy <- rbind(lower_bound, upper_bound)
      results <- gst1(r, na=K, inf=I, zbdy, theta=delta)                         # Function provided by Chris Jennison
      Pow = results[[2]]
    }
      Power[j] <<- Pow
    }))

   # Stop cluster
#  stopCluster(cl)
   Power = Power[Power != -99]
   return(list( Pow = Power, Counter = counter_zero_allocations))
}

#Power_condMVN(n=24, n_sim=1, K=3, RP="PBR", sfu="LDMPocock", sides = 1, alpha =0.025,  rb = 2, mti =3, p=2/3, delta=0, futility=FALSE, futility_binding=TRUE)
#Power_condMVN(n=24, n_sim=1, K=3, RP="PBR", sfu="LDMPocock", sides = 1, alpha =0.025,  rb = 2, mti =3, p=2/3, delta=0)


# Calculates power for 2 or 3 stages when given EV; Cov and boundaries. Calculation is performed manually, armitage formular not used.
MVN_calculation <- function(mean, Cov, K, lower_bound, upper_bound) {
  
  Power =  pmvnorm(lower = c(upper_bound[1]), upper = c(999), 
                   mean = mean[1], sigma = Cov[1,1], abseps = 1e-16)
  
  
  if (K>1) {
    Power = Power + pmvnorm(lower = c(lower_bound[1], upper_bound[2]), upper = c(upper_bound[1], 999), 
                            mean = mean[1:2], sigma = Cov[1:2,1:2], abseps = 1e-16)  
  }
  if (K==3) {
    
    Power = Power + pmvnorm(lower = c(lower_bound[1], lower_bound[2], upper_bound[3]), upper = c(upper_bound[1], upper_bound[2], 999), 
                            mean = mean[1:3], sigma = Cov[1:3,1:3], abseps = 1e-16)
  }
  
  return(Power)
}



# Checks if any stage only does not have allocations to both groups
check_subsequences <- function(current_seq, K) {
  zero_allocation_stage = FALSE
  n <- length(current_seq)
  subsequences <- split(current_seq, rep(1:K, each = n/K))  # Split the vector into K subsequences
  for (subseq in subsequences) {
    unique_values <- unique(subseq)
    if (length(unique_values) == 1) {  # If there is only one unique value in the subsequence
      #     cat("Subsequence", subseq, "has only", unique_values, "entries\n")
      zero_allocation_stage = TRUE
    }
  }
  return(zero_allocation_stage)
}

#### Calculation of Power for inverse normal combination test 
# only for K=2 and K=3 and only one-sided for alpha =.025 currently supported
Power_inverse_normal <- function(n, K, RP, n_sim, delta, sfu, futility=FALSE, futility_binding=FALSE) {
  Power = rep(999999, n_sim)   
  seq = Sequence_generation(n=n, K=K, RP = RP, n_sim= n_sim)
  counter_zero_allocations = 0                                                # Amount of allocations that have been skipped due to only allocations to one group
  Power <- unlist(lapply(1:n_sim, function(j) {                               # loop when no Cluster is in use
    k=n/K
    
    if (futility == FALSE && futility_binding== TRUE) {
      futility_binding = FALSE
      warning ("Futility_binding set to FALSE AS futility is FALSE in parameter settings for Power_inverse_normal.")
    }
    
    
    # Initalize vectors for calculation
    current_seq = seq@M[j,]
    if (check_subsequences(current_seq, K) == TRUE) {                       # Check if in any stage only allocations to one group
      counter_zero_allocations <<- counter_zero_allocations +1 
      Power[j] = -99
    } else {
      I = Summe = resp = numeric(K)
      n_A = n_B =  numeric(K+1)
     n_A[[1]] = n_B[[1]] = 0 

     for (i in (1:K)) {                                                          # For every stage
      subseq = current_seq[((i-1)*k+1):(i*k)]                                   # Allocations in Stage i 
      n_A[i+1] = sum(subseq)                                                   # total sample size in group A until stage i
      n_B[i+1] = (k)-n_A[[i+1]]                                               # total sample size in group B until stage i
      sigma = 1
      I[i] =  1 / ( sigma/ n_A[[i+1]]+ sigma /n_B[[i+1]] )                    # Information for each stage 

      if (futility_binding == FALSE) {
      testdesign = gsDesign(k=K, test.type = 1 , sfu = sfu, alpha= 0.025)  
      upper_bound = testdesign$upper$bound
      }
      if (futility_binding == TRUE) {
        design_types <- c("OF" = "OF", "Pocock" = "P")
        
        if (sfu %in% names(design_types)) {
          design <- getDesignGroupSequential(
            sided = 1, alpha = 0.025, typeOfDesign = design_types[[sfu]],
            futilityBounds = rep(0, K-1), bindingFutility = TRUE
          )
        }
        upper_bound <- design$criticalValues
      }
      if (futility == FALSE) {
      lower_bound=rep(-99, K)
      } else {
        lower_bound=rep(0, K)
        
      }
     }
    # Mean and Cov for inverse normal combination function
    mean <- sapply(1:K, function(k) sum(sqrt(I[1:k])) / sqrt(k) * delta)
    entries <- c(1, 1/sqrt(2), 1/sqrt(3), 
                 1/sqrt(2), 1, 2/sqrt(6), 
                 1/sqrt(3), 2/sqrt(6), 1)
    Cov <- matrix(entries, nrow = 3, ncol = 3, byrow = TRUE)
    Power[j] =MVN_calculation(mean=mean[1:K], Cov = Cov[1:K, 1:K], K=K, lower_bound = lower_bound, upper_bound=upper_bound) 
    }
  }))
  Power = Power[Power != -99]
  return(list(Pow = Power, Counter = counter_zero_allocations))
}

Power_inverse_normal(n=24, K=3, RP="PBR", n_sim=1, delta=0, sfu="OF", futility=TRUE, futility_binding=TRUE)
#Power_inverse_normal(n=24, K=3, RP="PBR", n_sim=10, delta=0, sfu="OF")


