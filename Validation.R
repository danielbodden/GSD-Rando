source("RandomizationGSD.R")


# Currently not needed for paper -> Use later for validation.
# Calculates the T1E probability based on Monte Carlo n_sim of the patients and randomization sequences
SimulateTeststatistics <- function(sfu, n, K, RP, n_sim) {
  if (sfu=="naive") {                                                        # naive Calculcation of boundaries of GSD without update of boundaries
    testdesign = gsDesign(k = K, test.type = 1, sfu = "Pocock")                 # Create GSD
    # TODO: Make sfu variable
    # Initialize vectors to store results 
    means1 = vector("list", K)
    means2 = vector("list", K)
    stDevs1 = vector("list", K)
    stDevs2 = vector("list", K)
    sequences = vector("list", n_sim)
    
    testdecision <- logical(n_sim)                                              # Initialize with the number of n_sim
    testdecision_comb <- logical()
    
    # Randomization object
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
    randseq = genSeq(randobj, n_sim) 
    
    
    for (i in 1:n_sim) {                                                        # For number of n_sim
      sequence =  randseq@M[i,]                                                 # Generate a randomization sequence
      patients <- rnorm(n)                                                      # Generate n N(0,1) distributed random variables
      
      # Initalize vectors for calculation
      n_group1=vector("numeric", K)                                             # Patients allocated to group 1 in each stage
      n_group2=vector("numeric", K)                                             # Patients allocated to group 2 in each stage
      means1 <- vector("numeric", K)                                            # cumulative means for group 1 in each stage
      means2 <- vector("numeric", K)                                            # cumulative means for group 2 in each stage
      I <- vector("numeric", K)                                                 # Information fraction for each stage
      t <- vector("numeric", K)                                                 # Test statistic for each stage
      continue <- vector("list", K)                                             # Continue / reject decision for each stage
      
      
      
      for (j in 1:K) {                                                          # for each stage
        n_group1[j]=sum(sequence[1:(n*j/K)])                                    # Calculate number of patients in stage j for group 1
        n_group2[j]=(n*j/K)-sum(sequence[1:(n*j/K)])                            # Calculate number of patients in stage j for group 2  
        
        if (n_group1[j]==0 | n_group2[j]==0) {                                  # if one group is empty, jump to next stage
          continue[j] = TRUE                                                    # TODO: For other designs, adjust boundaries if one gorup is empty 
        } else {
          means1[[j]] <-sum(patients[1:(n*j/K)]*sequence[1:(n*j/K)])/(n_group1[j])      # Mean of group 1 in stage j
          means2[[j]] <-sum(patients[1:(n*j/K)]*(1-sequence[1:(n*j/K)]))/(n_group2[j])  # mean of group 2 in stage j
          I[j] = 1 * 1/sqrt(1/n_group1[j]+1/n_group2[j])                                # Information stage j
          t[j] = I[j]*(means1[j]-means2[j])                                             # Test statistic stage j
          continue[j]= (t[j] < testdesign$upper$bound[j])                               # Decision for stage j
        }
      }
      testdecision[i] =all(unlist(continue))                                          # Overall decision: H0 not rejected if decision in all stages continue=TRUE
    }
    T1E <- (1 - mean(testdecision))                                      
    return(T1E)
  }
}
#SimulateTeststatistics(sfu="naive", RP="CR", n_sim=8, n=24,K=3)


