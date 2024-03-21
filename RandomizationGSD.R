###########################################################
# This code
# 1. Generates n_sim randomization sequences from a specific randomization procedure
# 2. Generates n_sim N(0,1) distritbuted variables (the simulated patients)
# 3. Calculate mean and std for every stage and both groups
# 4. Create group sequential design (with parameters for stages, stopping boundaries, etc.)
# 5. Outputs T1E probability (rejection probability of the trial)


source("GSD-allocation.R")
source("GSD-rand_procs.R")
library(randomizeR)
library(rpact)

set.seed(42)

# Brauch ich glaub ich nicht mehr.
#split_into_stages <- function(vector, K) {
#  n <- length(vector)  # Total number of elements in the vector
#  subgroups <- vector("list", K)  # Initialize list to store the subgroups
#  
#  for (i in 1:K) {
# Calculate the number of elements to include in the current subgroup
#    n_elements <- ceiling(i * n / K)
#    subgroups[[i]] <- vector[(i-1)*n/K:n_elements]
#  }
#  return(unlist(subgroups))
#}


simulateClinicalTrial <- function(n_sim, n_patients, K) {
  
  # Create group sequential design
  design <- getDesignGroupSequential(
    typeOfDesign=c("OF"),
    alpha = 0.05,
    sided = 2,
    informationRates = c(0.33,0.67,1)
  )
  #  print(design)
  design <- getDesignInverseNormal(kMax = 3)
  
  # Initialize vectors to store results
  means1 <- vector("list", K)
  means2 <- vector("list", K)
  stDevs1 <- vector("list", K)
  stDevs2 <- vector("list", K)
  sequences <- vector("list", n_sim)
  
  randobj = bsdPar(n_patients, mti = 3, groups = c("0", "1"))
  randseq = genSeq(randobj, n_sim, seed = 2) 
  decision <- list()
  for (i in 1:n_sim) {
    sequence =  randseq@M[i,]               # Generate randomisation sequence
    patients <- rnorm(n_patients)       # Generate N(0,1) distributed random variables
    
    #    print("toTest:")
    #    print(sequence)
    #    print(patients)
    # Split patients into groups based on randomisation sequence
    #    group1 <- patients[sequence == 1]
    #    group2 <- patients[sequence == 0]
    # Split groups into stages
    #    group1 <- split_into_stages(group1, K)
    #    group2 <- split_into_stages(group2, K)
    # Calculate means and standard deviations for each group
    
    n_group1=vector("numeric", K)
    n_group2=vector("numeric", K)
    for (j in 1:K) {
      n_group1[j]=sum(sequence[(1+n_patients*(j-1)/K):(n_patients*j/K)])
      n_group2[j]=(n_patients/K)-sum(sequence[(1+n_patients*(j-1)/K):(n_patients*j/K)])
      means1[j] <-sum(patients[(1+n_patients*(j-1)/K):(n_patients*j/K)]*sequence[(1+n_patients*(j-1)/K):(n_patients*j/K)])/(n_group1[j])
      means2[j] <-sum(patients[(1+n_patients*(j-1)/K):(n_patients*j/K)]*(1-sequence[(1+n_patients*(j-1)/K):(n_patients*j/K)]))/(n_group2[j])
      #    means1[[j]] <-mean(group1[j]) # ich denke wegen subgruppen macht es am meisten sinn den mean manuell zu berechnen, wie oben
      #    means2[[j]] <-mean(group2[j]) # ansonstne habe ich hier einen fehler drin....
      #   stDevs1[[j]] <-c(1,1,1)
      #    stDevs2[[j]] <-c(1,1,1)
    }
    #  print(means2)
    #  print(n_group2)
    dataset <- getDataset(means1=unlist(means1), means2=unlist(means2), stDevs1=c(1,1,1),stDevs2=c(1,1,1), n1=n_group1, n2=n_group2)
    results = getStageResults(design = design, dataInput = dataset)
    testaction <- getTestActions(results)
    testdecision <- all(unlist(testaction == c("continue", "continue", "accept"))) # checks if all elements are in this form
    decision <- c(decision, testdecision)
    # Calculate the number of incorrect sequences
    print(results)
    print(testdecision)
    print(testdecision)
  }
  
  percentage_false <- (1 - mean(unlist(decision)))
  return( percentage_false) 
  #  return(t1e)
}

result <- simulateClinicalTrial(n_sim = 1, n_patients = 24, K = 3)
print(result)




library(gsDesign)     # Adjusted significance levels
library(drcarlate)
############### Manual calculation





set.seed(1)
# Test statistic
# Calculate means and standard deviations for each group
simulateClinicalTrial_manual <- function(method, n_sim, n_patients, K) {
  if (method=="naive") {     # naive calculcation of boundaries of gsd without updates of boundaries
    # Create group sequential design
    testdesign = gsDesign(k = K, test.type = 1, sfu = "Pocock")
    
    # Initialize vectors to store results
    means1 <- vector("list", K)
    means2 <- vector("list", K)
    stDevs1 <- vector("list", K)
    stDevs2 <- vector("list", K)
    sequences <- vector("list", n_sim)
    
    testdecision <- logical(n_sim) # Initialize with the number of simulations
    testdecision_comb <- logical()
    
    # Randomization object
    randobj = rpbrPar(24, rb = 2, groups = c("0", "1"))
    randseq = genSeq(randobj, n_sim) 
    
    
    for (i in 1:n_sim) {
      sequence =  randseq@M[i,]               # Generate randomisation sequence
      patients <- rnorm(n_patients)       # Generate N(0,1) distributed random variables
      
      
      n_group1=vector("numeric", K)   # Patients allocated to group 1 in each stage
      n_group2=vector("numeric", K)   # Patients allocated to group 2 in each stage
      means1 <- vector("numeric", K)  # cumulative means for group 1 in each stage
      means2 <- vector("numeric", K)  # cumulative means for group 2 in each stage
      I <- vector("numeric", K)       # Information fraction for each stage
      t <- vector("numeric", K)       # Test statistic for each stage
      continue <- vector("list", K)   # Continue/Reject decision for each stage
      
      
      
      for (j in 1:K) {                # for each stage
        n_group1[j]=sum(sequence[1:(n_patients*j/K)])                                                 # Calculate number of patients in stage j
        n_group2[j]=(n_patients*j/K)-sum(sequence[1:(n_patients*j/K)])
        means1[[j]] <-sum(patients[1:(n_patients*j/K)]*sequence[1:(n_patients*j/K)])/(n_group1[j])
        means2[[j]] <-sum(patients[1:(n_patients*j/K)]*(1-sequence[1:(n_patients*j/K)]))/(n_group2[j])
        I[j] = 1 * 1/sqrt(1/n_group1[j]+1/n_group2[j])
        t[j] = I[j]*(means1[j]-means2[j])
        continue[j]= (t[j] < testdesign$upper$bound[j])
      }
      testdecision[i] =all(unlist(continue))
      
      
    }
    percentage_false <- (1 - mean(testdecision))
    print(percentage_false)
    
  } else { # hier habe ich das einmal für den inverse normal combination test probiert, wird aber aktuell nicht benötigt.
    
    testdesign = gsDesign(k = K, test.type = 1, sfu = "Pocock")
    
    testdesign_comb = gsDesign(k = K, test.type = 1, sfu = "Pocock", n.I=c(1/sqrt(2), 1))
    #  I_comb =  c(sqrt(1/sqrt(2)), 1)
    #  testdesign_comb=  gsDesign(k = K, test.type = 2, sfu = "Pocock", n.I = I_comb)
    #  print(testdesign_comb$upper$bound)
    #  print(testdesign$upper$bound)
    # Create group sequential design
    design <- getDesignGroupSequential(
      typeOfDesign=c("OF"),
      alpha = 0.05,
      sided = 2,
      informationRates = c(0.33,0.67,1)
    )
    #  print(design)
    #  design <- getDesignInverseNormal(kMax = 3)
    
    # Initialize vectors to store results
    means1 <- vector("list", K)
    means2 <- vector("list", K)
    stDevs1 <- vector("list", K)
    stDevs2 <- vector("list", K)
    sequences <- vector("list", n_sim)
    
    # randobj = bsdPar(24, mti = 3, groups = c("0", "1"))
    randobj = rpbrPar(24, rb = 2, groups = c("0", "1"))
    
    randseq = genSeq(randobj, n_sim) 
    testdecision <- logical()
    testdecision_comb <- logical()
    
    for (i in 1:n_sim) {
      sequence =  randseq@M[i,]               # Generate randomisation sequence
      patients <- rnorm(n_patients)       # Generate N(0,1) distributed random variables
      
      n_group1=vector("numeric", K)
      n_group2=vector("numeric", K)
      n_group1_comb=vector("numeric", K)
      n_group2_comb=vector("numeric", K)
      means1 <- vector("numeric", K)
      means2 <- vector("numeric", K)
      I <- vector("numeric", K)
      t <- vector("numeric", K)
      continue <- vector("list", K)
      p <- vector("numeric", K)
      t_combination <- vector("numeric", K)
      continue_combination <- vector("list", K)
      means1_comb <- vector("numeric", K)
      means2_comb <- vector("numeric", K)
      t_comb <- vector("numeric", K)
      
      
      
      for (j in 1:K) {
        n_group1[j]=sum(sequence[1:(n_patients*j/K)])
        n_group2[j]=(n_patients*j/K)-sum(sequence[1:(n_patients*j/K)])
        means1[[j]] <-sum(patients[1:(n_patients*j/K)]*sequence[1:(n_patients*j/K)])/(n_group1[j])
        means2[[j]] <-sum(patients[1:(n_patients*j/K)]*(1-sequence[1:(n_patients*j/K)]))/(n_group2[j])
        I[j] = 1 * 1/sqrt(1/n_group1[j]+1/n_group2[j])
        t[j] = I[j]*(means1[j]-means2[j])
        continue[j]= (t[j] < testdesign$upper$bound[j])
        
        
        # for inverse normal test
        n_group1_comb[j]=sum(sequence[(1+n_patients*(j-1)/K):(n_patients*j/K)])
        n_group2_comb[j]=(n_patients/K)-sum(sequence[(1+n_patients*(j-1)/K):(n_patients*j/K)])
        means1_comb[j] <-sum(patients[(1+n_patients*(j-1)/K):(n_patients*j/K)]*sequence[(1+n_patients*(j-1)/K):(n_patients*j/K)])/(n_group1_comb[j])
        means2_comb[j] <-sum(patients[(1+n_patients*(j-1)/K):(n_patients*j/K)]*(1-sequence[(1+n_patients*(j-1)/K):(n_patients*j/K)]))/(n_group2_comb[j])
        t_comb[j] = I[j]*(means1_comb[j]-means2_comb[j])
        p[j] = pnorm(q=t_comb[j], mean=0, sd=1, lower.tail=FALSE)
        
        
      }
      testdecision =c(testdecision, all(unlist(continue)))
      
      
      
      # for inverse normal test
      # p2 sollte aber nur die daten von p2 beinhalten
      t_combination[1]=qnorm((1-p[1]))
      t_combination[2] = (1/sqrt(2))*t_combination[1]+(1/sqrt(2))*qnorm((1-p[2]))
      #print("pvalue")
      #print(p)
      continue_combination[1]= (t_combination[1] < testdesign_comb$upper$spend[1]) ##### ??
      continue_combination[2]= (t_combination[2] < testdesign_comb$upper$bound[2])
      testdecision_comb =c(testdecision_comb, all(unlist(continue_combination)))
      
      
      
      
      #print("fehlerkontrolle")
      #print(t_combination)
      #print(testdecision)
      #print(testdecision_comb)
      
      #design <- getDesignInverseNormal(typeOfDesign="P", kMax=2)
      #dataset <- getDataset(means1=unlist(means1_comb), means2=unlist(means2_comb), stDevs1=c(1,1),stDevs2=c(1,1), n1=n_group1, n2=n_group2)
      #results = getStageResults(design = design, dataInput = dataset)
      #print(results)
    }
    percentage_false <- (1 - mean(testdecision))
    print(percentage_false)
    percentage_false_comb <- (1 - mean(testdecision_comb))
    print("Rejection prob comb test:")
    print(percentage_false_comb)
    #  print(continue_combination)
    
    
  }
}
simulateClinicalTrial_manual(method="naive", n_sim=10000, n_patients=24,K=2)



#### Calculation of T1E for inverse normal combination test
library(mvtnorm)      # Integral calculation

#mu= c(1.144,1.144)
#mu=c(0,0)
#cov =  matrix(c(1, sqrt(1/sqrt(2)), sqrt(1/sqrt(2)),1), nrow = 2, ncol = 2, byrow = TRUE) # nach Wassmer

#testdesign = gsDesign(k = 2, test.type = 1, sfu = "Pocock", n.I=c(1/sqrt(2),1)) # nach Buch von Wassmer muss Information anders gewählt werden
#print(testdesign)

#upper_bound = testdesign$upper$bound[1:2]
#mean_values = mu[1:2]
#cov_matrix = cov[1:2, 1:2]
#integral = pmvnorm(algorithm = Miwa(), lower = -Inf, upper = upper_bound, 
#                   mean = mean_values, sigma = cov_matrix, abseps = 1e-16)
#alpha = 1-integral
#alpha = round(alpha, digits=5)
#print(alpha)


###

n=24

#mu= c(1.144,1.144)
mu=c(0,0)
I1 = 1/(sqrt(1/(n/4)+1/(n/4))) 
I2 = 1/(sqrt(1/(n/2)+1/(n/2)))
cov =  matrix(c(1, I1*I2*(2/(n/2)), I1*I2*(2/(n/2)),1), nrow = 2, ncol = 2, byrow = TRUE)

testdesign = gsDesign(k = 2, test.type = 1, sfu = "Pocock")

#lower_bound = testdesign$lower$bound[1:2]
upper_bound = testdesign$upper$bound[1:2]
mean_values = mu[1:2]
cov_matrix = cov[1:2, 1:2]
integral = pmvnorm(algorithm = Miwa(), lower = -Inf, upper = upper_bound, 
                   mean = mean_values, sigma = cov_matrix, abseps = 1e-16)
beta = integral
beta = round(beta, digits=5)
print(beta)


