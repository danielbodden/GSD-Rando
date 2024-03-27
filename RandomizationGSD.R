###########################################################
# This code
# 1. Generates n_sim randomization sequences from a specific randomization procedure
# 2. Generates n_sim N(0,1) distritbuted variables (the simulated patients)
# 3. Calculate mean and std for every stage and both groups
# 4. Create group sequential design (with parameters for stages, stopping boundaries, etc.)
# 5. Outputs T1E probability (rejection probability of the trial)


#source("GSD-allocation.R")
#source("GSD-rand_procs.R")
source("standards.R")
library(randomizeR)
library(rpact)
library(gsDesign)     # Adjusted significance levels
library(drcarlate)


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
simulateClinicalTrial_manual(method="naive", n_sim=3, n_patients=24,K=2)



#### Calculation of T1E for inverse normal combination test
library(mvtnorm)      # Integral calculation

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


## Calculation of Power for 2-stage group sequential design via MVN distribution
# Example for n1=n2=20
delta=0.3
cov= matrix(c(1, 1/sqrt(2)*0.5, 1/sqrt(2)*0.5, 1), nrow =2, ncol = 2)
#cov= matrix(c(1, 1, 1, 1), nrow =2, ncol = 2)
mu=c(delta*sqrt(40)*0.5,delta*sqrt(80)*0.5)
#mu=c(0,0)
r=c(2.18, 2.18)
a=c(0, 2.18)
int1 =pmvnorm(mean=mu, sigma=cov, lower=c(-r[1], -Inf), upper=c(-a[1], -r[2]))
int2 = pmvnorm(mean=mu, sigma=cov, lower=c(-r[1], r[2]), upper=c(-a[1], Inf))
int3 =pmvnorm(mean=mu, sigma=cov, lower=c(a[1], -Inf), upper=c(r[1], -r[2]))
int4 =pmvnorm(mean=mu, sigma=cov, lower=c(a[1], r[2]), upper=c(r[1], Inf))
int5= 1-pmvnorm(mean=mu[1], sigma=cov[1,1], lower=c(-r[1]), upper=c(r[1]))
power2=int1+int2+int3+int4+int5
print(power2)


#für K=2 komme ich auf die gleiche Power, wie wenn ich Jennisons approach benutze.
#Für K=3 ist noch irgendwas falsch....

#ToDo: Do the same for K=3
delta=0.3
cov= matrix(c(1, 1/sqrt(2), 1/sqrt(3), 1/sqrt(2), 1, sqrt(2/3), 1/sqrt(3), sqrt(2/3), 1 ), nrow =3, ncol = 3)
print(cov)
mu=c(delta*sqrt(20),delta*sqrt(40))
mu=c(0,0,0)
r=c(2.29, 2.29,2.29)
a=c(0, 0, 2.29)
int1 =pmvnorm(mean=mu[1:2], sigma=cov[1:2,1:2], lower=c(-r[1], -Inf), upper=c(-a[1], -r[2]))
int2 = pmvnorm(mean=mu[1:2], sigma=cov[1:2,1:2], lower=c(-r[1], r[2]), upper=c(-a[1], Inf))
int3 =pmvnorm(mean=mu[1:2], sigma=cov[1:2,1:2], lower=c(a[1], -Inf), upper=c(r[1], -r[2]))
int4 =pmvnorm(mean=mu[1:2], sigma=cov[1:2,1:2], lower=c(a[1], r[2]), upper=c(r[1], Inf))
int5= 1-pmvnorm(mean=mu[1], sigma=cov[1,1], lower=c(-r[1]), upper=c(r[1]))

int6 =pmvnorm(mean=mu, sigma=cov, lower=c(-r[1], -r[2] -Inf), upper=c(-a[1], -a[2], -r[3]))
int7 = pmvnorm(mean=mu, sigma=cov, lower=c(-r[1],-r[2], r[3]), upper=c(-a[1], -a[2], Inf))

int8 =pmvnorm(mean=mu, sigma=cov, lower=c(a[1],a[2], -Inf), upper=c(r[1],r[2], -r[3]))
int9 =pmvnorm(mean=mu, sigma=cov, lower=c(a[1],a[2], r[3]), upper=c(r[1],r[2], Inf))

int10 =pmvnorm(mean=mu, sigma=cov, lower=c(-r[1], a[2] -Inf), upper=c(-a[1], r[2], -r[3]))
int11= pmvnorm(mean=mu, sigma=cov, lower=c(-r[1],a[2], r[3]), upper=c(-a[1], r[2], Inf))

int12 =pmvnorm(mean=mu, sigma=cov, lower=c(a[1],-r[2], -Inf), upper=c(r[1],-a[2], -r[3]))
int13 =pmvnorm(mean=mu, sigma=cov, lower=c(a[1],-r[2], r[3]), upper=c(r[1],-a[2], Inf))

power3= int1+int2+int3+int4+int5+int6+int7+int8+int9+int10+int11+int12+int13
print(power3)
print(cov[1:2,1:2])



###########################################################


# Function for power calculation
Power_calculation <- function(n, reps, K, randproc, sfu, rb = 4, mti =3, p=2/3, effect_size, information) {
  # error control
  if (!((n/K) %% 1 == 0)) {
    stop("The amount of stages is not divisible by the sample size.")
  }
  randobj <- switch(randproc,
                    "CR" = crPar(n, K = 2),
                    "RAR" = rarPar(n, K = 2, groups = c("0", "1")),
                    "BSD" = bsdPar(n, mti = mti, groups = c("0", "1")),
                    "EBC" = ebcPar(n, p, groups = c("0", "1")),
                    "CHEN" = chenPar(n, mti = mti, p = p, groups = c("0", "1")),
                    "PBR" = rpbrPar(n, rb = rb, groups = c("0", "1")),
                    "MP" = mpPar(n, mti = mti, ratio = c(1, 1)),
                    stop("Invalid randproc parameter.")
  )
  
  
  seq = genSeq(randobj, reps, seed = 42)         # Generates a randomization sequence from randobj. Second parameter for amount of sequences to be created.
  #      seq@M = matrix(c(1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0), nrow=1, byrow=TRUE)  #   Hardcoded sequence for testing
  Power = rep(999999999999, reps)        #
  
  r <- 120  # Assume this scalar somehow determines the mesh density

  # Hypothetical information fractions at each analysis
  # This might be proportional to the number of patients recruited by these times

  # Boundary values for stopping the trial at each interim analysis
  # Assuming standard normal boundaries for simplicity
  
  if (!(K==1)) {
    testdesign = gsDesign(k=K, test.type = 2 , sfu = sfu, alpha= 0.025)
    lower_bound =testdesign$lower$bound
    upper_bound = testdesign$upper$bound
    lower_bound =c(-8, -8)
    upper_bound = c(2.157,2.201)
  }
  zbdy <- rbind(lower_bound, upper_bound)
  print(zbdy)
  
  # Now you would call the function with these inputs
  results <- gst1(r, na=K, inf=information, zbdy, theta=effect_size)
  print("Power:")
  print(results)
  print(results[[2]]+results[[3]])
  
  # check with nquery if correct
  # then change up Cov according to given randomisation sequence
  # then allow to change boundaries
  # then what about power for inverse normal combination test? 
  
  
}

 
 Power_calculation(n=90, reps=1, randproc="CR", sfu="Pocock", K=2, effect_size=0.3, information= c((sqrt(40)/2)**2, (sqrt(80)/2)**2))
 
                   