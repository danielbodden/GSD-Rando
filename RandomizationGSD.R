source("standards.R")
library(randomizeR)
library(rpact)
library(gsDesign)    
library(drcarlate)
library(mvtnorm)    


# Calculates the T1E probability based on Monte Carlo simulations of the patients and randomization sequences
# Function currently not in use
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




## Calculation of Power manually for 2-stage group sequential design via MVN distribution
# Example for n1=n2=20
manual_calc_power_mvn <- function() {
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
# Leads to same results as using CJ implementation

# Same for K=3, some error here
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
}


#### Calculation of T1E for inverse normal combination test K=2
# hher einmal genau aufschreiben in Overleaf

n=24


delta= 0.3
mu=c(delta*sqrt(40)*0.5,delta*sqrt(80)*0.5)

I1 = 1/(sqrt(1/(n/4)+1/(n/4))) 
I2 = 1/(sqrt(1/(n/2)+1/(n/2)))
cov =  matrix(c(1, I1*I2*(2/(n/2)), I1*I2*(2/(n/2)),1), nrow = 2, ncol = 2, byrow = TRUE)

cov =  matrix(c(1,1/sqrt(2), 1/sqrt(2),1), nrow = 2, ncol = 2, byrow = TRUE)


testdesign = gsDesign(k = 2, test.type = 1, sfu = "Pocock")




#### Calculation of T1E for inverse normal combination test K=3

n=24

power_inverse_normal <- function(effect_size, n, K, sfu, onesided="yes", alpha=0.025) {
  mu=c(effect_size*sqrt(n/K)*0.5,effect_size*sqrt(n/K*2)*0.5, effect_size*sqrt(n/K*3)*0.5)
  
  cov =  matrix(c(1,1/sqrt(2), 1/sqrt(3), 1/sqrt(2),1, sqrt(2)/sqrt(3), 1/sqrt(3), sqrt(2)/sqrt(3), 1), nrow = 3, ncol = 3, byrow = TRUE)
  testdesign = gsDesign(k = 3, test.type = 1, sfu = "Pocock")
  
    upper_bound = testdesign$upper$bound[1:3]
  mean_values = mu[1:3]
  cov_matrix = cov[1:3, 1:3]
  integral = pmvnorm(algorithm = Miwa(), lower = -Inf, upper = upper_bound, 
                     mean = mean_values, sigma = cov_matrix, abseps = 1e-16)
  beta = integral
  beta = round(beta, digits=5)
  return(1-beta)
}

power_inverse_normal(effect_size=0.5,n=90, K=3, sfu="Pocock")
power_inverse_normal(effect_size=0,n=24, K=2, sfu="Pocock")



###########################################################
# Power calculation for LDM and naive Pocock and OF design and inverse normal combination test
# based on MVN distribution conditioned on randomization sequence
# calculated by monte carlo simulation over randomization sequences of a randomization procedure
# note: alpha is always one-sided!
Power_calculation <- function(n, reps, K, randproc, sfu, sides = 1, alpha =0.025,  rb = 4, mti =3, p=2/3, effect_size, futility="no") {
  # error control
  Power = rep(999999, reps)        #
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
    #    seq@M = matrix(c(0,0,0,0,0,0,1,0,0,0,0,0), nrow=1, byrow=TRUE)  #   Hardcoded sequence for testing
  for (j in 1:reps) {
    current_seq = seq@M[j,]

    n_A = n_B =  numeric(K+1)
    I = Summe = resp = numeric(K)
    n_A[[1]] = n_B[[1]] = 0 
    k=n/K
 
    for (i in (1:K)) {                      # Count the allocations to group A and B in each stage
      subseq = current_seq[((i-1)*k+1):(i*k)]                   # Allocations in Stage i 
      n_A[i+1] = n_A[i]+sum(subseq)                         # total sample size in group A until stage i
      n_B[i+1] = (i*k)-n_A[[i+1]]                         # total sample size in group B until stage i
      sigma = 1
      if(n_A[[i + 1]] == 0 | n_B[[i + 1]] == 0) {     # when there have been no allocations to one group, the information is 0 (due to checks in gsdesign function we set it very close to 0)
       I[i] = 0.000001 * i
      } else {
        I[i] =  1 / ( sigma/ n_A[[i+1]]+ sigma /n_B[[i+1]] )   # Information for each stage 
      }
    }
    r <- 120  # Assume this scalar somehow determines the mesh density
    if (!(K==1)) {
      testdesign = gsDesign(k=K, test.type = 2 , sfu = sfu, alpha= alpha, n.I=I)     # alpha is always one-sided in gsDesign, even if 2-sided test design is selected.
      lower_bound =testdesign$lower$bound
      upper_bound = testdesign$upper$bound
    }
    if (sides == 1) {            # for one-sided test set lover boundary to -infty
      lower_bound=rep(-99, K)
    }

    for (i in (1:K)) {                         # if there have been no allocations to one group jump to next stage
      if(n_A[[i+1]] == 0 | n_B[[i+1]] == 0) {
        lower_bound[i] <- -99
        upper_bound[i] <-  99
      }
    }
    if (futility == "yes") {                   # (binding) futility boundary 
      lower_bound=rep(0, K)
    }
 
    zbdy <- rbind(lower_bound, upper_bound)
    results <- gst1(r, na=K, inf=I, zbdy, theta=effect_size)  # implemented by CJ
    Power[j] =results[[2]]
  }
  return(mean(Power))
}
Power_calculation(n=24, reps=1, randproc="PBR", sfu=sfLDPocock, K=3, effect_size=0.5, sides=1)
Power_calculation(n=24, reps=1, randproc="PBR", sfu=sfLDOF, K=3, effect_size=0.5)

 
# Plots the power for LDM and naive Pocock and OF design for different randomization procedures
Plot_power <- function(n, reps, K, sfu, alpha=0.025, sides=1) {
  
  # Generate a sequence of effect_size values between 0 and 2
  effect_sizes <- seq(0, 2, by = 0.2)
  
  # Define sfu values and their corresponding plotting labels
  sfu_values <- list("Pocock" = sfu, "spendingPocock" = sfLDPocock)  # Assume sfu is a string or object, and sfldOF is an object
  
  # Function to calculate power for a given method and sfu object
  calculate_power_for_method <- function(method, sfu_object) {
    sapply(effect_sizes, function(es) Power_calculation(n = n, reps = reps, K = K, randproc = method, sfu = sfu_object, effect_size = es))
  }
  
  # Setup for calculation using both sfu and sfldOF objects, mapping them to labels
  methods_and_sfu_labels <- expand.grid(method = c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), label = names(sfu_values))
  
  # Add the calculation for power using power_inverse_normal
  calculate_power_inverse_normal <- function(effect_sizes) {
    sapply(effect_sizes, function(es) power_inverse_normal(effect_size = es, n = n, K = K))
  }
  
  
  data_inverse_normal <- data.frame(
    effect_size = effect_sizes,
    power = calculate_power_inverse_normal(effect_sizes),
    method = "Inverse Normal",  
    sfu_label = "Pocock" 
  )
  
  # Calculate power for each combination
  data_list <- lapply(1:nrow(methods_and_sfu_labels), function(i) {
    row <- methods_and_sfu_labels[i, ]
    method <- row$method
    sfu_label <- row$label
    sfu_object <- sfu_values[[sfu_label]]
    data.frame(effect_size = effect_sizes, power = calculate_power_for_method(method, sfu_object), method = method, sfu_label = sfu_label)
  })
  # Combine all data frames
  data_to_plot <- do.call(rbind, data_list)
  data_to_plot<- rbind(data_to_plot, data_inverse_normal)
  
  # Use a color-blind-friendly palette with distinct colors
  color_palette <- RColorBrewer::brewer.pal(length(unique(data_to_plot$method)), "Dark2")
  
  # Plot the results with improved aesthetics
  ggplot(data_to_plot, aes(x = effect_size, y = power, color = interaction(method, sfu_label), group = interaction(method, sfu_label))) +
    geom_line(size = 0.5) +
    geom_point(size = 1, alpha = 0.8, aes(shape = method)) +
#    scale_color_manual(values = color_palette) +
#    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5)) +
    labs(title = paste("Power by Effect size for ", sides, "-sided ", "Pocock w alpha=", alpha, " n=", n, ", K=", K) ,
         x = "Effect Size",
         y = "Power",
         color = "Method & SFU",
         shape = "Method") +
    theme_minimal() +
    theme(legend.position = "right")
}
Plot_power(n=24, reps=10, K=3, sfu="Pocock", sides=1, alpha=0.025)

