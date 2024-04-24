source("standards.R")
library(randomizeR)
library(rpact)
library(gsDesign) 
library(mvtnorm)    
library(openxlsx)
library(parallel)


# Calculates the T1E probability based on Monte Carlo simulations of the patients and randomization sequences
# currently only with Pocock
simulateClinicalTrial_manual <- function(method, n_sim, n, K, randproc) {
  if (method=="naive") {     # naive calculcation of boundaries of gsd without updates of boundaries
    # Create group sequential design
    testdesign = gsDesign(k = K, test.type = 1, sfu = "Pocock")
    n=n
    # Initialize vectors to store results
    means1 <- vector("list", K)
    means2 <- vector("list", K)
    stDevs1 <- vector("list", K)
    stDevs2 <- vector("list", K)
    sequences <- vector("list", n_sim)
    
    testdecision <- logical(n_sim) # Initialize with the number of simulations
    testdecision_comb <- logical()
    
    # Randomization object
    randobj <- switch(as.character(randproc),
                      "CR" = crPar(n, K = 2),
                      "RAR" = rarPar(n, K = 2, groups = c("0", "1")),
                      "BSD" = bsdPar(n, mti = mti, groups = c("0", "1")),
                      "EBC" = ebcPar(n, p, groups = c("0", "1")),
                      "CHEN" = chenPar(n, mti = mti, p = p, groups = c("0", "1")),
                      "PBR" = rpbrPar(n, rb = rb, groups = c("0", "1")),
                      "MP" = mpPar(n, mti = mti, ratio = c(1, 1)),
                      stop("Invalid randproc parameter.")
    )
    randseq = genSeq(randobj, n_sim) 
    
    
    for (i in 1:n_sim) {
      sequence =  randseq@M[i,]               # Generate randomisation sequence
      patients <- rnorm(n)       # Generate N(0,1) distributed random variables
      
      
      n_group1=vector("numeric", K)   # Patients allocated to group 1 in each stage
      n_group2=vector("numeric", K)   # Patients allocated to group 2 in each stage
      means1 <- vector("numeric", K)  # cumulative means for group 1 in each stage
      means2 <- vector("numeric", K)  # cumulative means for group 2 in each stage
      I <- vector("numeric", K)       # Information fraction for each stage
      t <- vector("numeric", K)       # Test statistic for each stage
      continue <- vector("list", K)   # Continue/Reject decision for each stage
      
      
      
      for (j in 1:K) {                # for each stage
        n_group1[j]=sum(sequence[1:(n*j/K)])                                                 # Calculate number of patients in stage j
        n_group2[j]=(n*j/K)-sum(sequence[1:(n*j/K)])
        
        if (n_group1[j]==0 | n_group2[j]==0) {
          continue[j] = TRUE
        } else {
        means1[[j]] <-sum(patients[1:(n*j/K)]*sequence[1:(n*j/K)])/(n_group1[j])
        means2[[j]] <-sum(patients[1:(n*j/K)]*(1-sequence[1:(n*j/K)]))/(n_group2[j])
        I[j] = 1 * 1/sqrt(1/n_group1[j]+1/n_group2[j])
        t[j] = I[j]*(means1[j]-means2[j])
        continue[j]= (t[j] < testdesign$upper$bound[j])
        }
       }
      testdecision[i] =all(unlist(continue))
      
      
    }
    percentage_false <- (1 - mean(testdecision))
    return(percentage_false)
  }
  #  else { # hier habe ich das einmal für den inverse normal combination test probiert, wird aber aktuell nicht benötigt.
  #   
  #   testdesign = gsDesign(k = K, test.type = 1, sfu = "Pocock")
  #   
  #   testdesign_comb = gsDesign(k = K, test.type = 1, sfu = "Pocock", n.I=c(1/sqrt(2), 1))
  #   #  I_comb =  c(sqrt(1/sqrt(2)), 1)
  #   #  testdesign_comb=  gsDesign(k = K, test.type = 2, sfu = "Pocock", n.I = I_comb)
  #   #  print(testdesign_comb$upper$bound)
  #   #  print(testdesign$upper$bound)
  #   # Create group sequential design
  #   design <- getDesignGroupSequential(
  #     typeOfDesign=c("OF"),
  #     alpha = 0.05,
  #     sided = 2,
  #     informationRates = c(0.33,0.67,1)
  #   )
  #   #  print(design)
  #   #  design <- getDesignInverseNormal(kMax = 3)
  #   
  #   # Initialize vectors to store results
  #   means1 <- vector("list", K)
  #   means2 <- vector("list", K)
  #   stDevs1 <- vector("list", K)
  #   stDevs2 <- vector("list", K)
  #   sequences <- vector("list", n_sim)
  #   
  #   # randobj = bsdPar(24, mti = 3, groups = c("0", "1"))
  #   randobj = rpbrPar(24, rb = 2, groups = c("0", "1"))
  #   
  #   randseq = genSeq(randobj, n_sim) 
  #   testdecision <- logical()
  #   testdecision_comb <- logical()
  #   
  #   for (i in 1:n_sim) {
  #     sequence =  randseq@M[i,]               # Generate randomisation sequence
  #     patients <- rnorm(n)       # Generate N(0,1) distributed random variables
  #     
  #     n_group1=vector("numeric", K)
  #     n_group2=vector("numeric", K)
  #     n_group1_comb=vector("numeric", K)
  #     n_group2_comb=vector("numeric", K)
  #     means1 <- vector("numeric", K)
  #     means2 <- vector("numeric", K)
  #     I <- vector("numeric", K)
  #     t <- vector("numeric", K)
  #     continue <- vector("list", K)
  #     p <- vector("numeric", K)
  #     t_combination <- vector("numeric", K)
  #     continue_combination <- vector("list", K)
  #     means1_comb <- vector("numeric", K)
  #     means2_comb <- vector("numeric", K)
  #     t_comb <- vector("numeric", K)
  #     
  #     
  #     
  #     for (j in 1:K) {
  #       n_group1[j]=sum(sequence[1:(n*j/K)])
  #       n_group2[j]=(n*j/K)-sum(sequence[1:(n*j/K)])
  #       means1[[j]] <-sum(patients[1:(n*j/K)]*sequence[1:(n*j/K)])/(n_group1[j])
  #       means2[[j]] <-sum(patients[1:(n*j/K)]*(1-sequence[1:(n*j/K)]))/(n_group2[j])
  #       I[j] = 1 * 1/sqrt(1/n_group1[j]+1/n_group2[j])
  #       t[j] = I[j]*(means1[j]-means2[j])
  #       continue[j]= (t[j] < testdesign$upper$bound[j])
  #       
  #       
  #       # for inverse normal test
  #       n_group1_comb[j]=sum(sequence[(1+n*(j-1)/K):(n*j/K)])
  #       n_group2_comb[j]=(n/K)-sum(sequence[(1+n*(j-1)/K):(n*j/K)])
  #       means1_comb[j] <-sum(patients[(1+n*(j-1)/K):(n*j/K)]*sequence[(1+n*(j-1)/K):(n*j/K)])/(n_group1_comb[j])
  #       means2_comb[j] <-sum(patients[(1+n*(j-1)/K):(n*j/K)]*(1-sequence[(1+n*(j-1)/K):(n*j/K)]))/(n_group2_comb[j])
  #       t_comb[j] = I[j]*(means1_comb[j]-means2_comb[j])
  #       p[j] = pnorm(q=t_comb[j], mean=0, sd=1, lower.tail=FALSE)
  #       
  #       
  #     }
  #     testdecision =c(testdecision, all(unlist(continue)))
  #     
  #     
  #     
  #     # for inverse normal test
  #     # p2 sollte aber nur die daten von p2 beinhalten
  #     t_combination[1]=qnorm((1-p[1]))
  #     t_combination[2] = (1/sqrt(2))*t_combination[1]+(1/sqrt(2))*qnorm((1-p[2]))
  #     #print("pvalue")
  #     #print(p)
  #     continue_combination[1]= (t_combination[1] < testdesign_comb$upper$spend[1]) ##### ??
  #     continue_combination[2]= (t_combination[2] < testdesign_comb$upper$bound[2])
  #     testdecision_comb =c(testdecision_comb, all(unlist(continue_combination)))
  #     
  #     
  #     
  #     
  #     #print("fehlerkontrolle")
  #     #print(t_combination)
  #     #print(testdecision)
  #     #print(testdecision_comb)
  #     
  #     #design <- getDesignInverseNormal(typeOfDesign="P", kMax=2)
  #     #dataset <- getDataset(means1=unlist(means1_comb), means2=unlist(means2_comb), stDevs1=c(1,1),stDevs2=c(1,1), n1=n_group1, n2=n_group2)
  #     #results = getStageResults(design = design, dataInput = dataset)
  #     #print(results)
  #   }
  #   percentage_false <- (1 - mean(testdecision))
  #   print(percentage_false)
  #   percentage_false_comb <- (1 - mean(testdecision_comb))
  #   print("Rejection prob comb test:")
  #   print(percentage_false_comb)
  #   #  print(continue_combination)
  #   
  #   
  # }
}
#simulateClinicalTrial_manual(method="naive", n_sim=8, n=24,K=3)


## Calculation of Power manually for 2-stage group sequential design via MVN distribution 
# Function currently not in use
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


#### Calculation of T1E for inverse normal combination test 
power_inverse_normal <- function(effect_size, n, K, sfu, onesided="yes", alpha=0.025) {
  mu <- numeric(K)
  cov <- matrix(0, K, K)
  for (i in 1:K) {
    mu[i] <- effect_size * sqrt(n / K * i) * 0.5
    for (j in 1:K) {
      if (i == j) {
        cov[i, j] <- 1
      } else {
        cov[i, j] <- sqrt(min(i, j)) / sqrt(max(i, j))
      }
    }
  }
  print(cov)
    testdesign = gsDesign(k = K, test.type = 1, sfu = sfu)
    upper_bound = testdesign$upper$bound[1:K]

  integral = pmvnorm(algorithm = Miwa(), lower = -Inf, upper = upper_bound, 
                     mean = mu, sigma = cov, abseps = 1e-16)
  beta = integral
  beta = round(beta, digits=5)
  return(1-beta)
}
#power_inverse_normal(effect_size=0.3,n=80, K=2, sfu="Pocock")


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
  randobj <- switch(as.character(randproc),
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

  
  # Set up parallel computing
#  cl <- makeCluster(detectCores() - 1) # Create a cluster using one less than the number of cores available
  cl <- makeCluster(4) # Create a cluster using one less than the number of cores available
  
  # Export needed variables and functions to the cluster
  clusterExport(cl, varlist=c("seq", "n", "K", "randproc", "sfu", "sides", "alpha", "rb", "mti", "p", "effect_size", "futility"), envir=environment())
  
  # Load required libraries on all cluster nodes. Also source the standards.R script.
  clusterEvalQ(cl, {
    source("standards.R") # Make sure the path to standards.R is correct and accessible from the cluster nodes
    library(randomizeR)
    library(rpact)
    library(gsDesign)
    library(mvtnorm)
    library(openxlsx)
  })
  

     sfu_name=deparse(substitute(sfu)) 
    
  
  
  # Execute the core of Power_calculation inw parallel
  Power <- unlist(parLapply(cl, 1:reps, function(j) {
# Power <- unlist(lapply(1:reps, function(j) {
      
 #   for (j in 1:reps) {
#  seq@M = matrix(c(1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1), nrow=1, byrow=TRUE)  #   Hardcoded sequence for testing

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
      if (is.function(sfu)) {   # TODO: hier Abfrage ob zweiseitig oder einseitig!
        testdesign = gsDesign(k=K, test.type = 2 , sfu = sfu, alpha= alpha, n.I=I)     # alpha is always one-sided in gsDesign, even if 2-sided test design is selected.
      } else  {
        if (sfu=="corPOC") {
          sfu = "Pocock"
          testdesign = gsDesign(k=K, test.type = 2 , sfu = sfu, alpha= alpha, n.I=I)  
        } else if (sfu=="corOF") {
          sfu="OF"
          testdesign = gsDesign(k=K, test.type = 2 , sfu = sfu, alpha= alpha, n.I=I)  
        } else {
        testdesign = gsDesign(k=K, test.type = 2 , sfu = sfu, alpha= alpha)     # alpha is always one-sided in gsDesign, even if 2-sided test design is selected.
        }
      }
      lower_bound =testdesign$lower$bound
      upper_bound = testdesign$upper$bound
    }
    if (sides == 1) {            # for one-sided test set lover boundary to -infty
      lower_bound=rep(-99, K)
    }

    for (i in (1:K)) {                         # if there have been no allocations to one group jump to next stage
      if(n_A[[i+1]] == 0 | n_B[[i+1]] == 0) {
        lower_bound[i] <- -99     # TODO: Grenzen neu berechnen basieren auf einer interim analyse weniger, wenn allokation 0 zu X
        upper_bound[i] <-  99     # mitlaufen lassen wie oft das vorkommt, bei inverse normal nochmal drüber nachdenken
      }
    }
    if (futility == "yes") {                   # (binding) futility boundary 
      lower_bound=rep(0, K)   # TODO obere Grenzen anpassen für binding futility
    }
 
    zbdy <- rbind(lower_bound, upper_bound)
    results <- gst1(r, na=K, inf=I, zbdy, theta=effect_size)  # implemented by CJ
    Power[j] =results[[2]]
  }))

   # Clean up
  stopCluster(cl)
 #write.xlsx(Power, file = paste("data/Power_n", n, "_K", K, "_randproc", randproc, "_sfu", sfu_name, "_reps", reps, "_delta", effect_size,  ".xlsx"), rowNames = FALSE)
  return(mean(Power))
}
#Power_calculation(n=24, reps=10, randproc="CR", sfu=sfLDOF, K=2, effect_size=0, sides=1)
#Power_calculation(n=24, reps=100, randproc="PBR", sfu=sfLDPocock, K=3, effect_size=0.5, sides=1)


# Plots the power for LDM and naive Pocock and OF design for different randomization procedures
Plot_power <- function(n, reps, K, sfu, alpha=0.025, sides=1) {
  
  # Generate a sequence of effect_size values between 0 and 2
  effect_sizes <- seq(0, 2, by = 0.2)
  
  # Define sfu values and their corresponding plotting labels
  sfu_values <- list("POC" = sfu, "spendingPocock" = sfLDPocock)  # Assume sfu is a string or object, and sfldOF is an object
  
  # Function to calculate power for a given method and sfu object
  calculate_power_for_method <- function(method, sfu_object) {
    sapply(effect_sizes, function(es) mean(Power_calculation(n = n, reps = reps, K = K, randproc = method, sfu = sfu_object, effect_size = es)))
  }
  
  # Setup for calculation using both sfu and sfldOF objects, mapping them to labels
  methods_and_sfu_labels <- expand.grid(method = c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), label = names(sfu_values))
  
  # Add the calculation for power using power_inverse_normal
  calculate_power_inverse_normal <- function(effect_sizes) {
    sapply(effect_sizes, function(es) power_inverse_normal(effect_size = es, n = n, K = K))
  }
  
  
#  data_inverse_normal <- data.frame(
#    effect_size = effect_sizes,
#    power = calculate_power_inverse_normal(effect_sizes),
#    method = "Inverse Normal",  
#    sfu_label = "Pocock" 
#  )
  
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
 # data_to_plot<- rbind(data_to_plot, data_inverse_normal)
  
  print(data_to_plot)
  
  write.xlsx(data_to_plot, file = paste("data_to_plotPOC", n, ".xlsx"), rowNames = FALSE)
  
  # Use a color-blind-friendly palette with distinct colors
#  color_palette <- RColorBrewer::brewer.pal(length(unique(data_to_plot$method)), "Dark2")
  
  # Plot the results with improved aesthetics
  ggplot(data_to_plot, aes(x = effect_size, y = power, color = interaction(method, sfu_label), group = interaction(method, sfu_label))) +
    geom_line(linewidth = 0.03) + # Thicker lines for better visibility
    geom_point(size = 0.01, alpha = 0.3) + # Larger, semi-transparent points
 #   scale_color_manual(values = color_palette) +
#    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5)) +
    labs(title = paste("Power by Effect size for ", sides, "-sided ", "Pocock w alpha=", alpha, " n=", n, ", K=", K) ,
         x = "Effect Size",
         y = "Power",
         color = "Method & SFU",
         shape = "Method") +
    theme_minimal() +
    theme(legend.position = "right")
}




# Calculates the power for different randomization procedures with the following boundaries:
# - Pocock / OF
# - LanDeMets alpha spending Pocock / OF
# - Pocock / Of with correction for actual information
# and saves it as an excel file.
save_to_excel <- function(n, reps, K, sides = 1, alpha = 0.025, rb = 4, mti = 3, p = 2/3, futility = "no") {
  effect_sizes <- seq(0, 2, by = 0.2)
  
  # Define sfu values and their corresponding plotting labels
  sfu_values <- list("POC" = "Pocock", "LDMPOC" = sfLDPocock, "corPOC" = "corPOC", "OF" = "OF", "LDMOF" = sfLDOF, "corOF" = "corOF" )

  # Function to calculate power for a given method and sfu object
  calculate_power_for_method <- function(method, sfu_object, es) {
    Power_calculation(n = n, reps = reps, K = K, randproc = method, sfu = sfu_object, effect_size = es)
  }

  # Create a new Excel workbook
  wb <- createWorkbook()
  
  # Process each effect size separately
  for (es in effect_sizes) {
    # Initialize a list to collect data frames for this effect size
    data_frames <- list()
    
    # Loop over each sfu value and method to calculate power
    for (sfu_label in names(sfu_values)) {
      for (method in c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN")) {
        sfu_object <- sfu_values[[sfu_label]]
        power_value <- calculate_power_for_method(method, sfu_object, es)
        # Create a data frame with columns in the correct order
        data_frames[[length(data_frames) + 1]] <- data.frame(
          Effect_size = es,
          Power = power_value,
          RP = method,
          GSD = sfu_label
        )
      }
    }
    
    # Combine all data frames into one for this effect size
    data_combined <- do.call(rbind, data_frames)
    
    # Add a worksheet for this effect size
    sheet_name <- paste("EffectSize", format(es, nsmall = 1), sep = "_")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data_combined)
  }
  # Save the workbook to a file
  filepath <- paste("data/ResultsPower_n", n, "_K", K, "_reps", reps, ".xlsx")  # Define the file path where you want to save the Excel file
  
  saveWorkbook(wb, filepath, overwrite = TRUE)
  print(paste("Workbook saved to", filepath))
}
#save_to_excel(n=24,reps=1,K=3)


library(readxl)
library(ggplot2)
library(dplyr)
plot_effect_size_flexible = function(file_path, rp_values = c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), gsd_values) {
  # List all sheets
  sheets <- excel_sheets(file_path)
  
  # Initialize a list to store combined data
  all_data <- list()
  
  # Loop through each RP and GSD combination
  for (rp in rp_values) {
    for (gsd in gsd_values) {
      # Initialize a list to store data for current RP and GSD combination
      data_list <- list()
      
      # Loop through each sheet
      for (sheet in sheets) {
        # Read the sheet
        data <- read_excel(file_path, sheet = sheet)
        
        # Filter for the specific RP and GSD values
        filtered_data <- filter(data, RP == rp, GSD == gsd)
        
        # Check if there is any data left after filtering
        if (nrow(filtered_data) > 0) {
          # Store the data for plotting with additional columns for RP and GSD
          filtered_data$RP <- rp  # Add RP column to identify groups in plot
          filtered_data$GSD <- gsd  # Add GSD column
          data_list[[sheet]] <- filtered_data[c("Effect_size", "Power", "RP", "GSD")]
        }
      }
      
      # Combine data for the current RP and GSD
      if (length(data_list) > 0) {
        all_data[[paste(rp, gsd)]] <- do.call(rbind, data_list)
      }
    }
  }
  
  # Combine all data into a single dataframe for plotting
  if (length(all_data) > 0) {
    plot_data <- do.call(rbind, all_data)
    
    # Generate the plot
    ggplot(plot_data, aes(x = Effect_size, y = Power, color = interaction(RP, GSD), group = interaction(RP, GSD))) +
      geom_point() +
      geom_line() +
      theme_minimal() +
      scale_color_discrete(name = "RP and GSD", labels = function(x) gsub("\\.", ", ", x)) +
      labs(title = "Effect Size and Power Relationships for Multiple Group Sequential Designs (GSDs) and Randomization Procedures (RPs)",
           x = "Effect Size",
           y = "Power")
  } else {
    print("No data met the filter criteria for the specified RP and GSD values.")
  }
}

# Example usage
gsd_values <- c("corPOC", "POC")  # Multiple GSD values
plot_effect_size_flexible("data/ResultsPower_n 24 _K 2 _reps 8000 .xlsx", rp_values=c("CR"), gsd_values = gsd_values)


sqrt(0.025*(1-0.025)/20000000)

test =gsDesign(k=10, test.type = 2 , sfu = "Pocock")  

test$upper$bound

test =gsDesign(k=10, test.type = 1 , sfu = "Pocock")  
test$upper$bound
