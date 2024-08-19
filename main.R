source("RandomizationGSD.R")


reps = 100000000
# Manual calculation for validation
#CR_manual = simulateClinicalTrial_manual(method="naive", n_sim=reps, n=24,K=2, randproc="CR")
#print(paste("T1E for manual calculation with CR:", CR_manual))
#RAR_manual = simulateClinicalTrial_manual(method="naive", n_sim=reps, n=24,K=2, randproc="RAR")
#print(paste("T1E for manual calculation with RAR:", RAR_manual))


# Calculate and save

#save_to_excel(n=24,reps=1,K=2)



save_to_excel(n=24,n_sim=2000,K=3)
