source("Plot_and_Save.R")

# Manual calculation for validation
#CR_manual = simulateClinicalTrial_manual(method="naive", n_sim=reps, n=24,K=2, randproc="CR")
#print(paste("T1E for manual calculation with CR:", CR_manual))
#RAR_manual = simulateClinicalTrial_manual(method="naive", n_sim=reps, n=24,K=2, randproc="RAR")
#print(paste("T1E for manual calculation with RAR:", RAR_manual))


# Calculate and save

#power_save_to_excel(n=120,n_sim=1000,K=3, futility=FALSE,futility_binding=FALSE)

T1E_save_to_excel(n=24,n_sim=100000,K=3, futility=FALSE,futility_binding=FALSE)

# Call the function to save the data and plot the line plot
#Power_diff_sample_sizes(n_sim =1000)

