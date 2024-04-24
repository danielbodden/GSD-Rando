source("RandomizationGSD.R")


reps = 100000000
# Manual calculation for validation
#CR_manual = simulateClinicalTrial_manual(method="naive", n_sim=reps, n=24,K=2, randproc="CR")
#print(paste("T1E for manual calculation with CR:", CR_manual))
#RAR_manual = simulateClinicalTrial_manual(method="naive", n_sim=reps, n=24,K=2, randproc="RAR")
#print(paste("T1E for manual calculation with RAR:", RAR_manual))


# Calculate and save

save_to_excel(n=24,reps=1,K=2)




save_to_excel(n=96,reps=8000,K=2)


# Über Validierung prüfen ob das genauer ist
# dann für höheres n
# dann für K=3 probieren auf wie viele Wiederholungen ich komme
# mal Tasks hochschrauben und schauen ob das schneller läuft