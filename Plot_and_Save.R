source("RandomizationGSD.R")
library(writexl)


# Calculates the power for different randomization procedures with the following boundaries:
# - Pocock / OF
# - LanDeMets alpha spending Pocock / OF
# - Pocock / Of with correction for actual information
# _ inverse normal combination function Pocock / OF boundaries
# and saves it as an excel file.
power_save_to_excel <- function(n, n_sim, K, sides = 1, alpha = 0.025, rb = 4, mti = 3, p = 2/3, futility = FALSE, futility_binding = FALSE) {
  deltas <- seq(0, 2, by = 0.2)           # Grid for effect sizes

  # List of group sequential designs used
 # gsd <- list("POC" = "Pocock", "LDMPOC" = "LDMPocock", "coRPOC" = "coRPOC", "OF" = "OF", "LDMOF" = "LDMOF", "corOF" = "corOF", "IVNOF" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")
  # only important
  gsd <- list("LDMPOC" = "LDMPocock", "LDMOF" = "LDMOF")
  gsd <- list("IVNOF" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")
  
  
  # Function to calculate power for a given sfu and sfu object
  calculate_power_for_sfu <- function(rp, gsd_object, es) {
    Power_list = Power_condMVN(n = n, n_sim = n_sim, K = K, RP = rp, sfu = gsd_object, delta = es, futility = futility, futility_binding = futility_binding)
    Power_mean = round(mean(Power_list$Pow), 4)
    return(list(Pow = Power_mean, Counter = Power_list$Counter))
    }
  calculate_power_for_inverse_normal <- function(rp, gsd_object, es) {
    Power_list = Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = rp, sfu=gsd_object, delta = es, futility = futility, futility_binding = futility_binding)
    Power_mean = round(mean(Power_list$Pow), 4)
    return(list(Pow = Power_mean, Counter = Power_list$Counter))
  }
  # Create a new Excel workbook
  wb <- createWorkbook()
  
  # Process each effect size separately
  for (es in deltas) {
    # Initialize a list to collect data frames for this effect size
    data_frames <- list()
    
    # Loop over each sfu value and sfu to calculate power
    for (gsd_label in names(gsd)) {
      for (RP in c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN")) {
        if (gsd_label == "IVNOF") {
          power_value <- calculate_power_for_inverse_normal(RP, gsd_object = "OF", es)
        } else if (gsd_label == "IVNPOC") {
          power_value <- calculate_power_for_inverse_normal(RP, gsd_object = "Pocock", es)
        } else {
          gsd_object <- gsd[[gsd_label]]
          power_value <- calculate_power_for_sfu(RP, gsd_object, es)
        }
        print(RP)
        
        # Create a data frame with columns in the correct order
        data_frames[[length(data_frames) + 1]] <- data.frame(
          delta = es,
          Power = power_value$Pow,
          RP = RP,
          GSD = gsd_label,
          skipped_sequences = power_value$Counter
        )
      }
      print(gsd_label)
    }
    # Combine all data frames into one for this effect size
    data_combined <- do.call(rbind, data_frames)
    
    # Add a worksheet for this effect size
    sheet_name <- paste("EffectSize", format(es, nsmall = 1), sep = "_")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data_combined)
  }
  # Save the workbook to a file
  futility_string <- as.character(futility)
  futility_binding_string <- as.character(futility_binding)
  
  filepath <- paste("data/ResultsPower n", n, "K", K, "n_sim", n_sim, "fut", futility_string, "futbind", futility_binding_string, ".xlsx")  
  saveWorkbook(wb, filepath, overwrite = TRUE)
  print(paste("Workbook saved to", filepath))
}#
#power_save_to_excel(n=6,n_sim=5,K=3)


# Calculates the T1E for different randomization procedures with the following boundaries 
# - Pocock / OF
# - LanDeMets alpha spending Pocock / OF
# - Pocock / Of with correction for actual information
# _ inverse normal combination function Pocock / OF boundaries
# and saves it as an Excel file.
# Additionally the amount of skipped randomization sequences is being saved, which happens when:
# For Pocock, OF and LDM: if no allocations to one group in the first stage
# For INV: if no allocations to one group in any stage
T1E_save_to_excel <- function(n, n_sim, K, sides = 1, alpha = 0.025, rb = 4, mti = 3, p = 2/3, futility = FALSE, futility_binding = FALSE) {
  RP_values <- c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN")
  
  # List of group sequential designs used
#  gsd <- list("POC" = "Pocock", "LDMPOC" = "LDMPocock", "coRPOC" = "coRPOC", "OF" = "OF", "LDMOF" = "LDMOF", "corOF" = "corOF", "IVNOF" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")
  # only some:
  gsd <- list("LDMPOC" = "LDMPocock", "LDMOF" = "LDMOF", "IVNOF" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")
  
  
  # Function to calculate T1E for a given randomization procedure (rp) and gsd (gsd_object)
  calculate_power_for_sfu <- function(rp, gsd_object, es) {
    Power_condMVN(n = n, n_sim = n_sim, K = K, RP = rp, sfu = gsd_object, delta = 1, futility = futility, futility_binding = futility_binding)
  }
  calculate_power_for_inverse_normal <- function(rp, gsd_object, es) {
    Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = rp, sfu=gsd_object, delta = 1, futility = futility, futility_binding = futility_binding)
  }

  wb <- createWorkbook()                              # Create a new Excel workbook
  addWorksheet(wb, sheetName = "Skipped_sequences")   # Add a worksheet for counting the skipped sequences
  counter_skipped <- data.frame(gsd = character(), RP = character(), Value = numeric())
  
  # Loop over each GSD to calculate T1E
  for (gsd_label in names(gsd)) {
    cat("Processing GSD:", gsd_label, "\n")
    
    # Loop over each randomisation procedure to calculate T1E
    for (RP in RP_values) {
      cat("Processing RP:", RP, "\n")
      if (gsd_label == "IVNOF") {
        power_value <- calculate_power_for_inverse_normal(RP, gsd_object = "OF")
      } else if (gsd_label == "IVNPOC") {
        power_value <- calculate_power_for_inverse_normal(RP, gsd_object = "Pocock")
      } else {
        gsd_object <- gsd[[gsd_label]]
        power_value <- calculate_power_for_sfu(RP, gsd_object)
      }
      
      # Create a data frame for the power values
      data_frame <- data.frame(
        Index = 1:length(power_value$Pow),
        gsd = gsd_label,
        RP = RP,
        Power = power_value$Pow
      )
      
      # Create a unique sheet name based on RP and gsd_label
      sheet_name <- paste(gsd_label, RP, sep = "_")
      
      # Add a new sheet to the workbook with the name of the current RP and gsd_label
      addWorksheet(wb, sheetName = sheet_name)
      
      # Add a row to counter_skipped for this iteration
      counter_skipped[nrow(counter_skipped) + 1, ] <- list(gsd=gsd_label, RP=RP, skipped=power_value$Counter)
      
      # Write the data to the current sheet
      writeData(wb, sheet = sheet_name, x = data_frame)
    }
    writeData(wb,sheet="Skipped_sequences",x=counter_skipped) # FIll worksheet with the amount of skipped sequences
  }
  # Save the workbook to a file
  futility_string <- as.character(futility)
  futility_binding_string <- as.character(futility_binding)
  
  # Save the workbook to the specified file
  filepath <- paste("data/Pow_RS n", n, " K", K, " n_sim", n_sim, "fut", futility_string, "futbind", futility_binding_string, ".xlsx")  

  saveWorkbook(wb, filepath, overwrite = TRUE)
  cat("Workbook saved to", filepath, "\n")
}


# Read in the saved excel file with T1E for the RPs and return a dataframe 
read_all_sheets_into_one_df <- function(file_path) {
  sheet_names <- excel_sheets(file_path)      # Get the sheet names
  combined_df <- data.frame()                   # Initialize an empty dataframe
  
  for (sheet in sheet_names[-1]) {                 # Loop through each sheet and read it into a dataframe
    df <- read_excel(file_path, sheet = sheet)
    colnames(df) <- c("Index", "gsd", "RP", "Power")      # Assign the specified column names
    combined_df <- rbind(combined_df, df)                # Combine the dataframes
  }
  return(combined_df)
}


create_boxplot <- function(file_path) {
  data <- read_all_sheets_into_one_df(file_path)        # Read data from the Excel file
  data <- subset(data, gsd == "POC" | gsd == "OF")      # Only for Pocock and OF
  # Ensure 'gsd' and 'RP' are treated as factors and ordered
  data$gsd <- factor(data$gsd, levels = unique(data$gsd))
  data$RP <- factor(data$RP, levels = unique(data$RP))
  
  # Create the interaction term with a specific order
  data$interaction_term <- with(data, interaction(gsd, RP, sep = ":", lex.order = TRUE))
  
  # Ensure the levels of interaction_term are ordered by gsd first, then by RP
  interaction_levels <- with(data, levels(interaction(gsd, RP, sep = ":", lex.order = TRUE)))
  data$interaction_term <- factor(data$interaction_term, levels = interaction_levels[order(as.numeric(data$gsd), as.numeric(data$RP))])
  
  # Round the Power values
  data$Power <- round(data$Power, 4)


  
  
  # Create the boxplots with the specified title
  p <- ggplot(data, aes(x = interaction_term, y = "Type I Error")) +
    geom_boxplot() +
    labs(x = "GSD and RP Combination", y = "Power") +
    ggtitle("T1E for GSD with n=24 and K=3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = levels(data$interaction_term)) 
 #   ylim(0.024, NA)  # Set lower limit of y-axis to 0.0225
  print(p)    # Display the plot
}

# Example usage
#file_path <- "data/T1E n 24  K 3  n_sim 1000 fut TRUE futbind TRUE .xlsx"
#create_boxplot(file_path)

create_violin_plot <- function(file_path) {
  data <- read_all_sheets_into_one_df(file_path)        # Read data from the Excel file
  data <- subset(data, gsd == "POC" | gsd == "OF")      # Only for Pocock and OF
  # Ensure 'gsd' and 'RP' are treated as factors and ordered
  data$gsd <- factor(data$gsd, levels = unique(data$gsd))
  data$RP <- factor(data$RP, levels = unique(data$RP))
  
  # Create the interaction term with a specific order
  data$interaction_term <- with(data, interaction(gsd, RP, sep = ":", lex.order = TRUE))
  
  # Ensure the levels of interaction_term are ordered by gsd first, then by RP
  interaction_levels <- with(data, levels(interaction(gsd, RP, sep = ":", lex.order = TRUE)))
  data$interaction_term <- factor(data$interaction_term, levels = interaction_levels[order(as.numeric(data$gsd), as.numeric(data$RP))])
  
  # Round the Power values
  data$Power <- round(data$Power, 4)
  
  # Create the violin plot with the specified title
  p <- ggplot(data, aes(x = interaction_term, y = Power, fill = interaction_term)) +
    geom_violin(width=1, alpha=.5) +
    stat_summary(fun = "mean",
                 geom = "crossbar", 
                 width = 0.2,
                 colour = "red") +
    labs(x = "GSD and RP Combination", y = "Type I Error") +
    ggtitle("T1E for GSD with n=24 and K=3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_minimal() + theme(legend.position = "none") +
    scale_x_discrete(limits = levels(data$interaction_term)) 
  #   ylim(0.024, NA)  # Set lower limit of y-axis to 0.0225
  print(p)    # Display the plot
}


# Example usage
#file_path <- "data/T1En24K3n_sim8000 .xlsx"
#create_violin_plot(file_path)


PlotPower = function(file_path, RP_values = c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), sfu) {
  sheets <- excel_sheets(file_path)         # List all sheets
  all_data <- list()                          # Initialize a list to store combined data
  
  for (rp in RP_values) {                     # Loop through each RP and GSD combination
    for (gsd in sfu) {
      data_list <- list()                     # Initialize a list to store data for current RP and GSD combination
      
      for (sheet in sheets) {                   # Loop through each sheet
        data <- read_excel(file_path, sheet = sheet)
        filtered_data <- filter(data, RP == rp, GSD == gsd)       # Filter for the specific RP and GSD values
        
        if (nrow(filtered_data) > 0) {              # Check if there is any data left after filtering
          # Store the data for plotting with additional columns for RP and GSD
          filtered_data$RP <- rp  # Add RP column to identify groups in plot
          filtered_data$GSD <- gsd  # Add GSD column
          data_list[[sheet]] <- filtered_data[c("delta", "Power", "RP", "GSD")]
        }
      }
      if (length(data_list) > 0) {        # Combine data for the current RP and GSD
        all_data[[paste(rp, gsd)]] <- do.call(rbind, data_list)
      }
    }
  }
  if (length(all_data) > 0) {                 # Combine all data into a single dataframe for plotting
    plot_data <- do.call(rbind, all_data)
    print(plot_data)
    
    # Generate the plot
    ggplot(plot_data, aes(x = delta, y = Power, color = interaction(RP, GSD), group = interaction(RP, GSD))) +
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


#sfu <- c("LDMOF")  # Multiple GSD values
#PlotPower("data/ResultsPower n 24 K 4 n_sim 100 fut FALSE futbind FALSE .xlsx", RP_values=c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), sfu = sfu)

#sfu <- c("Pocock")  # Multiple GSD values
#PlotPower("data/ResultsPower_n 24 _K 2 _reps 8000 .xlsx", RP_values=c("CR"), sfu = sfu)



#sfu <- c("IVNPOC")  # Multiple GSD values
#PlotPower("data/ResultsPower_n 24 _K 2 _n_sim 1000 .xlsx", RP_values=c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), sfu = sfu)



# IVNPOC has slighlty lower Power than Poc and LDMPOC for CR, but same power for PBR(4) as then no wrong information is used.
# I still have to describe coRPOC
# create Boxplots for RPs
# when creating th eboxplots, mention how many values are not shown (Removed 13 rows containing non-finite values (stat_boxplot()).)
# Problem with inverse  normal when a stage has only zero allocations -> in that case I remove the data from the patients completely.
# Other option would be to add it to another stage, but if it was in the last stage the test was alrady conducted.
# Recreate boxplots for n=5000 with new calculation



#####


# Diese Grafik wäre interessant für die Power, nicht für den T1E.

create_violin_plot_binding <- function() {
  data1 <- read_all_sheets_into_one_df("data/Pow_RS n 24  K 3  n_sim 1000 fut TRUE futbind TRUE .xlsx")        # Read data from the Excel file
  data1$fut <- "binding"
  data2 <- read_all_sheets_into_one_df("data/Pow_RS n 24  K 3  n_sim 1000 fut TRUE futbind FALSE .xlsx")        # Read data from the Excel file
  data2$fut <- "non_binding"
  data3 <- read_all_sheets_into_one_df("data/Pow_RS n 24  K 3  n_sim 1000 fut FALSE futbind FALSE .xlsx")        # Read data from the Excel file
  data3$fut <-"no_futility"
  data <- bind_rows(data1,data2,data3)
  data <- subset(data, gsd == "IVNOF")      # Only for Pocock and OF
  # Ensure 'gsd' and 'RP' are treated as factors and ordered
  data$gsd <- factor(data$gsd, levels = unique(data$gsd))
  data$RP <- factor(data$RP, levels = unique(data$RP))
  
  # Create the interaction term with a specific order
  data$interaction_term <- with(data, interaction(gsd, RP, sep = ":", lex.order = TRUE))
  
  # Ensure the levels of interaction_term are ordered by gsd first, then by RP
  interaction_levels <- with(data, levels(interaction(gsd, RP, sep = ":", lex.order = TRUE)))
  data$interaction_term <- factor(data$interaction_term, levels = interaction_levels[order(as.numeric(data$gsd), as.numeric(data$RP))])
  
  # Round the Power values
  data$Power <- round(data$Power, 4)
  
  data$fut <- factor(data$fut, levels=c("no_futility", "binding", "non_binding"))
  # Create the violin plot with the specified title
  p <- ggplot(data, aes(x = interaction_term, y = Power, fill = fut)) +
    geom_violin(width=1, alpha=.5) +
    stat_summary(aes(x = interaction_term, group = fut), fun = "mean",
                 geom = "crossbar", 
                 width = 0.2,
                 colour = "red", 
                 position = position_dodge(0.9)) +  
    labs(x = "GSD and RP Combination", y = "Power") +
    ggtitle("Power for GSD with n=24 and K=3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_minimal()  +
       coord_cartesian(ylim=c(0.58, 0.68)) 

  
  print(p)    # Display the plot
}


#create_violin_plot_binding()




# Different sample sizes

# Sample size on x axis
# Power on y axis
# for one effect size
# Function to save Power_condMVN results for n=12 to n=30 for both sfu = "LDMOF" and "IVNOF" and create a line plot
Power_diff_sample_sizes <- function(n_sim = 1, K = 3, delta = 1) {
  
  # Initialize an empty dataframe to store results for both "LDMOF" and "IVNOF"
  results_df <- data.frame(n = integer(), power = numeric(), sfu = character())
  
  for (n in seq(8, 30, by = 1)) {
    # Call Power_condMVN with sfu = "LDMOF"
    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = 2, RP = "PBR", sfu = "LDMOF", delta = delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_PBR_K2"))
    
    # Call Power_condMVN with sfu = "IVNOF"
    power_ivnof <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = 2, RP = "PBR", sfu="OF", delta=delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ivnof, sfu = "IVNOF_PBR_K2"))
  }
  # Loop over n from 12 to 30 in increments of 3 for both "LDMOF" and "IVNOF"
  for (n in seq(12, 30, by = 1)) {
    # Call Power_condMVN with sfu = "LDMOF"
    
    # Call Power_condMVN with sfu = "IVNOF"
#    power_ivnof <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "BSD", sfu="OF", delta=delta)$Pow)
#    results_df <- rbind(results_df, data.frame(n = n, power = power_ivnof, sfu = "IVNOF_BSD"))
    
    # Call Power_condMVN with sfu = "IVNOF"
 
    #   power_ivnof_RAR <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "RAR", sfu="OF", delta=delta)$Pow)
#    results_df <- rbind(results_df, data.frame(n = n, power = power_ivnof_RAR, sfu = "IVNOF_RAR"))
    
    # Call Power_condMVN with sfu = "LDMOF"
    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "PBR", sfu = "LDMOF", delta = delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_PBR_K3"))
    
    # Call Power_condMVN with sfu = "IVNOF"
    power_ivnof <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "PBR", sfu="OF", delta=delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ivnof, sfu = "IVNOF_PBR_K3"))

#    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "CHEN", sfu = "LDMOF", delta = delta)$Pow)
#    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_CHEN"))
    
    # Call Power_condMVN with sfu = "IVNOF"
#    power_ivnof <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "CHEN", sfu="OF", delta=delta)$Pow)
#    results_df <- rbind(results_df, data.frame(n = n, power = power_ivnof, sfu = "IVNOF_CHEN"))
    
    
#    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "RAR", sfu = "LDMOF", delta = delta)$Pow)
#    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_RAR"))
    
    # Call Power_condMVN with sfu = "IVNOF"
 #   power_ivnof <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "RAR", sfu="OF", delta=delta)$Pow)
#    results_df <- rbind(results_df, data.frame(n = n, power = power_ivnof, sfu = "IVNOF_RAR"))
    
    
#    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "BSD", sfu = "LDMOF", delta = delta)$Pow)
#    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_BSD"))
    
    # Call Power_condMVN with sfu = "IVNOF"
 #   power_ivnof <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "BSD", sfu="OF", delta=delta)$Pow)
#    results_df <- rbind(results_df, data.frame(n = n, power = power_ivnof, sfu = "IVNOF_BSD"))
  }
  
  for (n in seq(16, 30, by = 1)) {
    # Call Power_condMVN with sfu = "LDMOF"
    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = 4, RP = "PBR", sfu = "LDMOF", delta = delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_PBR_K4"))
    
    # Call Power_condMVN with sfu = "IVNOF"
    power_ivnof <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = 4, RP = "PBR", sfu="OF", delta=delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ivnof, sfu = "IVNOF_PBR_K4"))
  }
  # Save the results to an Excel file
  file_path <- paste0("data/Power_diff_sample_sizes_PBR_nsim_", n_sim, "_K_24", K, ".xlsx")
  
  write_xlsx(results_df, file_path)
  
  # Read back the Excel file
  read_data <- read_excel(file_path)
  
  # Plot the line plot for both sfu = "LDMOF" and "IVNOF"
  plot <- ggplot(read_data, aes(x = n, y = power, color = sfu, group = sfu)) +
    geom_line(size = 0.3) +        # Line connecting the points
    geom_point(size = 2) +       # Points at each value of n
    labs(x = "n", y = "Power", title = paste0("Power for different sample sizes for K=", K))  +
    theme_minimal() +
    scale_x_continuous(breaks = seq(8, 30, by = 1))   # Ensure breaks on x-axis correspond to the values of n
  
  file_name <- paste("data/plot_K", K, "_nSim", n_sim, ".pdf", sep = "")
  
  # Save the plot
  ggsave(filename = file_name, plot = plot)
}

#Power_diff_sample_sizes()

Power_diff_sample_sizes_plot <- function(data) {
  # Read back the Excel file
  read_data <- read_excel(data)
  
  # Plot the line plot for both sfu = "LDMOF" and "IVNOF"
  ggplot(read_data, aes(x = n, y = power, color = sfu, group = sfu)) +
    geom_line(size = 0.3) +        # Line connecting the points
    geom_point(size = 2) +       # Points at each value of n
    labs(x = "n", y = "Power", title = paste0("Power for different sample sizes")) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(12, 90, by = 3))   # Ensure breaks on x-axis correspond to the values of n
}


#Power_diff_sample_sizes_plot(data = "data/Power_diff_sample_sizes_nsim_1000_K_3.xlsx")
