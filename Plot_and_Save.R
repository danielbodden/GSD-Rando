source("RandomizationGSD.R")


# Calculates the power for different randomization procedures with the following boundaries:
# - Pocock / OF
# - LanDeMets alpha spending Pocock / OF
# - Pocock / Of with correction for actual information
# _ inverse normal combination function Pocock / OF boundaries
# and saves it as an excel file.
power_save_to_excel <- function(n, n_sim, K, sides = 1, alpha = 0.025, rb = 4, mti = 3, p = 2/3, futility = "no") {
  deltas <- seq(0, 2, by = 0.2)           # Grid for effect sizes

  # List of group sequential designs used
  gsd <- list("POC" = "Pocock", "LDMPOC" = sfLDPocock, "coRPOC" = "coRPOC", "OF" = "OF", "LDMOF" = sfLDOF, "corOF" = "corOF", "IVNOF" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")
  
  # Function to calculate power for a given sfu and sfu object
  calculate_power_for_sfu <- function(rp, gsd_object, es) {
    Power_list = Power_condMVN(n = n, n_sim = n_sim, K = K, RP = rp, sfu = gsd_object, delta = es)
    Power_mean = round(mean(Power_list$Pow), 4)
    return(list(Pow = Power_mean, Counter = Power_list$Counter))
    }
  calculate_power_for_inverse_normal <- function(rp, gsd_object, es) {
    Power_list = Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = rp, sfu=gsd_object, delta = es)
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
  filepath <- paste("data/ResultsPower_n", n, "_K", K, "_n_sim", n_sim, ".xlsx")  
  
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
T1E_save_to_excel <- function(n, n_sim, K, sides = 1, alpha = 0.025, rb = 4, mti = 3, p = 2/3, futility = "no") {
  RP_values <- c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN")
  
  # List of group sequential designs used
  gsd <- list("POC" = "Pocock", "LDMPOC" = sfLDPocock, "coRPOC" = "coRPOC", "OF" = "OF", "LDMOF" = sfLDOF, "corOF" = "corOF", "IVNOF" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")
  
  # Function to calculate T1E for a given randomization procedure (rp) and gsd (gsd_object)
  calculate_power_for_sfu <- function(rp, gsd_object, es) {
    Power_condMVN(n = n, n_sim = n_sim, K = K, RP = rp, sfu = gsd_object, delta = 0)
  }
  calculate_power_for_inverse_normal <- function(rp, gsd_object, es) {
    Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = rp, sfu=gsd_object, delta = 0)
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
  
  
  # Save the workbook to the specified file
  filepath <- paste("data/T1E n", n, " K", K, " n_sim", n_sim, ".xlsx")  
  saveWorkbook(wb, filepath, overwrite = TRUE)
  cat("Workbook saved to", filepath, "\n")
}


# Read in the saved excel file with T1E for the RPs and return a dataframe 
read_all_sheets_into_one_df <- function(file_path) {
  sheet_names <- excel_sheets(file_path)      # Get the sheet names
  combined_df <- data.frame()                   # Initialize an empty dataframe
  
  for (sheet in sheet_names) {                 # Loop through each sheet and read it into a dataframe
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
  p <- ggplot(data, aes(x = interaction_term, y = Power)) +
    geom_boxplot() +
    labs(x = "GSD and RP Combination", y = "Power") +
    ggtitle("T1E for GSD with n=24 and K=3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = levels(data$interaction_term)) 
 #   ylim(0.024, NA)  # Set lower limit of y-axis to 0.0225
  print(p)    # Display the plot
}

# Example usage
#file_path <- "data/T1E n 24  K 3  n_sim 1000 .xlsx"
#create_boxplot(file_path)


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


#sfu <- c("IVNOF")  # Multiple GSD values
#PlotPower("data/ResultsPower_n 24 _K 3 _n_sim 1000 .xlsx", RP_values=c("CR", "BSD", "CHEN", "PBR", "EBC", "RAR"), sfu = sfu)

#sfu <- c("Pocock")  # Multiple GSD values
#PlotPower("data/ResultsPower_n 24 _K 2 _reps 8000 .xlsx", RP_values=c("CR"), sfu = sfu)



#sfu <- c("IVNPOC")  # Multiple GSD values
#PlotPower("data/ResultsPower_n 24 _K 3 _n_sim 5000 .xlsx", RP_values=c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), sfu = sfu)



# IVNPOC has slighlty lower Power than Poc and LDMPOC for CR, but same power for PBR(4) as then no wrong information is used.
# I still have to describe coRPOC
# create Boxplots for RPs
# when creating th eboxplots, mention how many values are not shown (Removed 13 rows containing non-finite values (stat_boxplot()).)
# Problem with inverse  normal when a stage has only zero allocations -> in that case I remove the data from the patients completely.
# Other option would be to add it to another stage, but if it was in the last stage the test was alrady conducted.
# Recreate boxplots for n=5000 with new calculation

