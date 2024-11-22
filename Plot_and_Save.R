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
 #   deltas = c(1)
  # List of group sequential designs used
 # gsd <- list("POC" = "Pocock", "LDMPOC" = "LDMPocock", "coRPOC" = "coRPOC", "OF" = "OF", "LDMOF" = "LDMOF", "corOF" = "corOF", "INCT" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")
  # only important
  gsd <- list("LDMPOC" = "LDMPocock", "LDMOF" = "LDMOF", "INCT" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")

  
  # Function to calculate power for a given sfu and sfu object
  calculate_power_for_sfu <- function(rp, gsd_object, es) {
    Power_list = Power_condMVN(n = n, n_sim = n_sim, K = K, RP = rp, sfu = gsd_object, delta = es, futility = futility, futility_binding = futility_binding)
    Power_mean = round(mean(Power_list$Pow), 4)
    return(list(Pow = Power_mean, Counter = Power_list$Counter, stderr = Power_list$stderr))
    }
  calculate_power_for_inverse_normal <- function(rp, gsd_object, es) {
    Power_list = Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = rp, sfu=gsd_object, delta = es, futility = futility, futility_binding = futility_binding)
    Power_mean = round(mean(Power_list$Pow), 4)
    return(list(Pow = Power_mean, Counter = Power_list$Counter, stderr = Power_list$stderr))
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
        if (gsd_label == "INCT") {
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
          skipped_sequences = power_value$Counter,
          stderr=power_value$stderr
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
#  gsd <- list("POC" = "Pocock", "LDMPOC" = "LDMPocock", "coRPOC" = "coRPOC", "OF" = "OF", "LDMOF" = "LDMOF", "corOF" = "corOF", "INCT" = "inverse normal OF", "IVNPOC" = "inverse normal Pocock")
  # only some:
  gsd <- list("Pocock" = "Pocock", "OF" = "OF")
  
  
  # Function to calculate T1E for a given randomization procedure (rp) and gsd (gsd_object)
  calculate_power_for_sfu <- function(rp, gsd_object, es) {
    Power_condMVN(n = n, n_sim = n_sim, K = K, RP = rp, sfu = gsd_object, delta = 0, futility = futility, futility_binding = futility_binding)
  }
  calculate_power_for_inverse_normal <- function(rp, gsd_object, es) {
    Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = rp, sfu=gsd_object, delta = 0, futility = futility, futility_binding = futility_binding)
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
      if (gsd_label == "INCT") {
        power_value <- calculate_power_for_inverse_normal(RP, gsd_object = "OF")
      } else if (gsd_label == "IVNPOC") {
        power_value <- calculate_power_for_inverse_normal(RP, gsd_object = "Pocock")
      } else {
        gsd_object <- gsd[[gsd_label]]
        power_value <- calculate_power_for_sfu(RP, gsd_object)
      }
 #     power_value$Pow = round(mean(power_value$Pow),4)
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
#  filepath <- paste("data/Pow_RS n", n, " K", K, " n_sim", n_sim, "fut", futility_string, "futbind", futility_binding_string, "_delta1.4", ".xlsx")  
  filepath <- paste("data/T1E n", n, " K", K, " n_sim", n_sim, "fut", futility_string, "futbind", futility_binding_string, ".xlsx")  
  
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
#file_path <- "data/Pow_RS n 24  K 3  n_sim 16000 fut FALSE futbind FALSE .xlsx"
#create_boxplot(file_path)

create_violin_plot <- function(file_path) {
  data <- read_all_sheets_into_one_df(file_path)        # Read data from the Excel file
  data <- subset(data, gsd == "Pocock" | gsd == "OF")      # Only for Pocock and OF
  
  # Ensure 'gsd' and 'RP' are treated as factors and ordered
  data$gsd <- factor(data$gsd, levels = unique(data$gsd))
  data$RP <- factor(data$RP, levels = unique(data$RP))
  
  # Rename 'LDMPOC' and 'LDMOF' in the 'gsd' column
  data$gsd <- factor(data$gsd, levels = c("Pocock", "OF"), 
                     labels = c("Pocock", "O'Brien-Fleming"))

  # Mapping of randomization procedures to full names
  RP_full_names <- c(
    "CR" = "Complete\nRandomization",
    "PBR" = "Permuted Block\nRandomization (4)",
    "BSD" = "Big Stick\nDesign (3)",
    "RAR" = "Random\nAllocation Rule",
    "EBC" = "Efron's Biased\nCoin (2/3)",
    "CHEN" = "Chen's design\n(2/3, 3)"
  )
  
  
  # Replace the abbreviations in the RP column with the full names
  data$RP <- recode(data$RP, !!!RP_full_names)

  
  # Color mapping for the violin plot
  color_mapping <- c(
    "Complete\nRandomization" = "darkgray",
    "Permuted Block\nRandomization (4)" = "#33a02c",
    "Big Stick\nDesign (3)" = "#1f78b4",
    "Random\nAllocation Rule" = "#6a3d9a",
    "Efron's Biased\nCoin (2/3)" = "#ff7f00",
    "Chen's design\n(2/3, 3)" = "#e31a1c"
  )
  
  # Round the Power values
  data$Power <- round(data$Power, 4)
  
  # Create the violin plot with facet for each GSD
  p <- ggplot(data, aes(x = RP, y = Power, fill = RP)) +
    geom_violin(width = 1, alpha = 0.5) +
    stat_summary(fun = "mean",
                 geom = "crossbar", 
                 width = 0.2,
                 linewidth = 0.4,
                 colour = "red") +
    geom_hline(yintercept = 0.025,        # Add horizontal line at y = 0.025
               linetype = "dashed",       # Make the line dotted
               color = "darkgreen",       # Set the color to green
               size = 1) +                # Adjust the thickness of the line
    scale_fill_manual(values = color_mapping) +  # Apply the custom color mapping
    labs(x = "Randomization Procedure", y = "Type I error rate conditioned on the randomization sequence") +
    theme_minimal() + 
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
          axis.text.y = element_text(size = 14),                         # Increase y-axis text size
          axis.title.x = element_text(size = 14),                        # Increase x-axis title size
          axis.title.y = element_text(size = 14),                        # Increase y-axis title size
          strip.text = element_text(size = 14),                          # Increase facet title size
          plot.title = element_text(size = 20, face = "bold")) +         # Set plot title size
    facet_wrap(~ gsd, ncol = 1, scales = "free_y")    # Separate plots by 'gsd', placed vertically
#  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))   # Increase number of y-axis values
 #   scale_y_continuous(breaks = scales::pretty_breaks(n = 15), limits = c(0, 0.029))   # Increase number of y-axis values
    
    
  # Print the plot
  print(p)
}




# Example usage
#file_path <- "data/T1E n 24  K 3  n_sim 1e+05 fut FALSE futbind FALSE .xlsx"
#create_violin_plot(file_path)



PlotPower = function(file_path, RP_values = c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), sfu) {
  # Create a mapping of abbreviations to full names
  RP_full_names <- c(
    "CR" = "Complete Randomization",
    "PBR" = "Permuted Block Randomization (4)",
    "BSD" = "Big Stick Design (3)",
    "RAR" = "Random Allocation Rule",
    "EBC" = "Efron's Biased Coin (2/3)",
    "CHEN" = "Chen's design (2/3, 3)"
  )
  
  # Desired order for RP (first PBR, then CHEN, etc.)
  RP_order <- c("Permuted Block Randomization (4)", "Chen's design (2/3, 3)", "Big Stick Design (3)", 
                "Efron's Biased Coin (2/3)", "Random Allocation Rule", "Complete Randomization")
  
  sheets <- excel_sheets(file_path)         # List all sheets
  all_data <- list()                        # Initialize a list to store combined data
  
  for (rp in RP_values) {                   # Loop through each RP and GSD combination
    for (gsd in sfu) {
      data_list <- list()                   # Initialize a list to store data for current RP and GSD combination
      
      for (sheet in sheets) {               # Loop through each sheet
        data <- read_excel(file_path, sheet = sheet)
        filtered_data <- filter(data, RP == rp, GSD == gsd)  # Filter for the specific RP and GSD values
        
        if (nrow(filtered_data) > 0) {      # Check if there is any data left after filtering
          # Replace abbreviations with full names
          filtered_data$RP <- recode(filtered_data$RP, !!!RP_full_names)  # Replace RP with full name
          filtered_data$GSD <- gsd  # Add GSD column
          data_list[[sheet]] <- filtered_data[c("delta", "Power", "RP", "GSD")]
        }
      }
      
      if (length(data_list) > 0) {          # Combine data for the current RP and GSD
        all_data[[paste(rp, gsd)]] <- do.call(rbind, data_list)
      }
    }
  }
  
  if (length(all_data) > 0) {               # Combine all data into a single dataframe for plotting
    plot_data <- do.call(rbind, all_data)
    print(plot_data)
    
    # Convert RP to factor and set the levels in the desired order
    plot_data$RP <- factor(plot_data$RP, levels = RP_order)
    
    # Create a column to define linetypes based on GSD
    plot_data$linetype <- ifelse(plot_data$GSD == "LDM", "solid", "dashed")
    
    # Set colors for RP values
    color_mapping <- c(
      "Complete Randomization" = "darkgray",
      "Permuted Block Randomization (4)" = "#33a02c",
      "Big Stick Design (3)" = "#1f78b4",
      "Random Allocation Rule" = "#6a3d9a",
      "Efron's Biased Coin (2/3)" = "#ff7f00",
      "Chen's design (2/3, 3)" = "#e31a1c"
    )
    
    # Set shapes for RP values
    shape_mapping <- c(
      "Complete Randomization" = 25, 
      "Permuted Block Randomization (4)" = 23,
      "Big Stick Design (3)" = 21,
      "Random Allocation Rule" = 24,
      "Efron's Biased Coin (2/3)" = 20,
      "Chen's design (2/3, 3)" = 22
    )  
    
    # Generate the plot
    ggplot(plot_data, aes(x = delta, y = Power, color = RP, group = interaction(RP, GSD), linetype = linetype)) +
      geom_point(aes(shape = RP, fill = ifelse(GSD == "LDM", RP, NA)), size = 1.5) +  # Use shapes and fill only for LDM
      geom_line(size = 0.4) +
      scale_color_manual(values = color_mapping) +  # Set the colors for RPs
      scale_shape_manual(values = shape_mapping) +  # Map RP to different shapes
      scale_fill_manual(values = color_mapping, na.translate = FALSE) +  # Fill color for LDM points only
      scale_linetype_manual(values = c("solid", "dashed")) +  # Linetypes for GSD
      labs(x = "Effect Size", y = "Power", color = "Randomization Procedure", shape = "Randomization Procedure") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme_minimal() +
      theme(    legend.position = c(0.8, 0.2),         # Set the legend position inside the plot
                legend.background = element_rect(fill = "white", color = "black"),  # Add a background box for readability
                legend.title = element_text(size = 12),  # Adjust legend title size
                legend.text = element_text(size = 10),  # Adjust legend text size
                axis.title.x = element_text(size = 12),  # Increase x-axis label size
                axis.title.y = element_text(size = 12),  # Increase y-axis label size
                axis.text.x = element_text(size = 10),   # Increase x-axis tick size
                axis.text.y = element_text(size = 10)    # Increase y-axis tick size
      ) +
      guides(linetype = "none")  # Remove the linetype from the legend
  } else {
    print("No data met the filter criteria for the specified RP and GSD values.")
  }
}


#sfu <- c("INCT")  # Multiple GSD values
#PlotPower("data/ResultsPower n 120 K 3 n_sim 1000 fut FALSE futbind FALSE .xlsx", RP_values=c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), sfu = sfu)

#sfu <- c("Pocock")  # Multiple GSD values
#PlotPower("data/ResultsPower_n 24 _K 2 _reps 8000 .xlsx", RP_values=c("CR"), sfu = sfu)



#sfu <- c("LDMOF")  # Multiple GSD values
#PlotPower("data/ResultsPower n 24 K 3 n_sim 1000 fut FALSE futbind FALSE .xlsx", RP_values=c("CR", "PBR", "BSD", "RAR", "EBC", "CHEN"), sfu = sfu)



create_violin_plot_binding <- function() {
  data1 <- read_all_sheets_into_one_df("data/T1E n 24  K 3  n_sim 1000 fut TRUE futbind TRUE .xlsx")        # Read data from the Excel file
  data1$fut <- "Binding"
  data2 <- read_all_sheets_into_one_df("data/T1E n 24  K 3  n_sim 1000 fut TRUE futbind FALSE .xlsx")        # Read data from the Excel file
  data2$fut <- "Non binding"
  data3 <- read_all_sheets_into_one_df("data/T1E n 24  K 3  n_sim 8000 .xlsx")        # Read data from the Excel file
  data3$fut <-"No futility"
  data <- bind_rows(data1,data2,data3)
  data <- subset(data, gsd == "OF")      # Only for Pocock and OF
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
  
  data$fut <- factor(data$fut, levels=c("No futility", "Binding", "Non binding"))
  # Create the violin plot with the specified title
  p <- ggplot(data, aes(x = interaction_term, y = Power, fill = fut)) +
    geom_boxplot(width=1, alpha=.5) +
    geom_hline(yintercept = 0.025,        # Add horizontal line at y = 0.025
               linetype = "dashed",       # Make the line dotted
               color = "darkgreen",           # Set the color to green
               size = 1) +              # Adjust the thickness of the line
    labs(x = "GSD and RP Combination", y = "Type I Error") +
  #  ggtitle("Power for GSD with n=24 and K=3") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_minimal()  +
    theme(legend.position = "none") +
    theme(legend.position = c(0.3, 0.2),         # Set the legend position inside the plot (adjust coordinates as needed)
          legend.background = element_rect(fill = "white", color = "black"), # Add a background box for readability
          legend.title = element_text(size = 10),                           # Adjust legend title size
          legend.text = element_text(size = 8)) +                           # Adjust legend text size    ylim(0.023, 0.026)  
       coord_cartesian(ylim=c(0.023, 0.026)) 

  
  print(p)    # Display the plot
}


#create_violin_plot_binding()

# Different sample sizes

# Sample size on x axis
# Power on y axis
# for one effect size
# Function to save Power_condMVN results for n=12 to n=30 for both sfu = "LDMOF" and "INCT" and create a line plot


Power_diff_sample_sizes <- function(n_sim = 1, K = 3, delta = 1) {
  
  # Initialize an empty dataframe to store results for both "LDMOF" and "INCT"
  results_df <- data.frame(n = integer(), power = numeric(), sfu = character())

 # for (n in seq(12, 60, by = 6)) {

    # Call Power_condMVN with sfu = "INCT"
 #       power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "CR", sfu="OF", delta=delta)$Pow)
#        results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_CR"))
        
        # Call Power_condMVN with sfu = "INCT"
#        power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "PBR", sfu="OF", delta=delta)$Pow)
#        results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_PBR"))
        
        
        # Call Power_condMVN with sfu = "INCT"
#        power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "BSD", sfu="OF", delta=delta)$Pow)
#        results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_BSD"))
        
        
        # Call Power_condMVN with sfu = "INCT"
#        power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "RAR", sfu="OF", delta=delta)$Pow)
 #       results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_RAR"))
        
        
        # Call Power_condMVN with sfu = "INCT"
#        power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "EBC", sfu="OF", delta=delta)$Pow)
#        results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_EBC"))
        
        
        # Call Power_condMVN with sfu = "INCT"
 #       power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "CHEN", sfu="OF", delta=delta)$Pow)
  #      results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_CHEN"))
    
    # Call Power_condMVN with sfu = "INCT"


   #     power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "CR", sfu = "LDMOF", delta = delta)$Pow)
  #      results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_CR"))
    
        
   #     power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "PBR", sfu = "LDMOF", delta = delta)$Pow)
  #      results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_PBR"))
        
        
   #     power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "BSD", sfu = "LDMOF", delta = delta)$Pow)
  #      results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_BSD"))
        
        
  #      power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "EBC", sfu = "LDMOF", delta = delta)$Pow)
  #      results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_EBC"))
        
  #      power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "CHEN", sfu = "LDMOF", delta = delta)$Pow)
  #      results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_CHEN"))
        
#}
  for (n in seq(8, 30, by = 1)) {
    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = 2, RP = "PBR", sfu = "LDMOF", delta = delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_PBR_K2"))
    
    # Call Power_condMVN with sfu = "INCT"
    power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = 2, RP = "PBR", sfu="OF", delta=delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_PBR_K2"))
  }
  # Loop over n from 12 to 30 in increments of 3 for both "LDMOF" and "INCT"
  for (n in seq(12, 30, by = 1)) {

    # Call Power_condMVN with sfu = "LDMOF"
    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = K, RP = "PBR", sfu = "LDMOF", delta = delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_PBR_K3"))
    
    # Call Power_condMVN with sfu = "INCT"
    power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = K, RP = "PBR", sfu="OF", delta=delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_PBR_K3"))


  }
  
  for (n in seq(16, 30, by = 1)) {
    # Call Power_condMVN with sfu = "LDMOF"
    power_ldmof <- mean(Power_condMVN(n = n, n_sim = n_sim, K = 4, RP = "PBR", sfu = "LDMOF", delta = delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_ldmof, sfu = "LDMOF_PBR_K4"))
    
    # Call Power_condMVN with sfu = "INCT"
    power_INCT <- mean(Power_inverse_normal(n = n, n_sim = n_sim, K = 4, RP = "PBR", sfu="OF", delta=delta)$Pow)
    results_df <- rbind(results_df, data.frame(n = n, power = power_INCT, sfu = "INCT_PBR_K4"))
  }
  # Save the results to an Excel file
  file_path <- paste0("data/Power_diff_sample_sizes_PBR_nsim_", n_sim, "_K_24", K, "add_patients_div.xlsx")
 # file_path <- paste0("data/Power_diff_sample_sizes_nsim_", n_sim, "_K", K, ".xlsx")
  
  write_xlsx(results_df, file_path)
  
  # Read back the Excel file
  read_data <- read_excel(file_path)
  
  # Plot the line plot for both sfu = "LDMOF" and "INCT"
  plot <- ggplot(read_data, aes(x = n, y = power, color = sfu, group = sfu)) +
    geom_line(size = 0.3) +        # Line connecting the points
    geom_point(size = 2) +       # Points at each value of n
    labs(x = "n", y = "Power", title = paste0("Power for different sample sizes for K=", K, "patients_divided"))  +
    theme_minimal() +
    scale_x_continuous(breaks = seq(8, 30, by = 1))   # Ensure breaks on x-axis correspond to the values of n
  
  file_name <- paste("data/plot_K", K, "_nSim", n_sim, ".pdf", sep = "")
  
  # Save the plot
  ggsave(filename = file_name, plot = plot)
}


#Power_diff_sample_sizes()
Power_diff_sample_sizes_plot <- function(data) {
  # Read the Excel file
  read_data <- read_excel(data)
  print(read_data)
  
  # Filter data to only include rows where n <= 60
  read_data <- read_data[read_data$n <= 60, ]
  
  # Define new names for each category in sfu
  new_names <- c("IVNOF_BSD" = "INCT OF with BSD",
                 "IVNOF_CHEN" = "INCT OF with CHEN",
                 "IVNOF_PBR" = "INCT OF with PBR",
                 "IVNOF_RAR" = "INCT OF with RAR",
                 "LDMOF_BSD" = "LDM OF with BSD",
                 "LDMOF_CHEN" = "LDM OF with CHEN",
                 "LDMOF_PBR" = "LDM OF with PBR",
                 "LDMOF_RAR" = "LDM OF with RAR",
                 "IVNOF_PBR_K2" = "INCT OF with PBR for K=2",
                 "IVNOF_PBR_K3" = "INCT OF with PBR for K=3",
                 "IVNOF_PBR_K4" = "INCT OF with PBR for K=4",
                 "LDMOF_PBR_K2" = "LDM OF with PBR for K=2",
                 "LDMOF_PBR_K3" = "LDM OF with PBR for K=3",
                 "LDMOF_PBR_K4" = "LDM OF with PBR for K=4")
  
  # Create a column to distinguish INCT and LDM for line type and fill
  read_data$sfu_type <- ifelse(grepl("LDMOF", read_data$sfu), "Lan DeMets", "Inverse Normal Combination Test")
  
  # Extract the randomization method from 'sfu'
  read_data$randomization <- sub("^.*_(.*)$", "\\1", read_data$sfu)
  
  # Define full names for randomization procedures
  randomization_full_names <- c(
#    "CR" = "Complete Randomization",
#    "PBR" = "Permuted Block Randomization (4)",
#    "BSD" = "Big Stick Design (3)",
#    "RAR" = "Random Allocation Rule",
#    "EBC" = "Efron's Biased Coin (2/3)",
#    "CHEN" = "Chen's design (2/3, 3)"
    "K2" = "2",
    "K3" = "3",
    "K4" = "4"
  )
  
  # Replace the abbreviations in the randomization column with full names
  read_data$randomization <- recode(read_data$randomization, !!!randomization_full_names)
  
 # read_data$randomization <- factor(read_data$randomization, 
#                                    levels = c( "Permuted Block Randomization (4)", 
#                                                "Chen's design (2/3, 3)",
#                                               "Big Stick Design (3)", 
#                                               "Efron's Biased Coin (2/3)",
#                                               "Random Allocation Rule", 
#                                                "Complete Randomization"))
  
  # Set shapes for randomization methods
#  shape_mapping <- c(
#    "Complete Randomization" = 25, 
#    "Permuted Block Randomization (4)" = 23,
#    "Big Stick Design (3)" = 21,
#    "Random Allocation Rule" = 24,
#    "Efron's Biased Coin (2/3)" = 20,
#    "Chen's design (2/3, 3)" = 22
#  )
  
  # Improved, more visible colors for the lines and symbols
  color_mapping <- c(
    "Complete Randomization" = "darkgray",
    "Permuted Block Randomization (4)" = "#33a02c",
    "Big Stick Design (3)" = "#1f78b4",
    "Random Allocation Rule" = "#6a3d9a",
    "Efron's Biased Coin (2/3)" = "#ff7f00",
    "Chen's design (2/3, 3)" = "#e31a1c"
  )
  
  
  
  # Plot with shape, color, and linetype customizations
  ggplot(read_data, aes(x = n, y = power, group = sfu, shape = randomization, linetype = sfu_type)) +
    geom_line(aes(color = randomization), size = 0.5) +  # Line color depends on randomization method
    geom_point(aes(color = randomization, fill = randomization), size = 2) +  # Fill based on randomization (for LDM)
#    scale_shape_manual(values = shape_mapping) +  # Map randomization procedure to shapes
    scale_linetype_manual(values = c("Lan DeMets" = "solid", "Inverse Normal Combination Test" = "dashed")) +  # Linetypes for LDM and INCT
  #  scale_fill_manual(values = color_mapping) +  # Apply the same colors for fill
#    scale_color_manual(values = color_mapping) +  # Use the same colors for line and point color
    labs(x = "Maximum sample size", y = "Power", shape = "Number of Stages", linetype = "Group Sequential Design", color = "Number of Stages") +
    theme_minimal() +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(legend.position = c(0.85, 0.25),  # Adjust the legend position to be more central
          legend.background = element_rect(fill = "white", color = "black"),  # Add a background box for readability
          legend.title = element_text(size = 18),  # Adjust legend title size
          legend.text = element_text(size = 16),  # Adjust legend text size for better visibility
          legend.key.size = unit(0.9, "cm"),
          axis.title.x = element_text(size = 20),  # Increase x-axis label size
          axis.title.y = element_text(size = 20),  # Increase y-axis label size
          axis.text.x = element_text(size = 14),   # Increase x-axis tick size
          axis.text.y = element_text(size = 14)    # Increase y-axis tick size
          ) +  # Adjust size of legend keys for better readability
    scale_x_continuous(breaks = seq(8, 30, by = 1)) +  # Ensure breaks on x-axis correspond to the values of n
    guides(fill = "none")  # Remove the 'fill' legend if needed
  
}


#Power_diff_sample_sizes_plot(data = "data/Power_diff_sample_sizes_nsim_1000_K3.xlsx")


#Power_diff_sample_sizes_plot(data = "data/Power_diff_sample_sizes_PBR_nsim_1_K_243patients_divided.xlsx")
