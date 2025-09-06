# Function to diagnose sampling behavior
diagnose_sampling <- function(data, sample_frac) {
  
  # Get class information
  class_levels <- levels(data$outcome)
  class_counts <- table(data$outcome)
  
  cat("\nClass distribution in full data:\n")
  print(prop.table(class_counts))
  
  # Train single tree for inspection
  single_tree <- ranger(
    outcome ~ .,
    data = data,
    num.trees = 1,
    sample.fraction = sample_frac,
    replace = TRUE,
    probability = TRUE,
    keep.inbag = TRUE  # Keep track of which samples were used
  )
  
  # Check in-bag counts by class
  inbag <- single_tree$inbag.counts[[1]]
  
  cat("\nAverage times each class appears in bootstrap sample:\n")
  for(level in class_levels) {
    class_idx <- which(data$outcome == level)
    avg_count <- mean(inbag[class_idx])
    cat(level, ":", round(avg_count, 3), "\n")
  }
  
  # Predict on each class separately
  cat("\nPredictions by true class:\n")
  for(level in class_levels) {
    class_data <- data[data$outcome == level, ]
    if(nrow(class_data) > 0) {
      preds <- predict(single_tree, class_data)$predictions
      cat(level, "- Mean predicted prob for", class_levels[2], ":", 
          round(mean(preds[,2]), 3), "\n")
    }
  }
  
  return(single_tree)
}

# Run diagnosis
diag1 <- diagnose_sampling(small_train_data, c(frac_neg, frac_pos))
diag2 <- diagnose_sampling(small_train_data, rev(c(frac_neg, frac_pos)))
