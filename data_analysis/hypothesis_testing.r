# Q1: statistical testing to assess differential expression between non-tumor vs tumor tissue in cancer patients (normal vs cancerous)
# samples from same patient - adjacent tissue

#selected genes

genes <- c("TP53", "KRAS", "MYC", "MMP7", "EPCAM")

# load dataset
expression <- read.csv('data/processed/GSE44076_normalised.csv')

#dataset summary
head(expression)
dim(expression)
summary(expression)


metadata <- read.csv('data/raw/metadata_GSE44076.csv')

#metadata summary
head(metadata)
dim(metadata)
summary(metadata)

#join expression + metadata

data <- merge(
    expression,
    metadata,
    by = 'sample_id'

)
# check merged dataset integrity
dim(data)
head(data)
table(data$refinebio_disease)

#create subset groups
healthy <- subset(
    data,
    refinebio_disease == "mucosa"
)

normal <- subset(
    data, 
    refinebio_disease == "normal"
)

tumour <- subset(
    data,
    refinebio_disease == "tumor"
)

#verify new tables
nrow(healthy)
nrow(normal)
nrow(tumour)

#RQ1 - paired t-test for normal vs tumour
names(metadata)
head(metadata)

paired_norm_tumor <- merge(
    normal,
    tumour,
    by = "refinebio_subject",
    suffixes = c("_normal", "_tumour")
)

# verify dataset + column matches
nrow(paired_norm_tumor)
anyDuplicated(paired_norm_tumor$refinebio_subject)
head(paired_norm_tumor)

all(paired_norm_tumor$refinebio_age_normal == paired_norm_tumor$refinebio_age_normal)
all(paired_norm_tumor$refinebio_sex_normal == paired_norm_tumor$refinebio_sex_normal)

--------------------------------------------------------------------------
# HYPOTHESIS TESTING

# Genes to analyse
genes <- c("TP53", "KRAS", "MYC", "MMP7", "EPCAM")

----------------------------------------------------------------------
# RQ1 - NORMAL VS TUMOUR (PAIRED T-TESTS)

rq1_results <- data.frame()

# Open histogram image
png(
  "results/figures/RQ1_paired_histograms.png",
  width = 12,
  height = 8,
  units = "in",
  res = 200
)

par(mfrow = c(2, 3))

for (gene in genes) {

  # Calculate tumour - normal expression
  difference <- paired_norm_tumor[[paste0(gene, "_tumour")]] -
                paired_norm_tumor[[paste0(gene, "_normal")]]

  # Histogram
  hist(
    difference,
    main = gene,
    xlab = "Tumour - Normal",
    ylab = "Frequency",
    col = "lightblue"
  )

  abline(v = 0, col = "red", lty = 2)

  # Paired t-test on differences
  test <- t.test(difference, mu = 0)

  # Add results to table
  rq1_results <- rbind(
    rq1_results,
    data.frame(
      Gene = gene,
      N = length(difference),
      Mean = mean(difference),
      Median = median(difference),
      SD = sd(difference),
      Percent_Increased = mean(difference > 0) * 100,
      T = as.numeric(test$statistic),
      DF = as.numeric(test$parameter),
      P_Value = test$p.value,
      CI_Lower = test$conf.int[1],
      CI_Upper = test$conf.int[2]
    )
  )
}

dev.off()

# Adjust p-values for 5 genes
rq1_results$Adjusted_P <- p.adjust(
  rq1_results$P_Value,
  method = "BH"
)

print(rq1_results)

write.csv(
  rq1_results,
  "results/tables/RQ1_paired_ttests.csv",
  row.names = FALSE
)

# RQ1 - Q-Q PLOTS FOR PAIRED DIFFERENCES

png(
  "results/figures/RQ1_qq_plots.png",
  width = 12,
  height = 8,
  units = "in",
  res = 200
)

par(mfrow = c(2, 3))

for (gene in genes) {

  difference <- paired_norm_tumor[[paste0(gene, "_tumour")]] -
                paired_norm_tumor[[paste0(gene, "_normal")]]

  qqnorm(
    difference,
    main = paste(gene, "Q-Q Plot")
  )

  qqline(difference, col = "red")
}

dev.off()

---------------------------------------------------------------
# RQ1 PAIRED DIFFERENCE PLOTS


png(
  "results/figures/RQ1_paired_difference_plots.png",
  width = 12,
  height = 8,
  units = "in",
  res = 200
)

par(mfrow = c(2, 3))

for (gene in genes) {

  difference <- paired_norm_tumor[[paste0(gene, "_tumour")]] -
                paired_norm_tumor[[paste0(gene, "_normal")]]

  # Sort patients by expression difference
  difference <- sort(difference)

  # Colour increases green and decreases red
  point_colours <- ifelse(difference > 0, "darkgreen", "red")

  plot(
    difference,
    main = gene,
    xlab = "Patients (ordered by difference)",
    ylab = "Tumour - Normal expression",
    pch = 19,
    col = point_colours
  )

  # Zero = no expression difference
  abline(h = 0, col = "black", lty = 2)

  # Average expression difference
  abline(h = mean(difference), col = "blue", lwd = 2)
}

dev.off()

print(rq1_results)

#export table with results
write.csv(
  rq1_results,
  "results/tables/RQ1_paired_ttests.csv",
  row.names = FALSE
)
--------------------------------------------------------------
# RQ2 - HEALTHY VS NORMAL (INDEPENDENT WELCH T-TESTS)

rq2_results <- data.frame()

png(
  "results/figures/RQ2_boxplots.png",
  width = 12,
  height = 8,
  units = "in",
  res = 200
)

par(mfrow = c(2, 3))

for (gene in genes) {

  x <- healthy[[gene]]
  y <- normal[[gene]]

  # Side-by-side boxplots
  boxplot(
    x,
    y,
    names = c("Healthy", "Normal"),
    main = gene,
    xlab = "Tissue group",
    ylab = "Normalised expression",
    col = c("lightblue", "lightgreen")
  )

  # Independent Welch's t-test
  # Difference = Normal - Healthy
  test <- t.test(
    y, x,
    var.equal = FALSE
  )

  # Store results
  rq2_results <- rbind(
    rq2_results,
    data.frame(
      Gene = gene,
      N_Healthy = length(x),
      N_Normal = length(y),
      Mean_Healthy = mean(x),
      Mean_Normal = mean(y),
      Mean_Difference = mean(y) - mean(x),
      T = as.numeric(test$statistic),
      DF = as.numeric(test$parameter),
      P_Value = test$p.value,
      CI_Lower = test$conf.int[1],
      CI_Upper = test$conf.int[2]
    )
  )
}

dev.off()

print(rq2_results)

#export results table
write.csv(
  rq2_results,
  "results/tables/RQ2_healthy_vs_normal.csv",
  row.names = FALSE
)
----------------------------------------------------------------
# RQ3 - HEALTHY VS TUMOUR (INDEPENDENT WELCH T-TESTS)


rq3_results <- data.frame()

png(
  "results/figures/RQ3_boxplots.png",
  width = 12,
  height = 8,
  units = "in",
  res = 200
)

par(mfrow = c(2, 3))

for (gene in genes) {

  x <- healthy[[gene]]
  y <- tumour[[gene]]

  # Side-by-side boxplots
  boxplot(
    x,
    y,
    names = c("Healthy", "Tumour"),
    main = gene,
    xlab = "Tissue group",
    ylab = "Normalised expression",
    col = c("lightblue", "salmon")
  )

  # Independent Welch's t-test
  # Difference = Tumour - Healthy
  test <- t.test(
    y, x,
    var.equal = FALSE
  )

  rq3_results <- rbind(
    rq3_results,
    data.frame(
      Gene = gene,
      N_Healthy = length(x),
      N_Tumour = length(y),
      Mean_Healthy = mean(x),
      Mean_Tumour = mean(y),
      Mean_Difference = mean(y) - mean(x),
      T = as.numeric(test$statistic),
      DF = as.numeric(test$parameter),
      P_Value = test$p.value,
      CI_Lower = test$conf.int[1],
      CI_Upper = test$conf.int[2]
    )
  )
}

dev.off()

print(rq3_results)

#export tables
write.csv(
  rq3_results,
  "results/tables/RQ3_healthy_vs_tumour.csv",
  row.names = FALSE
)

---------------------------------------------------------------
# CHECK RESULTS


nrow(rq1_results)
nrow(rq2_results)
nrow(rq3_results)


