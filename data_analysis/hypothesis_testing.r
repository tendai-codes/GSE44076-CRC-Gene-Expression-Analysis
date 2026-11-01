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
# Paired ttests

# Genes to analyse
genes <- c("TP53", "KRAS", "MYC", "MMP7", "EPCAM")

# Empty results table
results <- data.frame()

# Open image file FIRST
png(
  "results/figures/RQ1_paired_histograms.png",
  width = 12,
  height = 8,
  units = "in",
  res = 200
)

# Arrange histograms in 2 rows and 3 columns
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

  # Paired t-test
  test <- t.test(difference, mu = 0)

  # Add results to table
  results <- rbind(
    results,
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

# Close image file
dev.off()


# Display results
print(results)

# Save results
write.csv(
  results,
  "results/tables/RQ1_paired_ttests.csv",
  row.names = FALSE
)
