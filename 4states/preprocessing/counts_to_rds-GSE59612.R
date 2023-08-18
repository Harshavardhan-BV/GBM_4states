source('./common_functions.R')

# Read the counts
df = read.delim("../Data/GSE59612_RAW/GSE59612_contrast-enhancing_glioma_FPKM.txt.gz")
# Check for non-numeric values and replace with 0
df[, -1] <- lapply(df[, -1], function(x) ifelse(is.na(as.numeric(x)), 0, as.numeric(x)))
# merge duplicate gene names
df = df %>%
  group_by(X) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE59612', 'contrast-enhancing')

# Read the counts
df = read.delim("../Data/GSE59612_RAW/GSE59612_non-enhancing_glioma_FPKM.txt.gz")
# merge duplicate gene names
df = df %>%
  group_by(X) %>%
  summarize_all(mean) %>% data.frame()
# Set the row names as gene
rownames(df) = df[,1]
df = df[,-1]
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'GSE59612', 'non-enhancing')