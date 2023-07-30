source('./common_functions.R')

# Read the counts
df = read.csv("../Data/GSE189577_RAW/GSM5704141_EarlyStage_Day1.csv.gz",row.names = 1)
# Convert to TPM
df = countToTpm(df)
# uppercase gene names
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE189577', 'GSM5704141', sc=T)

# Read the counts
df = read.csv("../Data/GSE189577_RAW/GSM5704142_EarlyStage_Day2.csv.gz",row.names = 1)
# Convert to TPM
df = countToTpm(df)
# uppercase gene names
rownames(df) = toupper(rownames(df))
# Save the counts matrix
exportCount(df, 'GSE189577', 'GSM5704142', sc=T)