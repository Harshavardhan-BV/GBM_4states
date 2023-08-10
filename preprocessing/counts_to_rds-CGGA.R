source('./common_functions.R')

# Read the counts
df = read.delim("../Data/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",row.names = 1)
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'CGGA', '325')

# Read the counts
df = read.delim("../Data/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",row.names = 1)
# Convert to TPM
df = FPKMToTPM(df)
# Save the counts matrix
exportCount(df, 'CGGA', '693')