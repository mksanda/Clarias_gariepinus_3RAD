# Load the packages into R
library(adegenet)
library(hierfstat)
library(pheatmap)

# load genepop data from Stacks populations run 
genind_obj <- read.genepop("catfish.snps.gen", ncode = 2) 

# Convert to a data frame
df <- genind2df(genind_obj)

write.csv(df, file = "fst_input.csv", row.names = FALSE)
df_fst <- read.csv("fst_input.csv", header = TRUE)

# Calculate F-statistics
fst_results <- wc(df_fst)
print(fst_results)

# Calculate pairwise Fst values
fst_matrix <- genet.dist(df_fst, method = "WC84")
print(fst_matrix)

# Plot heatmap
pheatmap(as.matrix(fst_matrix), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", display_numbers = TRUE)

#Fst confidence interval
perm_test <- boot.ppfst(genind_obj, nboot = 1000)
perm_test

