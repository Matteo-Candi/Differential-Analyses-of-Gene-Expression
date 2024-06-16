# 0 - Load packages --------------

# Install Bioconductor packages
BiocManager::install("BiocGenerics")  # Install BiocGenerics package
BiocManager::install("DESeq2")        # Install DESeq2 package

# Install GGally package from CRAN
install.packages("GGally")  # Install GGally package

# List of required packages
required_packages <- c("BiocGenerics", "DESeq2", "psych", "NetworkToolbox", "ggplot2",
                       "GGally", "sna", "network", "TCGAbiolinks", "SummarizedExperiment", "DT", "latex2exp", "gridExtra")

# Load required packages
lapply(required_packages, library, character.only = TRUE)  # Load all the required packages

# Set seed for reproducibility
seed <- 123
set.seed(seed)  # Set the seed to ensure reproducibility

# Define a vector of colors
colors <- c("#2EBFA5", "#F26419", "#FDCA40")  # Custom colors for plotting



# 1 - Load data and Preprocessing -----------------------------------------------------------

# Project ID
proj <- "TCGA-CHOL"
dir.create(file.path(proj))  # Create a directory for the project

# Look for all data linked to a "Primary Tumor" sample and store them in a dataframe
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")
GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")  # Download data
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata")  # Prepare data
rna.expr.data.C <- assay(rna.data.C)  # Extract assay data
genes.info.c <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))  # Get gene information

# Apply the same procedure to the "Normal tissue"
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")
GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")  # Download data
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata")  # Prepare data
rna.expr.data.N <- assay(rna.data.N)  # Extract assay data
genes.info.n <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))  # Get gene information

# Check if the gene information matches
all(na.omit(genes.info.n) == na.omit(genes.info.c))

# Query clinical data and save it to a file
clinical.query <- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
write.csv(clinical.query, file = file.path(proj, paste(proj, "_clinical_data.txt", sep = "")), row.names = FALSE, quote = FALSE)

# Display the table of cancer stages
table(clinical.query$ajcc_pathologic_stage)

# Investigate the distribution of age and cancer stage
ggplot(clinical.query[!is.na(clinical.query$age_at_index), ], aes(x = ajcc_pathologic_stage, y = age_at_index)) +
  geom_boxplot(fill = colors[3]) +
  labs(title = "Relationship between Age and Cancer Stage",
       x = "",
       y = "Age") +
  coord_cartesian(ylim = c(20, 90)) +
  theme_minimal() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  theme(plot.title = element_text(size = 12))

# View RNA expression data for normal tissue
View(rna.expr.data.N)

# Check the dimensions and uniqueness of the data
dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1, 12)))  # No duplicates for normal tissue

dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1, 12)))  # No duplicates for cancer tissue

# Total number of unique patients
length(unique(clinical.query$submitter_id))

# IDs of the "primal tumor" class patients
patients.C <- substr(colnames(rna.expr.data.C), 1, 12)

expr.C <- as.data.frame(rna.expr.data.C)
expr.N <- as.data.frame(rna.expr.data.N)

# Change the patients IDs to a shorter version
colnames(expr.C) <- substr(colnames(expr.C), 1, 12)
colnames(expr.N) <- substr(colnames(expr.N), 1, 12)

# For how many of them we have both a tumor and normal sample?
intersect(colnames(expr.N), colnames(expr.C))
length(intersect(colnames(expr.N), colnames(expr.C)))  # 8

# We have 35 samples of patients with cancer and 9 samples of patients without.
# We have just 8 samples of tissue from the same patients that lead to a group of 36 patients.

# Since we want to do a comparison between the sane and ill patients we need to get rid of the ones for which we don't have the "primal tumor" sample
idx <- match(setdiff(colnames(expr.N), colnames(expr.C)), colnames(expr.N))  # idx to remove
expr.N <- expr.N[, -c(idx)]

# And select those same patients from the "primal tumor" data
expr.C <- expr.C[, colnames(expr.N)]

length(intersect(colnames(expr.N), colnames(expr.C)))

# Check if the values stored are ok
# Check the actual counts
typeof(expr.C[1, 1])  # ok
any(is.na(expr.C))  # ok
any(is.nan(as.matrix(expr.C)))  # ok

typeof(expr.N[1, 1])  # ok
any(is.na(expr.N))  # ok
any(is.nan(as.matrix(expr.N)))  # ok

# Explanation for the normalization:
# Detailed explanation: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# Youtube explanation: https://www.youtube.com/watch?v=UFB993xufUU

# Merge the data and transform it in a dataframe
full.data <- cbind(expr.N, expr.C)
full.data <- data.frame(full.data)

# Create a new dataframe in which we specify which kind of tissue the information
# comes from "cancer" or "normal"
n_patients <- ncol(full.data)

metad <- rep("cancer", n_patients)
metad[1:(n_patients / 2)] <- "normal"
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[, 1] <- as.factor(metad[, 1])
full.data <- cbind(rownames(full.data), full.data)

# Transform all the data in a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = full.data, 
                              colData = metad, 
                              design = ~condition,
                              tidy = TRUE)

# Filtering: more than 10 counts on 90% of patients
limit.90 <- floor(n_patients * .9)  # 15

# Index of the rows that are at least 10 for 90% of the patients
keep <- rowSums(counts(dds) >= 10) >= limit.90
dds <- dds[keep, ]

# Number of genes we end up having
dim(counts(dds))

# Now let's normalize
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
sum(rowSums(normalized_counts == 0) == n_patients)  # no null rows

# Keep only the genes without zeros
normalized_counts <- normalized_counts[-which(rowSums(normalized_counts == 0) >= 1), ]

# Divide the cancerous and normal data
filtr.expr.n <- as.data.frame(normalized_counts[, 1:(n_patients / 2)])
filtr.expr.c <- as.data.frame(normalized_counts[, (n_patients / 2 + 1):n_patients])

# Cancerous sample names were added a ".1" in full.data because  
# they had the same names as the normal samples
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1, 12)




# 2 - Differentially Expressed Genes (DEGs) -------------------------------

# To Identify the Differently Expressed Genes we first need to compute the 
# FC (fold change)
fc <- log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n)) 
names(fc) <- rownames(filtr.expr.c)

# And the p-value using the t test
pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(filtr.expr.c[i,], filtr.expr.n[i,]))$p.value)

# And then we can adjust the values
pval.fc.fdr <- p.adjust(pval.fc, method = "fdr")

# Then we can collect both the fc and pvalues inside a dataframe
expr.table <- data.frame(cbind(fc, pval.fc.fdr))
expr.table[, 1] <- round(expr.table[, 1], 2)

# Now given the FCs and p-values we would like to identify which are the genes 
# that are different for both the measures

# So we can set a threshold for the FC and one for the p-value so that 
# we obtain hundreds of genes in the end
threshold.p <- 0.05
threshold.fc <- 1.2

# And then we select only the genes that have the absolute value of the fold change 
# over the threshold.fc and the p-value under threshold.p at the same time
# In this way we are sure that the genes are distributed differently depending 
# on the classes 
deg.genes <- rownames(expr.table[abs(expr.table$fc) >= threshold.fc & expr.table$pval.fc.fdr <= threshold.p,]) 
length(deg.genes)  # We have 3250 different genes

# Let's select just the DEG genes
filtr.expr.c.deg <- filtr.expr.c[deg.genes, ]
filtr.expr.n.deg <- filtr.expr.n[deg.genes, ]

# This difference can be visualized using a volcano plot
expr.table$diffexpressed <- "NO SIGNIFICANT"
# The genes with FC >= threshold.fc
expr.table$diffexpressed[expr.table$fc >= threshold.fc & expr.table$pval.fc.fdr <= threshold.p] <- "UP"
# The genes with FC <= -threshold.fc
expr.table$diffexpressed[expr.table$fc <= -threshold.fc & expr.table$pval.fc.fdr <= threshold.p] <- "DOWN"
expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)

line_col <- 'darkgrey'

ggplot(data = expr.table, aes(x = fc, y = -log10(pval.fc.fdr), col = diffexpressed)) +
  geom_point() +
  xlab(expression(log[2]~"(Fold Change)")) +
  ylab(expression(-log[10]~"(Adjusted p-value)")) +
  geom_hline(yintercept = -log10(threshold.p), col = line_col) +
  geom_vline(xintercept = c(-threshold.fc, threshold.fc), col = line_col) +
  labs(title = "", color = "") +
  coord_cartesian(ylim = c(-.5, 6.5), xlim = c(-8.5, 7)) +
  scale_color_manual(values = c("UP" = colors[1], "DOWN" = colors[2], "NO SIGNIFICANT" = "lightgrey")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))



# 3 - Co-expression networks ----------------------------------------------
## 3.1 - Computation  ----

# Now we would like to build the co-expression network, to do so we first need
# to compute the adjacency matrix

# We will define the matrix using the correlation between the genes so
# the edge between two genes will be present if the absolute value of the 
# correlation is over a threshold threshold.rho
threshold.rho <- 0.8 

# Before computing correlations, we log transform the data using log2(x + 1)
filtr.expr.c.deg <- log2(filtr.expr.c.deg + 1)
filtr.expr.n.deg <- log2(filtr.expr.n.deg + 1) 

# Create the CANCER NETWORK

# Correlation matrix for the cancer network - we will use Spearman correlation now
cor.mat.c <- corr.test(t(filtr.expr.c.deg), use = "pairwise", 
                       method = "spearman", adjust = "fdr", ci = FALSE)

# Matrix with the Spearman correlation coefficient
rho.c <- cor.mat.c$r
# Set the diagonal to zero (no self loops)
diag(rho.c) <- 0
# Binary matrix, keep the edge for big correlation
adj.mat.c <- abs(rho.c) > threshold.rho

# Same process for the NORMAL NETWORK

# Correlation for the normal network using Spearman
cor.mat.n <- corr.test(t(filtr.expr.n.deg), use = "pairwise", 
                       method = "spearman", adjust = "fdr", ci = FALSE)

rho.n <- cor.mat.n$r
diag(rho.n) <- 0

adj.mat.n <- abs(rho.n) > threshold.rho


## 3.2 - Analysis  ----

# Create the two networks given the adjacency matrices
net.n <- as.network(adj.mat.n, directed = FALSE)
net.c <- as.network(adj.mat.c, directed = FALSE)

# Check the number of nodes and edges

# Normal Network
network.size(net.n)  # Number of nodes: 3250
network.edgecount(net.n)  # Number of edges: 278221

# Cancer Network
network.size(net.c)  # Number of nodes: 3250
network.edgecount(net.c)  # Number of edges: 139842

# The number of edges in the normal network is higher, indicating that genes are more correlated in normal tissues.

# We look for a scale-free network in gene expression of cancer tissues because identifying key hub genes can reveal critical regulatory elements and potential therapeutic targets that drive cancer progression.
# To determine if a network is scale-free, we need to examine its degree distribution. In a scale-free network, a few nodes (hubs) have many connections, while most nodes have few. For a co-expression network, this means a small number of genes have numerous co-expression relationships, while most have only a few. If the degree distribution follows a power-law, the network is scale-free.

# First, for the NORMAL network:
# Extract the degree of each node
d.n <- sna::degree(net.n, gmode = 'graph')

# Remove the nodes with 0 degree
d.n <- d.n[d.n > 0]
names(d.n) <- network.vertex.names(net.n)

# Determine the 95% quantile to identify hubs
quantile_hubs <- 0.95
x_n <- quantile(d.n, quantile_hubs)

data.n <- data.frame(value = unlist(d.n))
row.names(head(sort(data.n, decreasing = TRUE), 10))  # Display the top 10 nodes with highest degree

# Plot the degree distribution
plot_1 <- ggplot(data.n, aes(x = value, fill = value > x_n)) +
  geom_histogram(binwidth = 15, color = "snow") +
  ggtitle("Normal Tissues") +
  xlab('Degree') + 
  ylab('Count') +
  scale_fill_manual(values = c("FALSE" = colors[1], "TRUE" = colors[3])) +  # Color based on hub threshold
  theme_minimal() +
  guides(fill = 'none') +
  theme(plot.title = element_text(size = 12), axis.text = element_text(size = 12))
show(plot_1)

# Number of hubs in the normal network
hubs.n <- d.n[d.n >= x_n]
length(hubs.n)  # Number of hubs: 71

# Check if the normal network follows a power-law distribution
d.n.table <- table(d.n)

# Convert the table to a data frame
d.n.fd <- data.frame(degree = as.numeric(names(d.n.table)),
                     count = as.numeric(d.n.table) / length(hubs.n))

# Plot the log-log degree distribution
plot_2 <- ggplot(d.n.fd, aes(x = log(degree), y = log(count))) +
  geom_point(fill = colors[1], color = colors[1], size = 2) +
  ggtitle("") +
  xlab('log (Degree)') + 
  ylab('log (Count)') +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))
show(plot_2)

# Now for the CANCER network:
d.c <- sna::degree(net.c, gmode = 'graph')

# Remove the nodes with 0 degree
d.c <- d.c[d.c > 0]
names(d.c) <- network.vertex.names(net.c)

# Determine the 95% quantile to identify hubs
x_c <- quantile(d.c, quantile_hubs)

data.c <- data.frame(value = unlist(d.c))
row.names(head(sort(data.c, decreasing = TRUE), 10))  # Display the top 10 nodes with highest degree

# Plot the degree distribution
plot_3 <- ggplot(data.c, aes(x = value, fill = value > x_c)) +
  geom_histogram(binwidth = 5, color = "snow") +
  ggtitle("Cancer Tissues") +
  xlab('Degree') + 
  ylab('Count') +
  scale_fill_manual(values = c("FALSE" = colors[2], "TRUE" = colors[3])) +  
  theme_minimal() +
  guides(fill = 'none') +
  theme(plot.title = element_text(size = 12)) 
show(plot_3)

# Number of hubs in the cancer network
hubs.c <- d.c[d.c >= x_c]
length(hubs.c)  # Number of hubs: 170

# Check if the cancer network follows a power-law distribution
d.c.table <- table(d.c)

# Convert the table to a data frame
d.c.fd <- data.frame(degree = as.numeric(names(d.c.table)),
                     count = as.numeric(d.c.table) / length(hubs.c))

# Plot the log-log degree distribution
plot_4 <- ggplot(d.c.fd, aes(x = log(degree), y = log(count))) +
  geom_point(fill = colors[2], color = colors[2], size = 2) +
  ggtitle("") +
  xlab('log (Degree)') + 
  ylab('log (Count)') +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))
show(plot_4)

# Arrange the plots in a grid
grid.arrange(plot_1, plot_3, plot_2, plot_4, nrow = 2)


## 3.3 Hubs - sub-network ----

# Compare the hubs sets

common <- intersect(names(hubs.c), names(hubs.n))
length(common)  # Number of common hubs: 1

# Identify hubs selectively characterizing each network
not_in_common_c <- names(hubs.c[!(names(hubs.c) %in% common)])
not_in_common_n <- names(hubs.n[!(names(hubs.n) %in% common)])

# Create the sub-network composed by the Hubs and their neighbors for Normal tissues
hubs.n.ids <- vector("integer", length(hubs.n))
for (i in 1:length(hubs.n)) {
  hubs.n.ids[i] <- match(names(hubs.n)[i], rownames(adj.mat.n))
}

hubs.n.neigh <- c()
for (f in hubs.n.ids) {
  hubs.n.neigh <- append(hubs.n.neigh, get.neighborhood(net.n, f))
}

hubs.n.neigh <- unique(hubs.n.neigh)
hubs.n.neigh.names <- rownames(adj.mat.n[hubs.n.neigh, ])
subnet.n <- unique(c(names(hubs.n), hubs.n.neigh.names))
hub.n.adj <- adj.mat.n[subnet.n, subnet.n]
names.hubs.n <- names(hubs.n)

# Save the adjacency matrix and hub names for the Normal tissues network
write.csv(hub.n.adj, file = "networks/hub_n_adj.csv", row.names = FALSE)
writeLines(names.hubs.n, "networks/hub_n_names.txt")

# Create the sub-network composed by the Hubs and their neighbors for Cancer tissues
hubs.c.ids <- vector("integer", length(hubs.c))
for (i in 1:length(hubs.c)) {
  hubs.c.ids[i] <- match(names(hubs.c)[i], rownames(adj.mat.c))
}

hubs.c.neigh <- c()
for (f in hubs.c.ids) {
  hubs.c.neigh <- append(hubs.c.neigh, get.neighborhood(net.c, f))
}

hubs.c.neigh <- unique(hubs.c.neigh)
hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh, ])
subnet.c <- unique(c(names(hubs.c), hubs.c.neigh.names))
hub.c.adj <- adj.mat.c[subnet.c, subnet.c]
names.hubs.c <- names(hubs.c)

# Save the adjacency matrix and hub names for the Cancer tissues network
write.csv(hub.c.adj, file = "networks/hub_c_adj.csv", row.names = FALSE)
writeLines(names.hubs.c, "networks/hub_c_names.txt")

# The plots have been made in the file network_visualization.py



# 4 - Differential Co-expression Network -----------------------------------

## 4.1 - Computation  ----

# Apply Fisher z-transformation to the correlation matrices
z.c <- 0.5*log((1 + rho.c)/(1- rho.c))
z.n <- 0.5*log((1 + rho.n)/(1- rho.n))

# z-scores
den <- sqrt(2/(n_patients - 3))
Z = (z.c - z.n)/den

# Create adjacency matrix
adj.mat.diff <- abs(Z) >= 3



##  4.2 - Analysis ----


# Create differential network
adj.mat.diff <- adj.mat.c != adj.mat.n  # Calculate differential adjacency matrix
net.diff <- as.network(adj.mat.diff, directed = FALSE)  # Create network from differential adjacency matrix

# Network properties
network.vertex.names(net.diff)[1:10]  # First 10 vertex names
network.size(net.diff)  # Number of nodes: 3250 
network.edgecount(net.diff)  # Number of edges: 512872
network.density(net.diff)  # Network density: 0.09714175

# Degree distribution
d.diff <- sna::degree(net.diff, gmode = 'graph')
d.diff <- d.diff[d.diff > 0]
names(d.diff) <- network.vertex.names(net.diff)

# Determine the 95% quantile to identify hubs
x_d <- quantile(d.diff, quantile_hubs)

data <- data.frame(value = unlist(d.diff))
plot_5 <- ggplot(data, aes(x = value, fill = value > x_d)) +
  geom_histogram(binwidth = 15, color = "snow") +
  ggtitle("Degree distribution - Differential Co-expressed Network") +
  xlab('Degree') + 
  ylab('Count') +
  scale_fill_manual(values = c("FALSE" = colors[3], "TRUE" = colors[4])) +  
  theme_minimal() +
  guides(fill = 'none') +
  theme(plot.title = element_text(size = 12))
show(plot_5)

# Identify hubs in the differential network
hubs.diff <- d.diff[d.diff >= x_d]

# Check if hubs.diff has common entries with not_in_common_n and not_in_common_c
length(intersect(names(hubs.diff), not_in_common_c))  # 0
length(intersect(names(hubs.diff), not_in_common_n))  # 8

# Power law: Check if the network is scale-free
d.diff.table <- table(d.diff)

# Convert the table to a data frame
d.diff.fd <- data.frame(degree = as.numeric(names(d.diff.table)),
                        count = as.numeric(d.diff.table) / length(hubs.diff))

x <- log(d.diff.fd$degree)
y <- log(d.diff.fd$count)

# Fit a linear model to log-log data
model <- lm(y[2:length(y)] ~ x[2:length(x)])
slope <- round(model$coefficients[2], 2)

plot_6 <- ggplot(d.diff.fd, aes(x = log(degree), y = log(count))) +
  geom_point(fill = colors[3], color = colors[3], size = 2) +
  ggtitle("Power Law Distribution of Degree - Differential Co-expressed Network") +
  xlab('log (Degree)') + 
  ylab('log (Count)') +
  coord_cartesian(ylim = c(-6, -1.5)) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 13))
show(plot_6)

# Plot the power law fit
plot(x[2:length(x)], model$fitted.values, type = 'l', ylim = c(-6, -1), xlab = 'logarithm of degree', ylab = 'logarithm of count', col = 'blue4', 
     main = 'Power Law - Differential Network')
points(log(d.diff.fd$degree), log(d.diff.fd$count), 
       xlab = "Degree", ylab = "Degree Count", 
       main = "Degree Distribution", 
       pch = 16, col = "gold")

text(x = 4.9, y = -2, labels = paste0("slope = ", slope), pos = 2, offset = 1, col = "blue4", cex = 1)

# Arrange the plots
grid.arrange(plot_5, plot_6, nrow = 2)



## 4.3 Hubs - sub-network ----

# Create the sub-network for the differential network

# Get the indices of the hubs in the differential network
hubs.diff.ids <- vector("integer", length(hubs.diff))
for (i in 1:length(hubs.diff)) {
  hubs.diff.ids[i] <- match(names(hubs.diff)[i], rownames(adj.mat.diff))
}

# Get the neighbors of the hubs in the differential network
hubs.diff.neigh <- c()
for (f in hubs.diff.ids) {
  hubs.diff.neigh <- append(hubs.diff.neigh, get.neighborhood(net.diff, f))
}

hubs.diff.neigh <- unique(hubs.diff.neigh)
hubs.diff.neigh.names <- rownames(adj.mat.diff[hubs.diff.neigh, ])
subnet.diff <- unique(c(names(hubs.diff), hubs.diff.neigh.names))
hub.diff.adj <- adj.mat.diff[subnet.diff, subnet.diff]
names.hubs.diff <- names(hubs.diff)

# Save the adjacency matrix and hub names for the differential network
write.csv(hub.diff.adj, file = "networks/hub_diff_adj.csv", row.names = FALSE)
writeLines(names.hubs.diff, "networks/hub_diff_names.txt")

# The plots have been made in the file network_visualization.py

# Compare the hubs obtained with the differential expressed network to those of the previous points
intersect(names(hubs.c), names(hubs.diff)) #1
intersect(names(hubs.n), names(hubs.diff)) #9




# 5 - Patient Similarity Network ------------------------------------------

# Full cancer data normalization and preparation

full.data.c <- as.data.frame(rna.expr.data.C)
colnames(full.data.c) <- substr(colnames(full.data.c), 1, 12)

metad <- rep("cancer", dim(full.data.c)[2])
metad <- data.frame(metad)
colnames(metad)[1] <- "condition"
metad[, 1] <- as.factor(metad[, 1])

full.data.c <- cbind(names = rownames(full.data.c), full.data.c)

dds.c <- DESeqDataSetFromMatrix(countData = full.data.c, 
                                colData = metad,
                                design = ~1,
                                tidy = TRUE)

dds.c <- estimateSizeFactors(dds.c)
normalized_counts.c <- counts(dds.c, normalized = TRUE)

# Keep only the genes without zeros
normalized_counts.c <- normalized_counts.c[-which(rowSums(normalized_counts.c == 0) >= 1), ]

filtr.expr.c.all <- as.data.frame(normalized_counts.c)

# Rename columns to shorter versions
colnames(filtr.expr.c.all) <- substr(colnames(filtr.expr.c.all), 1, 12)

# Correlation for the cancer network using Pearson correlation
cor.mat.cp <- corr.test(filtr.expr.c.all, use = "pairwise", 
                        method = "pearson", adjust = "fdr", ci = FALSE)

# Matrix with the Pearson correlation coefficient
rho.cp <- cor.mat.cp$r
# Set the diagonal to zero (no self-loops)
diag(rho.cp) <- 0

# Binary matrix, keep the edge for significant correlation
adj.mat.cp <- rho.cp * (abs(rho.cp) > 0.8)

rows_to_keep <- apply(adj.mat.cp, 1, function(row) any(row != 0))
cols_to_keep <- apply(adj.mat.cp, 2, function(col) any(col != 0))
adj.mat.cp <- adj.mat.cp[rows_to_keep, cols_to_keep]

# Create network from adjacency matrix
net.cp <- as.network(adj.mat.cp, directed = FALSE)

l.comp.p <- component.largest(net.cp, result = "graph")
l.comp.p <- adj.mat.cp[rownames(l.comp.p), rownames(l.comp.p)]

# Save the largest component of the cancer network
write.csv2(l.comp.p, "input-matrix.csv")

# Assume external community detection tool is run here
# python btc-community-c.py input-matrix-c.csv

# Save the adjacency matrix and cancer stage info
write.csv(adj.mat.cp, file = "networks/community_adj.csv", row.names = FALSE)

cancer.stage.info <- clinical.query[, c("ajcc_pathologic_stage", "submitter_id")]
selected.patients <- colnames(adj.mat.cp)
cancer.stage.info <- subset(cancer.stage.info, submitter_id %in% selected.patients)
colnames(cancer.stage.info) <- c('stage', 'name')

write.csv(cancer.stage.info, file = "networks/cancer_stage.csv", row.names = FALSE)

# The plots have been made in the file network_visualization.py


# 6 - OPTIONAL TASKS ----------------------


## 6.1 - Compute a different centrality index and check overlap with previous hubs -----
#previous hubs (degree-based):
hubs.c 
hubs.n 

bet_n <- betweenness(net.n)   #for normal tissue network
bet_c <- betweenness(net.c)   #for cancer tissue network

names(bet_n) <- rownames(filtr.expr.n.deg)  #assigning row names
names(bet_c) <- rownames(filtr.expr.c.deg)  #assigning row names
 
hubs.n.bet <- names(bet_n[bet_n > quantile(bet_n, 0.95)]) #top 5% hubs in normal tissue
hubs.c.bet <- names(bet_c[bet_c > quantile(bet_c, 0.95)]) #top 5% hubs in cancer tissue

# Normal tissue hubs in common
intersect(names(hubs.n), hubs.n.bet)   #1

# Cancer tissue hubs in common
intersect(names(hubs.c), hubs.c.bet)    #25


## 6.2 - Enrichment analysis on cancer nodes ---------------  

# Extract common names of the genes 
hubs.c_names <- genes.info.c[names(hubs.c), ]$gene_name


#not_in_common_c_names <- genes.info.c[not_in_common_c, ]$gene_name

library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()
if (websiteLive) head(dbs)
if (websiteLive) dbs$libraryName

dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
if (websiteLive) {  enriched <- enrichr(hubs.c_names, dbs) }

#Plotting the results

# Database GO 
if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 20, numChar = 45, y = "Count", orderBy = "P.value")}

# Database KEGG
if (websiteLive) {
  plotEnrich(enriched[[2]], showTerms = 20, numChar = 45, y = "Count", 
             orderBy = "P.value") }


## 6.3 - Survival Analysis ---------------  
install.packages("survival")
install.packages("survminer")

library(survival)
library(survminer)

survival.data <- clinical.query[, c("submitter_id", "days_to_death", "vital_status", "ajcc_pathologic_stage")]
colnames(survival.data) <- c('name', 'time', 'status', 'stage')
selected.patients <- network.vertex.names(net.final.p)

survival_patients <- survival.data$name
selected_patients <- selected.patients
common_patients <- intersect(survival_patients, selected_patients)
survival_data_filtered <- survival.data[survival.data$name %in% common_patients, ]

survival_data_filtered$status <- ifelse(survival_data_filtered$status == "Dead", 1, 0)

surv_object <- Surv(time = survival_data_filtered$time, event = survival_data_filtered$status)
fit <- survfit(surv_object ~ stage, data = survival_data_filtered)

ggsurvplot(fit, data = survival_data_filtered, pval = TRUE, risk.table = TRUE)

ggsurvplot(
  fit,                      # Fitted survival model
  data = survival_data_filtered ,            # Show p-value of the survival test
  risk.table = TRUE,        # Show risk table
  palette = c("#E41A1C", "#377EB8","green","black"),  # Optional: Customize color palette
  main = "Survival Curve",  # Optional: Main title for the plot
  xlab = "Time",   # Optional: X-axis label
  ylab = "Survival Probability",  # Optional: Y-axis label
  legend.title = "Stages:",
  legend.labs = c("Stage 1", "Stage 2", "Stage 3", "Stage 4"),
  ggtheme = theme_minimal()
)


