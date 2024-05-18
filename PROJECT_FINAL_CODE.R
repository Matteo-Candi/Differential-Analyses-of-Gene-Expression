# 0 - Load packages --------------

BiocManager::install("BiocGenerics")
BiocManager::install("DESeq2")
install.packages("GGally")

required_packages <- c("BiocGenerics", "DESeq2", "psych", "NetworkToolbox", "ggplot2",
                       "GGally", "sna", "network", "TCGAbiolinks", "SummarizedExperiment", "DT", "latex2exp", "gridExtra")

lapply(required_packages, library, character.only = TRUE)

seed <- 123
set.seed(seed)

colors <- c("#2EBFA5", "#F26419", "#FDCA40")


# 1 - Load data and Preprocessing -----------------------------------------------------------

# Project ID
proj <- "TCGA-CHOL"
dir.create(file.path(proj))

# Look for all data linked to a "Primary Tumor" sample and store them in a dataframe
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")
GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
rna.data.C <- GDCprepare(rna.query.C, directory = "GDCdata")
rna.expr.data.C <- assay(rna.data.C)
genes.info.c <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))


# Apply the same procedure to the "Normal tissue"
rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")
GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata" )
rna.expr.data.N <- assay(rna.data.N)
genes.info.n <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))



all(na.omit(genes.info.n) == na.omit(genes.info.c))
clinical.query<- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)

table(clinical.query$ajcc_pathologic_stage)

# Investigate the distribution of age and cancer stage of our data


ggplot(clinical.query[!is.na(clinical.query$age_at_index), ], aes(x = ajcc_pathologic_stage, y = age_at_index)) +
  geom_boxplot(fill = colors[3]) +
  labs(title = "Relationship between Age and Cancer Stage",
       x = "",
       y = "Age") +
  coord_cartesian(ylim = c(20, 90)) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 12))




boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical.query,
        col = "gold1", main = "Relationship between Age and Cancer Stage", xlab = "", ylab= "age", las=1 )





View(rna.expr.data.N)

# First let's look at the amount of data we have for each class (c and n)
dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) #there are no duplicates for normal tissue

dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) #there are no duplicates for cancer tissue 

# Total number of unique patients
length(unique(clinical.query$submitter_id))

# IDs of the "primal tumor" class patients
patients.C <- substr(colnames(rna.expr.data.C), 1,12)

#  [TO DELETE]
# # select index of the patients that are unique
# unique.patients.C <- names(which(table(patients.C) == 1))
# idx.unique.pats <- match(unique.patients.C, substr(colnames(rna.expr.data.C), 1,12))

expr.C <- as.data.frame(rna.expr.data.C)
expr.N <- as.data.frame(rna.expr.data.N)

# Change the patients IDs to a shorter version
colnames(expr.C) <- substr(colnames(expr.C), 1,12)
colnames(expr.N) <- substr(colnames(expr.N), 1,12)

# For how many of them we have both a tumor and normal sample?
intersect(colnames(expr.N), colnames(expr.C)) 
length(intersect(colnames(expr.N), colnames(expr.C))) # 8

# We have 35 sample of patients with cancer and 9 samples of patients without.We have just 8 samples of tissue from the same patients that lead to a group of 36 patients


# Since we want to do a comparison between the sane and ill patients we need to get rid of the one for which we don't have the "primal tumor" sample
idx <- match(setdiff(colnames(expr.N), colnames(expr.C)), colnames(expr.N)) #idx to remove
expr.N <- expr.N[,-c(idx)]

# And select those same patients from the "primal tumor" data
expr.C <- expr.C[, colnames(expr.N)]

length(intersect(colnames(expr.N), colnames(expr.C)))

# Now let's check if the values stored are ok

#let's check the actual counts
typeof(expr.C[1,1]) #ok
any(is.na(expr.C)) #ok
any(is.nan(as.matrix(expr.C))) #ok

typeof(expr.N[1,1]) #ok
any(is.na(expr.N)) #ok
any(is.nan(as.matrix(expr.N))) #ok


# Explanation for the normalization:
# detalied explanation: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# youtube explanation: https://www.youtube.com/watch?v=UFB993xufUU


# Merge the data and transform it in a dataframe
full.data <- cbind(expr.N, expr.C)
full.data <- data.frame(full.data)

# Create a new dataframe in which we specify which kind of tissue the information
# come from "cancer" or "normal"

n_patients <- ncol(full.data)

metad <- rep("cancer", n_patients)
metad[1:(n_patients / 2)] <- "normal"
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1])
full.data <- cbind(rownames(full.data), full.data)

# Transform all the data in a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=full.data, 
                              colData=metad, 
                              design= ~condition,
                              tidy=TRUE)


# Filtering: more than 10 counts on 90% of patients
limit.90 <- floor(n_patients * .9) #15

# Index of the rows that are at least 10 for 90% of the patients
keep <- rowSums(counts(dds) >= 10) >= limit.90
dds <- dds[keep,]

# Number of genes we end up having
dim(counts(dds))

# Now let's normalize
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == n_patients) #no null rows

# Keep only the genes w/o zeros
normalized_counts <- normalized_counts[-which(rowSums(normalized_counts == 0) >= 1),]

# Divide the cancerous and normal data
filtr.expr.n <- as.data.frame(normalized_counts[, 1:(n_patients/2)])
filtr.expr.c <- as.data.frame(normalized_counts[, (n_patients/2+1):n_patients])

#cancerous sample names were added a ".1" in full.data because  
#they had the same names as the normal samples
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1,12)


# 2 - Differentially Expressed Genes (DEGs) -------------------------------

# To Identify the Differently Expressed Genes we first need to compute the 
# FC (fold change)
fc <-  log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n)) 
names(fc) <- rownames(filtr.expr.c)

# And the p-value using the t test
pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(filtr.expr.c[i,], filtr.expr.n[i,] ))$p.value)
# And then we can adjust the values
pval.fc.fdr <- p.adjust(pval.fc, method="fdr")

# Then we can collect both the fc and pvalues inside a df
expr.table <- data.frame(cbind(fc, pval.fc.fdr))
expr.table[,1] <- round(expr.table[,1],2)


# Now given the FCs and p-values we would like to identify which are the genes 
# that are different for both the measures

# So we can set a threshold for the FC and one for the pvalue so that 
# we obtain hundreds of genes in the end
threshold.p <- 0.05
threshold.fc <- 1.2

# And then we select only the genes that have the absolute value of the fold change 
# over the threshold.fc and the p-value under threshold.p at the same time
# In this way we are sure that the genes are distributed differently depending 
# on the classes 

deg.genes <- rownames(expr.table[abs(expr.table$fc) >= threshold.fc & expr.table$pval.fc.fdr <= threshold.p,]) 
length(deg.genes) # We have 3250 different genes

# Let's select just the DEG genes
filtr.expr.c.deg <- filtr.expr.c[deg.genes,]
filtr.expr.n.deg <- filtr.expr.n[deg.genes,]

# This difference between the visualized using a volcano plot
expr.table$diffexpressed <- "NO SIGNIFICANT";
# the genes with FC >= threshold.fc
expr.table$diffexpressed[expr.table$fc >= threshold.fc & expr.table$pval.fc.fdr <= threshold.p] <- "UP"
# the genes with FC <= -threshold.fc
expr.table$diffexpressed[expr.table$fc <= -threshold.fc & expr.table$pval.fc.fdr <= threshold.p] <- "DOWN"
expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)

line_col <- 'darkgrey'

ggplot(data = expr.table, aes(x = fc, y = -log10(pval.fc.fdr), col = diffexpressed)) +
  geom_point() +
  xlab(expression(log[2]~"(Fold Change)")) +
  ylab(expression(-log[10]~"(Adjusted p-value)")) +
  geom_hline(yintercept = -log10(threshold.p), col = line_col) +
  geom_vline(xintercept = c(-threshold.fc, threshold.fc), col = line_col) +
  labs(title = "Volcano Plot", color = "") +
  coord_cartesian(ylim = c(-.5, 6.5), xlim=c(-8.5, 7)) +
  scale_color_manual(values = c("UP" = colors[1], "DOWN" = colors[2], "NO SIGNIFICANT" = "lightgrey")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_minimal() 


ggplot(data=expr.table, aes(x=fc, y=-log10(pval.fc.fdr), col=diffexpressed))+  
  geom_point() +
  xlab("log2(FC)") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=-log10(threshold.p), col=line_col)+
  geom_vline(xintercept=threshold.fc, col=line_col)+
  geom_vline(xintercept=-threshold.fc, col=line_col)



# 3 - Co-expression networks ----------------------------------------------
## 3.1 - Computation  ----

# Now we would like to build the co-expression network, to do so we first need
# to compute the adjacency matrix

# We will define the matrix using the correlation between the genes so
# the edge between two genes will be present if the absolute value of the 
# correlation is over a threshold threshold.rho

threshold.rho <- 0.8 

# Before computing correlations, we log transform the data using log2(x+1)

filtr.expr.c.deg <- log(filtr.expr.c.deg+1)
filtr.expr.n.deg <- log(filtr.expr.n.deg+1) 

# Create the CANCER NETWORK

# Correlation matrix for the cancer network - we will use Spearman correlation now
cor.mat.c <- corr.test(t(filtr.expr.c.deg), use = "pairwise", 
                       method = "spearman", adjust="fdr", ci=FALSE)

# matrix with the Spearman correlation coefficient
rho.c <- cor.mat.c$r
# Set the diagonal to zero (no self loops)
diag(rho.c) <- 0
# Binary matrix, keep the edge for big correlation
adj.mat.c <- abs(rho.c) > threshold.rho

# Same process for the NORMAL NETWORK

# Correlation for the normal network using Spearman
cor.mat.n <- corr.test(t(filtr.expr.n.deg), use = "pairwise", 
                       method = "spearman", adjust="fdr", ci=FALSE)

rho.n <- cor.mat.n$r
diag(rho.n) <- 0

adj.mat.n <- abs(rho.n) > threshold.rho


## 3.2 - Analysis  ----

# Create the two networks given the adjacency matrices
net.n <- as.network(adj.mat.n, directed = FALSE)
net.c <- as.network(adj.mat.c, directed = FALSE)

# Let's check the number of edges and nodes

# Normal Network
network.size(net.n) # 3250
network.edgecount(net.n) # 278221

# Cancer Network
network.size(net.c) # 3250
network.edgecount(net.c) # 139842


# The number of edge for a normal patient is higher than the other one meaning that genes ar more correlated to each other for normal patients.


# We look for a scale-free network in gene expression of cancer tissues because identifying key hub genes can reveal critical regulatory elements and potential therapeutic targets that drive cancer progression.
# To determine if a network is scale-free, we need to examine its degree distribution. In this case, a few nodes (hubs) have many connections, while most nodes have few. For a co-expression network, this means a small number of genes have numerous co-expression relationships, while most have only a few. If the degree distribution follows a power-law, the network is scale-free.




# First for the NORMAL graph: 
# extract the degree of each node
d.n <- sna::degree(net.n, gmode = 'graph')

# Remove the nodes with 0 degree.
d.n <- d.n[d.n>0]
names(d.n) <- network.vertex.names(net.n)

# Print the histogram of the degree together with a line for the 95% quantile
x_n <- quantile(d.n,0.95)

data <- data.frame(value = unlist(d.n))
plot_1 <- ggplot(data, aes(x = value, fill = value > x_n)) +
  geom_histogram(binwidth = 15, color = "snow") +
  ggtitle("Normal Tissues") +
  xlab('Degree') + 
  ylab('Count') +
  scale_fill_manual(values = c("FALSE" = colors[1], "TRUE" = colors[3])) +  # Specify colors based on condition
  # geom_segment(x = x_n, xend = x_n, yend = 450, y = 0, linetype = "dashed", color = colors[2], size = 0.9) +
  theme_minimal() +
  guides(fill = 'none') +
  theme(plot.title = element_text(size = 12))
show(plot_1)
  

hist(d.n,col = "gold", main = "Degree distribution - Normal Tissues", breaks = 50)
abline(v=x_n, col="blue4", lwd = 2, lty = 2)
box()
# It doesn't look scale free

# The number of hubs of the network is
hubs.n <- d.n[d.n>=x_n]
length(hubs.n) # 170

# Power law: Another way to check if a network is scale free is to verify if it
# follows the power law. i.e. there should be a link between the fraction of 
# nodes f(k) with degree k and k itself f(k) ⁓ k^(-γ) (usually 2< γ <3)
d.n.table <- table(d.n)

# Convert the table to a data frame
d.n.fd <- data.frame(degree = as.numeric(names(d.n.table)),
                        count = as.numeric(d.n.table)/length(hubs.n))

plot_2 <- ggplot(d.n.fd, aes(x = log(degree), y = log(count))) +
  geom_point(fill = colors[1], color= colors[1], size=2) +
  ggtitle("") +
  xlab('log (Degree)') + 
  ylab('log (Count)') +
  coord_cartesian(ylim = c(-6, -1)) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 12))
show(plot_2)

plot(log(d.n.fd$degree), log(d.n.fd$count), 
     xlab = "Degree", ylab = "Degree Count", 
     main = "Degree Distribution", 
     pch = 16, col = "gold")
# For sure not scale free


d.c <- sna::degree(net.c, gmode = 'graph')

# Remove the nodes with 0 degree.
d.c <- d.c[d.c>0]
names(d.c) <- network.vertex.names(net.c)

x_c <- quantile(d.c, 0.95)

data <- data.frame(value = unlist(d.c))
plot_3 <- ggplot(data, aes(x = value, fill = value > x_c + 4)) +
  geom_histogram(binwidth = 5, color = "snow") +
  ggtitle("Cancer Tissues") +
  xlab('Degree') + 
  ylab('Count') +
  scale_fill_manual(values = c("FALSE" = colors[2], "TRUE" = colors[3])) +  
  # geom_segment(x = x_n, xend = x_n, yend = 450, y = 0, linetype = "dashed", color = colors[2], size = 0.9) +
  theme_minimal() +
  guides(fill = 'none') +
  theme(plot.title = element_text(size = 12)) 
show(plot_3)

data <- data.frame(value = unlist(d.c))
ggplot(data, aes(x = value)) +
  geom_histogram(binwidth = 5, fill = colors[2], color = "snow") +
  ggtitle("Degree distribution - Cancer Tissues") +
  xlab('Degree') + 
  ylab('Count') +
  geom_segment(x=x_c, xend = x_c, yend=230, y=0, linetype = "dashed", color = colors[1], size = 1) +
  theme_minimal()
show(plot_3)

hist(d.c,col = "gold", main = "Degree distribution - cancer tissue", breaks = 50)
abline(v = x_c, col="blue4", lwd = 2, lty = 2)
box()
# It doesn't look scale free

# number of hubs
hubs.c <- d.c[d.c>=x_c]
length(hubs.c) # 170

# let's check the Power Law
d.c.table <- table(d.c)

# Convert the table to a data frame
d.c.fd <- data.frame(degree = as.numeric(names(d.c.table)),
                     count = as.numeric(d.c.table)/length(hubs.c))

ggplot(d.c.fd, aes(x = log(degree), y = log(count))) +
  geom_point(fill = colors[2], color= colors[2], size=2) +
  ggtitle("Power Law Distribution of Degree - Cancer Tissues") +
  xlab('log (Degree)') + 
  ylab('log (Count)')

plot_4 <- ggplot(d.c.fd, aes(x = log(degree), y = log(count))) +
  geom_point(fill = colors[2], color= colors[2], size=2) +
  ggtitle("") +
  xlab('log (Degree)') + 
  ylab('log (Count)') +
  coord_cartesian(ylim = c(-6, -.5)) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 12))
show(plot_4)

plot(log(d.c.fd$degree), log(d.c.fd$count), 
     xlab = "Degree", ylab = "Degree Count", 
     main = "Degree Distribution", 
     pch = 16, col = "gold")

# It is not scale free

grid.arrange(plot_1, plot_3, plot_2, plot_4, nrow = 2)


## 3.3 [TODO] hub - subnetwork ----


# Compare the hubs sets 

common <- intersect(names(hubs.c),names(hubs.n))
length(common) # 6 hubs in common

# let's identify those selectively characterizing the network - those we are REALLY interested in

not_in_common_c <- names(hubs.c[!(names(hubs.c) %in% common)])

not_in_common_n <- names(hubs.n[!(names(hubs.n) %in% common)])



# Create the hub subnetwork for the cancer graph (the one with a possible scale-free structure)

hubs.c.ids <- vector("integer",length(hubs.c))
for (i in 1:length(hubs.c)){hubs.c.ids[i] <- match(names(hubs.c)[i],rownames(adj.mat.c))}

# Identifying the neighborhood of the hubs

hubs.c.neigh <- c()
for (f in hubs.c.ids){
  hubs.c.neigh <- append(hubs.c.neigh, get.neighborhood(net.c, f))
}

hubs.c.neigh <- unique(hubs.c.neigh)
hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh,])

# Select only hubs and their neighbors
subnet <- unique(c(names(hubs.c), hubs.c.neigh.names))

# Creating the subnetwork
hub.c.adj <- adj.mat.c[subnet, subnet]

names.hubs <-names(hubs.c)

rownames(hub.c.adj)[1:length(hubs.c)] <- names.hubs
colnames(hub.c.adj)[1:length(hubs.c)] <- names.hubs


net.hub <- network(hub.c.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
network.density(net.hub)

sum(hub.c.adj > 0)

net.hub %v% "type" = ifelse(network.vertex.names(net.hub) %in% names.hubs,"hub", "non-hub")
net.hub %v% "color" = ifelse(net.hub %v% "type" == "non-hub", "deepskyblue3", "black")

ggnet2(net.hub,  color = "color",alpha = 0.9, size = 2,
       edge.color = "grey", edge.alpha = 0.15,  edge.size = 0.15,
       node.label = names.hubs, label.color = "black", label.size = 4)+
  guides(size = "none")






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
net.diff <- as.network(adj.mat.diff, directed = FALSE)

network.vertex.names(net.diff)[1:10]
network.size(net.diff) # 3250 
network.edgecount(net.diff) # 512872

network.density(net.diff) # 0.09714175


d.diff <- sna::degree(net.diff, gmode = 'graph')
d.diff <- d.diff[d.diff>0]
names(d.diff) <- network.vertex.names(net.diff)
sort(d.diff)

x_d <- quantile(d.diff,0.95)

data <- data.frame(value = unlist(d.diff))
plot_5 <- ggplot(data, aes(x = value, fill = value > x_d-1)) +
  geom_histogram(binwidth = 15, color = "snow") +
  ggtitle("Degree distribution - Differential Co-expressed Network") +
  xlab('Degree') + 
  ylab('Count') +
  scale_fill_manual(values = c("FALSE" = colors[3], "TRUE" = colors[4])) +  
  # geom_segment(x = x_n, xend = x_n, yend = 450, y = 0, linetype = "dashed", color = colors[2], size = 0.9) +
  theme_minimal() +
  guides(fill = 'none') +
  theme(plot.title = element_text(size = 12)) 
show(plot_5)


hist(d.diff,col = "gold", main = "Degree distribution - differential co-expressed network", xlab = 'Degree', breaks = 50)
abline(v=x_d, col="blue4", lwd = 2, lty = 2)
box()

hubs.diff <- d.diff[d.diff>=x_d]


# Check if hubs.diff has common entries with not_in_common_n and not_in_common_c

length(intersect(names(hubs.diff), not_in_common_c)) #7
length(intersect(names(hubs.diff), not_in_common_n)) #100

# Power law: Another way to check if a network is scale free is to verify if it
# follows the power law. i.e. there should be a link between the fraction of 
# nodes f(k) with degree k and k itself f(k) ⁓ k^(-γ) (usually 2< γ <3)
d.diff.table <- table(d.diff)

# Convert the table to a data frame
d.diff.fd <- data.frame(degree = as.numeric(names(d.diff.table)),
                     count = as.numeric(d.diff.table)/length(hubs.diff))


x <- log(d.diff.fd$degree)
y <- log(d.diff.fd$count)

model <- glm.fit(x[2:length(x)], y[2:length(y)])
slope <- round(model$coefficients, 2)

plot_6 <- ggplot(d.diff.fd, aes(x = log(degree), y = log(count))) +
  geom_point(fill = colors[3], color= colors[3], size=2) +
  ggtitle("Power Law Distribution of Degree - Differential Co-expressed Network") +
  xlab('log (Degree)') + 
  ylab('log (Count)') +
  coord_cartesian(ylim = c(-6, -1.5)) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 13))
show(plot_6)


plot(x[2:length(x)], model$fitted.values, type='l', ylim = c(-6, -1), xlab = 'logarithm of degree', ylab = 'logarith of count', col = 'blue4', 
     main = 'Power Law - Differential Network')
points(log(d.diff.fd$degree), log(d.diff.fd$count), 
     xlab = "Degree", ylab = "Degree Count", 
     main = "Degree Distribution", 
     pch = 16, col = "gold")

text(x = 4.9, y = -2, labels = paste0("slope = ", slope), pos = 2, offset = 1, col = "blue4", cex = 1)

# It doesn't look scale Free
grid.arrange(plot_5, plot_6, nrow = 2)


## 4.3 [TODO] hub - subnetwork ----

# Create the sub-network

hubs.diff.ids <- vector("integer",length(hubs.diff))
for (i in 1:length(hubs.diff)){hubs.diff.ids[i] <- match(names(hubs.diff)[i],rownames(adj.mat.diff))}
hubs.diff.ids

#identifying the neighborhood
hubs.diff.neigh <- c()
for (f in hubs.diff.ids){
  hubs.diff.neigh <- append(hubs.diff.neigh, get.neighborhood(net.diff, f))
}

hubs.diff.neigh <- unique(hubs.diff.neigh)
hubs.diff.neigh
hubs.diff.neigh.names <- rownames(adj.mat.diff[hubs.diff.neigh,])
subnet.diff <- unique(c(names(hubs.diff), hubs.diff.neigh.names))

#creating the subnetwork
hub.diff.adj <- adj.mat.diff[subnet.diff, subnet.diff]

names.hubs.diff <-names(hubs.diff)

rownames(hub.diff.adj)[1:length(hubs.diff)] <- names.hubs.diff
colnames(hub.diff.adj)[1:length(hubs.diff)] <- names.hubs.diff

net.hub.diff <- network(hub.diff.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
network.density(net.hub.diff) # 0.0987710

sum(hub.diff.adj > 0 )

node_degrees <- degree(net.hub.diff)

net.hub.diff %v% "type" = ifelse(network.vertex.names(net.hub.diff) %in% names.hubs.diff,"hub", "non-hub")
net.hub.diff %v% "color" = ifelse(net.hub.diff %v% "type" == "non-hub", "deepskyblue3", "tomato")

ggnet2(net.hub.diff,  color = "color",alpha = 0.9, size = 2,
       edge.color = "grey", edge.alpha = 0.5,  edge.size = 0.15)+
  guides(size = "none")

# Compare the hubs obtained with the differential expressed network to those of the previous points

intersect(names(hubs.c), names(hubs.diff))
intersect(names(hubs.n), names(hubs.diff))






# 5 - Patient Similarity Network ------------------------------------------

# Correlation for the cancer network using Pearson correlation
cor.mat.cp <- corr.test(filtr.expr.c, use = "pairwise", 
                       method = "pearson", adjust="fdr", ci=FALSE)

# matrix with the pearson correlation coefficient
rho.cp <- cor.mat.cp$r
# Set the diagonal to zero (we are not interested to the edge with a gene with itself)
diag(rho.cp) <- 0

# Binary matrix, keep the edge for big correlation
qval.cp <- cor.mat.cp$p
qval.cp[lower.tri(qval.cp)] <- t(qval.cp)[lower.tri(qval.cp)]
qval.cp

# correlation network of cancer samples

# Since the p-values are all 0 we can try to define the adjacency only looking at the correlation values
adj.mat.cp <- rho.cp * (abs(rho.cp) > 0.7)   

net.cp <- as.network(adj.mat.cp, directed = FALSE)


l.comp.p <- component.largest(net.cp, result = "graph")
l.comp.p <- adj.mat.cp[rownames(l.comp.p), rownames(l.comp.p)]

# Save the cancer matrix

write.csv2(l.comp.p, "input-matrix-c.csv")

# in the terminal: 

# python btc-community-c.py input-matrix-c.csv

# Read results 

comm.res.p <- read.csv2("output-c.txt", header = FALSE)
rownames(comm.res.p) <- rownames(l.comp.p)
sort(table(comm.res.p[,1]), decreasing = T)
length(table(comm.res.p[,1]))

net.final.p <- network(l.comp.p, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)


# Assign communities
net.final.p  %v% "community" <-  as.character(comm.res.p[,1])

ncom <- length(unique(net.final.p  %v% "community"))

pal <- sample(colors(distinct = T), ncom)
names(pal) <- 1:ncom

node_mapping <- data.frame(real_label = unique(net.final.p %v% "vertex.names"),
                           new_label = 1:length(unique(net.final.p %v% "vertex.names")))

# Update the network with new labels (to better identify patients we label them with numbers from 1 to 18)
net.final.p %v% "vertex.names" <- node_mapping$new_label


# Plot subnetwork

ggnet2(net.final.p, color = "community", palette =  pal, alpha = 1, 
       size = 5, edge.color = "grey", edge.alpha = 1, edge.size = 0.15, label = TRUE, label.size = 5)+
  guides(size = "none") 




# Correlation for the normal network
cor.mat.np <- corr.test(filtr.expr.n, use = "pairwise", 
                        method = "pearson", adjust="fdr", ci=FALSE)

# matrix with the pearson correlation coefficient
rho.np <- cor.mat.np$r
# Set the diagonal to zero (we are not interested to the edge with a gene with itself)
diag(rho.np) <- 0
# Binary matrix, keep the edge for big correlation

qval.np <- cor.mat.np$p
qval.np[lower.tri(qval.np)] <- t(qval.np)[lower.tri(qval.np)]

# correlation network of cancer samples
#adj.mat.np <- rho.np * (qval.np <= 1e-3)  
                                           
adj.mat.np <- rho.np * (abs(rho.np) > 0.7)


net.np <- as.network(adj.mat.np, directed = FALSE)


l.comp.pn <- component.largest(net.np, result = "graph") 
l.comp.pn <- adj.mat.np[rownames(l.comp.pn), rownames(l.comp.pn)]

write.csv2(l.comp.pn, "input-matrix-n.csv")

# in the terminal :
# python btc-community-n.py input-matrix-n.csv

comm.res.pn <- read.csv2("output-n.txt", header = FALSE)
rownames(comm.res.pn) <- rownames(l.comp.pn)
sort(table(comm.res.pn[,1]), decreasing = T)
length(table(comm.res.pn[,1]))

net.final.pn <- network(l.comp.pn, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)

net.final.pn  %v% "community" <-  as.character(comm.res.pn[,1])

ncom_n <- length(unique(net.final.pn  %v% "community"))

set.seed(123)
pal <- sample(colors(distinct = T), ncom_n)
names(pal) <- 1:ncom_n
pal

vertex_names <- as.character(net.final.pn %v% "vertex.names")
vertex_data <- data.frame(name = vertex_names, x = layout_nicely(net.final.pn)[,1], y = layout_nicely(net.final.pn)[,2])

net.final.pn %v% "vertex.names" <- node_mapping$new_label


ggnet2(net.final.pn, color = "community", palette =  pal, alpha = 1,
       size = 5, edge.color = "grey", edge.alpha = 1, edge.size = 0.15, label = TRUE, label.size = 5)+
  guides(size = "none")



# 6 - OPTIONAL TASKS ----------------------


## 6.1 - Compute a different centrality index and check overlap with previous hubs -----

hubs.c #previous hubs (degree-based)
hubs.n 


bet_n <- betweenness(net.n)
bet_c <- betweenness(net.c)

names(bet_n) <- rownames(filtr.expr.n.deg)
names(bet_c) <- rownames(filtr.expr.c.deg)

hubs.n.bet <- names(bet_n[bet_n > quantile(bet_n, 0.95)])
hubs.c.bet <- names(bet_c[bet_c > quantile(bet_c, 0.95)])

# Normal tissue hubs in common
intersect(names(hubs.n), hubs.n.bet)   #3

# Cancer tissue hubs in common
intersect(names(hubs.c), hubs.c.bet)    #3





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




## 6.3 - repeat the analysis with Spearman correlation --------------------------------------


## 6.4 - Co-expression networks ----

## 6.5 - Computation of the Adj matrices ----

# Now we would like to build the co-expression network, to do so we first need
# to compute the adjacency matrix

# We will define the matrix using the correlation between the genes so
# the edge between two genes will be present if the absolute value of the 
# correlation is over a threshold
threshold.rho <- 0.8

# Before computing correlations, we transform the data using log2(x+1)

filtr.expr.c.deg <- log2(filtr.expr.c.deg+1)

# Correlation for the cancer network
cor.mat.c <- corr.test(t(filtr.expr.c.deg), use = "pairwise", 
                       method = "spearman", adjust="fdr", ci=FALSE)

# matrix with the pearson correlation coefficient
rho.c <- cor.mat.c$r
# Set the diagonal to zero (we are not intrested to the edge with a gene with itself)
diag(rho.c) <- 0
# Binary matrix, keep the edge for big correlation
adj.mat.c <- 1 * (abs(rho.c) > threshold.rho)


# Same process for the normal network

filtr.expr.n.deg <- log2(filtr.expr.n.deg+1) 
# Correlation for the normal network
cor.mat.n <- corr.test(t(filtr.expr.n.deg), use = "pairwise", 
                       method = "spearman", adjust="fdr", ci=FALSE)

rho.n <- cor.mat.n$r
diag(rho.n) <- 0

adj.mat.n <- 1 * (abs(rho.n) > threshold.rho)


## 6.6 - Create the networks ----

# Create the normal network given the adjacency matrix
net.n <- as.network(adj.mat.n, directed = FALSE)

# Let's check the number of edges and nodes
network.size(net.n) 
network.edgecount(net.n) 

# and the density of the network
network.density(net.n) 

# Same for the cancer network
net.c <- as.network(adj.mat.c, directed = FALSE)

network.size(net.c) 
network.edgecount(net.c) 

network.density(net.c)


## 6.7 - Scale-free -----

# To check if a network is scale free we need to look at the degree distribution
# The network will be scale free if there are a few hubs and the other nodes have 
# just few connections

# First for the normal graph
# extract the degree of each node
d.n <- sna::degree(net.n, gmode = 'graph')
names(d.n) <- network.vertex.names(net.n)

# Print the histogram of the degree together with a line for the 95% quantile
x_n <- quantile(d.n[d.n>0],0.95)
hist(d.n,col = "lightblue", main = "Degree distribution - normal tissue")
abline(v=x_n, col="red")
# Not scale free


# The number of hubs the network has is
hubs.n <- d.n[d.n>=x_n]
length(hubs.n) # 26

# Power law: Another way to check if a network is scale free is to verify if it
# follows the power law. i.e. there should be a link between the fraction of 
# nodes f(k) with degree k and k itself f(k) ⁓ k^(-γ) (usually 2< γ <3)
d.n.table <- table(d.n)

# Convert the table to a data frame
d.n.fd <- data.frame(degree = as.numeric(names(d.n.table)),
                     count = as.numeric(d.n.table)/length(hubs.n))

plot(log(d.n.fd$degree), log(d.n.fd$count), 
     xlab = "Degree", ylab = "Degree Count", 
     main = "Degree Distribution", 
     pch = 16, col = "lightblue")
# Not scale free



# Same for the cancer graph
d.c <- sna::degree(net.c, gmode = 'graph')
names(d.c) <- network.vertex.names(net.c)

x_c <- quantile(d.c[d.c>0],0.95)
hist(d.c,col = "lightblue", main = "Degree distribution - cancer tissue - Spearman correlation")
abline(v = x_c, col="red")
# It looks better. More scale-free than before

# hubs
hubs.c <- d.c[d.c>=x_c]
length(hubs.c) # 17 ---> fewer hubs than before

# power law
d.c.table <- table(d.c)

# Convert the table to a data frame
d.c.fd <- data.frame(degree = as.numeric(names(d.c.table)),
                     count = as.numeric(d.c.table)/length(hubs.c))

plot(log(d.c.fd$degree), log(d.c.fd$count), 
     xlab = "Degree", ylab = "Degree Count", 
     main = "Degree Distribution - Spearman", 
     pch = 16, col = "lightblue")


# Scale free


## 6.8 - Compare the hubs sets -----

common <- intersect(names(hubs.c),names(hubs.n))
length(common)

# let's identify those selectively characterizing the network

not_in_common_c <- names(hubs.c[!(names(hubs.c) %in% common)])

not_in_common_n <- names(hubs.n[!(names(hubs.n) %in% common)])




hubs.c
hubs.c.ids <- vector("integer",length(hubs.c))
for (i in 1:length(hubs.c)){hubs.c.ids[i] <- match(names(hubs.c)[i],rownames(adj.mat.c))}
hubs.c.ids

#identifying the neighborhood
hubs.c.neigh <- c()
for (f in hubs.c.ids){
  hubs.c.neigh <- append(hubs.c.neigh, get.neighborhood(net.c, f))
}

hubs.c.neigh <- unique(hubs.c.neigh)
hubs.c.neigh
hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh,])
subnet <- unique(c(names(hubs.c), hubs.c.neigh.names))

#creating the subnetwork
hub.c.adj <- adj.mat.c[subnet, subnet]

names.hubs <-names(hubs.c)

rownames(hub.c.adj)[1:length(hubs.c)] <- names.hubs
colnames(hub.c.adj)[1:length(hubs.c)] <- names.hubs


net.hub <- network(hub.c.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
network.density(net.hub)

sum(hub.c.adj > 0 )

net.hub %v% "type" = ifelse(network.vertex.names(net.hub) %in% names.hubs,"hub", "non-hub")
net.hub %v% "color" = ifelse(net.hub %v% "type" == "non-hub", "deepskyblue3", "tomato")

ggnet2(net.hub,  color = "color",alpha = 0.9, size = 2, 
       edge.color = "grey", edge.alpha = 0.15,  edge.size = 0.15)+
  guides(size = "none")


# fewer hubs (the spearman correlation is lower than pearson's)


## 6.9 - Differential co-expression ------


z.c <- 0.5*log((1 + rho.c)/(1- rho.c))
z.n <- 0.5*log((1 + rho.n)/(1- rho.n))


den <- sqrt(2/(18 - 3))
Z = (z.c - z.n)/den

adj.mat.diff <- 1 * (abs(Z) > 3)


net.diff <- as.network(adj.mat.diff, directed = FALSE)

#checking to see if it is a scale free network --> look at degree distributions. in 
#scale free networks There is a small number of highly connected nodes, called hubs (tail of the distribution)
#Most nodes have few connections

d.diff <- sna::degree(net.diff, gmode = 'graph')
names(d.diff) <- network.vertex.names(net.diff)


x_n <- quantile(d.diff[d.diff>0],0.95)
hist(d.diff,col = "lightblue", main = "Degree distribution - differential network - Spearman", xlab = 'Degree')
abline(v=x_n, col="red")
text(x = x_n, y = par("usr")[4]*0.95, labels = "0.95 Quantile", pos = 2, offset = 1, col = "red", cex = 1)


hubs.diff <- d.diff[d.diff>=x_n]

d.diff.table <- table(d.diff)

# Convert the table to a data frame
d.diff.fd <- data.frame(degree = as.numeric(names(d.diff.table)),
                        count = as.numeric(d.diff.table)/length(hubs.diff))


x <- log(d.diff.fd$degree)
y <- log(d.diff.fd$count)

model <- glm.fit(x[2:length(x)], y[2:length(y)])
model$coefficients

plot(x[2:length(x)], model$fitted.values, type='l', ylim = c(-3, 1), xlab = 'logarithm of degree', ylab = 'logarith of count', col = 'red', 
     main = 'Power Law - Differential Network - Spearman')
points(log(d.diff.fd$degree), log(d.diff.fd$count), 
       xlab = "Degree", ylab = "Degree Count", 
       main = "Degree Distribution", 
       pch = 16, col = "lightblue")

text(x = 1.5, y = -1, labels = "slope = -0.61", pos = 2, offset = 1, col = "red", cex = 1)


# Scale Free


hubs.diff.ids <- vector("integer",length(hubs.diff))
for (i in 1:length(hubs.diff)){hubs.diff.ids[i] <- match(names(hubs.diff)[i],rownames(adj.mat.diff))}
hubs.diff.ids

#identifying the neighborhood
hubs.diff.neigh <- c()
for (f in hubs.diff.ids){
  hubs.diff.neigh <- append(hubs.diff.neigh, get.neighborhood(net.diff, f))
}

hubs.diff.neigh <- unique(hubs.diff.neigh)
hubs.diff.neigh.names <- rownames(adj.mat.diff[hubs.diff.neigh,])
subnet.diff <- unique(c(names(hubs.diff), hubs.diff.neigh.names))

#creating the subnetwork
hub.diff.adj <- adj.mat.diff[subnet.diff, subnet.diff]

names.hubs.diff <-names(hubs.diff)

rownames(hub.diff.adj)[1:length(hubs.diff)] <- names.hubs.diff
colnames(hub.diff.adj)[1:length(hubs.diff)] <- names.hubs.diff

net.hub.diff <- network(hub.diff.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

node_degrees <- degree(net.hub.diff)


net.hub.diff %v% "type" = ifelse(network.vertex.names(net.hub.diff) %in% names.hubs.diff,"hub", "non-hub")
net.hub.diff %v% "color" = ifelse(net.hub.diff %v% "type" == "non-hub", "deepskyblue3", "tomato")
net.hub.diff %v% "size" = ifelse(net.hub.diff %v% "type" == "non-hub", 1, 2)

ggnet2(net.hub.diff,  color = "color",alpha = 0.9, size = 2, 
       edge.color = "grey", edge.alpha = 0.5,  edge.size = 0.15)+guides(size = "none")


# Compare the hubs obtained with the differential expressed network to those of the previous points

intersect(names(hubs.c), names(hubs.diff))
intersect(names(hubs.n), names(hubs.diff))


# COMMUNITIES

# Correlation for the cancer network
cor.mat.cp <- corr.test(filtr.expr.c, use = "pairwise", 
                        method = "spearman", adjust="fdr", ci=FALSE)

# matrix with the pearson correlation coefficient
rho.cp <- cor.mat.cp$r
# Set the diagonal to zero (no self loops)
diag(rho.cp) <- 0
# Binary matrix, keep the edge for big correlation


# correlation network of cancer samples

adj.mat.cp <- rho.cp * (abs(rho.cp) > 0.85)   #we have to change the theshold 

net.cp <- as.network(adj.mat.cp, directed = FALSE)


l.comp.p <- component.largest(net.cp, result = "graph")
l.comp.p <- adj.mat.cp[rownames(l.comp.p), rownames(l.comp.p)]
#let's use the nodes name to index the weighted matrix

write.csv2(l.comp.p, "input-matrix-c.csv")

# in the terminal 
# python btc-community-c.py input-matrix-c.csv


comm.res.p <- read.csv2("output-c.txt", header = FALSE)
rownames(comm.res.p) <- rownames(l.comp.p)
sort(table(comm.res.p[,1]), decreasing = T)
length(table(comm.res.p[,1]))

net.final.p <- network(l.comp.p, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)

all(net.final.p %v% "vertex.names" == rownames(comm.res.p)) 
net.final.p  %v% "community" <-  as.character(comm.res.p[,1])

ncom <- length(unique(net.final.p  %v% "community"))

set.seed(123)
pal <- sample(colors(distinct = T), ncom)
names(pal) <- 1:ncom
pal

node_mapping <- data.frame(real_label = unique(net.final.p %v% "vertex.names"),
                           new_label = 1:length(unique(net.final.p %v% "vertex.names")))

# Update the network with new labels
net.final.p %v% "vertex.names" <- node_mapping$new_label


ggnet2(net.final.p, color = "community", palette =  pal, alpha = 1, 
       size = 5, edge.color = "grey", edge.alpha = 1, edge.size = 0.15, label = TRUE, label.size = 5)+
  guides(size = "none") 




# Correlation for the normal network
cor.mat.np <- corr.test(filtr.expr.n.deg, use = "pairwise", 
                        method = "spearman", adjust="fdr", ci=FALSE)

# matrix with the pearson correlation coefficient
rho.np <- cor.mat.np$r
# Set the diagonal to zero (no self loops)
diag(rho.np) <- 0
# Binary matrix, keep the edge for big correlation


adj.mat.np <- rho.np * (abs(rho.np) > 0.85)


net.np <- as.network(adj.mat.np, directed = FALSE)


l.comp.pn <- component.largest(net.np, result = "graph") 
l.comp.pn <- adj.mat.np[rownames(l.comp.pn), rownames(l.comp.pn)]


write.csv2(l.comp.pn, "input-matrix-n.csv")

# in the terminal 
# python btc-community-n.py input-matrix-n.csv

comm.res.pn <- read.csv2("output-n.txt", header = FALSE)
rownames(comm.res.pn) <- rownames(l.comp.pn)
sort(table(comm.res.pn[,1]), decreasing = T)
length(table(comm.res.pn[,1]))

net.final.pn <- network(l.comp.pn, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)

all(net.final.pn %v% "vertex.names" == rownames(comm.res.pn)) 
net.final.pn  %v% "community" <-  as.character(comm.res.pn[,1])

ncom_n <- length(unique(net.final.pn  %v% "community"))

set.seed(123)
pal <- sample(colors(distinct = T), ncom_n)
names(pal) <- 1:ncom
pal

vertex_names <- as.character(net.final.pn %v% "vertex.names")
vertex_data <- data.frame(name = vertex_names, x = layout_nicely(net.final.pn)[,1], y = layout_nicely(net.final.pn)[,2])

net.final.pn %v% "vertex.names" <- node_mapping$new_label


ggnet2(net.final.pn, color = "community", palette =  pal, alpha = 1,
       size = 5, edge.color = "grey", edge.alpha = 1, edge.size = 0.15, label = TRUE, label.size = 5)+
  guides(size = "none")




