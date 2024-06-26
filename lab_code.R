# 2/11/23

BiocManager::install("BiocGenerics")
BiocManager::install("DESeq2")
install.packages("GGally")

required_packages <- c("BiocGenerics", "DESeq2", "psych", "NetworkToolbox", "ggplot2",
                       "GGally", "sna", "network", "TCGAbiolinks", "SummarizedExperiment", "DT")

lapply(required_packages, library, character.only = TRUE)

#1: Downloading data from the TGCA -------
proj <- "TCGA-CHOL"
dir.create(file.path(proj), recursive = TRUE)


rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")

GDCdownload(query = rna.query.C, directory = "TCGA-CHOL/GDCdata", method = "api")
rna.data.C <- GDCprepare(rna.query.C, directory = "TCGA-CHOL/GDCdata")
rna.expr.data.C <- assay(rna.data.C)

View(BiocGenerics::as.data.frame(rowRanges(rna.data.C)))
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.C))

rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
rna.data.N <- GDCprepare(rna.query.N, directory = "GDCdata"  )
rna.expr.data.N <- assay(rna.data.N)
genes.info2 <- BiocGenerics::as.data.frame(rowRanges(rna.data.N))
all(na.omit(genes.info2) == na.omit(genes.info))

clinical.query<- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
#write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)

View(clinical.query)
table(clinical.query$ajcc_pathologic_stage)

boxplot(age_at_index ~ ajcc_pathologic_stage, data = clinical.query,
        col = "gold", main = "Title", xlab = "", ylab= "age", las=2 )


View(rna.expr.data.N)

dim(rna.expr.data.C)
dim(rna.expr.data.N)
length(unique(clinical.query$submitter_id))



#2: Data cleaning -----

ncol(rna.expr.data.C)
head(colnames(rna.expr.data.C))
head(substr(colnames(rna.expr.data.N), 1,12))

dim(rna.expr.data.N)
length(unique(substr(colnames(rna.expr.data.N), 1,12))) #no duplicates 
dim(rna.expr.data.C)
length(unique(substr(colnames(rna.expr.data.C), 1,12))) #duplicates!

patients.C <- substr(colnames(rna.expr.data.C), 1,12)
sort(table(patients.C)) 

unique.patients.C <- names(which(table(patients.C) == 1))
#let's get their index in the list of patients
idx.unique.pats <- match(unique.patients.C, substr(colnames(rna.expr.data.C), 1,12) )

#I'm only keeping patients with one sample
expr.C <- as.data.frame(rna.expr.data.C[,idx.unique.pats])
expr.N <- as.data.frame(rna.expr.data.N)

#let's rename patients in a shorter way
colnames(expr.C) <- substr(colnames(expr.C), 1,12)
unique(colnames(expr.C))

colnames(expr.N) <- substr(colnames(expr.N), 1,12)
unique(colnames(expr.N))

intersect(colnames(expr.N), colnames(expr.C)) #38 instead of 41
setdiff(colnames(expr.N), colnames(expr.C))
# 3 normal sample do not have a cancerous sample to compare to. Let's remove them
match(setdiff(colnames(expr.N), colnames(expr.C)), colnames(expr.N)) #idx to remove
expr.N <- expr.N[,-c(10,15,20)]

length(intersect(colnames(expr.N), colnames(expr.C)))
#normal samples are a subset of the cancer samples. great. 


#let's check the actual counts
typeof(expr.C[1,1]) #ok
any(is.na(expr.C)) #ok
any(is.nan(as.matrix(expr.C))) #ok

typeof(expr.N[1,1]) #ok
any(is.na(expr.N)) #ok
any(is.nan(as.matrix(expr.N))) #ok

#let's consider only patients for which we have both normal and cancer samples
expr.C <- expr.C[, colnames(expr.N)]

#let's assume we only want female patients
women <- (clinical.query[clinical.query$gender == "female", ])$submitter_id
head(women)
tail(women)
women <- women[1:216]

length(intersect(colnames(expr.C), women)) #19

expr.C <- expr.C[, intersect(colnames(expr.C), women)]
expr.N <- expr.N[, intersect(colnames(expr.C), women)]






#3: Normalizing data with Deseq2 ----- 

#detalied explanation: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#youtube explanation: https://www.youtube.com/watch?v=UFB993xufUU

all(rownames(expr.C) == rownames(expr.N))
full.data <- cbind(expr.N, expr.C)

dim(full.data)
full.data <- data.frame(full.data)

metad <- rep("cancer", 38)
metad[1:19] <- "normal"
metad
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1])
full.data <- cbind(rownames(full.data), full.data)

dds <- DESeqDataSetFromMatrix(countData=full.data, 
                              colData=metad, 
                              design= ~condition,
                              tidy=TRUE)

View(counts(dds))
dim(counts(dds))

# filtering: at least ten counts on 90% of patients? 
( 38*90 )/100
keep <- rowSums(counts(dds) >= 10) >= 34
dds <- dds[keep,]
dim(counts(dds))

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == 38) #no null rows


filtr.expr.n <- as.data.frame(normalized_counts[, 1:19])
filtr.expr.c <- as.data.frame(normalized_counts[, 20:38])
#cancerous sample names were added a ".1" in full.data because  
#they had the same names as the normal samples
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1,12)

#4: Gene selection ----

genes <- read.csv2("targets.csv", row.names = 1)
genes <- genes[,1]
head(genes) #gene symbols, not ensemble ids

#their enseble ids 
genes.info[ genes.info$gene_name %in% genes ,"gene_id"]
#not all of them might be among our genes!

genes.c <- intersect(rownames(filtr.expr.c), 
                     genes.info[ genes.info$gene_name %in% genes , "gene_id"]   ) 

genes.n <- intersect(rownames(filtr.expr.n),  
                     genes.info2[ genes.info2$gene_name %in% genes , "gene_id"]   )  

setdiff(genes.c, genes.n)

length(genes)
length(genes.c)
length(genes.n)

filtr.expr.n <- filtr.expr.n[genes.n, ]
filtr.expr.c <- filtr.expr.c[genes.c, ]

rownames(filtr.expr.n) <- genes.info2[genes.n, "gene_name"]
rownames(filtr.expr.c) <- genes.info[genes.c, "gene_name"]


#5: Differentially expressed genes (DEGs) (+ brief enrichment overview) ----- 

#what are DEGs?
fc <-  log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n) ) 
names(fc) <- rownames(filtr.expr.c)
head(fc)
#what is the fold change?

pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(filtr.expr.c[i,], filtr.expr.n[i,] ))$p.value)
pval.fc.fdr <- p.adjust(pval.fc, method="fdr")

expr.table <- data.frame(cbind(fc, pval.fc.fdr))
expr.table[,1] <- round(expr.table[,1],2)
#expr.table[,2] <- format(expr.table[,2], digits = 4) nice for printing but it converts to string

deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.5 & expr.table$pval.fc.fdr <=0.01,]) 
deg.genes

head(expr.table[deg.genes,], 10)
#write.table(expr.table[deg.genes,], file = "DEG.csv", sep = ";")

#volcano plot

expr.table$diffexpressed <- "NO";
expr.table$diffexpressed[expr.table$fc >= 1.5 & expr.table$pval.fc.fdr <= 0.01] <- "UP"
expr.table$diffexpressed[expr.table$fc <= -1.5 & expr.table$pval.fc.fdr <= 0.01] <- "DOWN"
head(expr.table)

expr.table$diffexpressed <- as.factor(expr.table$diffexpressed)
summary(expr.table$diffexpressed)

ggplot(data=expr.table, aes(x=fc, y=-log10(pval.fc.fdr), col=diffexpressed))+  
  geom_point() +
  xlab("fold change (log2)") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  geom_vline(xintercept=1.5, col="red")+
  geom_vline(xintercept=-1.5, col="red")


# print and enrichment 
cat(deg.genes , sep = "\n")








#6: Adjacency matrices of co-expression networks -----

#cancer network
cor.mat.c <- corr.test(t(filtr.expr.c), use = "pairwise", 
                       method = "spearman",adjust="fdr",ci=FALSE)

rho.c <- cor.mat.c$r
diag(rho.c) <- 0
qval.c <- cor.mat.c$p
#qvals are reported on the upper triangle only
qval.c[lower.tri(qval.c)] <- t(qval.c)[lower.tri(qval.c)]

adj.mat.c <- rho.c * (qval.c <= 0.05)

#normal network 
cor.mat.n <- corr.test(t(filtr.expr.n), use = "pairwise", 
                       method = "spearman",adjust="fdr",ci=FALSE)

rho.n <- cor.mat.n$r
diag(rho.n) <- 0
qval.n <- cor.mat.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]

adj.mat.n <- rho.n * (qval.n <= 0.05)

#7: Co-expression networks ----

#Cancer network 

net.c <- network(adj.mat.c, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.c)
network.size(net.c)
network.edgecount(net.c) 
#nrow(component.largest(net.c, result = "graph")) #1549
clustcoeff(adj.mat.c, weighted = FALSE)$CC

sum(adj.mat.c != 0)
#how many positive/negative correlations? 
sum(adj.mat.c > 0) 
sum(adj.mat.c < 0) 

degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c, decreasing = T)
head(degree.c,10)
sum(degree.c == 0) #unconnected nodes 

hist(degree.c)
x <- quantile(degree.c[degree.c>0],0.95) #how big is the degree of the most connected nodes?
x
hist(degree.c)
abline(v=x, col="red")

hubs.c <- degree.c[degree.c>=x]
names(hubs.c) #

net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c),"hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, "red", "blue"))

#coord.c <- gplot.layout.fruchtermanreingold(net.c, NULL)
#net.c %v% "x" = coord.c[, 1]
#net.c %v% "y" = coord.c[, 2]

ggnet2(net.c, color = "color", alpha = 0.7, size = 2,  #mode= c("x","y"),
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

#this is extremely dense... what if we change the pval threshold?
adj.mat.c <- rho.c * (qval.c <= 0.01)
#what if it's even lower?
adj.mat.c <- rho.c * (qval.c <= 1e-3)
adj.mat.c <- rho.c * (qval.c <= 1e-4) #too much?


#might be useful to look at negative and postive edges separately:
net.c1 <- network(adj.mat.c* (adj.mat.c > 0),
                  matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

ggnet2(net.c1, color = "deepskyblue3", alpha = 0.7, size = 2, 
       edge.color = "red", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 


net.c2 <- network(adj.mat.c* (adj.mat.c < 0),
                  matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

ggnet2(net.c2, color = "deepskyblue3", alpha = 0.7, size = 2, 
       edge.color = "blue", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

#for the negative one is too much!
adj.mat.c <- rho.c * (qval.c <= 5e-4)


#Normal network 

net.n <- network(adj.mat.n, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.n)
network.size(net.n)
network.edgecount(net.n)
clustcoeff(adj.mat.n, weighted = FALSE)$CC
#nrow(component.largest(net.n, result = "graph")) #1549

sum(adj.mat.n != 0)
#how many positive/negative correlations? 
sum(adj.mat.n > 0) 
sum(adj.mat.n < 0) 

degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = T)
head(degree.n,10)
sum(degree.n == 0) #unconnected nodes 

hist(degree.n)
y <- quantile(degree.n[degree.n>0],0.95) #how big is the degree of the most connected nodes?
y
hist(degree.n)
abline(v=y, col="red")

hubs.n <- degree.n[degree.n>=y]
names(hubs.n) #

net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n),"hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))

ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 

#this is extremely dense... what if we change the pval threshold?
adj.mat.n <- rho.c * (qval.c <= 1e-3) 
adj.mat.n <- rho.c * (qval.c <= 1e-4) #too much?


intersect(names(hubs.c), names(hubs.n))


#8: Plotting the hub subnetwork -----

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
head(rownames(hub.c.adj))
head(colnames(hub.c.adj))

net.hub <- network(hub.c.adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")
network.density(net.hub)

sum(hub.c.adj > 0 )
sum(hub.c.adj < 0)

net.hub %v% "type" = ifelse(network.vertex.names(net.hub) %in% names.hubs,"hub", "non-hub")
net.hub %v% "color" = ifelse(net.hub %v% "type" == "non-hub", "deepskyblue3", "tomato")
set.edge.attribute(net.hub, "ecolor", ifelse(net.hub %e% "weights" > 0, "red", "blue"))

ggnet2(net.hub,  color = "color",alpha = 0.9, size = 2, 
       edge.color = "ecolor", edge.alpha = 0.9,  edge.size = 0.15, 
       node.label = names.hubs, label.color = "black", label.size = 4)+
  guides(size = "none") 



