library(Seurat)
library(Ternary)
library(viridis)
library(ggpubr)

# LOAD pre-clustered all neuron data object

# Navigate to directory containing the R object THSN_Seurat_object.Rdata
load("AllNeuron_Seurat_object.Rdata") # 46,073 droplet dataset, pre-clustered
MyNeurons <- neuL_60d_110n_ng_vf13_CR50
total_droplets<-ncol(MyNeurons) #46,073 droplets

# Fig. 1D and S5: Visualize UMAPs for clusters, batches, and key marker genes

# Note that the AaegL5 gtf annotation used for the CellRanger alignments underlying these data uses
# "OrX", "IrX", "GrX" names for chemosensory receptors. Other genes are named using the associated
# AAEL number followed by the name of the Drosophila ortholog (*where orthologs include multiple
# mosquito genes that have the same best match in Drosophila). For example, "AAEL019818" was changed
# to "AAEL019818-nompC". See Methods

DimPlot(MyNeurons, pt.size=0.8, label=T) # Fig. 1D
DimPlot(MyNeurons, group.by="orig.ident", pt.size=0.8)
FeaturePlot(MyNeurons, feature="Orco", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(MyNeurons, feature="Ir25a", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(MyNeurons, feature="Ir8a", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(MyNeurons, feature="Ir76b", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(MyNeurons, feature="Ir93a", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
genes <- rownames(MyNeurons)
FeaturePlot(MyNeurons, features = genes[grep('nompC',genes)], pt.size=0.6, cols=rev(viridis(100, option='D')), order=T) #UMAP for gene of interest with "XXX" somewhere in name

# ----- Subsetting olfactory neurons only (eliminating thermo/hygrosensory, mechanosensory neurons and junk droplet clusters)
MyOlfactoryNeurons <- subset(MyNeurons, idents= c("41","39","29","54", "48", "45", "55", "1", "47"), invert = T)

# ----- Analyzing coreceptor expression in all OSN clusters
avgs <- AverageExpression(MyOlfactoryNeurons, slot='counts')$SCT # raw average expression
clusters <- as.vector(colnames(avgs))
co4 = avgs[c('Orco','Ir25a','Ir76b','Ir8a'),]
co4 = co4/apply(co4,1,max) #normalize values for each coreceptor by max value

# Fig. 1E. Histograms of coreceptor % maximal expression
par(mfrow=c(2,2))
hist(co4["Orco",], breaks=20, ylim=c(0,40))
hist(co4["Ir25a",], breaks=20, ylim=c(0,40))
hist(co4["Ir8a",], breaks=20, ylim=c(0,40))
hist(co4["Ir76b",], breaks=20, ylim=c(0,40))

# List of clusters in which key co-receptors called as expressed
orco <- colnames(co4)[co4["Orco",]>0.1]
Ir25a <- colnames(co4)[co4["Ir25a",]>0.1]
Ir8a <- colnames(co4)[co4["Ir8a",]>0.1]
Ir76b <-colnames(co4)[co4["Ir76b",]>0.1]

# Fig. 1F: Ternary plot of coreceptor % maximal expression
par(mfrow=c(1,1))
co3 <- co4[c('Orco','Ir76b','Ir8a'),]
co3<-co3[,apply(co3,2,max)>0.1] # Remove one cluster that does not express any of the 3 co-receptors (Amt cluster)
propss<-t(co3)
rownames(propss)
TernaryPlot(atip='Orco',btip='Ir76b',ctip='Ir8a',
            grid.lines = 4,
            grid.minor.lines = 0,
            grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), 
            axis.col = "black", ticks.col = "black",
            axis.rotate = FALSE)
AddToTernary(points,as.list(as.data.frame(t(propss)[c('Orco','Ir76b','Ir8a'),])), pch = 21, cex = 3.8)
AddToTernary(text,as.list(as.data.frame(t(propss)[c('Orco','Ir76b','Ir8a'),])), labels=rownames(propss),cex=1)

# ----- Analyzing expression of maximally expressed OrX or IrX in each OSN cluster
#avgs <- AverageExpression(MyOlfactoryNeurons, slot='counts')$SCT # raw average expression
avgs <- AverageExpression(MyOlfactoryNeurons, return.seurat=T)$SCT@data # log scaled average expression
clusters <- as.vector(colnames(avgs))

OrXlist <- rownames(avgs)[grep('^Or',rownames(avgs))]
OrXlist <- OrXlist[!(OrXlist %in% c("Orco"))]
length(OrXlist) # 111 ligand-specific ORs in the OSN dataset
OrXs <- avgs[OrXlist,]

IrXlist <- rownames(avgs)[grep('^Ir',rownames(avgs))]
IrXlist <- IrXlist[!(IrXlist %in% c("Ir25a","Ir8a","Ir76b","Ir93a"))]
length(IrXlist) # 90 ligand-specific IRs in the OSN dataset
IrXs <- avgs[IrXlist,]

maxOrX <- apply(OrXs,2,max) #list of max OrX values for each cluster
maxIrX <- apply(IrXs,2,max) #list of max IrX values for each cluster

# Order the clusters by decreasing maxOrX, down to a value of ~0.2 and then by increasing maxIrX...
OrIrXorder <- order(maxOrX,decreasing=T)[maxOrX[order(maxOrX,decreasing=T)]>0.2]
OrIrXorder <- c(OrIrXorder, order(maxIrX)[! order(maxIrX) %in% OrIrXorder])
OrIrXclusters <- clusters[OrIrXorder]

# Color clusters by which coreceptors they express
OrIrXcolors <- rep("black",length(OrIrXclusters))
OrIrXcolors[OrIrXclusters %in% orco] = "pink"
OrIrXcolors[OrIrXclusters %in% Ir8a] = "blue"
OrIrXcolors[OrIrXclusters %in% Ir76b] = "light blue"
OrIrXcolors[OrIrXclusters %in% c("46")] = "purple" #orco+Ir76b
OrIrXcolors[OrIrXclusters %in% c("42","52")] = "orange" #orco+Ir76b+Ir8a
OrIrXcolors[OrIrXclusters %in% c("35")] = "cyan" #Ir8a+Ir76b

# Fig. 1G: Plot highest OrX and highest IrX in each cluster
par(mfrow=c(2,2))
max = ceiling(max(c(maxOrX,maxIrX)))
barplot(maxOrX[OrIrXorder], ylim=c(0,max), col=OrIrXcolors)
barplot(maxIrX[OrIrXorder], ylim=c(0,max), col=OrIrXcolors) 
b=1
barplot(maxOrX[OrIrXorder], ylim=c(0,b), col=OrIrXcolors)
barplot(maxIrX[OrIrXorder], ylim=c(0,b), col=OrIrXcolors) 
abline(h=0.05)

