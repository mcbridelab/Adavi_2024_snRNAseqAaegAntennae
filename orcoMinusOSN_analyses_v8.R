library(Seurat)
library(viridis)
library(plotrix)
library(scCustomize)
library(scales)

# LOAD pre-clustered orco- OSN data object

# Navigate to directory containing the R object THSN_Seurat_object.Rdata
load("orcoMinusOSN_Seurat_object.Rdata")
MyOrcoMinusOSNs <- IR_Seurat
total_droplets<-ncol(MyOrcoMinusOSNs) # 12,243 droplet dataset
clustorder=c("0","7","13","2","10a","10b","1","17","6","12","15","3","16","11","5","4","9","14a","14b","8")

# Fig. 3A and S11: Visualize UMAPs for clusters, batches, and key marker genes

# Note that the AaegL5 gtf annotation used for the CellRanger alignments underlying these data uses
# "OrX", "IrX", "GrX" names for chemosensory receptors. Other genes are named using the associated
# AAEL number followed by the name of the Drosophila ortholog (*where orthologs include multiple
# mosquito genes that have the same best match in Drosophila). For example, "AAEL019818" was changed
# to "AAEL019818-nompC". See Methods

DimPlot(MyOrcoMinusOSNs, pt.size=0.5, label=T) # Fig. 3A
DimPlot(MyOrcoMinusOSNs, group.by="orig.ident", pt.size=0.5)
FeaturePlot(MyOrcoMinusOSNs, feature="Ir25a", pt.size=0.5, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(MyOrcoMinusOSNs, feature="Ir8a", pt.size=0.5, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(MyOrcoMinusOSNs, feature="Ir76b", pt.size=0.5, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(MyOrcoMinusOSNs, feature="Ir64a", pt.size=0.5, cols=rev(viridis(100, option='D')), order=T)
genes <- rownames(MyOrcoMinusOSNs)
FeaturePlot(MyOrcoMinusOSNs, features = genes[grep('Amt',genes)], pt.size=0.6, cols=rev(viridis(100, option='D')), order=T) #UMAP for gene of interest with "XXX" somewhere in name


# Fig 3B and S12A-B: Choosing and using log2FC' cutoff to call receptor expression

#------- OPTIONAL (used in Fig. 3B): Merging 4/9 and 14a/14b due to expression of same/similar complement of receptors
MyOrcoMinusOSNs<-RenameIdents(MyOrcoMinusOSNs,"4"="4_9", "9"="4_9", "14a"="14a_14b", "14b"="14a_14b")
clustorder=c("0","7","13","2","10a","10b","1","17","6","12","15","3","16","11","5","4_9","14a_14b","8")

#-------- Assembling list of tuning receptors detected in at least one droplet
receptorlist = rownames(MyOrcoMinusOSNs$SCT@data)[grep('^[OIG]r',rownames(MyOrcoMinusOSNs$SCT@data))] # Add all expressed ORs, GRs, IRs to list
receptorlist <- c(receptorlist, rownames(MyOrcoMinusOSNs$SCT@data)[grep('Amt',rownames(MyOrcoMinusOSNs$SCT@data))]) # Add Amt to list
receptorlist <- receptorlist[!rowSums(MyOrcoMinusOSNs$SCT@data[receptorlist,])==0] # Remove receptors that were not detected in any orco- OSN droplet

#-------- Calculating Log2FC' for each receptor gene in each cluster
MyMarkers <- FindAllMarkers(MyOrcoMinusOSNs, features=receptorlist, logfc.threshold=0, return.thresh=1, min.pct=0, test.use='t')
log2FC <- matrix(NA, nrow=length(receptorlist), ncol=length(levels(MyMarkers$cluster))) # create empty matrix
colnames(log2FC) <- levels(MyMarkers$cluster) # columns are clusters
rownames(log2FC) <- receptorlist # rows are receptors
for(i in 1:nrow(MyMarkers)){ # populate matrix with ave_log2FC values
  log2FC[MyMarkers$gene[i],as.character(MyMarkers$cluster[i])]<-MyMarkers$avg_log2FC[i]
}
log2FC[is.na(log2FC)] <- 0 #replacing NAs with 0s 
med.att <- apply(log2FC, 1, median) #calculating median Log2FC values (for each gene across clusters)
log2FCprime <- sweep(log2FC, 1, med.att) #subtract median Log2FC value from cluster-specific Log2FCs to get Log2FC'

#---------- Fig. S12A: Plotting log2FC` Histogram (includes values for each Ir/Gr/Or in each cluster)
par(mfrow=c(2,1))
hist(log2FCprime, xlim=c(-2,7), breaks=seq(-2,7,0.1))
abline(v=0.4, col='red') # add threshold value used for calling expression (see below)
hist(log2FCprime, xlim=c(-2,7), ylim=c(0,20), breaks=seq(-2,7,0.1))
abline(v=0.4, col='red') # add threshold value used for calling expression (see below)

#---------- Fig. 3B and S12B: Plotting heatmap of log2FC` values (with dots marking those expressed above cutoff of 0.4)
cutoffthresh = 0.4
displaythresh = cutoffthresh
displaythresh = 0.15 # used for Fig. S12B-C
abovedisplaythresh <- rownames(log2FCprime)[apply(log2FCprime,1,max)>displaythresh]
abovedisplaythresh <- abovedisplaythresh[!(abovedisplaythresh %in% c("Orco","Ir25a","Ir76b","Ir8a","Ir93a"))] # exclude co-receptors from log2FC' heatmap

log2FCP_Heatmap_tuning <-log2FCprime[abovedisplaythresh,clustorder]
receptororder <- abovedisplaythresh[order(apply(log2FCP_Heatmap_tuning,1,which.max))] # order receptors by  maximal expression in each sequential cluster
receptororder <- receptororder[c(grep('^[O]r',receptororder),grep('^[G]r',receptororder),grep('^[I]r',receptororder))] # pull Grs and Irs to front
log2FCP_Heatmap_tuning <-log2FCP_Heatmap_tuning[receptororder,]
log2FCP_Heatmap_tuning[log2FCP_Heatmap_tuning<0]<-0 #converting negative log2FC` values to zero
log2FCP_Heatmap <- log2FCP_Heatmap_tuning
color_palette=c(colorRampPalette(c(alpha("white", alpha = 0),'#FDE333'))(10),hcl.colors(450, palette = 'viridis', alpha = NULL, rev = T, fixup = TRUE))
par(mfrow=c(1,1))
image(t(log2FCP_Heatmap),col=color_palette,axes=F,zlim=c(0,round(max(log2FCP_Heatmap))), main='Log2FC`')
axis(3,at=seq(0,1,length=ncol(log2FCP_Heatmap)),las=2,labels=colnames(log2FCP_Heatmap),cex.axis=0.5,lwd=0)
axis(2,at=seq(0,1,length=nrow(log2FCP_Heatmap)),las=2,labels=rownames(log2FCP_Heatmap),cex.axis=0.5,lwd=0)
color.legend(0.05,-0.15,0.2,-0.17,rect.col=color_palette,legend=c(0,round(max(log2FCP_Heatmap))),align='rb')
for(i in 1:length(clustorder)){
  for(j in 1:length(receptororder)){
    if(log2FCP_Heatmap[receptororder[j],clustorder[i]]>cutoffthresh) text(seq(0,1,length=length(clustorder))[i],
                                                                          seq(0,1,length=nrow(log2FCP_Heatmap))[j],labels='.',cex=1, col = "black")
  }
}

# Fig S12C: Visualizing alternative log average expression cutoff for receptors

AveExp <- AverageExpression(MyOrcoMinusOSNs, features=receptorlist, return.seurat=T)
AveExp_tuning <- AveExp$SCT@data[receptororder,clustorder] # extracting log normalized average counts for same tuning receptors that met the log2FC display threshold used above
AveExp_co <- AveExp$SCT@data[c("Ir25a","Ir76b","Ir8a"),clustorder] # extracting log normalized average counts for co-receptors

#---------- Explore histogram of log average expression values for each receptor in each cluster to choose cutoff
par(mfrow=c(2,1))
hist(AveExp_tuning, breaks=seq(0,5,0.01), main='logAveExpression')
hist(AveExp_tuning, ylim=c(0,50), breaks=seq(0,5,0.01), xlim=c(0,1),main='logAveExpression')
abline(v=0.15)

#---------- Plotting heatmap of log average expression values (including co-receptors)
AveExpcutoffthresh = 0.15
AveExp_Heatmap <- rbind(AveExp_co,AveExp_tuning)
Corcols=hcl.colors(450, palette = 'rocket', alpha = NULL, rev = T, fixup = TRUE)
par(mfrow=c(1,1))
image(t(AveExp_Heatmap),col=Corcols,axes=F, zlim=c(0,ceiling(max(AveExp_Heatmap))), main='Log average expression')
axis(3,at=seq(0,1,length=ncol(AveExp_Heatmap)),las=2,labels=colnames(AveExp_Heatmap),cex.axis=0.5,lwd=0)
axis(2,at=seq(0,1,length=nrow(AveExp_Heatmap)),las=2,labels=rownames(AveExp_Heatmap),cex.axis=0.3,lwd=0)
color.legend(0.05,-0.15,0.2,-0.17,rect.col=Corcols,legend=c(0,ceiling(max(AveExp_Heatmap))),align='rb')
# add dots for receptors that exceed log ave expression threshold
for(i in 1:length(clustorder)){
  for(j in 1:length(receptororder)){
    y<-j+3 # accounts for space taken up by co-receptors at bottom of plot
    if(AveExp_Heatmap[receptororder[j],clustorder[i]]>AveExpcutoffthresh) text(seq(0,1,length=length(clustorder))[i],
                                                                    seq(0,1,length=nrow(AveExp_Heatmap))[y],labels='.',cex=1, col = "black")
  }
}


# Fig. S12D-E: Calculating and plotting pairwise Pearson's correlations

#---------- Refresh list of receptors with log2FC'>0.15 in at least one cluster and get count data for all droplets
myreceptors <- rownames(log2FCprime)[apply(log2FCprime,1,max)>displaythresh]
myreceptors <- myreceptors[!(myreceptors %in% c("Orco","Ir25a","Ir76b","Ir8a","Ir93a"))]
MyData <- subset(MyOrcoMinusOSNs,features=myreceptors)
MyMatrix <- t(as.matrix(MyData@assays$SCT@data))

#---------- Generate matrix of pairwise correlations and reorder the receptors
log2FCP_matrix <- log2FCprime[myreceptors,clustorder]
receptororder <- myreceptors[order(apply(log2FCP_matrix,1,which.max))] # order receptors by  maximal expression in each sequential cluster
receptororder <- receptororder[c(grep('^[G]r',receptororder),grep('^[I]r',receptororder))] # pull Grs and Irs to front

cormat <- cor(MyMatrix, method = "pearson", use = 'everything')
cormat <- cormat[receptororder,receptororder] # reorder receptors

#---------- Generate matrix showing which genes are called as coexpressed in at least one cluster
coexp = matrix(data=F, nrow=length(receptororder), ncol=length(receptororder))
colnames(coexp) <- receptororder
rownames(coexp) <- receptororder
for (i in 1:nrow(coexp)) { 
  for (j in 1:ncol(coexp)) {
    is_coexp = F
    for (k in 1:ncol(log2FCP_Heatmap)) { # check whether the two receptors BOTH have log2FC'>0.3 in any of k clusters
      if (sum(c(log2FCP_Heatmap[receptororder[i],k]>cutoffthresh,log2FCP_Heatmap[receptororder[j],k]>cutoffthresh))==2) { is_coexp = T }
    }
    if (is_coexp) { coexp[i,j] = T }
  }
}

#---------- Fig. S12E: Plot heatmap of pairwise correlations for all receptors with log2FC'>0.15 in at least one cluster
par(mfrow=c(1,1))
image(cormat,col=viridis(256, direction = -1, option = "B"),zlim=c(0,1),axes=F)
axis(1,at=seq(0,1,length=nrow(cormat)),las=2,labels=rownames(cormat),cex.axis=0.5,lwd=0)
axis(2,at=seq(0,1,length=ncol(cormat)),las=2,labels=colnames(cormat),cex.axis=0.5,lwd=0)
color.legend(0.01,0.99,0.1,1.0,rect.col=viridis(256, direction = -1, option = "B"),legend=c(0,1),align='rb')
# Add dots for squares that correspond to receptors called as co-expressed
for(i in 1:length(receptororder)){
  for(j in 1:length(receptororder)){
    if(coexp[receptororder[i],receptororder[j]]) text(seq(0,1,length=length(receptororder))[i],
                                                      seq(0,1,length=length(receptororder))[j],labels='.',cex=1, col = "black")
  }
}

#---------- Limit matrices to 'expressed' IRs (with log2FC'>0.4 in at least one cluster)
myexpressedIRs <- rownames(log2FCprime)[apply(log2FCprime,1,max)>cutoffthresh]
myexpressedIRs <- myexpressedIRs[!(myexpressedIRs %in% c("Orco","Ir25a","Ir76b","Ir8a","Ir93a"))]
myexpressedIRs <- myexpressedIRs[grep('^[I]r',myexpressedIRs)] #limit to IRs only
cormat2 <- cormat
cormat2 <- cormat2[rownames(cormat2) %in% myexpressedIRs,colnames(cormat2) %in% myexpressedIRs]
coexp2 <- coexp
coexp2 <- coexp2[rownames(coexp2) %in% myexpressedIRs,colnames(coexp2) %in% myexpressedIRs]

#---------- Fig. S12D: Plotting histogram of pairwise Pearson correlations for 'expressed' IRs
YEScoexpressed <- cormat2[upper.tri(cormat2, diag=F)][coexp2[upper.tri(cormat2, diag=F)]]
NOcoexpressed <- cormat2[upper.tri(cormat2, diag=F)][!coexp2[upper.tri(cormat2, diag=F)]]
par(mfrow=c(2,1))
hist(NOcoexpressed, xlim=c(-0.2,0.9), breaks=seq(-0.2,0.9,0.02))
hist(YEScoexpressed, xlim=c(-0.2,0.9), breaks=seq(-0.2,0.9,0.02))


# Fig. S10: Dotplots of marker expression

MyOrcoMinusOSNs_Marks<-FindAllMarkers(MyOrcoMinusOSNs,logfc.threshold=0.3, return.thresh=1,min.pct=0,test.use='t', only.pos = T)
MyOrcoMinusOSNs_markers <- Extract_Top_Markers(marker_dataframe = MyOrcoMinusOSNs_Marks, num_genes = 900, named_vector = FALSE,
                                     make_unique = TRUE)
Clustered_DotPlot(seurat_object = MyOrcoMinusOSNs, features = MyOrcoMinusOSNs_markers, exp_color_min=0, exp_color_max = 2,
                  colors_use_exp = viridis(50, direction = -1, option = "D"), plot_km_elbow=F, raster = T)

# Fig. 5A: Plot phylogenetic VS genome distance for pairs of receptors

IRcoexp <- coexp2 # use geneset that includes all IRs with log2FC'>0.4

#---------- Load genomic distances
# Use the following linux code to extract OR and IR start positions from the updated annotation file
# awk '$3 == "transcript" && $12 ~ /^\"[OI]r/ {split($12, a, "\""); print $1 "\t" $4 "\t" a[2]}' AaegyptiLVP_AGWG_ThreePrimeUTRextended_Adavi2024.gtf > ORIRpositions.txt
# Navigate to directory containing the new file
genomicPositions <- read.table("ORIRpositions.txt", header=T)
rownames(genomicPositions) <- genomicPositions$Gene
calculate_distance <- function(gene1, gene2) {
  if (gene1$Chr == gene2$Chr) {distance <- abs(gene1$Start - gene2$Start)} #if genes on same chromosome
  else {distance <- 600000000} #or else just return very large number
  return(distance)
}

#---------- Load pairwise phylogenetic distances from treefile
library(ggtree)
library(ape)
# Navigate to directory containing the new file
tree_Aedes_irs <- read.tree(file = "Aedes_IR_tree.txt") # Load IR tree
tree_Aedes_irs <- root(tree_Aedes_irs,outgroup='Ir25a') # Re-root tree
IRpairwise_phylo <- cophenetic.phylo(tree_Aedes_irs)

#---------- Populate matrices
IRgenodist <- matrix(nrow=nrow(IRcoexp), ncol=ncol(IRcoexp))
rownames(IRgenodist) <- rownames(IRcoexp)
colnames(IRgenodist) <- colnames(IRcoexp)
IRphylodist <- IRgenodist
IRcoexpcol <- matrix("black",nrow=nrow(IRcoexp), ncol=ncol(IRcoexp))
IRcoexpcol[IRcoexp] <- "green"
for (i in (1:nrow(IRcoexp))) {
  for (j in (1:ncol(IRcoexp))) {
    gene1 <- rownames(IRcoexp)[i]
    gene2 <- colnames(IRcoexp)[j]
    IRgenodist[i,j] <- calculate_distance(genomicPositions[gene1,],genomicPositions[gene2,])
    IRphylodist[i,j] <- IRpairwise_phylo[gene1,gene2]
  }
}

#---------- Fig. 5A: Plot genomic by phylogenetic distances (colored by whether genes are coexpressed)
par(mfrow=c(2,1))
MyGenodists = log10(IRgenodist[upper.tri(IRcoexp)])
for (i in 1:length(MyGenodists)) { if (MyGenodists[i]==log10(600000000)) MyGenodists[i]=jitter(9,factor=0.3) }
plot(MyGenodists,IRphylodist[upper.tri(IRcoexp)],xlab="Genomic distance (log10(bp))",ylab="Phylogenetic distance",col=IRcoexpcol[upper.tri(IRcoexp)])
plot(IRgenodist[upper.tri(IRcoexp)],IRphylodist[upper.tri(IRcoexp)],xlab="Genomic distance (bp)",ylab="Phylogenetic distance",col=IRcoexpcol[upper.tri(IRcoexp)])

#---------- Fig. 5A: Plot marginal densities for phylogenetic distances (separately for genes that are or are not coexpressed)
par(mfrow=c(2,1))
plot(density(IRphylodist[upper.tri(IRcoexp) & !IRcoexp]),xlim=c(0,0.14))
plot(density(IRphylodist[upper.tri(IRcoexp) & IRcoexp]),xlim=c(0,0.14),col="green")

