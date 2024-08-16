library(Seurat)
library(viridis)
library(plotrix)
library(scCustomize)
library(scales)
library(tidyr)

# LOAD pre-clustered thermo/hygrosensory neuron data object

# Navigate to directory containing the R object THSN_Seurat_object.Rdata
load("THSN_Seurat_object.Rdata") # Loading pre-clustered THSN data
THSN_Seurat <- THSN_12d_15n_15
total_droplets<-ncol(THSN_Seurat) #589 droplets
clustorder=c("0","1","4","8","11","6","5","3","10","7","2","9")

# Fig. 2B and S7: Visualize UMAPs for clusters, batches, and key marker genes

# Note that the AaegL5 gtf annotation used for the CellRanger alignments underlying these data uses
# "OrX", "IrX", "GrX" names for chemosensory receptors. Other genes are named using the associated
# AAEL number followed by the name of the Drosophila ortholog (*where orthologs include multiple
# mosquito genes that have the same best match in Drosophila). For example, "AAEL019818" was changed
# to "AAEL019818-nompC". See Methods.

DimPlot(THSN_Seurat, pt.size=0.8, label=T) # Fig. 2B
DimPlot(THSN_Seurat, group.by="orig.ident", pt.size=0.8) # by library
FeaturePlot(THSN_Seurat, feature="Ir25a", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(THSN_Seurat, feature="Ir93a", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(THSN_Seurat, feature="Ir76b", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(THSN_Seurat, feature="Ir21a", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(THSN_Seurat, feature="Ir40a", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)
FeaturePlot(THSN_Seurat, feature="Gr76", pt.size=0.6, cols=rev(viridis(100, option='D')), order=T)

# Fig. 2C: Dotplots of marker expression

MyTHSNs_Marks<-FindAllMarkers(THSN_Seurat,logfc.threshold=0.3, return.thresh=1,min.pct=0,test.use='t', only.pos = T)
Clustered_DotPlot(seurat_object=THSN_Seurat, features=MyTHSNs_Marks, exp_color_min=0, exp_color_max=2,
                  colors_use_exp=viridis(50, direction=-1, option="D"), plot_km_elbow=F, raster=T)

# Fig. S7D: Average expression heatmap for target receptors

TargetReceptors = rev(c("Ir40a","Ir21a","Gr76","Ir76b","Ir93a","Ir25a"))
clustorder <- reverse(c("0","1","4","8","11","6","5","3","10","7","2","9"))
AveExp <- AverageExpression(THSN_Seurat,features=TargetReceptors,return.seurat=T)
AveExp <- AveExp$SCT@data[,clustorder] # selecting scaled log average counts for target receptors
AveExp_Heatmap <- AveExp[,clustorder] # selecting log normalized average counts

par(mfrow=c(1,1))
Corcols=hcl.colors(450, palette = 'rocket', alpha = NULL, rev = T, fixup = TRUE)
max=ceiling(max(AveExp_Heatmap))
image(t(AveExp_Heatmap),col=Corcols,axes=F, zlim=c(0,max), main='Log average expression')
axis(3,at=seq(0,1,length=ncol(AveExp_Heatmap)),las=2,labels=colnames(AveExp_Heatmap),cex.axis=0.5,lwd=0)
axis(2,at=seq(0,1,length=nrow(AveExp_Heatmap)),las=2,labels=rownames(AveExp_Heatmap),cex.axis=0.3,lwd=0)
color.legend(0.05,-0.15,0.2,-0.17,rect.col=Corcols,legend=c(0,max),align='rb')

# Fig. S7E: log2FC' heatmap for target receptors

TargetReceptors = rev(c("Ir40a","Ir21a","Gr76","Ir76b"))
MyMarkers <-FindAllMarkers(THSN_Seurat,features=TargetReceptors,logfc.threshold=0, return.thresh=1, min.pct=0, test.use='t')
log2FC <- matrix(NA, nrow=length(TargetReceptors), ncol=length(levels(MyMarkers$cluster)))
colnames(log2FC) <- levels(MyMarkers$cluster)
rownames(log2FC) <- TargetReceptors
for(i in 1:nrow(MyMarkers)){
  log2FC[MyMarkers$gene[i],as.character(MyMarkers$cluster[i])]<-MyMarkers$avg_log2FC[i]
}
log2FC[is.na(log2FC)] <- 0 #replacing NAs with 0s 
med.att <- apply(log2FC, 1, median) #calculating median Log2FC values (for each gene across clusters)
log2FCprime <- sweep(log2FC, 1, med.att) #substract median Log2FC value from cluster-specific Log2FCs to get Log2FC'
log2FCP_Heatmap<-log2FCprime[TargetReceptors,clustorder]
log2FCP_Heatmap[log2FCP_Heatmap<0]<-0 #converting negative log2FC` values to zero

color_palette=c(colorRampPalette(c(alpha("white", alpha = 0),'#FDE333'))(10),hcl.colors(450, palette = 'viridis', alpha = NULL, rev = T, fixup = TRUE))
par(mfrow=c(1,1))
max = ceiling(max(log2FCP_Heatmap))
image(t(log2FCP_Heatmap),col=color_palette,axes=F,zlim=c(0,max), main='Log2FC`')
axis(3,at=seq(0,1,length=ncol(log2FCP_Heatmap)),las=2,labels=colnames(log2FCP_Heatmap),cex.axis=0.5,lwd=0)
axis(2,at=seq(0,1,length=nrow(log2FCP_Heatmap)),las=2,labels=rownames(log2FCP_Heatmap),cex.axis=0.3,lwd=0)
color.legend(0.05,-0.15,0.2,-0.17,rect.col=color_palette,legend=c(0,max),align='rb')

