library(pheatmap)
cols = colorRampPalette(c("blue", "white", "red"))(100)

marker_types = c('Contraction', 
    'Contraction', 
    'Contraction', 
    'Contraction', 
    'OxPhos', 
    'OxPhos', 
    'OxPhos', 
    'OxPhos', 
    'Calcium', 
    'Calcium', 
    'Calcium', 
    'Calcium', 
    'CellCycle',
    'CellCycle',
    'CellCycle',
    'CellCycle',
    'MYC',
    'MYC',
    'MYC',
    'MYC',
    'WNT/b-Catenin', 
    'WNT/b-Catenin', 
    'WNT/b-Catenin', 
    'WNT/b-Catenin', 
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'ECM',
    'Skeletal',
    'Skeletal',
    'Skeletal',
    'Skeletal',
    'Skeletal',
    'Skeletal',
    'Mesen/Fibro', 
    'Mesen/Fibro', 
    'Mesen/Fibro', 
    'Mesen/Fibro')

df = read.csv('Fig5_heatmap_erms.csv')
df = df[, -1]

markers = c('ACTN2', 
           'MYL4', 
           'TNNI1', 
           'TNNT2', 
           'ATP5J2',
           'COX6B1', 
           'NDUFA1', 
           'UQCR10', 
           'ATP2A1', 
           'ATP2B1', 
           'RYR1', 
           'RYR2', 
           'CCND2', 
           'CDK4', 
           'CENPF', 
           'MKI67', 
           'ENO1', 
           'MYC', 
           'NPM1', 
           'PTMA', 
           'CHD8', 
           'CYR61', 
           'GPC3', 
           'VCAN', 
           'LAMA2', 
           'LAMA4', 
           'LAMB1', 
           'LAMC1', 
           'FBLN5', 
           'SDC2', 
           'DCN', 
           'FN1', 
           'LGALS1', 
           'MFAP2', 
           'COL1A1', 
           'COL3A1', 
           'COL5A2', 
           'COL11A1', 
           'COL12A1', 
           'COL14A1', 
           'PCOLCE', 
           'FBLN2', 
           'LUM', 
           'SULF2', 
           'EBF2', 
           'MGP', 
           'OGN', 
           'OSR1', 
           'POSTN', 
           'TNMD', 
           'FAP', 
           'PDGFRA', 
           'THY1', 
           'VIM')

expr = df[, -c(1:2)]
expr[is.na(expr)] = 1e-6
colnames(expr) = markers
expr = t(expr)

annot = df[, 1:2]
colnames(annot) = c("CellState", "Tumor")
colnames(expr) = 1:28

annot2 = as.data.frame(marker_types)
rownames(annot2) = markers

expr = t(scale(t(expr)))

breaks = seq(-2, 2, length.out=101)

pdf('Fig5_erms.pdf', width=9, height=9.5)
## pheatmap(expr, cluster_rows = F, annotation_col=annot, annotation_row=annot2,
##          cluster_cols = F, scale='column') #annotation_col=as.data.frame(markers))
pheatmap(expr, cluster_rows = F, annotation_col=annot, annotation_row=annot2,
         cluster_cols = F, scale='none', color=cols, breaks=breaks) #annotation_col=as.data.frame(markers))
dev.off()

df = read.csv('Fig5_heatmap_arms.csv')
df = df[, -1]

markers = c('ACTN2', 
           'MYL4', 
           'TNNI1', 
           'TNNT2', 
           'ATP5J2',
           'COX6B1', 
           'NDUFA1', 
           'UQCR10', 
           'ATP2A1', 
           'ATP2B1', 
           'RYR1', 
           'RYR2', 
           'CCND2', 
           'CDK4', 
           'CENPF', 
           'MKI67', 
           'ENO1', 
           'MYC', 
           'NPM1', 
           'PTMA', 
           'CHD8', 
           'CYR61', 
           'GPC3', 
           'VCAN', 
           'LAMA2', 
           'LAMA4', 
           'LAMB1', 
           'LAMC1', 
           'FBLN5', 
           'SDC2', 
           'DCN', 
           'FN1', 
           'LGALS1', 
           'MFAP2', 
           'COL1A1', 
           'COL3A1', 
           'COL5A2', 
           'COL11A1', 
           'COL12A1', 
           'COL14A1', 
           'PCOLCE', 
           'FBLN2', 
           'LUM', 
           'SULF2', 
           'EBF2', 
           'MGP', 
           'OGN', 
           'OSR1', 
           'POSTN', 
           'TNMD', 
           'FAP', 
           'PDGFRA', 
           'THY1', 
           'VIM')

expr = df[, -c(1:2)]
expr[is.na(expr)] = 1e-6
colnames(expr) = markers
expr = t(expr)

annot = df[, 1:2]
colnames(annot) = c("CellState", "Tumor")
colnames(expr) = 1:16

expr = t(scale(t(expr)))
breaks = seq(-2, 2, length.out=101)

library(gplots)
library(RColorBrewer)

pdf('Fig5_arms.pdf', width=6, height=9.5)
## pheatmap(expr, cluster_rows = F, annotation_col=annot, annotation_row=annot2,
##          cluster_cols = F, scale='column') #annotation_col=as.data.frame(markers))
pheatmap(expr, cluster_rows = F, annotation_col=annot, annotation_row=annot2,
         cluster_cols = F, scale='none', color=cols, breaks=breaks) #annotation_col=as.data.frame(markers))
dev.off()
