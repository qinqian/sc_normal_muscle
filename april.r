library(Seurat)
library(foreach)
library(dplyr)

a = Sys.glob('/data/langenau/human_rms_pdxs/data_april/*Robj')

for (i in a) {
    load(i)
    print(i)
}

rm(a)
result = list()
rm(i)
for (i in  ls()) {
    if (grepl("AP", i)) {
        print(i)
        seu = UpdateSeuratObject(get(i))
        ## seu = as.loom(seu, filename=paste0(i, '.loom'), verbose=T)
        ## seu$close_all()
        result[[i]] = seu
        }
}

result = Reduce(merge, result)
dim(result@assays$RNA@counts)

## seu = as.loom(result, filename='normal.muscle.loom', verbose=T)
## seu$close_all()

result <- NormalizeData(object = result, normalization.method = "LogNormalize", scale.factor = 10000)
result <- FindVariableFeatures(object = result, selection.method = "vst",
                               mean.function = ExpMean, dispersion.function = LogVMR,
                               x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
result <- ScaleData(object = result)
highvar.genes <- head(VariableFeatures(object = result), 1000)
result <- RunPCA(object = result, features = highvar.genes,
                 npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                 reduction.name = "pca", reduction.key = "PC_", seed.use = 123)
result <- FindNeighbors(object = result, k.param = 20, dims = 1:20, reduction = "pca")
result <- FindClusters(object = result, reduction.type = "pca", dims.use = 1:20, 
                       algorithm = 1, 
                       resolution = c(0.4, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
                       random.seed = 123)
result <- RunUMAP(result, reduction.use = "pca", dims = 1:20, reduction.name = "umap")
result <- RunTSNE(result, reduction.use = "pca", dims.use = 1:20, reduction.name = "tsne")

Idents(result) = result$seurat_clusters = result$RNA_snn_res.0.8
allmarkers = FindAllMarkers(result)

source('../01_sc_rms/phaseA_explore_rms/DEGs_seurat3_sara.R')

allcluster <- names((result@meta.data$seurat_clusters) %>% table())

clusterde <- list()
for (i in allcluster) {
    print(i)
    print(allcluster[-(as.integer(i)+1)])
    de.up <- get_upregulated_genes_in_clusters(result, i, allcluster[-(as.integer(i)+1)])
    if (nrow(de.up) > 0) {
        de.up$cluster <- i
    }
    clusterde[[i]] <- de.up
}

seurat1.de <- do.call('rbind', clusterde) ## still too many genes
## again, filter by adjusted p value and fraction of cells expressing the genes
seurat1.de <- subset(seurat1.de, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1))
write.table(seurat1.de, file='normal_muscle_res0.8_diffgenes_orig.csv', sep='\t', quote=F, col.names=NA)

head(result@meta.data)

library(readr)
write.table(result@meta.data, file='normal_muscle_metadata.csv', sep='\t', quote=F, col.names=NA)
write_csv(as.data.frame(result@reductions$umap@cell.embeddings), 'normal_muscle_umap.csv')
write_csv(allmarkers, 'normal_muscle_res0.8_diffgenes.csv')

erms.topsign = c("ITM2A", "FABP5", "EMILIN1", "RRBP1", "CCND1", "EMP3", "EIF4EBP1", "TSTA3", "ADA", "HMGA2")
arms.topsign = c("HSPB2", "MYL4", "PIPOX", "TNNT2", "MYOG", "ENO3", "NRN1", "GYPC", "TSPAN3", "TFF3")
library(patchwork)
library(ggplot2)
result@meta.data$age = apply(matrix(result@meta.data$orig.ident), 1, function(x) {paste(unlist(strsplit(x, '_'))[2:3], collapse='_')})
result@meta.data$age = factor(result@meta.data$age, levels=c('Emb_Wk06.0', 'Emb_Wk06.5', 'Emb_Wk07.0', 'Emb_Wk07.25', 'Emb_Wk07.75', 'Fet_Wk12', 'Fet_Wk14', 'Fet_Wk17', 'Fet_Wk18', 'Juv_Yr07', 'Juv_Yr11', 'Adt_Yr34', 'Adt_Yr42'))

pdf('erms_arms_markers.pdf', width=15, height=5)
p5 <- DotPlot(result, features=c(erms.topsign, arms.topsign), cols = c("lightgrey", "red"), group.by='RNA_snn_res.0.8') + xlab('Top RMS markers') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(plot.margin = unit(c(3,1,1,2), "cm"))
## p6 <- DotPlot(result, features=arms.topsign, cols = c("lightgrey", "red"), group.by='RNA_snn_res.0.8') + xlab('ARMS top markers')
p3 <- DimPlot(result, reduction='umap',
              group.by='age', label = F)+theme(plot.margin = unit(c(1,0.5,0.5,1), "cm"))
result@meta.data$age = factor(result@meta.data$age, levels=rev(c('Emb_Wk06.0', 'Emb_Wk06.5', 'Emb_Wk07.0', 'Emb_Wk07.25', 'Emb_Wk07.75', 'Fet_Wk12', 'Fet_Wk14', 'Fet_Wk17', 'Fet_Wk18', 'Juv_Yr07', 'Juv_Yr11', 'Adt_Yr34', 'Adt_Yr42')))
p1 <- DotPlot(result, features=c(erms.topsign, arms.topsign), cols = c("lightgrey", "red"), group.by='age') + xlab('Top RMS markers')+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.margin = unit(c(1,0.5,0.5,1), "cm"))
## p2 <- DotPlot(result, features=arms.topsign, cols = c("lightgrey", "red"), group.by='age') + xlab('ARMS top markers')
## p4 <- FeaturePlot(result, features='PAX7')
## p4 <- DimPlot(result, group.by='RNA_snn_res.0.8')
## print(p5+p6+p1 + p2 + p3 + p4 + plot_layout(ncol=2))
print(p3+p1+plot_layout(ncol=2))
dev.off()

pdf('erms_umap.pdf', width=18, height=12)
p4 <- FeaturePlot(result, features=erms.topsign)
print(p4)
dev.off()

pdf('arms_umap.pdf', width=18, height=12)
p4 <- FeaturePlot(result, features=arms.topsign)
print(p4)
dev.off()

a=result[erms.topsign, ]@assays$RNA
result@meta.data[, 'ERMS core signature'] = colMeans(as.matrix(a@data))
a=result[arms.topsign, ]@assays$RNA
result@meta.data[, 'ARMS core signature'] = colMeans(as.matrix(a@data))

png('core_score_umap.png', width=3600, height=1350, res=300)
p4 <- FeaturePlot(result, features=c("ARMS core signature", 'TNNT2', 'MYOG', 'MYL4',
                                     'ERMS core signature', 'ITM2A', 'CCND1', 'HMGA2'), ncol=4)
print(p4)
dev.off()

rms.core = read.table('../01_sc_rms/phaseA_explore_rms/RMS_core_t_test_pval0.05_fold1.5.xls')
fc = apply(rms.core, 1, function(x) {
    mean(log2(x[1:4])) - mean(log2(x[5:11]))
})

pdf('erms_all_umap.pdf', width=36, height=150)
p4 <- FeaturePlot(result, features=rownames(rms.core)[fc<0])
print(p4)
dev.off()
pdf('arms_all_umap.pdf', width=36, height=150)
p4 <- FeaturePlot(result, features=rownames(rms.core)[fc>0])
print(p4)
dev.off()

MAST111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')
MAST139 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')
MAST39 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST39.rds')
RH74 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_RH74.rds')
MAST95 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST95.rds')
MAST85 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST85.rds')
MSK82489 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MSK82489.rds')
MAST35 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST35.rds')

## new added samples
MSK72117_10cells = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds')
MAST118 = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/MAST118_seurat-object.rds')
MSK74711 = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK74711_seurat-object.rds')
MSK72117_10cells@meta.data$orig.ident = gsub('20191031_', '', MSK72117_10cells@meta.data$orig.ident)
MSK74711@meta.data$orig.ident = gsub('20191031_', '', MSK74711@meta.data$orig.ident)

for (i in list(MAST111, MAST139, MAST39, RH74, MAST95, MAST85, MSK82489, MAST35,
               MSK72117_10cells, MAST118, MSK74711)) {
    seu = as.loom(i, filename=paste0(i@meta.data$orig.ident[1], '.loom'), verbose=T)
    seu$close_all()
}

all.obj = Reduce('merge', list(MAST111, MAST139, MAST39, RH74, MAST95, MAST85, MSK82489, MAST35,
                               MSK72117_10cells, MAST118, MSK74711))

primary.obj <- Reduce("merge", list(
                                readRDS('../01_sc_rms/phaseA_explore_rms/20082_recluster2_tumor_only.rds'),
                                readRDS('../01_sc_rms/figures/20696_hg19_tumoronly_res0.8_umap.rds'),
                                readRDS('../01_sc_rms/figures/21202_hg19_premrna_tumoronly_res0.8_umap.rds'),
                                readRDS('../01_sc_rms/figures/29806_hg19_premrna_tumoronly_res0.8_umap.rds')))

primary.obj@meta.data$orig.ident <- gsub('_hg19_premrna', '', primary.obj@meta.data$orig.ident)

library(stringr)
result.all = list()
rm(i)
for (i in  ls()) {
    if (grepl("AP", i)) {
        print(i)
        seu = UpdateSeuratObject(get(i))
        result.all[[i]] = seu
    }
}

result.all = Reduce(merge, result.all)
result.all = merge(result.all, primary.obj)
result.all = merge(result.all, all.obj)

result.all <- NormalizeData(object = result.all, normalization.method = "LogNormalize", scale.factor = 10000)
result.all <- FindVariableFeatures(object = result.all, selection.method = "vst",
                               mean.function = ExpMean, dispersion.function = LogVMR,
                               x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
result.all <- ScaleData(object = result.all)
highvar.genes <- head(VariableFeatures(object = result.all), 1000)
result.all <- RunPCA(object = result.all, features = highvar.genes,
                 npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                 reduction.name = "pca", reduction.key = "PC_", seed.use = 123)
result.all <- FindNeighbors(object = result.all, k.param = 20, dims = 1:20, reduction = "pca")
result.all <- FindClusters(object = result.all, reduction.type = "pca", dims.use = 1:20, 
                       algorithm = 1, 
                       resolution = c(0.4, 0.8, 1, 1.2),
                       random.seed = 123)
result.all <- RunUMAP(result.all, reduction.use = "pca", dims = 1:20, reduction.name = "umap")
## result.all <- RunTSNE(result.all, reduction.use = "pca", dims.use = 1:20, reduction.name = "tsne")

library(ggplot2)
pdf('tumor_with_normal.pdf', width=24, height=18)
p3 <- DimPlot(result.all, reduction='umap',
              group.by='orig.ident', label = TRUE) + theme(legend.position='none') 
print(p3)
dev.off()

integration <- function(x, y, label, method='seurat') {
    seurat.pseudo.list = list(x, y)
    for (i in 1:length(seurat.pseudo.list)) {
        print(1)
        seurat.pseudo.list[[i]] <- NormalizeData(seurat.pseudo.list[[i]], verbose = FALSE)
        seurat.pseudo.list[[i]] <- FindVariableFeatures(seurat.pseudo.list[[i]],
                                                        selection.method = "vst", 
                                                        nfeatures = 2000, verbose = FALSE)
    }
    names(seurat.pseudo.list) <- label
    if (method == 'seurat') {
        anchors <- FindIntegrationAnchors(object.list = seurat.pseudo.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        DefaultAssay(integrated) <- "integrated"
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
    } else if (method == 'conos') {
        seurat.pseudo.list <- lapply(seurat.pseudo.list, function(x) {
            ScaleData(x) %>% RunPCA(verbose=F)
        })
        pseudo.con <- Conos$new(seurat.pseudo.list)
        pseudo.con$buildGraph(k=15, k.self=5, space="PCA", ncomps=30, n.odgenes=2000, matching.method='mNN',
                              metric = 'angular', score.component.variance=T, verbose=T)
        pseudo.con$findCommunities()
        pseudo.con$embedGraph()
        integrated = as.Seurat(pseudo.con)
    }
    integrated
}

result.mast139 = integration(result, MAST139, c("Normal", "MAST139"))
result.mast139$tumor = result.mast139$orig.ident == 'MAST139'
pdf('integration_normal_mast139.pdf', width=18, height=8)
p1=DimPlot(result.mast139, reduction='umap', 
          group.by='tumor', label = F) # + theme(legend.position='none')
p2=DimPlot(result.mast139, reduction='umap', 
           group.by='orig.ident', label = F) # + theme(legend.position='none')
print(p1+p2)
dev.off()

result.mast118 = integration(result, MAST118, c("Normal", "MAST118"))
result.mast118$tumor = result.mast118$orig.ident == 'MAST118'
pdf('integration_normal_mast118.pdf', width=18, height=8)
p1=DimPlot(result.mast118, reduction='umap', 
          group.by='tumor', label = F) # + theme(legend.position='none')
p2=DimPlot(result.mast118, reduction='umap', 
           group.by='orig.ident', label = F) # + theme(legend.position='none')
print(p1+p2)
dev.off()


primary.list = list(readRDS('../01_sc_rms/phaseA_explore_rms/20082_recluster2_tumor_only.rds'),
                    readRDS('../01_sc_rms/figures/20696_hg19_tumoronly_res0.8_umap.rds'),
                    readRDS('../01_sc_rms/figures/21202_hg19_premrna_tumoronly_res0.8_umap.rds'),
                    readRDS('../01_sc_rms/figures/29806_hg19_premrna_tumoronly_res0.8_umap.rds'))

for (i in primary.list) {
    result.mast118 = integration(result, i, c("Normal", i@meta.data$orig.ident[1]))
    result.mast118$tumor = result.mast118$orig.ident == i@meta.data$orig.ident[1]
    pdf(paste0('integration_normal_', i@meta.data$orig.ident[1], '.pdf'), width=18, height=8)
    p1=DimPlot(result.mast118, reduction='umap', 
               group.by='tumor', label = F) # + theme(legend.position='none')
    p2=DimPlot(result.mast118, reduction='umap', 
               group.by='orig.ident', label = F) # + theme(legend.position='none')
    print(p1+p2)
    dev.off()
}


for (i in list(MAST111, MAST139, MAST39, RH74, MAST95, MAST85, MSK82489, MAST35,
               MSK72117_10cells, MAST118, MSK74711))
{
    result.mast118 = integration(result, i, c("Normal", i@meta.data$orig.ident[1]))
    result.mast118$tumor = result.mast118$orig.ident == i@meta.data$orig.ident[1]
    pdf(paste0('integration_normal_', i@meta.data$orig.ident[1], '.pdf'), width=18, height=8)
    p1=DimPlot(result.mast118, reduction='umap', 
               group.by='tumor', label = F) # + theme(legend.position='none')
    p2=DimPlot(result.mast118, reduction='umap', 
               group.by='orig.ident', label = F) # + theme(legend.position='none')
    print(p1+p2)
    dev.off()
}

library(CytoTRACE)
test = list(as.matrix(result@assays$RNA@counts), as.matrix(MAST139@assays$RNA@counts),
            as.matrix(MAST118@assays$RNA@counts))
test.out = iCytoTRACE(test)

cols = read.delim('../01_sc_rms/final_annotations/Final_clusters.txt')
rownames(cols) = cols[,1]
cols = cols[,-1]

mast139.muscle = rownames(MAST139@meta.data)[MAST139@meta.data$RNA_snn_res.0.8 == which(as.vector(cols['MAST139', ])=='Muscle')-1]
mast118.muscle = rownames(MAST118@meta.data)[MAST118@meta.data$RNA_snn_res.0.8 == which(as.vector(cols['MAST118', ])=='Muscle')-1]

meta = rbind(result@meta.data[,3,drop=F], MAST139@meta.data[,1,drop=F], MAST118@meta.data[,1,drop=F])
test.out = cbind(test.out$CytoTRACE, meta)
colnames(test.out) = c("CytoTRACE", "Ident")

normal.muscle = test.out[grepl('AP', rownames(test.out)), ]
mast139.muscle = test.out[rownames(test.out) %in% mast139.muscle, ]
mast118.muscle = test.out[rownames(test.out) %in% mast118.muscle, ]
muscles=  rbind(normal.muscle, mast139.muscle, mast118.muscle)

library(ggplot2)
ggplot(muscles, aes(Ident, CytoTRACE)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90))
ggsave('cyto_mast139_mast118_norm_muscle.pdf')
ggplot(test.out, aes(Ident, CytoTRACE)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90))
ggsave('cyto_mast139_mast118_norm.pdf')

## Normal muscle markers
pdf('muscle_biomarker_all_umap.pdf', width=16, height=7)
p4 <- FeaturePlot(result, features=c("DES", "MKI67"))
print(p4)
dev.off()

## other paper biomarkers
pdf('other_paper_biomarkers_all_umap.pdf', width=36, height=150)

dev.off()
