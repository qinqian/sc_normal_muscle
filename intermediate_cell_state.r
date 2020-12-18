library(CytoTRACE)
library(patchwork)
library(ggplot2)

pdxs = Sys.glob('../01_sc_rms/data/seurat_obj/*rds')[1:10]
cols = read.delim('../01_sc_rms/final_annotations/Final_clusters.txt', sep='\t', row.names=1,
                  header=T,
                  check.names=F, stringsAsFactors=F)
pdxs = c(pdxs, '../01_sc_rms/results/seurat_sara/20191031_MSK74711_seurat-object.rds',
         '../01_sc_rms/results/seurat_sara/MAST118_seurat-object.rds',
         '../01_sc_rms/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds',
         '../01_sc_rms/results/seurat_sara/MAST139_1cells_seurat-object.rds')
pdxs = lapply(pdxs, readRDS)

fish <- lapply(c('../01_sc_rms/results/seurat_v6/Tumor21_recluster1.8.rds',
                 '../01_sc_rms/results/seurat_v6/Tumor22_recluster1.8.rds',
                 '../01_sc_rms/results/seurat_intersect_velocity/Tumor24_seu.rds'), readRDS)

mast139 = Sys.glob('../01_sc_rms/phaseA_explore_rms/*_cytotraceCytoTRACE_plot_table.txt')[7:8]
test = read.table(mast139[2])
test2 = read.table(mast139[1])

pdf('Rplots.pdf', width=18, height=12)
p3 = ggplot(data=test2) + geom_point(aes(x=Component1, y=Component2, color=Phenotype))
p4 = ggplot(data=test2) + geom_point(aes(x=Component1, y=Component2, color=CytoTRACE))
p1 = ggplot(data=test) + geom_point(aes(x=Component1, y=Component2, color=Phenotype))
p2 = ggplot(data=test) + geom_point(aes(x=Component1, y=Component2, color=CytoTRACE))
print(p1 + p2 + p3 + p4 + plot_layout(ncol=2))
dev.off()
