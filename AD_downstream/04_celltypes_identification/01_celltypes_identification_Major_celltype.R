# 01_celltypes_identification_Major_celltype
rm(list = ls())
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyr)
setwd("~/AD_downstream/04_celltypes_identification/")

states=readRDS("~/AD_downstream/states_mixsfind.rds")

Idents(states) <- "states_nn_alg1_new"
table(states@active.ident)
markers <- FindAllMarkers(
  object          = states,
  assay           = "totalRNA",
  test.use        = "wilcox",
  only.pos        = TRUE,
  logfc.threshold = 0.1,
  min.pct         = 0.25
)   
all.markers = markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox_alg1_sub.csv", row.names = F)
write.csv(top10, "top10_diff_genes_wilcox_alg1_sub.csv", row.names = F)

new.cluster.ids <- c("0" = "Non_Neuron",
                     "1" = "Neuron",
                     "2" = "Neuron",
                     "3" = "Non_Neuron",
                     "4" = "Non_Neuron",
                     "5" = "Non_Neuron",
                     "6" = "Neuron",
                     "7" = "Non_Neuron",
                     "8" = "Non_Neuron",#
                     "9" = "Non_Neuron",
                     "10" = "Mix",
                     "11" = "Neuron",
                     "12" = "Neuron",
                     "13" = "Neuron",
                     "14" = "Neuron",
                     "15" = "Non_Neuron",
                     "16" = "Non_Neuron",
                     "17" = "Non_Neuron",
                     "18" = "Non_Neuron",
                     "19" = "Neuron",
                     "20" = "Non_Neuron",
                     "21" = "Non_Neuron",
                     "22" = "Neuron",
                     "23" = "Non_Neuron",
                     "24" = "Neuron",
                     "25" = "Non_Neuron")
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label1<- states@active.ident
Idents(states) <- "states_nn_alg1_label1"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()

Idents(states) <- "states_nn_alg1_new"
new.cluster.ids <- c( "0" = "OLG",
                      "1" = "TEPN",
                      "2" = "TEPN",
                      "3" = "VAS",
                      "4" = "AC",
                      "5" = "OLG",
                      "6" = "TEPN",
                      "7" = "OLG",
                      "8" = "AC",#
                      "9" = "MLG",
                      "10" = "Mix",
                      "11" = "TEPN",
                      "12" = "TEPN",
                      "13" = "TEPN",
                      "14" = "INH",
                      "15" = "AC",
                      "16" = "CHOR/EPEN",
                      "17" = "VAS",
                      "18" = "AC",
                      "19" = "TEPN",
                      "20" = "CHOR/EPEN",
                      "21" = "OPC",
                      "22" = "CHO/PEP",
                      "23" = "VAS",
                      "24" = "DE/MEN",
                      "25" = "VAS")
states <- RenameIdents(states, new.cluster.ids)                        
states$states_nn_alg1_label2<- states@active.ident
Idents(states) <- "states_nn_alg1_label2"
table(states@active.ident)
DimPlot(states,reduction = "states.umap",split.by = "type",label = T)+ NoLegend()
DimPlot(states,reduction = "states.umap",label = T)+ NoLegend()
table(states$type, Idents(states))
table(states$group, Idents(states))

saveRDS(states, file = "states_Major_celltype.rds")

