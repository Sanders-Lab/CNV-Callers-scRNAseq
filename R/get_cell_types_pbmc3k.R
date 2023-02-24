library(SeuratData)
library(tidyverse)
library(Seurat)
library(SeuratDisk)

# one time affair
InstallData("pbmc3k")

# load the pbmc3k data
data("pbmc3k")
pbmc3k

median(pbmc3k@meta.data$nCount_RNA)
median(pbmc3k@meta.data$nFeature_RNA)


reference <- LoadH5Seurat("/fast/groups/ag_sanders/work/tools/pbmc_mapping_ref_seurat/pbmc_multimodal.h5seurat")


# The reference was normalized using SCTransform(),
# so we use the same approach to normalize the query here.
pbmc3k <- SCTransform(pbmc3k, verbose = FALSE)


# We then find anchors between reference and query.
# As described in the manuscript, we used a precomputed 
# supervised PCA (spca) transformation for this example. 
# We recommend the use of supervised PCA for CITE-seq datasets, 
# and demonstrate how to compute this transformation on the next 
# tab of this vignette. However, you can also use a standard PCA transformation.)
anchors <- FindTransferAnchors(reference = reference,
                               query = pbmc3k,
                               normalization.method = "SCT",
                               reference.reduction = "spca",
                               dims = 1:50)



pbmc3k <- MapQuery(anchorset = anchors,
                   query = pbmc3k,
                   reference = reference,
                   refdata = list(celltype.l1 = "celltype.l1",
                                  celltype.l2 = "celltype.l2",
                                  predicted_ADT = "ADT"
                                  ),
                   reference.reduction = "spca",
                   reduction.model = "wnn.umap")




# Idents(pbmc3k) <- 'predicted.celltype.l2'
Idents(pbmc3k) <- 'predicted.celltype.l1'
Idents(pbmc3k)


annot <- pbmc3k@meta.data$predicted.celltype.l1
pbmc3k@meta.data$predicted.celltype.l1 %>% is.na() %>% table()

# checking the freq of the annotations
pbmc3k@meta.data$seurat_annotations %>% is.na()
pbmc3k@meta.data


# using the seurat annotations for the cell group
pbmc3k@meta.data$cell <- rownames(pbmc3k@meta.data)
# annot <- pbmc3k@meta.data$seurat_annotations
cell <- pbmc3k@meta.data$cell

sample_annot <- as.data.frame(cell)
sample_annot$group <- annot
sample_annot %>% head
sample_annot %>% nrow
sample_annot$group %>% table()


write.table(sample_annot,
            gzfile('../proc/pbmc3k_cell_group_annot.tsv.gz'),
            row.names = F)
