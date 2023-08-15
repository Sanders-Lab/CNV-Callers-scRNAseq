library(data.table)
library(GenomicRanges)
library(biomaRt)


raw_segs <- fread("../proc/ct_loc.tsv")
raw_segs

raw_gro <- makeGRangesFromDataFrame(raw_segs,
                                    seqnames.field = "chrom",
                                    start.field = "start_loc",
                                    end.field = "end_loc",
                                    keep.extra.columns = T
)
raw_gro


cont_ranges <- reduce(raw_gro)
cont_dt <- as.data.table(cont_ranges)
cont_dt

unified_ranges <- data.table(seqnames = "chr6",
                             start = cont_dt[1, start],
                             end = cont_dt[6, end]
)
unified_ranges[, width := (end - start)]
unified_ranges[, strand := "*"]

unified_ranges
unified_gro <- makeGRangesFromDataFrame(unified_ranges)
unified_gro


# reading in rna counts of tall
raw_counts <- fread("../proc/pbmc3k_tall_ref_merged_counts.tsv.gz")
raw_counts[, 1:5]
gene_names <- raw_counts[, V1]
gene_names


# get the loc of the genes
hmart <- useMart(
                 biomart = "ENSEMBL_MART_ENSEMBL",
                 dataset = "hsapiens_gene_ensembl"
)

gene_pos <- biomaRt::getBM(attributes = c(
                                          "ensembl_gene_id",
                                          "hgnc_symbol",
                                          "chromosome_name",
                                          "start_position",
                                          "end_position"
                                          ),
                           filters = "hgnc_symbol",
                           values = gene_names,
                           mart = hmart,
                           useCache = F
)
gene_pos_dt <- as.data.table(gene_pos)
gene_pos_dt


# adding chr to chrom names
gene_pos_dt[, chromosome_name := paste0("chr", chromosome_name)]
gene_pos_dt


# finding out ref_hom segments
sv_calls <- fread("/fast/groups/ag_sanders/work/data/strand_seq/external/scTRIP_Final_SV_calls_Nov2018/TALL03-DEA5.strict.filtered.txt")
sv_calls

sv_gro <- makeGRangesFromDataFrame(sv_calls,
                                   seqnames.field = "chrom",
                                   keep.extra.columns = T
)
sv_gro
sv_calls[sv_call_name == "ref_hom"]


ideo_dt <- fread("../../digital_karyotype/proc/ideogram_scaffold.tsv")
ideo_dt


ideo_gro <- makeGRangesFromDataFrame(ideo_dt,
                                    seqnames.field = "chrom",
                                    start.field = "start_loc",
                                    end.field = "end_loc",
                                    keep.extra.columns = T
)
ideo_gro


# main logic to get the ref regions
dis_sv <- c(ideo_gro, sv_gro) |>
    disjoin() 
dis_sv
sv_overlaps <- subsetByOverlaps(dis_sv, sv_gro)
sv_overlaps

ref_regions <- setdiff(dis_sv, sv_overlaps)
findOverlaps(ref_regions, sv_overlaps)

ref_regions_dt <- as.data.table(ref_regions)
ref_regions_dt[, unique(seqnames)]


ref_regions_dt
fwrite(ref_regions_dt,
       "../proc/ref_regions_tall.tsv.gz",
       sep = "\t"
)

ref_ct_comb <- rbind(cont_dt,
      ref_regions_dt)

ref_ct_comb <- rbind(unified_ranges,
      ref_regions_dt)



################################################################################

# taking only the genes which overlap with the segments
gene_pos_gro <- makeGRangesFromDataFrame(gene_pos_dt,
                                         seqnames.field = "chromosome_name",
                                         start.field = "start_position",
                                         end.field = "end_position",
                                         keep.extra.columns = T
)
gene_pos_gro

ref_ct_gro <- makeGRangesFromDataFrame(ref_ct_comb)
ref_ct_gro


segs_gro <- subsetByOverlaps(gene_pos_gro, ref_ct_gro)
gene_segs_dt <- as.data.table(segs_gro)
gene_segs_dt

chr6 <- gene_segs_dt[seqnames == "chr6"]
chr6



raw_counts[, 1:5]
filt_counts <- raw_counts[V1 %in% gene_segs_dt[, hgnc_symbol]]

filt_counts <- raw_counts[V1 %in% chr6[, hgnc_symbol]]
filt_counts[, 1:5]

filt_counts_df <- as.data.frame(filt_counts)
rownames(filt_counts_df) <- filt_counts[, V1]
filt_counts_df$V1 <- NULL
rownames(filt_counts_df)
filt_counts_df[1:5, 1:5]

write.table(filt_counts_df,
      file = "../proc/ct_segs_ref_comb_filt_counts.tsv.gz",
      sep = "\t"
)

write.table(filt_counts_df,
      file = "../proc/ct_segs_ref_comb_chr6_filt_counts.tsv.gz",
      sep = "\t"
)

# for numbat, removing the pbmc cells

filt_counts_numbat <- filt_counts_df[,grep("tall_*", colnames(filt_counts_df))]
filt_counts_numbat[1:5,1:5]
dim(filt_counts_numbat)

write.table(filt_counts_numbat,
      file = "../proc/ct_segs_ref_comb_filt_counts_numbat.tsv.gz",
      sep = "\t"
)



write.table(filt_counts_numbat,
      file = "../proc/ct_segs_ref_comb_chr6_filt_counts_numbat.tsv.gz",
      sep = "\t"
)

read.table( "../proc/ct_segs_filt_counts.tsv.gz")[1:5, 1:5]



# checking the locations of segments and ref
source("../../digital_karyotype/R/utils.R", chdir = T)


ct_seg_genes <- as.data.table(subsetByOverlaps(gene_pos_gro, cont_ranges))
ct_seg_genes <- as.data.table(subsetByOverlaps(gene_pos_gro, unified_gro))
ct_seg_genes

gene_segs_dt[start %in% ct_seg_genes[, start],
             `:=`(cell_name = "filtered_matrix",
                    sv_state = "ct")]
gene_segs_dt[!start %in% ct_seg_genes[, start],
             `:=`(cell_name = "filtered_matrix",
                    sv_state = "ref")]
gene_segs_dt


setnames(gene_segs_dt,
         c("seqnames", "start", "end"),
         c("chrom", "start_loc", "end_loc")
)
gene_segs_dt
gene_segs_dt[sv_state == "ct"]

colors <- c("red", "green")
names(colors) <- c("ct", "ref")



str(Plot_Digital_Karyotype)
Plot_Digital_Karyotype(cell_sv_dt = gene_segs_dt,
                       plot_both_haplotypes = F,
                       kar_label = "CT seg and ref segs",
                       color_pal = colors,
                       save_digital_karyotype = T,
                       plot_dir = "filtered_counts_segs"

)

chr6_plot <- gene_segs_dt[chrom == "chr6"]
chr6_plot[, cell_name := "chr6"]

Plot_Digital_Karyotype(cell_sv_dt = chr6_plot,
                       plot_both_haplotypes = F,
                       kar_label = "CT seg and ref segs",
                       color_pal = colors,
                       save_digital_karyotype = T,
                       plot_dir = "filtered_counts_segs"

)


# making the bed file for numbat regenotyping
ideo_dt <- fread("../../digital_karyotype/proc/segmented_ideo.tsv.gz")
ideo_dt


ideo_gro <- makeGRangesFromDataFrame(ideo_dt,
                                    seqnames.field = "chrom",
                                    start.field = "start_loc",
                                    end.field = "end_loc",
                                    keep.extra.columns = T
)
ideo_gro

nb_seg_input <- as.data.table(disjoin(c(ideo_gro, unified_gro)))


setnames(nb_seg_input,
         c("seqnames", "start", "end"),
         c("CHROM", "seg_start", "seg_end")
)

nb_seg_input[, `:=`(width = NULL, strand = NULL)]
nb_seg_input

nb_seg_input[!seg_start %in% unified_ranges[, start],
             cnv_state := "neu"]
nb_seg_input[seg_start %in% unified_ranges[, start],
             cnv_state := "amp"]

# will manually change the segs
nb_seg_input[, seg := "1a"]
nb_seg_input


setcolorder(nb_seg_input,
            c("CHROM", "seg")
)
nb_seg_input


fwrite(nb_seg_input,
       "../proc/ct_segs_regenotype_numbat.tsv.gz",
       sep = "\t"
)


