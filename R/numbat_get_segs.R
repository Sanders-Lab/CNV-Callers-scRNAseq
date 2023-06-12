library(data.table)
library(stringr)
library(magrittr)


joint_post <- fread("../outputs/tnbc1/numbat_pbmc3k_epithelium_tnbc1/joint_post_2.tsv")
joint_post

# joint_post[, .N, by = cell]
# joint_post[cnv_state == "amp", .N, by = cell]


nb_calls <- joint_post[, .(
    seqnames = CHROM,
    start = seg_start,
    end = seg_end,
    cell_name = cell,
    cnv_state = cnv_state_map,
    n_genes
)]
nb_calls


# replacing neu with ""
nb_calls[, cnv_state := str_replace(cnv_state, "neu", "")]
nb_calls


# changing the cell_name to match other tools
nb_calls[, cell_name := gsub(
    "^",
    "tnbc1_",
    cell_name
)]
nb_calls[, cell_name := gsub(
    "-1",
    "",
    cell_name
)]
nb_calls


# adding chr to seqnames
nb_calls[, seqnames := paste0("chr", seqnames)]
nb_calls



fwrite(
    x = nb_calls,
    file = "../proc/numbat_tnbc_cna_segs.tsv.gz",
    sep = "\t"
)



nb_calls[cnv_state_map == "amp", .N, by = cell_name]
nb_calls[cnv_state_map == "del", .N, by = cell_name]



nb_calls$cell %>%
    unique() %>%
    length()

grep("amp", nb_calls$cnv_state) %>%
    length()
nb_calls[, count()]

nb_calls[cnv_state == "amp", .N, by = cell][, unique(N)]

nb_calls[, unique(cnv_state)]

nb_calls[, length(grep("amp", cnv_state)), by = cell][, unique(V1)]

nb_calls[, length(grep("loh", cnv_state)), by = cell][, unique(V1)]
nb_calls[, length(grep("del", cnv_state)), by = cell][, unique(V1)]
