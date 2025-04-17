#!/usr/bin/env Rscript
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
s <- readDNAStringSet("/mnt/ceph/454_data/Pisum_assembly_ver_2/assembly/230509_release_3/Cameor_v2_release_3.fasta")

# smallest version
gr_nano <- GRanges(seqnames = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
                          ranges = IRanges(start = 1, end = c(5e6, 5e6, 3e6, 3e6, 1e6))
)
s_nano <- getSeq(s, gr_nano)
names(s_nano) <- seqnames(gr_nano)

gr_micro <- GRanges(seqnames = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
                          ranges = IRanges(start = 1, end = c(1e7, 1e7, 1e7, 5e6, 5e6))
)

s_micro <- getSeq(s, gr_micro)
names(s_micro) <- seqnames(gr_micro)


gr_smallest <- GRanges(seqnames = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
                          ranges = IRanges(start = 1, end = c(3e7, 2e7, 1e7, 1e7, 1e7))
)

s_smallest <- getSeq(s, gr_smallest)
names(s_smallest) <- seqnames(gr_smallest)

gr_mid <- GRanges(seqnames = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
                          ranges = IRanges(start = 1, end = c(9e7, 8e7, 6e7, 5e7, 4e7))
)
s_mid <- getSeq(s, gr_mid)
names(s_mid) <- seqnames(gr_mid)

gr_large <- GRanges(seqnames = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
                          ranges = IRanges(start = 1, end = c(1.5e8, 1.2e8, 1e8, 9e7, 8e7))
)
s_large <- getSeq(s, gr_large)
names(s_large) <- seqnames(gr_large)

# Save the sequences to a file

out_dir <- "/mnt/ceph/454_data/workshop/2025/example_analysis/input_data"

writeXStringSet(s_smallest, filepath = file.path(out_dir, "tiny_pea.fasta"), format = "fasta")
writeXStringSet(s_mid, filepath = file.path(out_dir, "lite_pea.fasta"), format = "fasta")
writeXStringSet(s_large, filepath = file.path(out_dir, "chunky_pea.fasta"), format = "fasta")
writeXStringSet(s_micro, filepath = file.path(out_dir, "micro_pea.fasta"), format = "fasta")
writeXStringSet(s_nano, filepath = file.path(out_dir, "nano_pea.fasta"), format = "fasta")