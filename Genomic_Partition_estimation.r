library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


df <- read.delim("input_table.tsv", stringsAsFactors = FALSE)
query_gr <- makeGRangesFromDataFrame(df, seqnames.field="distal_chr",
start.field="distal_start", end.field="distal_end")
n_points <- 500 #can be more for better precision
steps <- seq(0, 1, length.out = n_points)
results <- lapply(1:length(query_gr), function(i) {
start <- start(query_gr[i])
end <- end(query_gr[i])
chr <- as.character(seqnames(query_gr[i]))
# Generate nucleotide position to be tested
pos <- round(start + (end - start) * steps)
pts_gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
# Test hierarchical grouping
is_p <- countOverlaps(pts_gr, proms) > 0
is_e <- (countOverlaps(pts_gr, exons) > 0) & !is_p
is_i <- (countOverlaps(pts_gr, introns) > 0) & !is_p & !is_e
# Ratio over total length computaion
w <- width(query_gr[i])
bp_p <- (sum(is_p) / n_points) * w
bp_e <- (sum(is_e) / n_points) * w
bp_i <- (sum(is_i) / n_points) * w
bp_inter <- w - (bp_p + bp_e + bp_i)
# Partitioning
vec <- c("Promoter"=bp_p, "Exon"=bp_e, "Intron"=bp_i, "Distal Intergenic"=bp_inter)
return(c(bp_p, bp_e, bp_i, bp_inter, names(vec)[which.max(vec)]))
})
# 3. Assembly
res_mat <- do.call(rbind, results)
df$bp_Promoter <- as.numeric(res_mat[,1])
df$bp_Exon <- as.numeric(res_mat[,2])
df$bp_Intron <- as.numeric(res_mat[,3])
df$bp_Intergenic <- as.numeric(res_mat[,4])
df$Genomic_Partition <- res_mat[,5]
# 4. Save
df$genomic_partition_distal <- assign_partition(query_gr, promoters_gr, exons_gr, introns_gr)
with(df, ftable(Genomic_Partition ~ genomic_partition_distal))

write.table(df[,-c(31,32,33,34)], file = "annotations_combined_all_WITH_PARTITION_2_columns.tsv", sep="\t", quote=FALSE, row.names=FALSE)
