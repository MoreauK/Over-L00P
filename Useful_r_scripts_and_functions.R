
# This script gives the key functions to produce the input table for OveR-L00P.
# The First function:  simplify_maser_table, take as input the SE table produce by maser package from rMATS output: 
#
# library(maser)
# RMATS_output <- maser("./", c("cond1", "cond2"), ftype = "JC")
# res_filtered <- filterByCoverage(RMATS_output, avg_reads = 10)
# res_filtered_top <- topEvents(res_filtered, fdr = 0.05, deltaPSI = 0.1)
# SE <- summary(res_filtered_top, type="SE") # This can be directly used by simplify_maser_table

library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)

simplify_maser_table <- function(df, overlap_threshold = 0.5) {
  df <- df %>%
    separate(exon_target, into = c("start", "end"), sep = "-", convert = TRUE)
  
  df$IncLevelDifference <- as.numeric(df$IncLevelDifference)
  
  gr <- GRanges(seqnames = df$Chr,
                ranges = IRanges(start = df$start, end = df$end),
                strand = df$Strand,
                geneSymbol = df$geneSymbol)
  
  hits <- findOverlaps(gr, gr, ignore.strand = FALSE)
  clusters <- split(subjectHits(hits), queryHits(hits))
  
  grouped <- list()
  used <- rep(FALSE, length(gr))
  
  for (i in seq_along(gr)) {
    if (used[i]) next
    curr_group <- clusters[[as.character(i)]]
    group_members <- c(i)
    
    for (j in curr_group) {
      if (i == j || used[j]) next
      ov_start <- max(start(gr[i]), start(gr[j]))
      ov_end <- min(end(gr[i]), end(gr[j]))
      ov_len <- max(0, ov_end - ov_start + 1)
      len_i <- width(gr[i])
      len_j <- width(gr[j])
      perc_overlap <- ov_len / min(len_i, len_j)
      
      if (perc_overlap >= overlap_threshold &&
          df$geneSymbol[i] == df$geneSymbol[j] &&
          df$Strand[i] == df$Strand[j] &&
          df$Chr[i] == df$Chr[j]) {
        group_members <- c(group_members, j)
        used[j] <- TRUE
      }
    }
    used[i] <- TRUE
    grouped[[length(grouped) + 1]] <- unique(group_members)
  }
  
  aggregate_rows <- function(indices) {
    subset_df <- df[indices, ]
    inc_diff <- subset_df$IncLevelDifference
    regulation <- if (all(inc_diff < 0)) {
      "Excluded"
    } else if (all(inc_diff > 0)) {
      "Included"
    } else {
      "Ambiguous"
    }
    
  
    exon_target_concat <- paste(paste0(subset_df$start, "-", subset_df$end), collapse = ",")
    
    data.frame(
      exon_range = paste0(min(subset_df$start), "-", max(subset_df$end)),
      geneSymbol = unique(subset_df$geneSymbol),
      PValue = paste(subset_df$PValue, collapse = ","),
      FDR = paste(subset_df$FDR, collapse = ","),
      min_FDR = min(subset_df$FDR),
      IncLevelDifference = paste(subset_df$IncLevelDifference, collapse = ","),
      max_dPSI = ifelse(inc_diff > 0, max(inc_diff), min(inc_diff)),
      mean_dPSI = mean(inc_diff),
      Regulation = regulation,
      PSI_1 = paste(subset_df$PSI_1, collapse = ","),
      PSI_2 = paste(subset_df$PSI_2, collapse = ","),
      Chr = unique(subset_df$Chr),
      Strand = unique(subset_df$Strand),
      exon_upstream = paste(subset_df$exon_upstream, collapse = ","),
      exon_downstream = paste(subset_df$exon_downstream, collapse = ","),
      exon_target = exon_target_concat,
      stringsAsFactors = FALSE
    )
  }

  simplified_df <- do.call(rbind, lapply(grouped, aggregate_rows))
  return(simplified_df)
}

# The following funciton is just a quick way to rename and replace a column 
library(dplyr)
library(tidyr)
Split_rename_column_exon_range2 <- function(x) {
x <- within(x, exon_range<-data.frame(do.call('rbind', strsplit(as.character(exon_range), '-', fixed=TRUE))))
write.table(x, file="x", quote=F, row.names=F, sep="\t")
x <- read.table("x", header=TRUE, sep="\t")
colnames(x)[3] <- "start"
colnames(x)[4] <- "end"
return(x)
}

# This is the pivotal point, this function will generate Regoins of Interests, by adding flanking introns to you list of exons.
process_single_exon_parallel <- function(i, my_exons_gr_global, all_gff_exons_global, gene_name_col_in_gff_global) {
# This function is executed independently by each parallel worker/core.
# It receives an index 'i' and retrieves the specific target exon from the global GRanges object.
current_exon_gr_single <- my_exons_gr_global[i] # Access the specific exon by its index
current_gene_name <- mcols(current_exon_gr_single)$gene_name
current_strand_char <- as.character(strand(current_exon_gr_single))
current_chr_char <- as.character(seqnames(current_exon_gr_single))
# Initialize the results vector for this exon (defaults to NA if no flanking introns are mapped)
res <- c(upstream_intron_start = NA_integer_,
upstream_intron_end = NA_integer_,
downstream_intron_start = NA_integer_,
downstream_intron_end = NA_integer_)
# Filter global GFF exons to match the current chromosome, strand, and gene name (coerced to characters)
gene_specific_gff_exons <- all_gff_exons_global[
as.character(seqnames(all_gff_exons_global)) == current_chr_char &
as.character(strand(all_gff_exons_global)) == current_strand_char &
as.character(mcols(all_gff_exons_global)[[gene_name_col_in_gff_global]]) == as.character(current_gene_name)
]
if (length(gene_specific_gff_exons) == 0) {
warning_msg <- paste0(
"Exon #", i, " (", as.character(current_exon_gr_single), ") : ",
"No exon found for this gene in the given GFF '", current_gene_name, "' ",
"on ", current_chr_char, ":", current_strand_char, ". This exon will be ignored.\n"
)
gff_genes_on_chr_strand <- unique(as.character(mcols(all_gff_exons_global)[[gene_name_col_in_gff_global]][
as.character(seqnames(all_gff_exons_global)) == current_chr_char &
as.character(strand(all_gff_exons_global)) == current_strand_char
]))
if (length(gff_genes_on_chr_strand) > 0) {
warning_msg <- paste0(warning_msg,
" Gene found in the GFF for ", current_chr_char, ":", current_strand_char, " : ",
paste(head(gff_genes_on_chr_strand, 10), collapse = ", "),
ifelse(length(gff_genes_on_chr_strand) > 10, ", ...", ""), "\n"
)
} else {
warning_msg <- paste0(warning_msg,
"  WARNING : No exon found for this gene in the given GFF for ", current_chr_char, ":", current_strand_char, " at all. \n"
)
}
warning(warning_msg)
return(res)
}
gene_specific_gff_exons_sorted <- sort(gene_specific_gff_exons, ignore.strand = FALSE)
upstream_flanking_exon <- GRanges()
downstream_flanking_exon <- GRanges()
if (current_strand_char == "+" || current_strand_char == "-") {
overlapping_gff_exons_idx <- GenomicRanges::findOverlaps(current_exon_gr_single, gene_specific_gff_exons_sorted, type = "any", ignore.strand = FALSE)
is_overlapping <- rep(FALSE, length(gene_specific_gff_exons_sorted))
if(length(overlapping_gff_exons_idx) > 0) {
is_overlapping[subjectHits(overlapping_gff_exons_idx)] <- TRUE
}
non_overlapping_gff_exons <- gene_specific_gff_exons_sorted[!is_overlapping]
if (length(non_overlapping_gff_exons) > 0) {
if (current_strand_char == "+") {
potential_upstream_candidates <- non_overlapping_gff_exons[end(non_overlapping_gff_exons) < start(current_exon_gr_single)]
potential_downstream_candidates <- non_overlapping_gff_exons[start(non_overlapping_gff_exons) > end(current_exon_gr_single)]
if (length(potential_upstream_candidates) > 0) {
upstream_flanking_exon <- potential_upstream_candidates[which.max(end(potential_upstream_candidates))]
}
if (length(potential_downstream_candidates) > 0) {
downstream_flanking_exon <- potential_downstream_candidates[which.min(start(potential_downstream_candidates))]
}
} else { # current_strand_char == "-"
potential_upstream_candidates <- non_overlapping_gff_exons[start(non_overlapping_gff_exons) > end(current_exon_gr_single)]
potential_downstream_candidates <- non_overlapping_gff_exons[end(non_overlapping_gff_exons) < start(current_exon_gr_single)]
if (length(potential_upstream_candidates) > 0) {
upstream_flanking_exon <- potential_upstream_candidates[which.min(start(potential_upstream_candidates))]
}
if (length(potential_downstream_candidates) > 0) {
downstream_flanking_exon <- potential_downstream_candidates[which.max(end(potential_downstream_candidates))]
}
}
}
} else {
warning(paste0("Unkown strand for exon : ", as.character(current_exon_gr_single), ". Unable to determine flanking introns."))
return(res)
}
# Calculate genomic coordinates for the flanking introns based on strand directionality
if (length(upstream_flanking_exon) > 0) {
if (current_strand_char == "+") {
res["upstream_intron_start"] <- end(upstream_flanking_exon) + 1
res["upstream_intron_end"] <- start(current_exon_gr_single) - 1
} else { # brin '-'
res["upstream_intron_start"] <- end(current_exon_gr_single) + 1
res["upstream_intron_end"] <- start(upstream_flanking_exon) - 1
}
if (res["upstream_intron_start"] > res["upstream_intron_end"]) {
res["upstream_intron_start"] <- NA_integer_
res["upstream_intron_end"] <- NA_integer_
}
}
if (length(downstream_flanking_exon) > 0) {
if (current_strand_char == "+") {
res["downstream_intron_start"] <- end(current_exon_gr_single) + 1
res["downstream_intron_end"] <- start(downstream_flanking_exon) - 1
} else { # brin '-'
res["downstream_intron_start"] <- end(downstream_flanking_exon) + 1
res["downstream_intron_end"] <- start(current_exon_gr_single) - 1
}
if (res["downstream_intron_start"] > res["downstream_intron_end"]) {
res["downstream_intron_start"] <- NA_integer_
res["downstream_intron_end"] <- NA_integer_
}
}
return(res)
}


coord_prog_AvT <- prog_AvT_SE_simplified3[,c(8,9,1,5,4,12,2)]
coord_prog_AvT <- Split_rename_column_exon_range2(coord_prog_AvT)
my_exons_df <- coord_prog_AvT
colnames(my_exons_df)[8] <- "gene_name"

library(GenomicRanges)
library(rtracklayer) # For importing and handling GFF/GTF annotation files
library(dplyr)       # For data manipulation (optional but useful)
library(BiocParallel) # For parallel computing execution across multiple cores

my_exons_gr <- makeGRangesFromDataFrame(my_exons_df, keep.extra.columns = TRUE)
seqlevelsStyle(my_exons_gr) <- "UCSC"
common_seqlevels <- unique(c(seqlevels(gff_annotation), seqlevels(my_exons_gr)))
results_list <- bplapply(seq_along(my_exons_gr), # Iterate across all exon indices in parallel
FUN = process_single_exon_parallel,
my_exons_gr_global = my_exons_gr, # Pass the GRanges object explicitly as a global reference to workers
all_gff_exons_global = all_gff_exons,
gene_name_col_in_gff_global = gene_name_col_in_gff,
BPPARAM = param)
results_matrix <- do.call(rbind, results_list)
my_exons_df[, c("upstream_intron_start", "upstream_intron_end",
"downstream_intron_start", "downstream_intron_end")] <- results_matrix
View(my_exons_df)
my_exons_df$roi <- ifelse(my_exons_df$Strand == "+", paste(my_exons_df$upstream_intron_start, my_exons_df$downstream_intron_end, sep="-"), paste(my_exons_df$downstream_intron_start, my_exons_df$upstream_intron_end, sep="-"))
write.table(my_exons_df, file="~/Bureau/NEw_KO2/AvT_exons_with_flanking_introns", quote=F, row.names = F, sep="\t")

## ROI is now created, next step: cross with MicroC bedpe files (a list of bedpe files)

find_roi_loop_overlaps_multi <- function(roi_df, dir_list) {
    all_results <- list()
    
    for (cond in names(dir_list)) {
        dir_df <- dir_list[[cond]] %>% as.data.frame()
        
        # create loop_id if non exists
        if (!"loop_id" %in% colnames(dir_df)) {
            dir_df$loop_id <- paste0("loop_", cond, "_", seq_len(nrow(dir_df)))
        }
        
        # safe coercions & name assumptions pour DIR table
        # colonnes attendues : Chr1, start1, end1, chr2, start2, end2, loop_id
        if (!all(c("Chr1","start1","end1","chr2","start2","end2","loop_id") %in% colnames(dir_df))) {
            stop(paste("dir_df for", cond, "is missing expected columns:",
                       paste(setdiff(c("Chr1","start1","end1","chr2","start2","end2","loop_id"),
                                     colnames(dir_df)), collapse=", ")))
        }
        
        # --- 1. prep ROIs ---
        roi_df_tmp <- roi_df %>%
            mutate(
                roi_chr = as.character(Chr),
                roi_start = as.integer(gsub("-.*", "", as.character(roi))),
                roi_end   = as.integer(gsub(".*-", "", as.character(roi)))
            )
        
        gr_rois <- GRanges(
            seqnames = roi_df_tmp$roi_chr,
            ranges = IRanges(start = roi_df_tmp$roi_start, end = roi_df_tmp$roi_end),
            geneSymbol = roi_df_tmp$gene_name,
            roi = roi_df_tmp$roi,
            Regulation = roi_df_tmp$Regulation
        )
        
        # --- 2. prep anchors ---
        gr_anchor1 <- GRanges(
            seqnames = as.character(dir_df$Chr1),
            ranges = IRanges(start = as.integer(dir_df$start1), end = as.integer(dir_df$end1)),
            anchor_side = "anchor1",
            loop_id = as.character(dir_df$loop_id)
        )
        
        gr_anchor2 <- GRanges(
            seqnames = as.character(dir_df$chr2),
            ranges = IRanges(start = as.integer(dir_df$start2), end = as.integer(dir_df$end2)),
            anchor_side = "anchor2",
            loop_id = as.character(dir_df$loop_id)
        )
        
        gr_anchors <- c(gr_anchor1, gr_anchor2)
        
        # --- 3. Overlaps ---
        hits <- findOverlaps(gr_rois, gr_anchors)
        if (length(hits) == 0) next
        
        res_df <- data.frame(
            roi_idx = queryHits(hits),
            anchor_idx = subjectHits(hits),
            stringsAsFactors = FALSE
        ) %>%
            mutate(
                geneSymbol = mcols(gr_rois)$geneSymbol[roi_idx],
                roi = mcols(gr_rois)$roi[roi_idx],
                roi_chr = as.character(seqnames(gr_rois)[roi_idx]),
                roi_start = start(gr_rois)[roi_idx],
                roi_end = end(gr_rois)[roi_idx],
                Regulation = mcols(gr_rois)$Regulation[roi_idx],
                anchor_chr = as.character(seqnames(gr_anchors)[anchor_idx]),
                anchor_start = start(gr_anchors)[anchor_idx],
                anchor_end = end(gr_anchors)[anchor_idx],
                anchor_side = mcols(gr_anchors)$anchor_side[anchor_idx],
                loop_id = mcols(gr_anchors)$loop_id[anchor_idx],
                Condition = cond,
                stringsAsFactors = FALSE
            )
        
        # --- 4. add DIR ---
        loop_info <- dir_df %>%
            dplyr::select(loop_id, Chr1, start1, end1, chr2, start2, end2)
        
        res_df <- res_df %>%
            left_join(loop_info, by = "loop_id") %>%
            rowwise() %>%
            mutate(
                distal_chr = ifelse(anchor_side == "anchor1", chr2, Chr1),
                distal_start = ifelse(anchor_side == "anchor1", start2, start1),
                distal_end = ifelse(anchor_side == "anchor1", end2, end1)
            ) %>%
            ungroup()
        
        # --- 5. Final Format  ---
        final_df <- res_df %>%
            dplyr::select(
                geneSymbol, roi, roi_chr, roi_start, roi_end,
                Regulation, Condition,
                anchor_side, anchor_chr, anchor_start, anchor_end,
                distal_chr, distal_start, distal_end, loop_id
            )
        
        all_results[[cond]] <- final_df
    }
    
    # Combine conditions
    bind_rows(all_results)
}

#Once you have you list of DIRs connected to your ROIs, it's time to look for chromatin features in these DIRs:


annotate_overlap2 <- function(region_df, bed_list, tolerance_bp = 0) {
    stopifnot(
        all(c("distal_chr", "distal_start", "distal_end") %in% colnames(region_df))
    )
    region_gr <- GRanges(
        seqnames = region_df$distal_chr,
        ranges = IRanges(
            start = region_df$distal_start,
            end   = region_df$distal_end
        )
    )
    for (i in seq_along(bed_list)) {
        bed <- bed_list[[i]]
        feat_name <- names(bed_list)[i]
        bed_gr <- GRanges(
            seqnames = bed[[1]],
            ranges = IRanges(start = bed[[2]], end = bed[[3]])
        )
        hits <- findOverlaps(region_gr, bed_gr, maxgap = tolerance_bp)
        overlap_flag <- logical(length(region_gr))
        overlap_flag[unique(queryHits(hits))] <- TRUE
        region_df[[paste0("overlap_", feat_name)]] <- overlap_flag
    }
    region_df
}


dir_list <- list(
MCF10A   = DIR_MCF10A,
MCF10AT1 = DIR_MCF10AT1,
MCF10CA1 = DIR_MCF10CA1
) # DIR_MCF10A and co are the bedpe files loaded as table in R environment.

prog_AvT_ROI <- read.table("AvT_exons_with_flanking_introns", header = T)
prog_AvT_ROI_loops_all_new <- find_roi_loop_overlaps_multi(prog_AvT_ROI, dir_list)
prog_AvT_ROI_loops_all_new$roi_id <- paste(prog_AvT_ROI_loops_all_new$roi_chr, prog_AvT_ROI_loops_all_new$roi_start, prog_AvT_ROI_loops_all_new$roi_end, sep = "_")
DIR_prog_for_alluvial_AvT <- prog_AvT_ROI_loops_all_new[,c(16,12,13,14,7,1,6)]
DIR_prog_for_alluvial_AvT <- unique(DIR_prog_for_alluvial_AvT)
DIR_prog_for_alluvial_AvT$Context <- "AvT"

library(readr)
library(dplyr)
library(purrr)
bed_files3 <- list.files(path = "/home/kevin/Bureau/NEw_KO2/Chip_progression/ChIP_to_overloop/DIFF/", pattern = "*.bed", full.names = TRUE)
bed_list3 <- map(bed_files3, ~ read_tsv(.x, col_names = FALSE))
names(bed_list3) <- tools::file_path_sans_ext(basename(bed_files3))
df2 <- annotate_overlap2(data.frame(df), bed_list3) # df is the rbind of several objects like DIR_prog_for_alluvial_AvT

