library(bedtoolsr)
# library(UpSetR)
library(stringr) 
library(reshape) 
library(ggplot2) 
library(tidyverse) 
library(RColorBrewer) 
library(ggpubr) 
library(ggsci) 

args <- commandArgs(TRUE)

big_table_file <- args[1]
filters_file <- args[2]
library <- args[3]
plots_path <- args[4]
classified_peaks <- args[5]
summary <- args[6]
error_prone <- args[7]
truth_set <- args[8]

# Read input and define columns
canonical_chrs <- c(
  "chr1","chr2","chr3","chr4","chr5","chr6",
  "chr7","chr8","chr9","chr10","chr11","chr12",
  "chr13","chr14","chr15","chr16","chr17","chr18",
  "chr19","chr20","chr21","chr22","chrX","chrY"
)

big_table <- read.table(big_table_file, sep="\t", header=TRUE)

colnames(big_table) <- c(
  "chrm","start","end","peak","shape","strand",
  "num_reads","num_gmotif_reads","num_polyA_reads",
  "peak_width","gmotif_percent","polyA_percent",
  "num_unique_reads","max_usp_distance","num_templates","num_endpoints",
  "template_ratio","start_end_ratio","unique_read_ratio","RPM","TPM",
  "nearest_peak","RPM_nearest","nearest_RPM_ratio",
  "satellite","segdup","homopolymer",
  "nearest_KR","nearest_KNR","nearest_KR_dist","nearest_KNR_dist"
)

# Ignore peaks in non-canonical chromosomes
big_table <- big_table[(big_table$chrm %in% canonical_chrs),]

#### BENCHMARKING ####
# Intersect peaks with truth set
true_peak_IDs = list()
if (truth_set != "-NA-") {
  big_bed <- big_table[, c("chrm","start","end","peak","RPM","strand")]

  big_bed[big_bed$strand == "+", "start"] <- big_bed[big_bed$strand == "+", "end"] # 3' end
  big_bed[big_bed$strand == "-", "end"] <- big_bed[big_bed$strand == "-", "start"]

  true_bed <- read_tsv(
    truth_set, 
    col_names=c("chrm","start","end","name","score","strand"), 
    col_types="ciicdc"
  )

  true_bed_slop <- bt.slop(i=true_bed, g="hg38", b=50)
  true_peak_list <- bt.intersect(a=big_bed, b=true_bed_slop, wo=TRUE, S=TRUE)[, c("V4","V10")]
  colnames(true_peak_list) <- c("peak","true_insertion_ID")
  true_peak_IDs = true_peak_list$peak
  
  # Add true_insertion_ID column
  big_table <- merge(big_table, true_peak_list, by="peak", all.x=TRUE) 
  big_table$true_insertion_ID[is.na(big_table$true_insertion_ID)] <- "-NA-"
}

#### CLASSIFICATION AND FILTERING ####
summary_file <- file(summary, open='a')
cat(paste0("########## Filtering results ", Sys.time(), " ##########\n"), file=summary_file, sep="")
cat(paste0("Total peaks: ", nrow(big_table)), file=summary_file, sep="\n")

# Annotate peaks overlapping centromere or SegDup and peaks without Gmotif
cat("\nANNOTATIONS", file=summary_file, sep="\n")
cat(
  paste0("Total peaks overlapping ALR/Alpha satellites (centromeres): ", 
    nrow(big_table[grepl("ALR/Alpha", big_table$satellite),])), 
  file=summary_file, 
  sep="\n"
)
cat(
  paste0("Total peaks overlapping SegDups: ", 
    nrow(big_table[big_table$segdup != '.',])), 
  file=summary_file, 
  sep="\n"
)
cat(
  paste0("Total peaks overlapping homopolymers: ", 
    nrow(big_table[big_table$homopolymer != '.',])), 
  file=summary_file, 
  sep="\n"
)
cat(
  paste0("Total peaks with Gmotif % less than 50%: ", 
    nrow(big_table[big_table$gmotif_percent < 0.5,])), 
  file=summary_file, 
  sep="\n"
)

# Initialize columns for filtering
big_table$classification <- "Candidate"
big_table$filter <- "-NA-"
big_table$candidate <- TRUE
big_table$FP <- FALSE

# Filters are successive (FP peak will be labeled by the first filter it fails)
cat("\nFILTERS", file=summary_file, sep="\n")

# Global filters (all peaks)
if (library != "bulk") {
  cat("\n\tAll peaks", file=summary_file, sep="\n")

  # In PTA libraries, remove smaller peaks nearby larger ones (amplification artifacts)
  big_table$filter[!(big_table$FP) & (big_table$nearest_RPM_ratio < 1)] <- "nearby_peak"
  big_table$FP[big_table$filter == "nearby_peak"] <- TRUE
  cat(
    paste0("\t\tNearby larger peak: ", 
      nrow(big_table[big_table$filter == "nearby_peak",])), 
    file=summary_file, 
    sep="\n"
  )
}

# Label known reference (KR) peaks - priority
# cat("\n\tKR/KNR peaks", file=summary_file, sep="\n")
big_table$classification[(grepl("L1HS", big_table$nearest_KR) | 
                          grepl("L1Hs", big_table$nearest_KR)) &
                          !(big_table$FP)] <- "KR"
big_table$filter[big_table$classification == "KR"] <- "PASS"
big_table$filter[(big_table$classification == "KR") & 
                 ((big_table$peak_width < 150) | (big_table$RPM < 100))] <- "lowConf"
# cat(
#   paste0("\t\tLow-confidence KRs: ", 
#     nrow(big_table[(big_table$classification == "KR") & (big_table$filter == "lowConf"),])), 
#   file=summary_file, 
#   sep="\n"
# )

# Label Non-specific peaks
big_table$classification[grepl("L1PA2|L1PA3|L1PA4|L1PA5", big_table$nearest_KR) & 
                         (big_table$nearest_KR_dist <= big_table$nearest_KNR_dist) &
                         !(big_table$classification == "KR") & 
                         !(big_table$FP)] <- "Non-specific"
big_table$filter[big_table$classification == "Non-specific"] <- "PASS"
big_table$filter[(big_table$classification == "Non-specific") & 
                 ((big_table$peak_width < 150) | (big_table$RPM < 100))] <- "lowConf" 

# Label known non-reference (KNR) peaks
# Also, artifically remove KNR labels from peaks in truth set for benchmarking
big_table$classification[(grepl("LINE1", big_table$nearest_KNR) | 
                         grepl("L1", big_table$nearest_KNR)) & 
                         (big_table$nearest_KNR_dist < big_table$nearest_KR_dist) &
                         !(big_table$peak %in% true_peak_IDs) &
                         !(big_table$classification == "KR") &
                         !(big_table$FP)] <- "KNR"
big_table$filter[big_table$classification == "KNR"] <- "PASS"
big_table$filter[(big_table$classification == "KNR") & 
                 ((big_table$peak_width < 150) | 
                 (big_table$unique_read_ratio < 0.1) | 
                 (big_table$polyA_percent == 0) | 
                 (big_table$RPM < 50))] <- "lowConf"
# big_table$filter[(big_table$classification == "KNR") & 
#                  ((big_table$peak_width < 150) & 
#                  (big_table$unique_read_ratio < 0.1) & 
#                  (big_table$polyA_percent == 0) & 
#                  (big_table$RPM < 50))] <- "noise"
# cat(
#   paste0("\t\tLow-confidence KNRs: ", 
#     nrow(big_table[(big_table$classification == "KNR") & (big_table$filter == "lowConf"),])), 
#   file=summary_file, 
#   sep="\n"
# ) 

# Separate category for noise near known elements
if (library == "single") {RPM_threshold <- 5} else {RPM_threshold <- 1}
big_table$filter[(big_table$classification != "Candidate") & (big_table$RPM < RPM_threshold)] <- "noise"

# Filter remaining peaks for FPs
big_table$candidate[big_table$classification != "Candidate"] <- FALSE
cat("\n\tCandidate peaks", file=summary_file, sep="\n")
cat(paste0("\t\tTotal peaks: ", nrow(big_table[big_table$candidate,])), file=summary_file, sep="\n")

if (library != "bulk") {
  # Require multi-template support for PTA libraries
  if (library == "single") {
    big_table$filter[(big_table$candidate) & !(big_table$FP) & (big_table$TPM < 0.25)] <- "templates"
  } else if (library == "micro") {
    big_table$filter[(big_table$candidate) & !(big_table$FP) & (big_table$num_templates < 2)] <- "templates"
  }

  big_table$FP[big_table$filter == "templates"] <- TRUE
  cat(
    paste0("\t\tToo few templates: ", 
      nrow(big_table[big_table$filter == "templates",])), 
    file=summary_file, 
    sep="\n"
  )

  big_table$filter[(big_table$candidate) & !(big_table$FP) & (big_table$template_ratio < 0.2)] <- "templateRatio"
  big_table$FP[big_table$filter == "templateRatio"] <- TRUE
  cat(
    paste0("\t\tTemplate ratio: ",
      nrow(big_table[big_table$filter == "templateRatio",])), 
    file=summary_file, 
    sep="\n"
  )
}

big_table$filter[(big_table$candidate) & !(big_table$FP) & (big_table$unique_read_ratio < 0.1)] <- "uniqueReads"
big_table$FP[big_table$filter == "uniqueReads"] <- TRUE
cat(
  paste0("\t\tFew unique reads: ", 
    nrow(big_table[big_table$filter == "uniqueReads",])), 
  file=summary_file, 
  sep="\n"
)

if (library == "single") {
  big_table$filter[(big_table$candidate) & !(big_table$FP) & (big_table$polyA_percent < 20)] <- "polyApercent"
} else if (library == "micro") {
  big_table$filter[(big_table$candidate) & !(big_table$FP) & (big_table$polyA_percent < 10)] <- "polyApercent"
} else {
  big_table$filter[(big_table$candidate) & !(big_table$FP) & (big_table$polyA_percent == 0)] <- "polyApercent"
}

big_table$FP[big_table$filter == "polyApercent"] <- TRUE
cat(
  paste0("\t\tNo polyA containing reads (junction spanning): ", 
    nrow(big_table[big_table$filter == "polyApercent",])), 
  file=summary_file, 
  sep="\n"
)

big_table$filter[(big_table$candidate) & !(big_table$FP) & (big_table$RPM < RPM_threshold)] <-"RPM"
big_table$FP[big_table$filter == "RPM"] <- TRUE
cat(
  paste0("\t\tLow RPM: ", 
    nrow(big_table[(big_table$filter == "RPM"),])), 
  file=summary_file, 
  sep="\n"
)

big_table$classification[big_table$FP] <- "FP"
big_table$candidate[big_table$classification != "Candidate"] <- FALSE
big_table$filter[big_table$candidate] <- "PASS"

# Label candidate peaks as UNK/SOM
if (library == "bulk") {
  big_table$classification[big_table$classification == "Candidate" & big_table$RPM >= 100 & big_table$TPM >= 4] <- "UNK"
  big_table$classification[big_table$classification == "Candidate" & ((big_table$num_templates > 1) & (big_table$max_usp_distance >= 10))] <- "SOM_clonal"
  big_table$classification[big_table$classification == "Candidate" & ((big_table$num_templates == 1) | (big_table$max_usp_distance < 10))] <- "SOM_private"
} else if (library == "micro") {
  big_table$classification[big_table$classification == "Candidate" & big_table$RPM >= 100 & big_table$TPM >= 4] <- "UNK"
  big_table$classification[big_table$classification == "Candidate"] <- "SOM"
} else {
  big_table$classification[big_table$classification == "Candidate"] <- "SOM" # need multiple cells to distinguish
}

# Filter out SOM peaks in error-prone regions if regions are provided
if (error_prone != "-NA-") {
  big_bed <- big_table[, c("chrm","start","end","peak","RPM","strand")]
  error_bed <- read_tsv(
    error_prone, 
    col_select=c(1,2,3), 
    col_names=c("chrm","start","end"), 
    col_types="cii"
  )
    
  error_peak_list <- bt.intersect(a=big_bed, b=error_bed, wo=TRUE)[, "V4"]
  
  big_table$filter[(big_table$peak %in% error_peak_list) & 
                          (grepl("SOM", big_table$classification)) & 
                          !(big_table$FP)] <- "errorProne"
  big_table$FP[big_table$filter == "errorProne"] <- TRUE
  big_table$classification[big_table$FP] <- "FP"
  cat(
    paste0("\t\tOverlapping error-prone region: ", 
      nrow(big_table[big_table$filter == "errorProne",])), 
    file=summary_file, 
    sep="\n"
  )
}

big_table$candidate <- NULL # drop column
big_table$FP <- NULL 

cat(
  paste0("\nTotal KR: ", 
    nrow(big_table[big_table$classification == "KR",])), 
  file=summary_file, 
  sep="\n"
)
cat(
  paste0("Total KNR: ", 
    nrow(big_table[big_table$classification == "KNR",])), 
    file=summary_file, 
    sep="\n"
)
cat(
  paste0("Total Non-specific: ", 
    nrow(big_table[big_table$classification == "Non-specific",])), 
  file=summary_file, 
  sep="\n"
)
cat(
  paste0("Total UNK: ", 
    nrow(big_table[big_table$classification == "UNK",])), 
  file=summary_file, 
  sep="\n"
)
if (library == "bulk") {
  cat(
    paste0("Total SOM_clonal: ", 
      nrow(big_table[big_table$classification == "SOM_clonal",])), 
    file=summary_file, 
    sep="\n"
  )
  cat(
    paste0("Total SOM_private: ", 
      nrow(big_table[big_table$classification == "SOM_private",])), 
    file=summary_file, 
    sep="\n"
  )
} else {
  cat(
    paste0("Total SOM: ", 
      nrow(big_table[big_table$classification == "SOM",])), file=summary_file, 
    sep="\n"
  )
}
cat(
  paste0("Total FP: ", 
    nrow(big_table[big_table$classification == "FP",])), 
    file=summary_file, 
    sep="\n"
  )

#### Plotting ####
filters <- ggplot(big_table[big_table$classification == "FP",], aes(x=filter, group=classification, color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~classification) +
  ggtitle("Number of FP Peaks per Filter") + 
  ylab("Peak Count") + 
  xlab("Filter")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 

pdf(plots_path)
filters
dev.off()

# Candidate table
bed_cols = c(
  "chrm","start","end","peak","RPM","strand",
  "classification","filter"
)

if (truth_set != "-NA-") {
  peak_table <- big_table[, c(bed_cols, "true_insertion_ID")]
} else {
  peak_table <- big_table[, bed_cols]
}

peak_table <- peak_table %>% arrange(chrm, start)
colnames(peak_table)[1] <- "#chrm" # for bedtools compatibility
write.table(
  peak_table, 
  classified_peaks, 
  sep="\t", 
  row.names=FALSE, 
  col.names=TRUE, 
  quote=FALSE, 
  na="-NA-"
)
write.table(
  big_table[, lapply(big_table, class) != "list"], 
  filters_file, 
  sep="\t", 
  row.names=FALSE, 
  col.names=TRUE, 
  quote=FALSE, 
  na="-NA-"
)
