library(bedtoolsr)
library(UpSetR)
library(stringr) 
library(reshape) 
library(ggplot2) 
library(tidyverse) 
library(RColorBrewer) 
library(ggpubr) 
library(ggsci) 

args <- commandArgs(TRUE)

big_table_file <- args[1]
filter_reasons_file <- args[2]
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
  "reads","Ntag","gmotif","polyA","breakpoint","bp_chimera_len","bp_chimera_ratio",
  "repeatmasker","evrony","homopolymers","i1kgp","gnomad","nyuwa",
  "xtea","hgsvc3","melt_lra","ont","bamreads","peakreads","uniqreads",
  "usp","unique_read_ratio","RPM","nearest_peak",
  "SegDups","max_distance"
)

# Define more columns
big_table$peak_width <- big_table$end - big_table$start
big_table$gmotif_percent <- big_table$gmotif / big_table$reads

big_table <- merge(
  big_table, 
  big_table[, c("RPM","peak")], 
  by.x="nearest_peak", 
  by.y="peak", 
  suffixes = c("","_nearest"), 
  all.x=TRUE
)
big_table$nearest_RPM_ratio <- big_table$RPM / big_table$RPM_nearest
big_table$nearest_RPM_ratio[is.na(big_table$nearest_RPM_ratio)] <- 1
big_table$RPM_log <- log(big_table$RPM)

big_table$number_templates <- as.numeric(
  unlist(lapply(strsplit(big_table$usp, ";"), `[[`, 1))
)
big_table$TPM <- (big_table$number_templates / big_table$bamreads) * 1000000
big_table$template_ratio <- NA
template_list <- lapply(
  strsplit(big_table$usp[big_table$number_templates > 1], ";"), `[[`, 2)
big_table$template_ratio[big_table$number_templates > 1] <- 
  as.numeric(lapply(strsplit(unlist(template_list), ","), `[[`, 2)) / 
  as.numeric(lapply(strsplit(unlist(template_list), ","), `[[`, 1))
big_table$polyA_percent <- big_table$polyA / big_table$reads

# Ignore peaks in non-canonical chromosomes
big_table <- big_table[(big_table$chrm %in% canonical_chrs),]

# If benchmarking, intersect peaks with truth set
true_peak_IDs = list()
if (truth_set != "-NA-") {
  big_bed <- big_table[, c("chrm","start","end","peak","RPM","strand")]
  true_bed <- read_tsv(
    truth_set, 
    col_names=c("chrm","start","end","name","score","strand"), 
    col_types="ciicdc"
  )

  true_bed_slop <- bt.slop(i=true_bed, g="hg38", b=50)
  true_peak_list <- bt.intersect(
    a=big_bed, 
    b=true_bed_slop, 
    wo=TRUE,
    S=TRUE)[, c("V4","V10")]
  colnames(true_peak_list) <- c("peak","true_insertion_ID")
  true_peak_IDs = true_peak_list$peak
  
  # Add true_insertion_ID column
  big_table <- merge(big_table, true_peak_list, by="peak", all.x=TRUE) 
  big_table$true_insertion_ID[is.na(big_table$true_insertion_ID)] <- "-NA-"
}

#### Classification and Filtering ####
summary_file <- file(summary, open='a')
cat(paste0("#################### Filtering results ", Sys.time(), " ####################\n"), file=summary_file, sep="")
cat(paste0("Total peaks: ", nrow(big_table)), file=summary_file, sep="\n")

# Annotate peaks overlapping centromere or SegDup and peaks without Gmotif
cat("\nANNOTATIONS", file=summary_file, sep="\n")
cat(
  paste0("Total peaks overlapping ALR/Alpha satellites (centromeres): ", 
    nrow(big_table[grepl("ALR/Alpha", big_table$repeatmasker),])), 
  file=summary_file, 
  sep="\n"
)
cat(
  paste0("Total peaks overlapping SegDups: ", 
    nrow(big_table[big_table$SegDups != '.',])), 
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
big_table$filter_reason <- "-NA-"
big_table$candidate <- TRUE
big_table$FP <- FALSE

# Filters are successive
# An FP peak will be labeled by the first filter that catches it
if (library != "bulk") {
  cat("\nFILTERS", file=summary_file, sep="\n")
  cat("\n\tAll peaks", file=summary_file, sep="\n")

  # In PTA libraries, remove smaller peaks nearby larger ones
  # that are due to incorrect amplification 
  big_table$filter_reason[!(big_table$FP) & (big_table$nearest_RPM_ratio < 1)] <- "nearby_peak"
  big_table$FP[big_table$filter_reason == "nearby_peak"] <- TRUE
  cat(
    paste0("\t\tNearby larger peak: ", 
      nrow(big_table[big_table$filter_reason == "nearby_peak",])), 
    file=summary_file, 
    sep="\n"
  )
}

# Label known reference (KR) peaks
cat("\n\tKR/KNR peaks", file=summary_file, sep="\n")
big_table$classification[!(big_table$FP) & 
                         (grepl("L1HS", big_table$repeatmasker) | 
                          grepl("L1Hs", big_table$evrony))] <- "KR"
big_table$filter_reason[big_table$classification == "KR" & big_table$RPM < 50] <- "KR_low_RPM" 
cat(
  paste0("\t\tKR RPM less than 50: ", 
    nrow(big_table[big_table$filter_reason == "KR_low_RPM",])), 
  file=summary_file, 
  sep="\n"
)

# Label known non-reference (KNR) peaks
# Also, artifically remove KNR labels from peaks in truth set for benchmarking
big_table$classification[(grepl("INS:ME:LINE1", big_table$i1kgp) | 
                         grepl("INS:ME:LINE1", big_table$gnomad) | 
                         grepl("INS:ME:LINE1", big_table$nyuwa) | 
                         grepl("INS:ME:LINE1", big_table$xtea) | 
                         grepl("LINE/L1", big_table$hgsvc3) | 
                         grepl("LINE1", big_table$melt_lra) | 
                         grepl("L1-INS", big_table$ont)) & 
                         !(big_table$peak %in% true_peak_IDs) &
                         !(big_table$FP) &
                         !(big_table$classification == "KR")] <- "KNR"
big_table$filter_reason[big_table$classification == "KNR" & big_table$RPM < 50] <- "KNR_low_RPM"
cat(
  paste0("\t\tKNR RPM less than 50: ", 
    nrow(big_table[big_table$filter_reason == "KNR_low_RPM",])), 
  file=summary_file, 
  sep="\n"
)

# Label Off-target peaks
big_table$classification[grepl("L1PA2|L1PA3|L1PA4|L1PA5|L1PA6|L1PA7", big_table$repeatmasker) & 
                        !(big_table$classification == "KR") & 
                        !(big_table$classification == "KNR") &
                        !(big_table$FP)] <- "Off-target"
big_table$filter_reason[big_table$classification == "Off-target"] <- "Off-target"

# Filter remaining peaks for FPs
big_table$candidate[big_table$classification != "Candidate"] <- FALSE
cat("\n\tCandidate peaks", file=summary_file, sep="\n")

# big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & 
#                         ((big_table$bp_chimera_len > 12) | (big_table$bp_chimera_ratio >= 0.8))] <- "chimera"
# big_table$FP[big_table$filter_reason == "chimera"] <- TRUE
# cat(
#   paste0("\t\tChimera: ", 
#     nrow(big_table[big_table$filter_reason == "chimera",])), 
#   file=summary_file, 
#   sep="\n"
# )

if (library != "bulk") {
  # Require multi-template support for PTA libraries
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$number_templates < 2)] <- "templates"
  big_table$FP[big_table$filter_reason == "templates"] <- TRUE
  cat(
    paste0("\t\tToo few templates: ", 
      nrow(big_table[big_table$filter_reason == "templates",])), 
    file=summary_file, 
    sep="\n"
  )

  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$template_ratio < 0.25)] <- "template_ratio"
  big_table$FP[big_table$filter_reason == "template_ratio"] <- TRUE
  cat(
    paste0("\t\tTemplate ratio: ",
      nrow(big_table[big_table$filter_reason == "template_ratio",])), 
    file=summary_file, 
    sep="\n"
  )
}

big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$unique_read_ratio < 0.1)] <- "read_ratio"
big_table$FP[big_table$filter_reason == "read_ratio"] <- TRUE
cat(
  paste0("\t\tFew unique reads: ", 
    nrow(big_table[big_table$filter_reason == "read_ratio",])), 
  file=summary_file, 
  sep="\n"
)

if (library == "single") {
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$polyA_percent < 0.2)] <- "polyApercent"
} else if (library == "micro") {
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$polyA_percent < 0.1)] <- "polyApercent"
} else {
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$polyA_percent == 0)] <- "polyApercent"
}

big_table$FP[big_table$filter_reason == "polyApercent"] <- TRUE
cat(
  paste0("\t\tNo polyA containing reads (junction spanning): ", 
    nrow(big_table[big_table$filter_reason == "polyApercent",])), 
  file=summary_file, 
  sep="\n"
)

# big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$breakpoint < 0.8)] <- "breakpoint"
# big_table$FP[big_table$filter_reason == "breakpoint"] <- TRUE
# cat(
#   paste0("\t\tInconsistent breakpoints: ", 
#     nrow(big_table[big_table$filter_reason == "breakpoint",])), 
#   file=summary_file, 
#   sep="\n"
# )

if (library == "single") {
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$RPM < 5)] <-"RPM"
} else {
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$RPM < 1)] <-"RPM"
}

big_table$FP[big_table$filter_reason == "RPM"] <- TRUE
cat(
  paste0("\t\tLow RPM: ", 
    nrow(big_table[(big_table$filter_reason == "RPM"),])), 
  file=summary_file, 
  sep="\n"
)

if (library == "single") {
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$peak_width < 150)] <- "peak_width"
  big_table$FP[big_table$filter_reason == "peak_width"] <- TRUE
  cat(
    paste0("\t\tPeak width: ",
      nrow(big_table[big_table$filter_reason == "peak_width",])), 
    file=summary_file, 
    sep="\n"
  )
}

big_table$classification[big_table$FP] <- "FP"
big_table$candidate[big_table$classification != "Candidate"] <- FALSE

# Label candidate peaks as UNK/SOM
if (library == "bulk") {
  big_table$classification[big_table$classification == "Candidate" & big_table$RPM >= 100 & big_table$TPM >= 4] <- "UNK"
  big_table$classification[big_table$classification == "Candidate" & big_table$number_templates > 1] <- "SOM_clonal"
  big_table$classification[big_table$classification == "Candidate" & big_table$number_templates == 1] <- "SOM_private"
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
  
  big_table$filter_reason[(big_table$peak %in% error_peak_list) & 
                          (grepl("SOM", big_table$classification)) & 
                          !(big_table$FP)] <- "error-prone"
  big_table$FP[big_table$filter_reason == "error-prone"] <- TRUE
  big_table$classification[big_table$FP] <- "FP"
  cat(
    paste0("\t\tOverlapping error-prone region: ", 
      nrow(big_table[big_table$filter_reason == "error-prone",])), 
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
  paste0("Total Off-target: ", 
    nrow(big_table[big_table$classification == "Off-target",])), 
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
filter_reasons <- ggplot(big_table[big_table$classification == "FP",], aes(x=filter_reason, group=classification, color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~classification) +
  ggtitle("Number of FP Peaks per Filter") + 
  ylab("Peak Count") + 
  xlab("Filter")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 

pdf(plots_path)
filter_reasons
dev.off()

# Candidate table
bed_cols = c(
  "chrm","start","end","peak","RPM","strand",
  "classification","filter_reason"
)

if (truth_set != "-NA-") {
  peak_table <- big_table[, c(bed_cols, "true_insertion_ID")]
} else {
  peak_table <- big_table[, bed_cols]
}

peak_table <- peak_table %>% arrange(chrm, start)
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
  filter_reasons_file, 
  sep="\t", 
  row.names=FALSE, 
  col.names=TRUE, 
  quote=FALSE, 
  na="-NA-"
)
