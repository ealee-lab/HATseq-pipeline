library(bedtoolsr)
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

#### INPUT ####
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

# Define filtering conditions
wide_peaks <- big_table$peak_width >= 150

# In PTA libraries, flag smaller peaks nearby larger ones (amplification artifacts)
if (library == "single") {
  solo_peaks <- big_table$nearest_RPM_ratio == 1
} else if (library == "micro") {
  solo_peaks <- big_table$nearest_RPM_ratio == 1
} else {
  solo_peaks <- big_table$nearest_RPM_ratio > 0 # all peaks
}

# Require that peaks have multiple templates
if (library == "single") {
  template_peaks <- big_table$TPM >= 0.25
} else if (library == "micro") {
  template_peaks <- big_table$num_templates >= 2
} else {
  template_peaks <- big_table$num_templates >= 1 # all peaks
}

# In PTA libraries, flag peaks with uneven amplification of templates
if (library == "single") {
  even_peaks <- big_table$template_ratio >= 0.2
} else if (library == "micro") {
  even_peaks <- big_table$template_ratio >= 0.2
} else {
  even_peaks <- big_table$template_ratio > 0 # all peaks
}

# Require that at least 10% of peak reads are unique
unique_peaks <- big_table$unique_read_ratio >= 0.1

# Require that peaks have polyA reads
if (library == "single") {
  polyA_peaks <- big_table$polyA_percent >= 20
} else if (library == "micro") {
  polyA_peaks <- big_table$polyA_percent >= 10
} else {
  polyA_peaks <- big_table$polyA_percent > 0
}

# Require that peaks have a minimum signal level
highRPM_peaks <- big_table$RPM >= 100
midRPM_peaks <- big_table$RPM >= 50

if (library == "single") {
  lowRPM_peaks <- big_table$RPM >= 5
} else {
  lowRPM_peaks <- big_table$RPM >= 1
}

# Filter peaks by set conditions
# Filters are successive (FP peak will be labeled by the first filter it fails)
cat("\nFILTERS", file=summary_file, sep="\n")

# Global filters (all peaks)
if (library != "bulk") {
  cat("\n\tAll peaks", file=summary_file, sep="\n")

  big_table$filter[!(big_table$FP) & !(solo_peaks)] <- "nearbyPeak"
  big_table$FP[big_table$filter == "nearbyPeak"] <- TRUE
  cat(
    paste0("\t\tNearby larger peak: ", 
      nrow(big_table[big_table$filter == "nearbyPeak",])), 
    file=summary_file, 
    sep="\n"
  )
}

# Label known reference (KR) peaks - priority
big_table$classification[(grepl("L1HS", big_table$nearest_KR) | 
                          grepl("L1Hs", big_table$nearest_KR)) &
                          !(big_table$FP)] <- "KR"
big_table$filter[big_table$classification == "KR"] <- "lowConf"
big_table$filter[(big_table$classification == "KR") & (wide_peaks) & (highRPM_peaks)] <- "PASS"

# Label L1PA (Non-specific) peaks
big_table$classification[grepl("L1PA2|L1PA3|L1PA4|L1PA5", big_table$nearest_KR) & 
                         (big_table$nearest_KR_dist <= big_table$nearest_KNR_dist) &
                         !(big_table$classification == "KR") & 
                         !(big_table$FP)] <- "Non-specific"
big_table$filter[big_table$classification == "Non-specific"] <- "lowConf"
big_table$filter[(big_table$classification == "Non-specific") & (wide_peaks) & (highRPM_peaks)] <- "PASS"

# Label known non-reference (KNR) peaks
# Also, artifically remove KNR labels from peaks in truth set for benchmarking
big_table$classification[(grepl("LINE1", big_table$nearest_KNR) | 
                         grepl("L1", big_table$nearest_KNR)) & 
                         !(big_table$peak %in% true_peak_IDs) &
                         !(big_table$classification == "KR") &
                         !(big_table$classification == "Non-specific") &
                         !(big_table$FP)] <- "KNR"
big_table$filter[big_table$classification == "KNR"] <- "lowConf"
big_table$filter[(big_table$classification == "KNR") &
                 (wide_peaks) &
                 (unique_peaks) &
                 (polyA_peaks) &
                 (midRPM_peaks)] <- "PASS"

# Label noise near known elements
big_table$filter[(big_table$classification != "Candidate") & !(lowRPM_peaks)] <- "noise"

# Filter remaining peaks for false positives
big_table$candidate[big_table$classification != "Candidate"] <- FALSE
cat("\n\tCandidate peaks", file=summary_file, sep="\n")
cat(paste0("\t\tTotal peaks: ", nrow(big_table[big_table$candidate,])), file=summary_file, sep="\n")

big_table$filter[(big_table$candidate) & !(big_table$FP) & !(template_peaks)] <- "templates"
big_table$FP[big_table$filter == "templates"] <- TRUE
cat(
  paste0("\t\tToo few templates: ", 
    nrow(big_table[big_table$filter == "templates",])), 
  file=summary_file, 
  sep="\n"
)

big_table$filter[(big_table$candidate) & !(big_table$FP) & !(even_peaks)] <- "templateRatio"
big_table$FP[big_table$filter == "templateRatio"] <- TRUE
cat(
  paste0("\t\tTemplate ratio: ",
    nrow(big_table[big_table$filter == "templateRatio",])), 
  file=summary_file, 
  sep="\n"
)

big_table$filter[(big_table$candidate) & !(big_table$FP) & !(unique_peaks)] <- "uniqueReads"
big_table$FP[big_table$filter == "uniqueReads"] <- TRUE
cat(
  paste0("\t\tFew unique reads: ", 
    nrow(big_table[big_table$filter == "uniqueReads",])), 
  file=summary_file, 
  sep="\n"
)

big_table$filter[(big_table$candidate) & !(big_table$FP) & !(polyA_peaks)] <- "polyApercent"
big_table$FP[big_table$filter == "polyApercent"] <- TRUE
cat(
  paste0("\t\tNo polyA containing reads (junction spanning): ", 
    nrow(big_table[big_table$filter == "polyApercent",])), 
  file=summary_file, 
  sep="\n"
)

big_table$filter[(big_table$candidate) & !(big_table$FP) & !(lowRPM_peaks)] <-"RPM"
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

# Label peaks passing filters as UNK/SOM
if (library == "bulk") {
  big_table$classification[big_table$candidate & (midRPM_peaks) & (wide_peaks)] <- "UNK"
  # big_table$classification[big_table$classification == "Candidate" & big_table$RPM >= 100 & big_table$TPM >= 4] <- "UNK"
  big_table$classification[big_table$candidate & ((big_table$num_templates > 1) & (big_table$max_usp_distance >= 10))] <- "SOM_clonal"
  big_table$classification[big_table$candidate & ((big_table$num_templates == 1) | (big_table$max_usp_distance < 10))] <- "SOM_private"
} else if (library == "micro") {
  big_table$classification[big_table$candidate & big_table$RPM >= 100 & big_table$TPM >= 4] <- "UNK"
  big_table$classification[big_table$candidate] <- "SOM"
} else {
  big_table$classification[big_table$candidate] <- "SOM" # need multiple cells to distinguish
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

#### PLOTTING ####
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

#### OUTPUT ####
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
