library(bedtoolsr)
library(UpSetR) 
library(stringr) 
library(reshape) 
library(ggplot2) 
library(tidyverse) 
library(RColorBrewer) 
library(ggpubr) 
library(ggsci) 

args<-commandArgs(TRUE)

big_table_file <- args[1]
library <- args[2]
plots_path <- args[3]
classified_peaks <- args[4]
summary <- args[5]
error_prone <- args[6]
truth_set <- args[7]
genome <- args[8]
big_table_filter_annotation_file <- paste0(sub("_classified_peaks.bed$", "", classified_peaks),
                                               "_big_table_filter_reasons.tsv")

# Read input and define columns
canonical_chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                    "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                    "chr20","chr21","chr22","chrX","chrY")

big_table <- read.table(big_table_file, sep="\t")
colnames(big_table) <- c("chrm","start","end","peak","shape","strand",
                         "reads","Ntag","polyA","gmotif","chimera","misaligned",
                         "repeatmasker","evrony","homopolymers",
                         "i1kgp","gnomad","nyuwa","xtea","hgsvc3","melt_lra","ont",
                         "bamreads","peakreads","usp","max_usp_diff","median_usp_diff","mean_usp_diff",
                         "RPM","nearest_peak","SegDups","max_distance")

big_table$gmotif_percent <- big_table$gmotif / big_table$reads
big_table <- merge(big_table, big_table[, c("RPM","peak")], by.x="nearest_peak", by.y="peak", suffixes = c("","_nearest"), all.x=TRUE)
big_table$nearest_RPM_ratio <- big_table$RPM / big_table$RPM_nearest
big_table$nearest_RPM_ratio[is.na(big_table$nearest_RPM_ratio)] <- 1
big_table$number_templates <- unlist(lapply(strsplit(big_table$usp, ";"), `[[`,1))
big_table$number_templates <- as.numeric(big_table$number_templates)
big_table$template_ratio <- NA
big_table$template_ratio[big_table$number_templates > 1] <- as.numeric(lapply(strsplit(unlist(lapply(strsplit(big_table$usp[big_table$number_templates >1], ";"), `[[`, 2)), ","), `[[`, 2)) / as.numeric(lapply(strsplit(unlist(lapply(strsplit(big_table$usp[big_table$number_templates >1] , ";"), `[[`, 2)), ","), `[[`, 1))
# big_table$most_duplicates <- as.numeric(lapply(strsplit(unlist(lapply(strsplit(big_table$usp, ";"), `[[`, 2)), ","), `[[`, 1)) 
# big_table$breadth_ratio <- as.numeric(big_table$reads) / as.numeric(big_table$number_templates)
# big_table$breadth <- as.numeric(big_table$end) - as.numeric(big_table$start)
# big_table$peak_shape <- strsplit(unlist(lapply(strsplit(big_table$shape, ";"), `[[`, 4)), ",")
# big_table$max_depth <- unlist(lapply(strsplit(big_table$shape, ";"), `[[`, 2))
big_table$polyA_percent <- big_table$polyA / big_table$reads
# big_table$segdups_annotation <- FALSE
# big_table$segdups_annotation[big_table$SegDups != '.'] <- TRUE
# big_table$chimera[big_table$chimera == ""] <- "chimera"
# big_table$misalignment[big_table$misaligned == ""] <- "misalignment"
# big_table$polymer_annotation <- FALSE
# big_table$polymer_annotation[big_table$homopolymers != '.'] <- TRUE
big_table$RPM_log <- log(big_table$RPM)


#### Classification and Filtering ####
summary_file <- file(summary, open='a')
cat(paste0("#################### Filtering results ", Sys.time(), "####################\n"), file=summary_file, sep="")
cat(paste0("Total peaks: ", nrow(big_table)), file=summary_file, sep="\n")

# Annotate peaks overlapping centromere or SegDup and peaks with Gmotif
cat("\nANNOTATIONS", file=summary_file, sep="\n")
cat(paste0("Total peaks overlapping ALR/Alpha satellites (centromeres): ", nrow(big_table[grepl("ALR/Alpha", big_table$repeatmasker),])), file=summary_file, sep="\n")
cat(paste0("Total peaks overlapping SegDups: ", nrow(big_table[big_table$SegDups != '.',])), file=summary_file, sep="\n")
cat(paste0("Total peaks with Gmotif % less than 50%: ", nrow(big_table[big_table$gmotif_percent < 0.5,])), file=summary_file, sep="\n")

# Initialize columns
big_table$classification <- "Candidate"
big_table$filter_reason <- "-NA-"
big_table$candidate <- TRUE
big_table$FP <- FALSE

# Label KR peaks
cat("\nKR/KNR FILTERS", file=summary_file, sep="\n")
big_table$classification[(grepl("L1HS", big_table$repeatmasker) | grepl("L1Hs", big_table$evrony))] <- "KR"
big_table$filter_reason[big_table$classification == "KR" & big_table$RPM < 50] <- "KR_low_RPM" 
big_table$FP[big_table$filter_reason == "KR_low_RPM"] <- TRUE
cat(paste0("KR RPM less than 50: ", nrow(big_table[big_table$filter_reason == "KR_low_RPM",])), file=summary_file, sep="\n")

# Intersect peaks with truth set if benchmarking
true_peak_IDs = list()
if (truth_set != "-NA-") {
  big_bed <- big_table[, c("chrm","start","end","peak","RPM","strand")]
  true_bed <- read_tsv(truth_set, col_names=c("chrm","start","end","name","score","strand"), col_types="ciicdc")

  true_bed_slop <- bt.slop(i=true_bed, g=genome, b=50)
  true_peak_list <- bt.intersect(a=big_bed, b=true_bed_slop, wo=TRUE, S=TRUE)[, c("V4","V10")]
  colnames(true_peak_list) <- c("peak","true_insertion_ID")
  true_peak_IDs = true_peak_list$peak

  big_table <- merge(big_table, true_peak_list, by="peak", all.x=TRUE) # adds true_peak_ID column
  big_table$true_insertion_ID[is.na(big_table$true_insertion_ID)] <- "-NA-"
}

# Label KNR peaks
# Also, artifically remove KNR labels from peaks in truth set if benchmarking
big_table$classification[(grepl("INS:ME:LINE1", big_table$i1kgp) | 
                         grepl("INS:ME:LINE1", big_table$gnomad) | 
                         grepl("INS:ME:LINE1", big_table$nyuwa) | 
                         grepl("INS:ME:LINE1", big_table$xtea) | 
                         grepl("LINE/L1", big_table$hgsvc3) | 
                         grepl("LINE1", big_table$melt_lra) | 
                         grepl("L1-INS", big_table$ont)) & 
                        !(big_table$peak %in% true_peak_IDs) &
                        !(big_table$classification == "KR")] <- "KNR"
big_table$filter_reason[big_table$classification == "KNR" & big_table$RPM < 50] <- "KNR_low_RPM"
big_table$FP[big_table$filter_reason == "KNR_low_RPM"] <- TRUE
cat(paste0("KNR RPM less than 50: ", nrow(big_table[big_table$filter_reason == "KNR_low_RPM",])), file=summary_file, sep="\n")

# Label Off-target peaks
big_table$classification[grepl("L1PA2|L1PA3|L1PA4|L1PA5|L1PA6|L1PA7", big_table$repeatmasker) & 
                        !(big_table$classification == "KR" | big_table$classification == "KNR")] <- "Off-target"
big_table$filter_reason[big_table$classification == "Off-target"] <- "Off-target"


# Filter candidate peaks for FPs
big_table$candidate[big_table$classification != "Candidate"] <- FALSE

# Filters are successive, meaning that a FP peak will be labeled by the first filter that catches it
# Filter peaks in non-canonical chromosomes
cat("\nCANDIDATE FILTERS", file=summary_file, sep="\n")
big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & !(big_table$chrm %in% canonical_chrs)] <- "chrms"
big_table$FP[big_table$filter_reason == "chrms"] <- TRUE
cat(paste0("Non-canonical chroms: ", nrow(big_table[big_table$filter_reason == "chrms",])), file=summary_file, sep="\n")

if (library != "bulk") {
  # In PTA libraries, remove smaller peaks nearby larger ones that are due to incorrect amplification during PTA 
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$nearest_RPM_ratio < 0.2)] <- "nearby_peak"
  big_table$FP[big_table$filter_reason == "nearby_peak"] <- TRUE
  cat(paste0("Nearby larger peak (less than 20% ratio): ", nrow(big_table[big_table$filter_reason == "nearby_peak",])), file=summary_file, sep="\n")

  # In single-cell libraries, remove peaks with too few original fragments
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$number_templates <= 2)] <- "num_templates"
  big_table$FP[big_table$filter_reason == "num_templates"] <- TRUE
  cat(paste0("Number of templates less than 3: ", nrow(big_table[big_table$filter_reason == "num_templates",])), file=summary_file, sep="\n")

  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$template_ratio <= 0.25)] <- "template_ratio"
  big_table$FP[big_table$filter_reason == "template_ratio"] <- TRUE
  cat(paste0("Template ratio: ", nrow(big_table[big_table$filter_reason == "template_ratio",])), file=summary_file, sep="\n")
}

big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$chimera == "chimera")] <- "chimera"
big_table$FP[big_table$filter_reason == "chimera"] <- TRUE
cat(paste0("Chimera: ", nrow(big_table[big_table$filter_reason == "chimera",])), file=summary_file, sep="\n")

big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$misaligned == "misalignment")] <- "misaligned_reads"
big_table$FP[big_table$filter_reason == "misaligned_reads"] <- TRUE
cat(paste0("Misaligned: ", nrow(big_table[big_table$filter_reason == "misaligned_reads",])), file=summary_file, sep="\n")

big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$polyA_percent == 0)] <- "polyApercent"
big_table$FP[big_table$filter_reason == "polyApercent"] <- TRUE
cat(paste0("No polyA containing reads (junction spanning): ", nrow(big_table[big_table$filter_reason == "polyApercent",])), file=summary_file, sep="\n")

if (library == "bulk") {
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$RPM <= 1)] <- "RPM"
  big_table$FP[big_table$filter_reason == "RPM"] <- TRUE
  cat(paste0("RPM less than or equal 1: ", nrow(big_table[(big_table$filter_reason == "RPM"),])), file=summary_file, sep="\n")
} else {
  big_table$filter_reason[(big_table$candidate) & !(big_table$FP) & (big_table$RPM <= 5)] <-"RPM"
  big_table$FP[big_table$filter_reason == "RPM"] <- TRUE
  cat(paste0("RPM less than or equal 5: ", nrow(big_table[big_table$filter_reason == "RPM",])), file=summary_file, sep="\n")
}

big_table$classification[big_table$FP] <- "FP"

# Label candidate peaks as UNK/SOM
if (library == "bulk") {
  big_table$classification[big_table$classification == "Candidate" & big_table$number_templates == 1] <- "SOM_private"
  big_table$classification[big_table$classification == "Candidate" & big_table$number_templates > 1] <- "SOM_clonal"
  big_table$classification[big_table$classification == "SOM_clonal" & big_table$RPM >= 100] <- "UNK" # multiple templates and RPM >= 100
} else if (library == "single") {
  big_table$classification[big_table$classification == "Candidate"] <- "UNK" # need multiple tissues to distinguish
} else if (library == "micro") { # TODO: Examine criteria
  big_table$classification[big_table$classification == "Candidate" & big_table$RPM >= 5] <- "UNK"
  big_table$classification[big_table$classification == "Candidate" & big_table$number_templates > 1] <- "SOM_clonal"
  big_table$classification[big_table$classification == "Candidate" & big_table$number_templates == 1] <- "SOM_private"
}

# Filter out SOM peaks in error-prone regions if regions are provided
if (error_prone != "-NA-") {
  big_bed <- big_table[, c("chrm","start","end","peak","RPM","strand")]
  error_bed <- read_tsv(error_prone, col_select=c(1,2,3), col_names=c("chrm","start","end"), col_types="cii")
    
  error_peak_list <- bt.intersect(a=big_bed, b=error_bed, wo=TRUE)[, "V4"] # TODO: debug

  # big_table$error-prone <- FALSE
  # big_table$error-prone[big_table$peak %in% error_peak_list] <- TRUE
  
  big_table$filter_reason[(big_table$peak %in% error_peak_list) & 
                          (big_table$classification == "SOM_private" | big_table$classification == "SOM_clonal") & 
                          !(big_table$FP)] <- "error-prone"
  big_table$FP[big_table$filter_reason == "error-prone"] <- TRUE
  big_table$classification[big_table$FP] <- "FP"
  cat(paste0("Overlapping error-prone region: ", nrow(big_table[big_table$filter_reason == "error-prone",])), file=summary_file, sep="\n")
}

big_table$FP <- NULL # drop column

cat(paste0("\nTotal KR: ", nrow(big_table[big_table$classification == "KR",])), file=summary_file, sep="\n")
cat(paste0("Total KNR: ", nrow(big_table[big_table$classification == "KNR",])), file=summary_file, sep="\n")

# big_table$classification[(big_table$classification == "KR") & (big_table$gmotif_percent >= 0.5)] <- "KR_gmotif_containing"
# big_table$classification[(big_table$classification == "KR") & (big_table$gmotif_percent < 0.5)] <- "KR_gmotif_missing"
# big_table$classification[(big_table$classification == "KNR") & (big_table$gmotif_percent >= 0.5)] <- "KNR_gmotif_containing"
# big_table$classification[(big_table$classification == "KNR") & (big_table$gmotif_percent < 0.5)] <- "KNR_gmotif_missing"

cat(paste0("Total Off-target: ", nrow(big_table[big_table$classification == "Off-target",])), file=summary_file, sep="\n")
cat(paste0("Total UNK: ", nrow(big_table[big_table$classification == "UNK",])), file=summary_file, sep="\n")
if (library == "bulk") {
  cat(paste0("Total SOM_clonal: ", nrow(big_table[big_table$classification == "SOM_clonal",])), file=summary_file, sep="\n")
  cat(paste0("Total SOM_private: ", nrow(big_table[big_table$classification == "SOM_private",])), file=summary_file, sep="\n")
}
cat(paste0("Total FP: ", nrow(big_table[big_table$classification == "FP",])), file=summary_file, sep="\n")


# TODO: Update plotting to work with new classification order
#### Plotting (pre-filter) ####
# filter_chrm <- ggplot(big_table, aes(x=(chrm %in% canonical_chrs), group=classification, color=classification, fill=classification)) +
#   geom_bar() +
#   scale_fill_nejm() +
#   scale_color_nejm()+
#   facet_wrap(~classification) +
#   ggtitle("Distribution of Peaks in Canonical Chromosomes") + 
#   ylab("Peak Count") + 
#   xlab("Mapped to Canonical Chromosome")

# filter_PTA_artifact  <- ggplot(big_table, aes(x=nearest_RPM_ratio, group=classification, color=classification, fill=classification)) +
#   geom_histogram() +
#   scale_fill_nejm() +
#   scale_color_nejm() +
#   geom_vline(xintercept = 0.2, color = "red", linetype = "dashed", size = 1) +
#   annotate("text", x = .55, y = .5*nrow(big_table), label = "Threshold \n 0.2%", color = "red", angle = 0, vjust = -0.5) +
#   facet_wrap(~classification) +
#   ggtitle("Distribution of Ratio of Artifact Peak to Nearest Peak") + 
#   ylab("Peak Count") + 
#   xlab("Peak RPM / Largest Nearby Peak RPM")

filter_templates <- ggplot(big_table, aes(x=classification, y=number_templates, color=classification, fill=classification)) +
  geom_violin() +
  geom_jitter(height=0, width=0.1) + 
  scale_fill_nejm() +
  scale_color_nejm() +
#  geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
#  annotate("text", x = 70, y = .5*nrow(big_table), label = "Threshold \n 3", color = "red", angle = 0, vjust = -0.5) +
#  facet_wrap(~classification) +
  ggtitle("Distribution of Number of Templates per Peak") + 
  ylab("Peak Count") + 
  xlab("Number of Templates")

filter_RPM <- ggplot(big_table, aes(x=classification, y=RPM_log, color=classification, fill=classification)) +
  geom_violin() +
  geom_jitter(height=0, width=0.1) + 
  scale_fill_nejm() +
  scale_color_nejm() +
#  facet_wrap(~classification, scales="free") +
  ggtitle("Distribution of RPM") + 
  ylab("Peak Count") + 
  xlab("ln(RPM)")

# filter_gmotif <-ggplot(big_table, aes(x = classification, y=gmotif_percent, group=classification, color=classification, fill=classification)) +
#   geom_violin() +
#   scale_fill_nejm() +
#   scale_color_nejm() +
#   geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
#   annotate("text", x=0.65,  y = 0.75, label = "Threshold \n 80%", color = "red", angle = 0, vjust = -0.5) +
#   #facet_wrap(~classification) +
#   ggtitle("Distribution of Gmotif Containing Read Support per Peak") + 
#   ylab("Percent of Reads Containing Gmotif") + 
#   xlab("Classification")

# filter_template_ratio <-ggplot(big_table[!is.na(big_table$template_ratio),], aes(x = classification, y=template_ratio, group=classification, color=classification, fill=classification)) +
#   geom_violin() +
#   scale_fill_nejm() +
#   scale_color_nejm() +
#   geom_hline(yintercept = 0.25, color = "red", linetype = "dashed", size = 1) +
#   annotate("text", x=2.5,  y = 0.25, label = "Threshold \n 25%", color = "red", angle = 0, vjust = -0.5) +
#   #facet_wrap(~classification) +
#   ggtitle("Distribution of Ratio of Two Highest Duplicate Depth Templates per Peak") + 
#   ylab("Ratio of Two Highest Duplicate Depth Templates per Peak") + 
#   xlab("Classification")

# filter_chimera <- ggplot(big_table, aes(x=chimera, group=classification , color=classification, fill=classification)) +
#   geom_bar() +
#   scale_fill_nejm() +
#   scale_color_nejm()+
#   facet_wrap(~classification) +
#   ggtitle("Distribution of Chimeric Peaks") + 
#   ylab("Peak Count") + 
#   xlab("Chimera")

# filter_misalignment <- ggplot(big_table, aes(x=misaligned, group=classification, color=classification, fill=classification)) + 
#     geom_bar()+
#     scale_fill_nejm()+
#     scale_color_nejm()+
#     facet_wrap(~classification) + 
#     ggtitle("Distribution of Misaligned reads Peaks") + 
#     ylab("Peak Count")+
#     xlab("Misalignment")

#  filter_polyA_spanning <-ggplot(big_table, aes(x = classification, y=polyA_percent, group=classification, color=classification, fill=classification)) +
#    geom_violin() +
#    scale_fill_nejm() +
#    scale_color_nejm() +
#    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
#    annotate("text", x=2.5,  y = 0.1, label = "Threshold \n 0%", color = "red", angle = 0, vjust = -0.5) +
#    #facet_wrap(~classification) +
#    ggtitle("Distribution of Percent of Reads Spanning the Insertion Junction per Peak") + 
#    ylab("Percent of Reads Spanning Insertion Junction") + 
#    xlab("Classification")

# filter_satellite <-ggplot(big_table, aes(x = (grepl("ALR/Alpha", repeatmasker)), group=classification, color=classification, fill=classification)) +
#   geom_bar() +
#   scale_fill_nejm() +
#   scale_color_nejm() +
#   facet_wrap(~classification) +
#   ggtitle("Distribution of Peaks in Satellite Regions") + 
#   ylab("Peak Count") + 
#   xlab("In Satellite Region")

# filter_segdup <-ggplot(big_table, aes(x = (SegDups != '.'), group=classification, color=classification, fill=classification)) +
#   geom_bar() +
#   scale_fill_nejm() +
#   scale_color_nejm() +
#   facet_wrap(~classification) +
#   ggtitle("Distribution of Peaks in Segmental Duplications") + 
#   ylab("Peak Count") + 
#   xlab("In Segmental Duplication")


#### Plotting (post-filter) ####
# filter_reasons <- ggplot(big_table[big_table$classification == "FP",], aes(x=filter_reason, group=classification, color=classification, fill=classification)) +
#   geom_bar() +
#   scale_fill_nejm() +
#   scale_color_nejm() +
#   facet_wrap(~classification) +
#   ggtitle("Distribution of Filtering FP Peaks") + 
#   ylab("Peak Count") + 
#   xlab("Filter")+
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 

# filter_reasons_KR <- ggplot(big_table[big_table$KR,], aes(x=filter_reason, group=classification, color=classification, fill=classification)) +
#   geom_bar() +
#   scale_fill_nejm() +
#   scale_color_nejm() +
#   facet_wrap(~classification)  + 
#   ggtitle("Distribution of Filtering KR Peaks") + 
#   ylab("Peak Count") + 
#   xlab("Filter")+
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# filter_reasons_KNR <- ggplot(big_table[big_table$KNR,], aes(x=filter_reason, group=classification, color=classification, fill=classification)) +
#   geom_bar() +
#   scale_fill_nejm() +
#   scale_color_nejm() +
#   facet_wrap(~classification) +
#   ggtitle("Distribution of Filtering KNR Peaks") + 
#   ylab("Peak Count") + 
#   xlab("Filter") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 

pdf(plots_path)
# filter_chrm
# filter_PTA_artifact
filter_templates
filter_RPM
# filter_gmotif
# filter_template_ratio
# filter_chimera
# filter_misalignment
# filter_polyA_spanning
# filter_satellite
# filter_segdup
# filter_reasons
# filter_reasons_KR
# filter_reasons_KNR
dev.off()

# Candidate table
bed_cols = c("chrm","start","end","peak","RPM","strand","classification","filter_reason")

if (truth_set != "-NA-") {
  peak_table <- big_table[, c(bed_cols, "true_insertion_ID")]
} else {
  peak_table <- big_table[, bed_cols]
}

write.table(peak_table, classified_peaks, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, na="-NA-")
write.table(big_table[, lapply(big_table, class) != "list"], big_table_filter_annotation_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, na="-NA-")
