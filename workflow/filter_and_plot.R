library(bedtoolsr)
library(UpSetR) 
library(stringr) 
library(reshape) 
library(ggplot2) 
library(tidyverse) 
library(RColorBrewer) 
library(ggpubr) 
library(ggsci) 
#library(circlize)

args<-commandArgs(TRUE)

big_table_file <- args[1]
library <- args[2]
plots_path <- args[3]
filtered_peaks <- args[4]
summary <- args[5]
error_prone <- args[6]
truth_set <- args[7]
genome <- args[8]

# File not tracked by Snakemake
big_table_filter_annotation_file <- paste0(sub("_filtered_peaks.tsv$", "", filtered_peaks), "_big_table_filter_reasons.tsv")

# Read input and define columns
canonical_chrs <- c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
big_table <- read.table(big_table_file, sep="\t")
colnames(big_table) <- c("chrm","start","end","peak","shape","strand","reads","Ntag","polyA","gmotif","chimera","misaligned","repeatmasker","evrony","homopolymers","gnomad","i1gp","nyuwa","xTea","bamreads","peakreads","usp" ,"max_usp_diff", "median_usp_diff","mean_usp_diff","RPM","nearest_peak","SegDups","max_distance")
#colnames(big_table) <- c("chrm","start","end","peak","shape","strand","reads","Ntag","polyA","gmotif","s3p","s5p","h5p","h3p","chimera","repeatmasker","evrony","homopolymers","gnomad","i1gp","nyuwa","bamreads","peakreads","usp","RPM","custompeakname","nearest_peak","SegDups","max_distance","chrgts","startgts","endgts","chrgts2","startgts2","sam1","sam2","sam3","sam4","sam5","sam6","something","something2","something3","level","number","overlap")

big_table$gmotif_percent <- big_table$gmotif / big_table$reads
# big_table$map_to_source = lapply(strsplit(paste(big_table$s3p, big_table$h3p, sep=","),","), unique)
# big_table$map_to_target = lapply(strsplit(paste(big_table$s5p, big_table$h5p, sep=","),","), unique)
big_table <- merge(big_table,big_table[,c("RPM","peak")],by.x="nearest_peak",by.y="peak",suffixes = c("","_nearest"),all.x=TRUE)
big_table$nearest_RPM_ratio <- big_table$RPM / big_table$RPM_nearest
big_table$nearest_RPM_ratio[is.na(big_table$nearest_RPM_ratio)] <- 1
big_table$number_templates <- unlist(lapply(strsplit(big_table$usp, ";"), `[[`,1))
big_table$number_templates <- as.numeric(big_table$number_templates)
big_table$template_ratio <- NA
big_table$template_ratio[big_table$number_templates > 1] <- as.numeric(lapply(strsplit(unlist(lapply(strsplit(big_table$usp[big_table$number_templates >1] , ";"), `[[`,2)),","),`[[`,2)) / as.numeric(lapply(strsplit(unlist(lapply(strsplit(big_table$usp[big_table$number_templates >1] , ";"), `[[`,2)),","),`[[`,1))
big_table$most_duplicates <- as.numeric(lapply(strsplit(unlist(lapply(strsplit(big_table$usp , ";"), `[[`,2)),","),`[[`,1)) 
big_table$breadth_ratio <- as.numeric(big_table$reads) / as.numeric(big_table$number_templates)
big_table$breadth <- as.numeric(big_table$end) - as.numeric(big_table$start)
big_table$peak_shape <- strsplit(unlist(lapply(strsplit(big_table$shape,";"),`[[`,4)),",")
big_table$max_depth <- unlist(lapply(strsplit(big_table$shape,";"),`[[`,2))
big_table$polyA_percent <- big_table$polyA / big_table$reads
big_table$segdups_annotation <- FALSE
big_table$segdups_annotation[big_table$SegDups != '.'] <- TRUE
big_table$chimera[big_table$chimera == ""] <- "chimera"
big_table$misalignment[big_table$misaligned == ""] <- "misalignment"
big_table$polymer_annotation <- FALSE
big_table$polymer_annotation[big_table$homopolymers != '.'] <- TRUE
big_table$RPM_log <- log(big_table$RPM)

if (error_prone != "None") { 
  big_table$error_prone <- FALSE # add optional column
}

####

# 
# library(pdc)
# #library(rowr)
# 
# require(plyr) # requires plyr for rbind.fill()
# library(rlist)
# #test_shapes <- as.data.frame(t(as.numeric(unlist(big_table$peak_shape[464]))))
# #test_shapes <- rbind.fill(test_shapes, as.data.frame(t(as.numeric(unlist(big_table$peak_shape[465])))), fill=NA )
# #as.matrix(as.numeric(unlist(big_table$peak_shape[464:500])))
# clust_shapes <- c()
# for (i in 1001:1100){
#   
#   if(is.null(nrow(clust_shapes)) ){
#     if(big_table$strand[i] == "-"){
#       clust_shapes <- as.data.frame(t(as.numeric(list.reverse(unlist(big_table$peak_shape[i])))))
#     }else{
#       clust_shapes <- as.data.frame(t(as.numeric(unlist(big_table$peak_shape[i]))))
#     }
#     
#   }else{
#     if(big_table$strand[i] == "-"){
#       clust_shapes <- rbind.fill(clust_shapes, as.data.frame(t(as.numeric(list.reverse(unlist(big_table$peak_shape[i]))))))
#     }else{
#       clust_shapes <- rbind.fill(clust_shapes, as.data.frame(t(as.numeric(unlist(big_table$peak_shape[i])))))
#     }
#     
#   }
# }
# clust_shapes_t <- t(clust_shapes)
# test_clust <- pdclust(clust_shapes_t, m=2, t=1)
# plot(test_clust)
# #library(mixtools)
# # 
# # #Fit two Normals
# #test_table <- big_table[,c(40,47)]
# #test = normalmixEM(test_table$RPM_log,  k=2,verb=TRUE)
# 
# # plot(wait1, density=TRUE, loglik=FALSE)
# # big_table$RPM_distribution[big_table$KR] <- wait1$posterior[,which.max(wait1$mu)] > wait1$posterior[,which.min(wait1$mu)]
# # 
# # 
# # wait2 = normalmixEM(big_table$RPM[big_table$KNR],  k=2,verb=TRUE)
# # plot(wait2, density=TRUE, loglik=FALSE)
# # big_table$RPM_distribution[big_table$KNR] <- wait2$posterior[,which.max(wait2$mu)] > wait2$posterior[,which.min(wait2$mu)]
# # 
# # 
# # wait3 = normalmixEM(big_table$RPM[!(big_table$KNR | big_table$Off_target_amplification | big_table$KR)],  k=3,verb=TRUE)
# # plot(wait3, density=TRUE, loglik=FALSE)
# # big_table$RPM_distribution[big_table$KR] <- wait3$posterior[,which.max(wait3$mu)] > wait3$posterior[,which.min(wait3$mu)]
# # 
# # 
# # wait4 = normalmixEM(big_table$RPM[big_table$Off_target_amplification],  k=2,verb=TRUE)
# # plot(wait4, density=TRUE, loglik=FALSE)
# # big_table$RPM_distribution[big_table$KR] <- wait4$posterior[,which.max(wait4$mu)] > wait4$posterior[,which.min(wait4$mu)]
# # 


big_table$FP <- FALSE

big_table$KR <- FALSE
big_table$KR[(grepl("L1HS", big_table$repeatmasker) | grepl("L1Hs", big_table$evrony)) & !big_table$FP  ] <- TRUE

big_table$Off_target_amplification <- FALSE
big_table$Off_target_amplification[(grepl("L1PA2|L1PA3|L1PA4|L1PA5", big_table$repeatmasker) ) & !(big_table$KR) ] <- TRUE

big_table$KNR <- FALSE
big_table$KNR[(grepl("LINE1", big_table$gnomad) | grepl("LINE1", big_table$i1gp) | grepl("LINE1", big_table$nyuwa) | grepl("LINE1",big_table$xTea)) & !(big_table$KR | big_table$FP)  ] <- TRUE

# If benchmarking, artifically remove KNR labels from peaks in truth set
if (truth_set != "None") {
  big_bed <- big_table[, c("chrm", "start", "end", "peak", "RPM", "strand")]
  true_bed <- read_tsv(truth_set, col_names=c("chrm", "start", "end", "name", "score", "strand"), col_types="ciicdc")

  true_bed_slop <- bt.slop(i=true_bed, g=genome, b=50)
  true_peak_list <- bt.intersect(a=big_bed, b=true_bed_slop, wo=TRUE, S=TRUE)[, c("V4", "V10")]
  colnames(true_peak_list) <- c("peak", "true_insertion_ID")

  big_table$KNR[big_table$peak %in% true_peak_list$peak] <- FALSE 
}

big_table$UNK <- FALSE
big_table$SOM_clonal <- FALSE
big_table$SOM_private <- FALSE

big_table$classification  <- "Candidate"
big_table$classification[big_table$FP] <- "FP"
big_table$classification[big_table$KNR] <-"KNR"
big_table$classification[big_table$KR] <- "KR"
#big_table$classification[big_table$UNK] <-"UNK"
#big_table$classification[big_table$SOM_clonal] <- "Clonal_Somatic"
#big_table$classification[big_table$SOM_private] <- "Private_Somatic"
big_table$classification[big_table$Off_target_amplification] <- "Off_target"

#### Plotting ####

filter_chrm <- ggplot(big_table, aes(x=(chrm %in% canonical_chrs),group=classification , color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm()+
  facet_wrap(~classification) +
  ggtitle("Distribution of Peaks in Canonical Chromosomes") + 
  ylab("Peak Count") + 
  xlab("Mapped to Canonical Chromosome")
filter_PTA_artifact  <- ggplot(big_table, aes(x=nearest_RPM_ratio, group=classification, color=classification, fill=classification)) +
  geom_histogram() +
  scale_fill_nejm() +
  scale_color_nejm() +
  geom_vline(xintercept = 0.2, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = .55, y = .5*nrow(big_table), label = "Threshold \n 0.2%", color = "red", angle = 0, vjust = -0.5) +
  facet_wrap(~classification) +
  ggtitle("Distribution of Ratio of Artifact Peak to Nearest Peak") + 
  ylab("Peak Count") + 
  xlab("Peak RPM / Largest Nearby Peak RPM")
filter_templates <-ggplot(big_table, aes(x=number_templates, group=classification, color=classification, fill=classification)) +
  geom_histogram(binwidth=3) +
  scale_fill_nejm() +
  scale_color_nejm() +
  geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 70, y = .5*nrow(big_table), label = "Threshold \n 3", color = "red", angle = 0, vjust = -0.5) +
  facet_wrap(~classification) +
  ggtitle("Distribution of Number of Templates per Peak") + 
  ylab("Peak Count") + 
  xlab("Number of Templates")
filter_RPM<-ggplot(big_table, aes(x=RPM_log, group=classification, color=classification, fill=classification)) +
  geom_histogram() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~classification) +
  ggtitle("Distribution of RPM") + 
  ylab("Peak Count") + 
  xlab("ln(RPM)")
filter_gmotif <-ggplot(big_table, aes(x = classification, y=gmotif_percent, group=classification, color=classification, fill=classification)) +
  geom_violin() +
  scale_fill_nejm() +
  scale_color_nejm() +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x=0.65,  y = 0.75, label = "Threshold \n 80%", color = "red", angle = 0, vjust = -0.5) +
  #facet_wrap(~classification) +
  ggtitle("Distribution of Gmotif Containing Read Support per Peak") + 
  ylab("Percent of Reads Containing Gmotif") + 
  xlab("Classification")
filter_template_ratio <-ggplot(big_table[!is.na(big_table$template_ratio),], aes(x = classification, y=template_ratio, group=classification, color=classification, fill=classification)) +
  geom_violin() +
  scale_fill_nejm() +
  scale_color_nejm() +
  geom_hline(yintercept = 0.25, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x=2.5,  y = 0.25, label = "Threshold \n 25%", color = "red", angle = 0, vjust = -0.5) +
  #facet_wrap(~classification) +
  ggtitle("Distribution of Ratio of Two Highest Duplicate Depth Templates per Peak") + 
  ylab("Ratio of Two Highest Duplicate Depth Templates per Peak") + 
  xlab("Classification")
filter_chimera <- ggplot(big_table, aes(x=chimera,group=classification , color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm()+
  facet_wrap(~classification) +
  ggtitle("Distribution of Chimeric Peaks") + 
  ylab("Peak Count") + 
  xlab("Chimera")
filter_misalignment <- ggplot(big_table, aes(x=misaligned,group=classification,color=classification, fill=classification)) + 
    geom_bar()+
    scale_fill_nejm()+
    scale_color_nejm()+
    facet_wrap(~classification) + 
    ggtitle("Distribution of Misaligned reads Peaks") + 
    ylab("Peak Count")+
    xlab("Misalignment")
 filter_polyA_spanning <-ggplot(big_table, aes(x = classification, y=polyA_percent, group=classification, color=classification, fill=classification)) +
   geom_violin() +
   scale_fill_nejm() +
   scale_color_nejm() +
   geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
   annotate("text", x=2.5,  y = 0.1, label = "Threshold \n 0%", color = "red", angle = 0, vjust = -0.5) +
   #facet_wrap(~classification) +
   ggtitle("Distribution of Percent of Reads Spanning the Insertion Junction per Peak") + 
   ylab("Percent of Reads Spanning Insertion Junction") + 
   xlab("Classification")
filter_satellite <-ggplot(big_table, aes(x = (grepl("ALR/Alpha", repeatmasker)), group=classification, color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~classification) +
  ggtitle("Distribution of Peaks in Satellite Regions") + 
  ylab("Peak Count") + 
  xlab("In Satellite Region")
filter_segdup <-ggplot(big_table, aes(x = (SegDups != '.'), group=classification, color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~classification) +
  ggtitle("Distribution of Peaks in Segmental Duplications") + 
  ylab("Peak Count") + 
  xlab("In Segmental Duplication")
####

#### Printing number of each filter ####
summary_file <- file(summary, open='a')
cat(paste0("\n#################### Filtering results ",Sys.time(),"####################\n"),file=summary_file,sep="")
big_table$filter_reason <- "-NA-"  
cat(paste0("Total peaks: " , nrow(big_table[!big_table$FP,])),file=summary_file,sep="\n")

# Annotate candidate peaks with G-motif and SegDup
cat("\nANNOTATIONS",file=summary_file,sep="\n")
cat(paste0("Peaks with Gmotif percent less than 50% and not KR or KNR overlapping: ", 
 nrow(big_table[big_table$gmotif_percent < 0.5 & !(big_table$KNR | big_table$KR),])),file=summary_file,sep="\n")
cat(paste0("Peaks overlapping SegDups: ", nrow(big_table[big_table$SegDups != '.',])),file=summary_file,sep="\n")

#Filter peaks in non-canonical chromosomes
cat("\nFILTERS",file=summary_file,sep="\n")
big_table$filter_reason[!(big_table$chrm %in% canonical_chrs) & !(big_table$FP)] <- "chrms"
big_table$FP[!(big_table$chrm %in% canonical_chrs) & !(big_table$FP)] <- TRUE
cat(paste0("Non-canonical chroms: " , nrow(big_table[big_table$filter_reason == "chrms",] )),file=summary_file,sep="\n")

#In PTA libraries, remove smaller peaks nearby larger ones that are due to incorrect amplification during PTA 
if (library != "bulk"){
  big_table$filter_reason[big_table$nearest_RPM_ratio < 0.2 & !(big_table$FP)] <- "nearby_peak"
  big_table$FP[big_table$nearest_RPM_ratio < 0.2 & !(big_table$FP)]<- TRUE
  cat(paste0("Nearby larger peak (less than 20% ratio): " , nrow(big_table[big_table$filter_reason == "nearby_peak",] )),file=summary_file,sep="\n")
}

#In Single Cell libraries remove peaks with two few original fragments
if (library != "bulk") {
  big_table$filter_reason[big_table$number_templates <= 2 & !(big_table$FP)] <- "num_templates"
  big_table$FP[big_table$number_templates <= 2 & !(big_table$FP)]<- TRUE
  cat(paste0("Number of templates less than 3: " , nrow(big_table[big_table$filter_reason == "num_templates",] )),file=summary_file,sep="\n")
}


# if(library != "single"){
#   big_table$filter_reason[(big_table$RPM <= 5 & (big_table$KR | big_table$KNR) & !(big_table$FP)) ] <-"RPM" 
#   big_table$FP[(big_table$RPM <= 5 & (big_table$KR | big_table$KNR) & !(big_table$FP))  ]<- TRUE  #| (big_table$RPM <= 1 & !(big_table$KR | big_table$KNR)& !(big_table$FP))
#   cat(paste0("RPM less than or equal 5 for KR and KNR overlap: " , nrow(big_table[(big_table$filter_reason == "RPM" ) ,] )),file=summary_file,sep="\n")
# }else {
#   big_table$filter_reason[(big_table$RPM <= 5 & (big_table$KR ) & !(big_table$FP)) | (big_table$RPM <= 2 & (big_table$KNR) & !(big_table$FP)) | (big_table$RPM <= 10 & !(big_table$KR | big_table$KNR) & !(big_table$FP)) ] <-"RPM" #<- "RPM"
#   big_table$FP[(big_table$RPM <= 5 & (big_table$KR ) & !(big_table$FP)) | (big_table$RPM <= 2 & (big_table$KNR) & !(big_table$FP)) | (big_table$RPM <= 10 & !(big_table$KR | big_table$KNR) & !(big_table$FP)) ]<- TRUE
#   cat(paste0("RPM less than or equal 5 for KR overlap and less than 2 for KNR overlap and less than or equal 10 all other: " , nrow(big_table[big_table$filter_reason == "RPM" ,] )),file=summary_file,sep="\n")
# }


if(library != "bulk"){
  big_table$filter_reason[big_table$template_ratio <= 0.25 & !(big_table$FP)] <- "template_ratio"
  big_table$FP[big_table$template_ratio <= 0.25 & !(big_table$FP)]<- TRUE
  cat(paste0("Template ratio: " , nrow(big_table[big_table$filter_reason == "template_ratio",] )),file=summary_file,sep="\n")
}

big_table$filter_reason[((big_table$chimera == "chimera") & !(big_table$KR) & !(big_table$FP))] <- "chimera"
big_table$FP[((big_table$chimera == "chimera") & !(big_table$KR) & !(big_table$FP))]<- TRUE
cat(paste0("Chimera: " , nrow(big_table[big_table$filter_reason == "chimera",] )),file=summary_file,sep="\n")

big_table$filter_reason[((big_table$misaligned == "misalignment") & !(big_table$FP) )] <- "misaligned_reads"
big_table$FP[((big_table$misaligned == "misalignment") & !(big_table$FP) )]<- TRUE
cat(paste0("Misaligned: " , nrow(big_table[big_table$filter_reason == "misaligned_reads",] )),file=summary_file,sep="\n")

big_table$filter_reason[(big_table$polyA_percent == 0) & !(big_table$FP)] <- "polyApercent"
big_table$FP[(big_table$polyA_percent == 0) & !(big_table$FP)]<- TRUE
cat(paste0("No polyA containing reads (junction spanning): " , nrow(big_table[big_table$filter_reason == "polyApercent",] )),file=summary_file,sep="\n")

big_table$filter_reason[grepl("ALR/Alpha",big_table$repeatmasker) & !(big_table$FP)] <- "ALR/Alpha"
big_table$FP[grepl("ALR/Alpha",big_table$repeatmasker) & !(big_table$FP)] <- TRUE
cat(paste0("ALR/Alpha Satellite overlapping: " , nrow(big_table[big_table$filter_reason == "ALR/Alpha",] )),file=summary_file,sep="\n")

if(library == "bulk"){
  big_table$filter_reason[(big_table$RPM <= 1) & !(big_table$FP)] <- "RPM" #& (big_table$KR | big_table$KNR) & !(big_table$FP))
  big_table$FP[(big_table$RPM <= 1) & !(big_table$FP)] <- TRUE  #| (big_table$RPM <= 1 & !(big_table$KR | big_table$KNR)& !(big_table$FP))
  cat(paste0("RPM less than or equal 1: " , nrow(big_table[(big_table$filter_reason == "RPM" ) ,] )),file=summary_file,sep="\n")
}else {
  big_table$filter_reason[(big_table$RPM <= 5) ] <-"RPM" # <- "RPM" & (big_table$KR ) & !(big_table$FP)) | (big_table$RPM <= 2 & (big_table$KNR) & !(big_table$FP)) | (big_table$RPM <= 10 & !(big_table$KR | big_table$KNR) & !(big_table$FP))
  big_table$FP[(big_table$RPM <= 5)]<- TRUE
  cat(paste0("RPM less than or equal 5: " , nrow(big_table[big_table$filter_reason == "RPM" ,] )),file=summary_file,sep="\n")
}

# Filter out peaks in error_prone regions if regions are provided
if (error_prone != "None") {
  big_bed <- big_table[, c("chrm", "start", "end", "peak", "RPM", "strand")]
  error_bed <- read_tsv(error_prone, col_select=c(1,2,3), col_names=c("chrm", "start", "end"), col_types="cii")
    
  error_peak_list <- bt.intersect(a=big_bed, b=error_bed, wa=TRUE)[, "V4"]
    
  if (truth_set != "None") {
    error_peak_list <- error_peak_list[!error_peak_list %in% true_peak_list$peak]
  }

  big_table$error_prone[big_table$peak %in% error_peak_list] <- TRUE
  big_table$filter_reason[(big_table$peak %in% error_peak_list) & !(big_table$FP)] <- "error-prone"
  big_table$FP[(big_table$peak %in% error_peak_list) & !(big_table$FP)] <- TRUE
  cat(paste0("Overlapping error-prone region: ", nrow(big_table[big_table$filter_reason == "error-prone",])), file=summary_file, sep="\n")
}

# Add truth set info to table if provided
if (truth_set != "None") {
  big_table <- merge(big_table, true_peak_list, by="peak", all.x=TRUE)
  big_table[is.na(big_table)] <- "-NA-"
}

# #testing clustering of each group after all other filters have been applied 
# library(mixtools)
#  
#  #Fit two Normals
# test_table <- big_table[,c(40,47)]
# 
# wait1 = normalmixEM(big_table$RPM_log[big_table$classification == "KR"],  k=3,verb=TRUE)
# boot.comp(big_table$RPM_log[big_table$classification == "KR"],max.comp=4, mix.type = c("normalmix"), B=100)
# plot(wait1, density=TRUE, loglik=FALSE)
# big_table$RPM_distribution[big_table$KR] <- wait1$posterior[,which.max(wait1$mu)] > wait1$posterior[,which.min(wait1$mu)]  #just to find the index of the higher distribution in k=2 situation: which.min(wait1$mu)
# # 
# # # 
# wait2 = normalmixEM(big_table$RPM[big_table$KNR],  k=2,verb=TRUE)
# plot(wait2, density=TRUE, loglik=FALSE)
# big_table$RPM_distribution[big_table$KNR] <- wait2$posterior[,which.max(wait2$mu)] > wait2$posterior[,which.min(wait2$mu)]
# # # 
# # # 
#  wait3 = normalmixEM(big_table$RPM[!(big_table$KNR | big_table$Off_target_amplification | big_table$KR)],  k=3,verb=TRUE)
#  plot(wait3, density=TRUE, loglik=FALSE)
#  big_table$RPM_distribution[big_table$KR] <- wait3$posterior[,which.max(wait3$mu)] > wait3$posterior[,which.min(wait3$mu)]
# # # 
# # # 
#  wait4 = normalmixEM(big_table$RPM[big_table$Off_target_amplification],  k=2,verb=TRUE)
#  plot(wait4, density=TRUE, loglik=FALSE)
#  big_table$RPM_distribution[big_table$KR] <- wait4$posterior[,which.max(wait4$mu)] > wait4$posterior[,which.min(wait4$mu)]
# # # 


cat(paste0("Total peaks after filtering: " , nrow(big_table[!big_table$FP,])),file=summary_file,sep="\n")
big_table$classification[big_table$FP] <- "FP"
if(library == "single"){
  big_table$classification[big_table$classification == "Candidate"] <- "UNK"
}else {  
  if(library == "bulk"){  
    big_table$classification[big_table$classification == "Candidate" & big_table$number_templates == 1 ] <- "SOM_private"
    big_table$classification[big_table$classification == "Candidate" & big_table$number_templates > 1 ] <- "SOM_clonal"
    big_table$classification[big_table$classification == "SOM_clonal" & big_table$RPM >= 5 ] <- "UNK" # UNK must have multiple templates and RPM >= 5
  }else{
    big_table$classification[big_table$classification == "Candidate" & big_table$RPM >= 5 ] <- "UNK"
    big_table$classification[big_table$classification == "Candidate" & big_table$number_templates > 1 ] <- "SOM_clonal"
    big_table$classification[big_table$classification == "Candidate" & big_table$number_templates == 1 ] <- "SOM_private"
  }
}
# }else{ #microbulk case 
#   big_table$classification[big_table$classification == "Candidate" & big_table$RPM ] 
# }
big_table$classification[(big_table$classification == "KR") & (big_table$gmotif_percent >= 0.5)  ] <- "KR_gmotif_containing"
big_table$classification[(big_table$classification == "KR") & (big_table$gmotif_percent < 0.5)  ] <- "KR_gmotif_missing"
big_table$classification[(big_table$classification == "KNR") & (big_table$gmotif_percent >= 0.5)  ] <- "KNR_gmotif_containing"
big_table$classification[(big_table$classification == "KNR") & (big_table$gmotif_percent < 0.5)  ] <- "KNR_gmotif_missing"


big_table$classification <- as.factor(big_table$classification)

cat(paste0("\nTotal KR: " , nrow(big_table[grepl("KR",big_table$classification),])),file=summary_file,sep="\n")
cat(paste0("Total KNR: " , nrow(big_table[grepl("KNR",big_table$classification),])),file=summary_file,sep="\n")
cat(paste0("Total UNK: " , nrow(big_table[big_table$classification == "UNK",])),file=summary_file,sep="\n")
if(library == "bulk"){
  cat(paste0("Total SOM_clonal: " , nrow(big_table[big_table$classification == "SOM_clonal",])),file=summary_file,sep="\n")
  cat(paste0("Total SOM_private: " , nrow(big_table[big_table$classification == "SOM_private",])),file=summary_file,sep="\n")
  
}
cat(paste0("Total Off target: " , nrow(big_table[big_table$classification == "Off_target",])),file=summary_file,sep="\n")
cat(paste0("Total FP: " , nrow(big_table[big_table$classification == "FP",])),file=summary_file,sep="\n")



filter_reasons <- ggplot(big_table[big_table$classification == "FP",], aes(x=filter_reason, group=classification, color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~classification) +
  ggtitle("Distribution of Filtering FP Peaks") + 
  ylab("Peak Count") + 
  xlab("Filter")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 

filter_reasons_KR <- ggplot(big_table[big_table$KR,], aes(x=filter_reason, group=classification, color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~classification)  + 
  ggtitle("Distribution of Filtering KR Peaks") + 
  ylab("Peak Count") + 
  xlab("Filter")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

filter_reasons_KNR <- ggplot(big_table[big_table$KNR,], aes(x=filter_reason, group=classification, color=classification, fill=classification)) +
  geom_bar() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~classification) +
  ggtitle("Distribution of Filtering KNR Peaks") + 
  ylab("Peak Count") + 
  xlab("Filter") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 


all_plots <- c(filter_chrm, filter_PTA_artifact, filter_templates, filter_RPM,filter_gmotif,filter_template_ratio,filter_chimera, filter_misalignment, filter_polyA_spanning, filter_satellite,filter_segdup)#
pdf(plots_path)
filter_chrm
filter_PTA_artifact
filter_templates
filter_RPM
filter_gmotif
filter_template_ratio
filter_chimera
filter_misalignment
filter_polyA_spanning
filter_satellite
filter_segdup
filter_reasons
filter_reasons_KR
filter_reasons_KNR
dev.off()

# 
# #Transduction support
# transduction_table <- big_table[!(big_table$classification %in% c("FP","Off_target")) ,c("peak","custompeakname","map_to_source","map_to_target","classification")]
# transduction_table$transduction_support <- NA
# for (row in (1:nrow(transduction_table))){
#   #print("##################################################")
#   source_list <- unlist(transduction_table$map_to_source[row]) #[unlist(transduction_table$map_to_source[row]) != c("", transduction_table$custompeakname[row])]
#   source_list <- source_list[!(source_list %in%  c("", transduction_table$custompeakname[row])) ]
#   #transduction_table$map_to_target[row] <- unlist(transduction_table$map_to_target[row])[unlist(transduction_table$map_to_target[row]) != c("", transduction_table$custompeakname[row])]
#   if (length(source_list) > 0){
#     #print(source_list)
#     for (custompeak in source_list){
#       #print(transduction_table$map_to_target[transduction_table$custompeakname == custompeak])
#       if (transduction_table$custompeakname[row] %in% unlist(transduction_table$map_to_target[transduction_table$custompeakname == custompeak])){
#         if (!is.na(transduction_table$transduction_support[row])){
#           transduction_table$transduction_support[row] <- paste(transduction_table$transduction_support[row], custompeak, sep=", ")
#         }else{
#           transduction_table$transduction_support[row] <- custompeak
#         }
#       }
#     }
#   }
# }
# transduction_table <- transduction_table[!is.na(transduction_table$transduction_support) & !grepl(",",transduction_table$transduction_support),c("peak","transduction_support")]
# transduction_bed <- read.table(transduction_bed_file, sep="\t",header=FALSE,col.names = c("chr","start","end","custompeakname","elements"))
# transduction_table <- merge(transduction_table,transduction_bed,by.x="transduction_support",by.y="custompeakname")
# colnames(transduction_table) <- c("transduction_support","peak","source_chr","source_start","source_end","source_elements")
# 
# 

#Filtered table
filtered_table <- big_table[big_table$classification != "FP", c("chrm","start","end","peak","classification","strand","RPM","shape","usp","gmotif_percent","polyA_percent","repeatmasker","evrony","homopolymers","gnomad","i1gp","nyuwa","xTea")]
#filtered_table <- merge(filtered_table, transduction_table, by="peak", all.x=TRUE)
# 
# filtered_table$transduction_support[grepl("KR",filtered_table$classification)] <- NA
# filtered_table$source_chr[grepl("KR",filtered_table$classification )] <- NA
# filtered_table$source_start[grepl("KR",filtered_table$classification)] <- NA
# filtered_table$source_end[grepl("KR",filtered_table$classification)] <- NA
# filtered_table$source_elements[grepl("KR",filtered_table$classification)] <- NA

write.table(filtered_table,filtered_peaks,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,na = "-NA-")
write.table(big_table[,lapply(big_table, class) != "list"] , big_table_filter_annotation_file ,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE,na = "-NA-")

# 
# 
# 
# circos_table <- filtered_table[!is.na(filtered_table$source_chr), c("chrm","start","end","source_chr","source_start","source_end","classification")]
# #circos.initializeWithIdeogram(species = "hg38")
# colnames(circos_table) <- c("chrom1","start1","end1","chrom2","start2","end2","classification")
# circos_table_polymorphic <- circos_table[grepl("KR|KNR",circos_table$classification),]
# circos_table_novel <- circos_table[circos_table$classification %in% c("UNK","SOM"),]
# 
# png(paste0(plots_path,"_circos_plot_ref_polymorphic.png"), width = 1000, height = 1000)
# par(mar = c(1, 1, 4, 1))
# circos.initializeWithIdeogram(species = "hg38")
# 
# for (i in 1:nrow(circos_table_polymorphic)) {
#   circos.link(
#     sector.index1 = circos_table_polymorphic$chrom1[i],
#     point1 = c(circos_table_polymorphic$start1[i], circos_table_polymorphic$end1[i]),
#     sector.index2 = circos_table_polymorphic$chrom2[i],
#     point2 = c(circos_table_polymorphic$start2[i], circos_table_polymorphic$end2[i]),
#     col = "blue",
#     directional = -1
#   )
# }
# title("Reference and Polymorphic Insertions with Transductions", cex.main = 1.5)
# 
# circos.clear()
# #circos.clear()
# dev.off()
# png(paste0(plots_path,"_circos_plot_novel.png"), width = 1000, height = 1000)
# par(mar = c(1, 1, 4, 1))
# 
# circos.initializeWithIdeogram(species = "hg38")
# 
# for (i in 1:nrow(circos_table_novel)) {
#   circos.link(
#     sector.index1 = circos_table_novel$chrom1[i],
#     point1 = c(circos_table_novel$start1[i], circos_table_novel$end1[i]),
#     sector.index2 = circos_table_novel$chrom2[i],
#     point2 = c(circos_table_novel$start2[i], circos_table_novel$end2[i]),
#     col = "red",
#     directional = -1
#   )
# }
# title("Novel Insertions with Transductions", cex.main = 1.5)
# 
# circos.clear()
# #circos.clear()
# #dev.off()
# dev.off()

