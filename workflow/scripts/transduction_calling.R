library(UpSetR)
library(stringr)
library(reshape)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(circlize)

args<-commandArgs(TRUE)

big_file <- args[1]
plots_path <- args[2]
transduction_peaks <- args[3]
transduction_bed_file <- args[4]
classified_peaks <- args[5]

canonical_chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                    "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                    "chr20","chr21","chr22","chrX","chrY")

big_table <- read.table(big_file, sep="\t", skip=1, header=FALSE)
colnames(big_table) <- c("peak","read","s3p","s5p","h5p","h3p","custompeakname")
big_table$map_to_source = lapply(strsplit(paste(big_table$s3p, big_table$h3p, sep=","),","), unique)
big_table$map_to_target = lapply(strsplit(paste(big_table$s5p, big_table$h5p, sep=","),","), unique)

# Transduction support
transduction_table <- big_table[c("peak","custompeakname","map_to_source","map_to_target")]
transduction_table$transduction_support <- NA
for (row in (1:nrow(transduction_table))) {
  source_list <- unlist(transduction_table$map_to_source[row])
  source_list <- source_list[!(source_list %in%  c("", transduction_table$custompeakname[row]))]

  if (length(source_list) > 0) {
    for (custompeak in source_list) {
      if (transduction_table$custompeakname[row] %in% unlist(transduction_table$map_to_target[transduction_table$custompeakname == custompeak])) {
        if (!is.na(transduction_table$transduction_support[row])) {
          transduction_table$transduction_support[row] <- paste(transduction_table$transduction_support[row], custompeak, sep=",")
        } else {
          transduction_table$transduction_support[row] <- custompeak
        }
      }
    }
  }
}

transduction_table <- transduction_table[!is.na(transduction_table$transduction_support) & !grepl(",",transduction_table$transduction_support), c("peak","transduction_support")]
transduction_bed <- read.table(transduction_bed_file, sep="\t", header=FALSE, col.names = c("chr","start","end","custompeakname","elements"))
transduction_table <- merge(transduction_table, transduction_bed, by.x="transduction_support", by.y="custompeakname")
colnames(transduction_table) <- c("transduction_support","peak","source_chr","source_start","source_end","source_elements")

filtered_table <- read.table(classified_peaks, sep="\t", header=TRUE)
filtered_table <- filtered_table[,1:7]
colnames(filtered_table) <- c("chrm","start","end","peak","RPM","strand","classification")
filtered_table <- filtered_table[filtered_table$classification != "FP",]
filtered_table <- merge(filtered_table, transduction_table, by="peak", all.x=TRUE)

filtered_table$transduction_support[grepl("KR",filtered_table$classification)] <- NA
filtered_table$source_chr[grepl("KR",filtered_table$classification)] <- NA
filtered_table$source_start[grepl("KR",filtered_table$classification)] <- NA
filtered_table$source_end[grepl("KR",filtered_table$classification)] <- NA
filtered_table$source_elements[grepl("KR",filtered_table$classification)] <- NA
circos_table <- filtered_table[!is.na(filtered_table$source_chr), c("chrm","start","end","source_chr","source_start","source_end","classification")]
write.table(circos_table, transduction_peaks, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, na="-NA-")

colnames(circos_table) <- c("chrom1","start1","end1","chrom2","start2","end2","classification")
circos_table <- circos_table[(circos_table$chrom1 %in% canonical_chrs) & (circos_table$chrom2 %in% canonical_chrs),]
circos_table_polymorphic <- circos_table[grepl("KR|KNR",circos_table$classification),]
circos_table_novel <- circos_table[circos_table$classification %in% c("UNK","SOM","SOM_clonal","SOM_private"),]

pdf(plots_path)

par(mar = c(1, 1, 4, 1))
circos.initializeWithIdeogram(species = "hg38")

for (i in 1:nrow(circos_table_polymorphic)) {
  circos.link(
    sector.index1 = circos_table_polymorphic$chrom1[i],
    point1 = c(circos_table_polymorphic$start1[i], circos_table_polymorphic$end1[i]),
    sector.index2 = circos_table_polymorphic$chrom2[i],
    point2 = c(circos_table_polymorphic$start2[i], circos_table_polymorphic$end2[i]),
    col = "blue",
    directional = -1
  )
}
title("Polymorphic Insertions with Transductions", cex.main = 1.5)

circos.clear()

par(mar = c(1, 1, 4, 1))
circos.initializeWithIdeogram(species = "hg38")

for (i in 1:nrow(circos_table_novel)) {
  circos.link(
    sector.index1 = circos_table_novel$chrom1[i],
    point1 = c(circos_table_novel$start1[i], circos_table_novel$end1[i]),
    sector.index2 = circos_table_novel$chrom2[i],
    point2 = c(circos_table_novel$start2[i], circos_table_novel$end2[i]),
    col = "red",
    directional = -1
  )
}
title("Novel Insertions with Transductions", cex.main = 1.5)

circos.clear()
dev.off()
