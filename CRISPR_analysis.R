##Loading Packages
library(MAGeCKFlute)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(magrittr)
library(dplyr)

##Gene Summary outputs from MAGeCK pipeline
gdata.het <- read.delim("BRCA_WT.het_gene_summary.txt", check.names = F)
gdata.homo <- read.delim("BRCA_WT.homo_gene_summary.txt", check.names = F)
gdata.het_homo <- read.delim("BRCA.hetvshomo.gene_summary.txt", check.names = F)
gdata.day10_all <- read.delim("BRCA_day10.gene_summary.txt")

#sdata <- read.delim("~/Documents/BRCA_Screen/MageckResults_correctDesign/BRCA_day10.sgrna_summary.txt", check.names = F)


# # Master table from day 10
# gdata <- read.delim("~/Documents/BRCA_Screen/MageckResults_correctDesign/BRCA_day10.gene_summary.txt", check.names = F)

# Read beta scores
gdata_beta.het = ReadBeta("BRCA_WT.het_gene_summary.txt")
gdata_beta.homo = ReadBeta("BRCA_WT.homo_gene_summary.txt")
gdata_beta.het_homo = ReadBeta("BRCA.hetvshomo.gene_summary.txt")
gdata_beta.day10_all = ReadBeta("BRCA_day10.gene_summary.txt")

## Control beta scores for Het
ctrlname = "day10_WT"
treatname = "day10_het"
gdata_cc = NormalizeBeta(gdata_beta.het, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_cc)
p <- DensityView(gdata_cc, samples=c(ctrlname, treatname))
MAView(gdata_cc, ctrlname, treatname)

gdata_cc$Control = rowMeans(gdata_cc[,ctrlname, drop = FALSE])
gdata_cc$Treatment = rowMeans(gdata_cc[,treatname, drop = FALSE])
gdata_cc <- OmitCommonEssential(gdata_cc, symbol = "Gene")

p1 = ScatterView(gdata_cc, "day10_WT", "day10_het", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, 
                 label = "Gene", label.top = T, top = 10)
print(p1)

gdata_cc$Diff_het = gdata_cc$Treatment - gdata_cc$Control
gdata_cc$Rank_Het = rank(gdata_cc$Diff)
p1 = ScatterView(gdata_cc, x = "Diff_het", y = "Rank_Het", label = "Gene", 
                 top = 5, model = "rank")
print(p1)

rankdata = gdata_cc$Treatment - gdata_cc$Control
names(rankdata) = gdata_cc$Gene
RankView(rankdata, width = 5, height = 5)
RankView(rankdata, genelist = c("AMOTL2", "PAWR", "MYH9", "NF1", "YAP1",
                                "MDM2", "PTEN"), width = 5, height = 5, top = 5)

# Get stats with beta score from top hits
het_table <- rbind(gdata.het[gdata.het$Gene %in% gdata_cc[order(gdata_cc$Rank, decreasing = T), ]$Gene, ],
                   gdata.het[gdata.het$Gene %in% gdata_cc[order(gdata_cc$Rank, decreasing = F), ]$Gene, ])

# prep gdata_cc table to merge
aux <- gdata_cc[, c('Gene', 'Diff_het', 'Rank_Het')]
#het_table <- merge(het_table, aux, by = 'Gene')
het_table <- merge(gdata.het, aux, by = 'Gene')

#                           loose threshold for ORA pathway analysis of het hits                       #
het_table <- het_table[het_table$`day10_het|fdr` <= 0.05, ]
het_pos_l <- het_table[het_table$Diff_het >= 0, ]$Gene
het_neg_l <- het_table[het_table$Diff_het <= 0, ]$Gene

#                                     ORA                             #

#plotting function using ggplot
plot.go.ora <- function(ora.obj, nCategory = NULL, title = NULL) {
  ggplot(ora.obj, showCategory = nCategory, 
         aes(FoldEnrichment, reorder(Description, FoldEnrichment))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_viridis_c(guide=guide_colorbar(reverse=T)) +
    scale_size_continuous(range=c(2, 10)) +
    theme_minimal() + 
    xlab("FoldEnrichment") +
    ylab(NULL) + 
    ggtitle(title)
} 
#
parse_ratio <- function(ratio) {
  ratio <- sub("^\\s*", "", as.character(ratio))
  ratio <- sub("\\s*$", "", ratio)
  numerator <- as.numeric(sub("/\\d+$", "", ratio))
  denominator <- as.numeric(sub("^\\d+/", "", ratio))
  return(numerator/denominator)
}

ora.het.pos <- enrichGO(gene =  het_pos_l,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'SYMBOL',
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.1)
ora.het.pos <- mutate(ora.het.pos, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ora.het.pos <- mutate(ora.het.pos, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

ora.het.neg <- enrichGO(gene =  het_neg_l,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'SYMBOL',
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.1)
ora.het.neg <- mutate(ora.het.neg, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ora.het.neg <- mutate(ora.het.neg, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

plot.go.ora(ora.het.pos, nCategory = 20, title = 'WT vs Het positive')
plot.go.ora(ora.het.neg, nCategory = 20, title = 'WT vs Het negative')



# Now we have rank by pos and neg selected gene - next step filter by pvalue in het
het_table <- het_table[het_table$`day10_het|fdr` <= 0.05 & abs(het_table$Diff_het) > 0.2, ]
# Filtered table sorted by positive hits het
het_table <-het_table[order(het_table$Rank_Het, decreasing = T), ]

write.table(het_table,  sep = "\t",file ='WT_vs_Het_CRISPR_Hits_Jan2023.tsv', 
            col.names=TRUE, row.names=F)


## Control beta scores for Homo
ctrlname = "day10_WT"
treatname = "day10_homo"
gdata_cc2 = NormalizeBeta(gdata_beta.homo, samples=c(ctrlname, treatname), method="cell_cycle")

gdata_cc2$Control = rowMeans(gdata_cc2[,ctrlname, drop = FALSE])
gdata_cc2$Treatment = rowMeans(gdata_cc2[,treatname, drop = FALSE])
gdata_cc2 <- OmitCommonEssential(gdata_cc2, symbol = "Gene")

p1 = ScatterView(gdata_cc2, "day10_WT", "day10_homo", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, 
                 label = "Gene", label.top = T, top = 10)
print(p1)

gdata_cc2$Diff_Homo = gdata_cc2$Treatment - gdata_cc2$Control
gdata_cc2$Rank_Homo = rank(gdata_cc2$Diff)
p1 = ScatterView(gdata_cc2, x = "Diff_Homo", y = "Rank_Homo", label = "Gene", 
                 top = 5, model = "rank")
print(p1)

rankdata = gdata_cc$Treatment - gdata_cc$Control
names(rankdata) = gdata_cc$Gene
RankView(rankdata, width = 5, height = 5)
RankView(rankdata, genelist = c("AMOTL2", "PAWR", "MYH9", "NF1", "YAP1",
                                "MDM2", "PTEN"), width = 5, height = 5, top = 5)


# Get stats with beta score from top hits
homo_table <- rbind(gdata.homo[gdata.homo$Gene %in% gdata_cc2[order(gdata_cc2$Rank, decreasing = T), ]$Gene, ],
                    gdata.homo[gdata.homo$Gene %in% gdata_cc2[order(gdata_cc2$Rank, decreasing = F), ]$Gene, ])

# prep gdata_cc table to merge
aux <- gdata_cc2[, c('Gene', 'Diff_Homo', 'Rank_Homo')]
#homo_table <- merge(homo_table, aux, by = 'Gene')
homo_table <- merge(gdata.homo, aux, by = 'Gene')

# loose threshold for ORA pathway analysis of homo hits #
homo_table <- homo_table[homo_table$`day10_homo|fdr` <= 0.05, ]
homo_pos_l <- homo_table[homo_table$Diff_Homo >= 0, ]$Gene
homo_neg_l <- homo_table[homo_table$Diff_Homo <= 0, ]$Gene
#
#                   ORA                    #
ora.homo.pos <- enrichGO(gene =  homo_pos_l,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.1)
ora.homo.pos <- mutate(ora.homo.pos, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ora.homo.pos <- mutate(ora.homo.pos, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))


ora.homo.neg <- enrichGO(gene =  homo_neg_l,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.1)
ora.homo.neg <- mutate(ora.homo.neg, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ora.homo.neg <- mutate(ora.homo.neg, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

plot.go.ora(ora.homo.pos, nCategory = 20, title = 'WT vs Homo positive')
plot.go.ora(ora.homo.neg, nCategory = 20, title = 'WT vs Homo negative')

#

# Now we have rank by pos and neg selected gene - next step filter by pvalue in homo
homo_table <- homo_table[homo_table$`day10_homo|fdr` <= 0.05 & abs(homo_table$Diff_Homo) > 0.2, ]
# Filtered table sorted by positive hits homo
homo_table <-homo_table[order(homo_table$Rank, decreasing = T), ]
##                          ##

write.table(homo_table,  sep = "\t",file ='WT_vs_Homo_CRISPR_Hits_Jan2023.tsv', 
            col.names=TRUE, row.names=F)

## Combined table with het and homo pos and neg selected genes
het_homo_table <- merge(het_table, homo_table[, c('Gene', 'Diff_Homo', 'Rank_Homo')], by = 'Gene')

write.table(het_homo_table,  sep = "\t",file ='Het_Homo_CRISPR_Hits_Jan2023.tsv', 
            col.names=TRUE, row.names=F)


### Check protein data for those hits
brca.prot.all <- read.csv('ProteinLevelQuant_msstats_DeepikaTotalProtein.csv')
brca.prot.all <-
  brca.prot.all %>%
  mutate(Gene.Name.s. = strsplit(as.character(Gene.Name.s.), " ")) %>%
  unnest(Gene.Name.s.)

brca.prot.all <- as.data.frame(brca.prot.all[brca.prot.all$Gene.Name.s. %in% het_homo_table$Gene, ])
#zcore
brca.prot.all[, 4:9] <- t(scale(t(brca.prot.all[, 4:9])))
#logFC Het and Homo vs WT
brca.prot.all$FC_Het <- (brca.prot.all$HeteroZ_24H + brca.prot.all$HeteroZ_48H) / (brca.prot.all$WT_24H + brca.prot.all$WT_48H)
brca.prot.all$FC_Homo <- (brca.prot.all$HomoZ_24H + brca.prot.all$HomoZ_48H) / (brca.prot.all$WT_24H + brca.prot.all$WT_48H)


### Find overlap with any DEG in comparisons of BRCA1 rnaseq 
deg.list <- readRDS("DEG_List.rds")

lapply(deg.list, function(x) intersect(het_homo_table$Gene, x))



### RNASeq
                                                                 ##
ft282_brca1_ocr <- readRDS("FT282-BRCA1_OCR_vs_Controls_Jan2023.rds")
common_genes <- ft282_brca1_ocr[ft282_brca1_ocr$Gene_ID %in% het_homo_table$Gene, ]

## Het hits non-homo
#het_only <- het_table[het_table$Gene %!in% het_homo_table$Gene, ]
het_only <- het_table[!(het_table$Gene %in% het_homo_table$Gene), ]
het_only <- het_only[order(het_only$Rank_Het, decreasing = T), ]


## Homo hits non-het
#homo_only <- homo_table[homo_table$Gene %!in% het_homo_table$Gene, ]
homo_only <- homo_table[!(homo_table$Gene %in% het_homo_table$Gene), ]
homo_only <- homo_only[order(homo_only$Rank_Homo, decreasing = T), ]

write.table(het_only,  sep = "\t",file ='Het_ONLY_CRISPR_Hits_Jan2023.tsv', 
            col.names=TRUE, row.names=F)
write.table(homo_only,  sep = "\t",file ='Homo_ONLY_CRISPR_Hits_Jan2023.tsv', 
            col.names=TRUE, row.names=F)

common_genes <- ft282_brca1_ocr[ft282_brca1_ocr$Gene_ID %in% het_homo_table$Gene, ]
common_genes.1 <- ft282_brca1_ocr[ft282_brca1_ocr$Gene_ID %in% het_table$Gene, ]
common_genes.2 <- ft282_brca1_ocr[ft282_brca1_ocr$Gene_ID %in% homo_table$Gene, ]
common_genes.3 <- ft282_brca1_ocr[ft282_brca1_ocr$Gene_ID %in% het_only$Gene, ]
common_genes.4 <- ft282_brca1_ocr[ft282_brca1_ocr$Gene_ID %in% homo_only$Gene, ]

## Control beta scores for Het vs Homo
ctrlname = "day10_het"
treatname = "day10_homo"
gdata_cc3 = NormalizeBeta(gdata_beta.het_homo, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_cc3)
p <- DensityView(gdata_cc3, samples=c(ctrlname, treatname))
MAView(gdata_cc3, ctrlname, treatname)

gdata_cc3$Control = rowMeans(gdata_cc3[,ctrlname, drop = FALSE])
gdata_cc3$Treatment = rowMeans(gdata_cc3[,treatname, drop = FALSE])
gdata_cc3 <- OmitCommonEssential(gdata_cc3, symbol = "Gene")

p1 = ScatterView(gdata_cc3, "day10_het", "day10_homo", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, 
                 label = "Gene", label.top = T, top = 10)
print(p1)

gdata_cc3$Diff_Combined = gdata_cc3$Treatment - gdata_cc3$Control
gdata_cc3$Rank_Combined = rank(gdata_cc3$Diff)
p1 = ScatterView(gdata_cc3, x = "Diff_Combined", y = "Rank_Combined", label = "Gene", 
                 top = 5, model = "rank")
print(p1)

rankdata = gdata_cc3$Treatment - gdata_cc3$Control
names(rankdata) = gdata_cc3$Gene
RankView(rankdata, width = 5, height = 5, top = 5)
RankView(rankdata, genelist = c("AMOTL2", "PAWR", "MYH9", "NF1", "YAP1",
                                "MDM2", "PTEN"), width = 5, height = 5, top = 5)


# Get stats with beta score from top hits
homo_vs_het_table <- rbind(gdata.het_homo[gdata.het_homo$Gene %in% gdata_cc3[order(gdata_cc3$Rank, decreasing = T), ]$Gene, ],
                           gdata.het_homo[gdata.het_homo$Gene %in% gdata_cc3[order(gdata_cc3$Rank, decreasing = F), ]$Gene, ])

# prep gdata_cc3 table to merge
aux <- gdata_cc3[, c('Gene', 'Diff_Combined', 'Rank_Combined')]
#het_table <- merge(het_table, aux, by = 'Gene')
homo_vs_het_table <- merge(gdata.het_homo, aux, by = 'Gene')

#                           loose threshold for ORA pathway analysis of homo hits                       #
homo_vs_het_table <- homo_vs_het_table[homo_vs_het_table$`day10_homo|fdr` <= 0.05, ]
het_homo_pos_l <- homo_vs_het_table[homo_vs_het_table$Diff_Combined >= 0, ]$Gene
het_homo_neg_l <- homo_vs_het_table[homo_vs_het_table$Diff_Combined <= 0, ]$Gene

#                                     ORA                             #

ora.homo.het.pos <- enrichGO(gene =  het_homo_pos_l,
                             OrgDb         = org.Hs.eg.db,
                             keyType       = 'SYMBOL',
                             ont           = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.1)
ora.homo.het.pos <- mutate(ora.homo.het.pos, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ora.homo.het.pos <- mutate(ora.homo.het.pos, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

ora.homo.het.neg <- enrichGO(gene =  het_homo_neg_l,
                             OrgDb         = org.Hs.eg.db,
                             keyType       = 'SYMBOL',
                             ont           = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.1)
ora.homo.het.neg <- mutate(ora.homo.het.neg, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ora.homo.het.neg <- mutate(ora.homo.het.neg, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

plot.go.ora(ora.homo.het.pos, nCategory = 20, title = 'Het vs Homo positive')
plot.go.ora(ora.homo.het.neg, nCategory = 20, title = 'Het vs Homo negative')



# Now we have rank by pos and neg selected gene - next step filter by pvalue in het
homo_vs_het_table <- homo_vs_het_table[homo_vs_het_table$`day10_homo|fdr` <= 0.05 & abs(homo_vs_het_table$Diff_Combined) > 0.2, ]
# Filtered table sorted by positive hits het
homo_vs_het_table <- homo_vs_het_table[order(homo_vs_het_table$Rank_Combined, decreasing = T), ]

write.table(homo_vs_het_table,  sep = "\t",file ='Het_vs_Homo_CRISPR_Hits_Jan2023.tsv', 
            col.names=TRUE, row.names=F)

combined_gene_list <- intersect(het_homo_table$Gene, homo_vs_het_table$Gene)
combined_gene_list

common_genes.5 <- ft282_brca1_ocr[ft282_brca1_ocr$Gene_ID %in% homo_vs_het_table$Gene, ]


## Control beta scores for day10 WT_Het vs WT_Homo
ctrlname = "WT+plasmidvsHet"
treatname = "WT+plasmidvsHomo"
gdata_cc4 = NormalizeBeta(gdata_beta.day10_all, samples=c(ctrlname, treatname), method="cell_cycle")
head(gdata_cc4)
p <- DensityView(gdata_cc4, samples=c(ctrlname, treatname))
MAView(gdata_cc4, ctrlname, treatname)

gdata_cc4$Control = rowMeans(gdata_cc4[,ctrlname, drop = FALSE])
gdata_cc4$Treatment = rowMeans(gdata_cc4[,treatname, drop = FALSE])
gdata_cc4 <- OmitCommonEssential(gdata_cc4, symbol = "Gene")

p1 = ScatterView(gdata_cc4, "WT+plasmidvsHet", "WT+plasmidvsHomo", groups = c("top", "bottom"), auto_cut_diag = TRUE, display_cut = TRUE, 
                 label = "Gene", label.top = T, top = 10)
print(p1)

gdata_cc4$Diff_Combined = gdata_cc4$Treatment - gdata_cc4$Control
gdata_cc4$Rank_Combined = rank(gdata_cc4$Diff)
p1 = ScatterView(gdata_cc4, x = "Diff_Combined", y = "Rank_Combined", label = "Gene", 
                 top = 5, model = "rank")
print(p1)

rankdata = gdata_cc4$Treatment - gdata_cc4$Control
names(rankdata) = gdata_cc4$Gene
RankView(rankdata, width = 5, height = 5, top = 5)
RankView(rankdata, genelist = c("AMOTL2", "PAWR", "MYH9", "NF1", "YAP1",
                                "MDM2", "PTEN"), width = 5, height = 5, top = 5)


# Get stats with beta score from top hits
day10_table <- rbind(gdata.day10_all[gdata.day10_all$Gene %in% gdata_cc4[order(gdata_cc4$Rank, decreasing = T), ]$Gene, ],
                     gdata.day10_all[gdata.day10_all$Gene %in% gdata_cc4[order(gdata_cc4$Rank, decreasing = F), ]$Gene, ])

# prep gdata_cc4 table to merge
aux <- gdata_cc4[, c('Gene', 'Diff_Combined', 'Rank_Combined')]
#het_table <- merge(het_table, aux, by = 'Gene')
day10_table <- merge(gdata.day10_all, aux, by = 'Gene')

#                           loose threshold for ORA pathway analysis of homo hits                       #
day10_table <- day10_table[day10_table$WT.plasmidvsHomo.fdr <= 0.05, ]
het_homo_pos_l <- day10_table[day10_table$Diff_Combined >= 0, ]$Gene
het_homo_neg_l <- day10_table[day10_table$Diff_Combined <= 0, ]$Gene

#                                     ORA                             #

ora.homo.het.pos <- enrichGO(gene =  het_homo_pos_l,
                             OrgDb         = org.Hs.eg.db,
                             keyType       = 'SYMBOL',
                             ont           = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.1)
ora.homo.het.pos <- mutate(ora.homo.het.pos, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ora.homo.het.pos <- mutate(ora.homo.het.pos, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

ora.homo.het.neg <- enrichGO(gene =  het_homo_neg_l,
                             OrgDb         = org.Hs.eg.db,
                             keyType       = 'SYMBOL',
                             ont           = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.1)
ora.homo.het.neg <- mutate(ora.homo.het.neg, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ora.homo.het.neg <- mutate(ora.homo.het.neg, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

plot.go.ora(ora.homo.het.pos, nCategory = 20, title = 'WT_Het vs WT_Homo positive')
plot.go.ora(ora.homo.het.neg, nCategory = 20, title = 'WT_Het vs WT_Homo negative')



# Now we have rank by pos and neg selected gene - next step filter by pvalue in het
day10_table <- day10_table[day10_table$WT.plasmidvsHomo.fdr <= 0.05 & abs(day10_table$Diff_Combined) > 0.2, ]

# Filtered table sorted by positive hits het
day10_table <- day10_table[order(day10_table$Rank_Combined, decreasing = T), ]

write.table(day10_table,  sep = "\t",file ='WT.Het_vs_WT.Homo_CRISPR_Hits_Jan2023.tsv', 
            col.names=TRUE, row.names=F)


common_genes.6 <- ft282_brca1_ocr[ft282_brca1_ocr$Gene_ID %in% day10_table$Gene, ]



list.1 <- read.delim("WT.Het_vs_WT.Homo_CRISPR_Hits_Jan2023.tsv")
list.2 <- read.delim("Het_vs_Homo_CRISPR_Hits_Jan2023.tsv")

combined_gene_list <- intersect(list.1$Gene, list.2$Gene)
combined_gene_list