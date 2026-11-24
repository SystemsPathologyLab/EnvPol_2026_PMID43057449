# Environmental Pollution 2026
# RNA-seq Analysis
# PMID: 43057449
#
# Description:
# This script generates figures presented in the
# Environmental Pollution manuscript.
#
# Input:
# data/RNAseq


library(DESeq2)
library(DT)
library(dplyr)
library(tibble)
library(ggplot2)
library(plotly)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(scales)
library(fgsea)
library(ComplexHeatmap)

# Load summer RNA-seq count matrix


summer_read_count <- read.table(
  "data/RNAseq/GSE319410_summer_raw_counts.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

summer_meta <- read_excel("data/RNAseq/META.xlsx", sheet = 1, col_names = TRUE)

summer_read_count <- summer_read_count[
  rowSums(summer_read_count >= 15) >= 1,
]

summer_read_count <- summer_read_count[
  !grepl("Mir", rownames(summer_read_count), ignore.case = TRUE),
]

summer_read_count <- as.matrix(summer_read_count)

dim(summer_read_count)

# DEG analysis using DESeq2

ddsFullCountTable_summer <- DESeqDataSetFromMatrix(
  countData = round(summer_read_count),
  colData = summer_meta,
  design = ~combination)

dds_summer<-DESeq(ddsFullCountTable_summer)
rld_summer <- rlogTransformation(dds_summer)
expression_summer<-assay(rld_summer)



summer_sham_pm_v_summer_sham_fa37<-results(dds_summer,contrast=c("combination", "summer_sham_PM", "summer_sham_FA"),pAdjustMethod = "BH")

summer_iri_pm_v_summer_iri_fa37<-results(dds_summer,contrast=c("combination", "summer_IRI_PM", "summer_IRI_FA"),pAdjustMethod = "BH")



################ Supplementary Figure S1b ################

summer_sham_pm_v_summer_sham_fa37 <- as.data.frame(summer_sham_pm_v_summer_sham_fa37)
summer_sham_pm_v_summer_sham_fa37$neglog <- -log10(summer_sham_pm_v_summer_sham_fa37$pvalue)
summer_sham_pm_v_summer_sham_fa37$GeneName <- rownames(summer_sham_pm_v_summer_sham_fa37)

summer_sham_pm_v_summer_sham_fa37.all <- summer_sham_pm_v_summer_sham_fa37


set.seed(123)

pathways.hallmark <- gmtPathways(
  "data/RNAseq/mh.all.v2023.2.Mm.symbols.gmt"
)


PMSHAM_v_FASHAM <- summer_sham_pm_v_summer_sham_fa37.all


PMSHAM_v_FASHAM <- PMSHAM_v_FASHAM[,-c(1, 3:7)]
colnames(PMSHAM_v_FASHAM) <- c("logFC", "genes")
PMSHAM_v_FASHAM <- PMSHAM_v_FASHAM[, c("genes", "logFC")]
PMSHAM_v_FASHAM <- PMSHAM_v_FASHAM[order(PMSHAM_v_FASHAM$logFC, decreasing = TRUE),]
ranks <- deframe(PMSHAM_v_FASHAM)


fgseaRes.sham <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)


fgseaResTidy.sham <- fgseaRes.sham %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  dplyr::select(-ES, -nMoreExtreme) %>%
  arrange(pval)


fgseaResTidy.sham_significant <- fgseaResTidy.sham %>% 
  filter(pval < 0.05)  

ggplot(fgseaResTidy.sham_significant, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +  
  scale_fill_manual(
    values = c("TRUE" = "red", "FALSE" = "blue"),  
    name = "Enrichment Direction",
    labels = c("Increased", "Decreased") 
  ) +
  coord_flip() +
  labs(
    x = "Pathway", 
    y = "Normalized Enrichment Score (NES)",
    title = "Significant Hallmark Pathways (p < 0.05)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))



# Summer_iri_PM v Summer_IRI_FA



summer_iri_pm_v_summer_iri_fa37 <- as.data.frame(summer_iri_pm_v_summer_iri_fa37)
summer_iri_pm_v_summer_iri_fa37$neglog <- -log10( summer_iri_pm_v_summer_iri_fa37$pvalue)
summer_iri_pm_v_summer_iri_fa37$GeneName <- rownames(summer_iri_pm_v_summer_iri_fa37)
summer_iri_pm_v_summer_iri_fa37.all <- summer_iri_pm_v_summer_iri_fa37


summer_iri_pm_v_summer_iri_fa37.sig<-subset(summer_iri_pm_v_summer_iri_fa37,summer_iri_pm_v_summer_iri_fa37$pvalue<0.05)
summer_iri_pm_v_summer_iri_fa37.nonsig <- subset(summer_iri_pm_v_summer_iri_fa37, pvalue >= 0.05)

summer_iri_pm_v_summer_iri_fa37.sig.up<-subset(summer_iri_pm_v_summer_iri_fa37.sig,summer_iri_pm_v_summer_iri_fa37.sig$log2FoldChange>=0)
summer_iri_pm_v_summer_iri_fa37.sig.down<-subset(summer_iri_pm_v_summer_iri_fa37.sig,summer_iri_pm_v_summer_iri_fa37.sig$log2FoldChange<=0)

summer_iri_meta <- summer_meta[summer_meta$injury == "IRI", ]
summer_iri_expression <- expression_summer[, summer_iri_meta$sample]
sig_genes <- rownames(summer_iri_pm_v_summer_iri_fa37.sig)
summer_iri_pm_v_summer_iri_fa37_sig_normalised_expression <- summer_iri_expression[rownames(summer_iri_expression) %in% sig_genes, ]

par(mai = c(1.4, 1, 1, 1))

num_up_regulated <- length(summer_iri_pm_v_summer_iri_fa37.sig.up$GeneName)
num_down_regulated <- length(summer_iri_pm_v_summer_iri_fa37.sig.down$GeneName)

summer_iri_pm_v_summer_iri_fa37 <- summer_iri_pm_v_summer_iri_fa37[order(summer_iri_pm_v_summer_iri_fa37$pvalue), ]
top_genes <- head(summer_iri_pm_v_summer_iri_fa37, 20)  



############ Figure  1g #################

## HALLMARK

set.seed(123)

pathways.hallmark <- gmtPathways("data/RNAseq/mh.all.v2023.2.Mm.symbols.gmt")

PMIRI_v_FAIRI <- summer_iri_pm_v_summer_iri_fa37.all

# create ranks
PMIRI_v_FAIRI <- PMIRI_v_FAIRI[,-c(1, 3:7)]
colnames(PMIRI_v_FAIRI) <- c("logFC", "genes")
PMIRI_v_FAIRI <- PMIRI_v_FAIRI[, c("genes", "logFC")]
PMIRI_v_FAIRI <- PMIRI_v_FAIRI[order(PMIRI_v_FAIRI$logFC, decreasing = TRUE),]
ranks <- deframe(PMIRI_v_FAIRI)


fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)


fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  dplyr::select(-ES, -nMoreExtreme) %>% 
  arrange(pval) 

fgseaResTidy %>%
  DT::datatable(
    rownames = FALSE,
    extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf'),
      pageLength = 20,
      scrollX = TRUE  
    )
  ) %>%
  formatSignif(columns = c('pval', 'padj'), digits = 5) %>%
  formatStyle(
    'pval',
    target = 'row',
    backgroundColor = styleInterval(
      c(10^(-6), 0.05),
      colorRampPalette(c("darkred", "white"))(3)
    ),
    color = styleInterval(
      10^(-6),
      c("white", "black")
    )
  )


fgseaResTidy_significant <- fgseaResTidy %>% 
  filter(pval < 0.05)  

ggplot(fgseaResTidy_significant, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +  
  scale_fill_manual(
    values = c("TRUE" = "red", "FALSE" = "blue"),  
    name = "Enrichment Direction",
    labels = c("Increased", "Decreased") 
  ) +
  coord_flip() +
  labs(
    x = "Pathway", 
    y = "Normalized Enrichment Score (NES)",
    title = "Significant Hallmark Pathways (p < 0.05)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))




############### Figure 3a ################


summer_iri_pm_v_summer_iri_fa37_up_BP <- enrichGO (gene = rownames(summer_iri_pm_v_summer_iri_fa37.sig.up),
                                                   keyType = "SYMBOL",
                                                   OrgDb = org.Mm.eg.db,
                                                   ont = "BP", 
                                                   pAdjustMethod = "BH",
                                                   qvalueCutoff = 0.05)
x_summer_iri_pm_v_summer_iri_fa37_up_BP <- summer_iri_pm_v_summer_iri_fa37_up_BP@result

x_summer_iri_pm_v_summer_iri_fa37_up_BP2 <- x_summer_iri_pm_v_summer_iri_fa37_up_BP[,-c(3, 4, 7)]




selected_pathways <- c(
  "response to hypoxia",
  "response to decreased oxygen levels",
  "cellular response to hypoxia",
  "cellular response to oxygen levels",
  "cellular response to decreased oxygen levels"
)

df <- as.data.frame(summer_iri_pm_v_summer_iri_fa37_up_BP)

df_sel <- df %>%
  filter(Description %in% selected_pathways)

ggplot(df_sel, aes(
  x = p.adjust,
  y = reorder(Description, p.adjust)
)) +
  geom_point(aes(size = Count), color = "red") +
  labs(
    x = "Adjusted P-value",
    y = NULL,
    title = "Upregulated GO-BP Hypoxia Pathways"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 12)
  )



############### Figure 3f ################


selected_pathways <- c("glucose catabolic process",
                       "glycolytic process",
                       "glucose catabolic process to pyruvate",
                       "glycolytic process through glucose-6-phosphate",
                       "glycolytic process through fructose-6-phosphate",
                       "glucose metabolic process")

df <- as.data.frame(summer_iri_pm_v_summer_iri_fa37_up_BP)

df_sel <- df %>%
  filter(Description %in% selected_pathways) %>%
  mutate(Description = factor(Description, 
                              levels = Description[order(p.adjust, decreasing = FALSE)]))


ggplot(df_sel, aes(x = p.adjust, y = Description)) +
  geom_point(aes(size = Count), color = "red") +
  labs(
    x = "Adjusted P-value",
    y = NULL,
    title = "Upregulated GO-BP Glycolysis Pathways"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 12)
  )


############# Supplementary S2a ############


my_palette <- colorRampPalette(c("blue", "white", "red2"))(n = 200)
color_function <- function(x) {my_palette(length(x))}
breaks <- seq(-4, 4, length.out = 200)
levels <- unique(summer_iri_meta$combination)
colors <- rainbow(length(levels))
color_list <- setNames(colors, levels)


ihypoxia_mouse_kidney <- read_excel("data/RNAseq/ihypoxia_mouse_kidney.xlsx", col_names = TRUE)

common_ihypoixa_iricompr <- intersect(ihypoxia_mouse_kidney$GeneName, rownames(summer_iri_pm_v_summer_iri_fa37.sig))
print(common_ihypoixa_iricompr)

ihypoxia_iricompr_deg <- summer_iri_pm_v_summer_iri_fa37.sig %>% 
  filter(row.names(summer_iri_pm_v_summer_iri_fa37.sig) %in% common_ihypoixa_iricompr)


ihypoxia_iricompr_deg_up <- ihypoxia_iricompr_deg %>% 
  filter(log2FoldChange > 0)


num_up_regulated <- length(ihypoxia_iricompr_deg_up$GeneName)

ihypoxia_iricompr_deg_down <- ihypoxia_iricompr_deg %>% 
  filter(log2FoldChange < 0)

num_down_regulated<- length(ihypoxia_iricompr_deg_down$GeneName)



summer_common_ihypoixa_iricompr_expression <- summer_iri_expression[common_ihypoixa_iricompr, ]

summer_common_ihypoixa_iricompr_expression  <- as.matrix(summer_common_ihypoixa_iricompr_expression )
summer_common_ihypoixa_iricompr_expression_t <- t(apply(summer_common_ihypoixa_iricompr_expression, 1, scale))
colnames(summer_common_ihypoixa_iricompr_expression_t) <- colnames(summer_common_ihypoixa_iricompr_expression)


log2fc_values <- summer_iri_pm_v_summer_iri_fa37.sig[common_ihypoixa_iricompr, "log2FoldChange"]
up_down_status <- ifelse(log2fc_values > 0, "Upregulated", "Downregulated")





Heatmap(
  summer_common_ihypoixa_iricompr_expression_t,
  col = my_palette,
  column_split = summer_iri_meta$combination,
  column_title_gp = gpar(fontsize = 12),
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  row_dend_reorder = TRUE,
  cluster_row_slices = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  row_dend_side = "left",
  row_dend_width = unit(12, "mm"),
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  show_row_dend = TRUE,
  cluster_columns = FALSE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  column_dend_side = "bottom",
  column_dend_height = unit(12, "mm"),
  show_column_dend = TRUE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Expression",
    at = c(-4, -2, 0, 2, 4),
    labels = c("-4", "-2", "0", "2", "4"),
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 12)
  ),
  top_annotation = HeatmapAnnotation(
    group = summer_iri_meta$combination,
    show_annotation_name = FALSE,
    show_legend = FALSE,
    col = list(group = color_list)
  )
)





############ Supplementary Figure S3a ##############

set.seed(123)
pathways.fao <- gmtPathways( "data/RNAseq/GOBP_FATTY_ACID_BETA_OXIDATION.v2024.1.Mm.gmt")

head(pathways.fao)
head(ranks, 20)

fgseaRes <- fgsea(pathways = pathways.fao, stats = ranks) 




geneSetName <- "GOBP_FATTY_ACID_BETA_OXIDATION"
geneSet <- pathways.fao[[geneSetName]]


gseaResults <- fgseaRes[pathway == geneSetName, .(pval, padj, NES)]  


plotEnrichment(geneSet, ranks) +
  labs(
    title = paste("GSEA Plot:", geneSetName),
    subtitle = sprintf("NES = %.2f, p = %.4f (adj. p = %.4f)", 
                       gseaResults$NES, 
                       gseaResults$pval, 
                       gseaResults$padj)
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )





############ Supplementary Figure S3e ##############

summer_iri_pm_v_summer_iri_fa37_up_CC <- enrichGO (gene = rownames(summer_iri_pm_v_summer_iri_fa37.sig.up),
                                                   keyType = "SYMBOL",
                                                   OrgDb = org.Mm.eg.db,
                                                   ont = "CC", 
                                                   pAdjustMethod = "BH",
                                                   qvalueCutoff = 0.05)
x_summer_iri_pm_v_summer_iri_fa37_up_CC <- summer_iri_pm_v_summer_iri_fa37_up_CC@result

x_summer_iri_pm_v_summer_iri_fa37_up_CC2 <- x_summer_iri_pm_v_summer_iri_fa37_up_CC[,-c(3, 4, 7)]





ggplot(summer_iri_pm_v_summer_iri_fa37_up_CC, 
       aes(x = p.adjust, 
           y = Description),
       showCategory = 20) +
  geom_point(aes(size = Count), color = "red") +
  labs(
    x = "Adjusted P-value",
    y = NULL,
    title = "Upregulated GO-CC Pathways"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 12)
  )



############ Figure 4b ##############

set.seed(123)

pathways.lcp <- gmtPathways( "data/RNAseq/GOBP_LIPID_CATABOLIC_PROCESS.v2024.1.Mm.gmt")



fgseaRes <- fgsea(pathways = pathways.lcp, stats = ranks) 


geneSetName <- "GOBP_LIPID_CATABOLIC_PROCESS"
geneSet <- pathways.lcp[[geneSetName]]


gseaResults <- fgseaRes[pathway == geneSetName, .(pval, padj, NES)]  


plotEnrichment(geneSet, ranks) +
  labs(
    title = paste("GSEA Plot:", geneSetName),
    subtitle = sprintf("NES = %.2f, p = %.4f (adj. p = %.4f)", 
                       gseaResults$NES, 
                       gseaResults$pval, 
                       gseaResults$padj)
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )



############### Figure 4f ################

num_up_regulated <- length(summer_iri_pm_v_summer_iri_fa37.sig.up$GeneName)
num_down_regulated <- length(summer_iri_pm_v_summer_iri_fa37.sig.down$GeneName)

summer_iri_pm_v_summer_iri_fa37_plot <- ggplot(summer_iri_pm_v_summer_iri_fa37,
                                               aes(x=log2FoldChange, y=neglog)) + 
  geom_point(data=summer_iri_pm_v_summer_iri_fa37.sig.up, 
             aes(x=log2FoldChange, y=neglog), 
             color="darkorange3", 
             size=1.2) +
  geom_point(data=summer_iri_pm_v_summer_iri_fa37.sig.down, 
             aes(x=log2FoldChange, y=neglog), 
             color="mediumorchid3", 
             size=1.2) +
  geom_point(data=summer_iri_pm_v_summer_iri_fa37.nonsig, 
             aes(x=log2FoldChange, y=neglog), 
             color="grey", 
             size=1.2) +
  geom_text_repel(data=top_genes, 
                  aes(x=log2FoldChange, y=neglog, label=GeneName), 
                  size=3, 
                  box.padding=0.35, 
                  point.padding=0.5, 
                  segment.color='grey50') +
  xlim(-5, 5) +   ylim(0, 5) + 
  scale_color_manual(values=c("darkorange3", "mediumorchid3", "grey")) +
  labs(x="Expression (log2FC)", 
       y="-log10(p-value)", 
       title="Summer_IRI_PM_v_Summer_IRI_FA DEG") +
  theme_minimal() + 
  guides(color=FALSE) +  
  annotate("text", 
           x=Inf, 
           y=Inf, 
           label=paste("Up:", num_up_regulated), 
           hjust=1.1, 
           vjust=2, 
           color="darkorange3", 
           size=4) +
  annotate("text", 
           x=-Inf, 
           y=Inf, 
           label=paste("Down:", num_down_regulated), 
           hjust=-0.1, 
           vjust=2, 
           color="mediumorchid3", 
           size=4)

print(summer_iri_pm_v_summer_iri_fa37_plot)


############### Figure 5a ################

pathways.senescence <- read_excel("data/RNAseq/SenMayo_geneset.xlsx", sheet = 2)
pathways.senescence <- pathways.senescence[,-c(2:4)]
colnames(pathways.senescence) <- c("SENESCENCE_SenMayo")
pathways.senescence1<-as.list(pathways.senescence)
pathways<-pathways.senescence1
head(pathways.senescence1, 20)


fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()


geneSetName <- "SENESCENCE_SenMayo"
geneSet <- pathways.senescence[[geneSetName]]

gseaResults <- fgseaRes[pathway == geneSetName, .(pval, padj, NES)]  

plotEnrichment(geneSet, ranks) +
  labs(
    title = paste("GSEA Plot:", geneSetName)
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )


############### supplementary Figure 5a ##############

pathways.rtos <- gmtPathways("data/RNAseq/GOBP_RESPONSE_TO_OXIDATIVE_STRESS.v2026.1.Mm.gmt")

fgseaRes <- fgsea(pathways = pathways.rtos, stats = ranks) 

geneSetName <- "GOBP_RESPONSE_TO_OXIDATIVE_STRESS"
geneSet <- pathways.rtos[[geneSetName]]


gseaResults <- fgseaRes[pathway == geneSetName, .(pval, padj, NES)]  


plotEnrichment(geneSet, ranks) +
  labs(
    title = paste("GSEA Plot:", geneSetName),
    subtitle = sprintf("NES = %.2f, p = %.4f (adj. p = %.4f)", 
                       gseaResults$NES, 
                       gseaResults$pval, 
                       gseaResults$padj)
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )




################ Figure 7b ##############

library(STRINGdb)
library(visNetwork)
library(igraph)
library(dplyr)
library(biomaRt)
library(ggraph)
library(tidygraph)
library(readxl)
library(tidyr)




summer_iri_genes_for_drug_screening_deg <- read_excel(
  "data/RNAseq/summer_iri_genes_for drug_screening.xlsx", 
  sheet = "DEG"
) %>% as.data.frame()


download.file(
  url = "http://stringdb-downloads.org/download/protein.aliases.v12.0/10090.protein.aliases.v12.0.txt.gz",
  destfile = "10090.protein.aliases.v12.0.txt.gz",
  method = "curl",
  extra = "-k"
)

download.file(
  url = "http://stringdb-downloads.org/download/protein.info.v12.0/10090.protein.info.v12.0.txt.gz",
  destfile = "10090.protein.info.v12.0.txt.gz",
  method = "curl",
  extra = "-k"
)

download.file(
  url = "http://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz",
  destfile = "10090.protein.links.v12.0.txt.gz",
  method = "curl",
  extra = "-k"
)


string_db <- STRINGdb$new(
  version = "12.0",
  species = 10090,  # Mouse species code
  score_threshold = 600, # Minimum combined score for interactions
  input_directory = getwd()
)



mapped_genes <- string_db$map(
  summer_iri_genes_for_drug_screening_deg, 
  "GeneName", 
  removeUnmappedRows = TRUE
)

ppi_plot <- string_db$plot_network(mapped_genes$STRING_id) 

ppi_data <- string_db$get_interactions(mapped_genes$STRING_id)

connected_genes <- unique(c(ppi_data$from, ppi_data$to))
connected_mapped_genes <- mapped_genes[mapped_genes$STRING_id %in% connected_genes, ]
ppi_plot <- string_db$plot_network(connected_mapped_genes$STRING_id)


ppi_data <- string_db$get_interactions(connected_mapped_genes$STRING_id)
ppi_graph <- graph_from_data_frame(
  ppi_data[, c("from", "to", "combined_score")],
  directed = FALSE
)

degree_values <- degree(ppi_graph)
betweenness_values <- betweenness(ppi_graph)
closeness_values <- closeness(ppi_graph)


hub_genes <- data.frame(
  STRING_id = names(degree_values),
  Degree = degree_values,
  Betweenness = betweenness_values,
  Closeness = closeness_values
) %>%
  merge(connected_mapped_genes[, c("STRING_id", "GeneName")], by = "STRING_id")


hub_genes$Degree_z <- scale(hub_genes$Degree)
hub_genes$Betweenness_z <- scale(hub_genes$Betweenness)
hub_genes$Closeness_z <- scale(hub_genes$Closeness)


hub_genes$WeightedScore <- (
  hub_genes$Degree_z * 0.5 +
    hub_genes$Betweenness_z * 0.3 +
    hub_genes$Closeness_z * 0.2
)


top_hub_genes <- hub_genes %>%
  arrange(desc(WeightedScore)) %>%
  head(50)


drug_data <- read.table(
  "data/RNAseq/repurposing_drugs_20200324.txt", 
  header = TRUE, 
  sep = "\t",         
  fill = TRUE,        
  quote = "",       
  comment.char = ""
) %>%
  separate_rows(target, sep = "\\|") %>%
  filter(!is.na(target), target != "")


hub_drug_targets <- drug_data %>%
  filter(target %in% top_hub_genes$GeneName) %>%
  group_by(target) %>%
  summarise(
    Drug_Count = n(),
    Drugs = paste(unique(pert_iname), collapse = "|")
  ) %>%
  arrange(desc(Drug_Count))


final_result <- merge(
  top_hub_genes,
  hub_drug_targets,
  by.x = "GeneName",
  by.y = "target",
  all.x = TRUE
) %>%
  filter(!is.na(Drug_Count)) %>%  
  arrange(desc(Degree))


print(final_result)

