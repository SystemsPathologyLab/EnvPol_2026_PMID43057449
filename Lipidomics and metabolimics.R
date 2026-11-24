# Environmental Pollution 2026
# lIPIDOMICS and METABOLICS Analysis
# PMID: 43057449
#
# Description:
# This script  generates figures of lipidomics and metbaolommics presented in the
# Environmental Pollution manuscript.
#
# Input:
# data/LIPIDOMICS_METABOLICS




library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)

lipid_data <- read_excel("data/LIPIDOMICS_METABOLICS/Statistical_Data_lipidomics.xlsx")
metabolmic <-  read_excel("data/LIPIDOMICS_METABOLICS/Statistical_Data_metabolics.xlsx")
meta <- read_excel("data/LIPIDOMICS_METABOLICS/meta.xlsx")
metabolimc_pathway <- read_excel("data/LIPIDOMICS_METABOLICS/pathway.xlsx")



df <- lipid_data %>% 
  mutate(
    DEG = p_value < 0.05,
    regulation = case_when(
      DEG & Log2_fold_change > 0 ~ "Up",
      DEG & Log2_fold_change < 0 ~ "Down",
      TRUE ~ "NotSig"
    ),
    negLog10P = -log10(p_value),
    label = ifelse(regulation %in% c("Up", "Down"), name, NA)
  )

top30 <- df %>%
  filter(DEG) %>%
  arrange(p_value) %>%
  slice(1:20)


up_n   <- sum(df$regulation == "Up")
down_n <- sum(df$regulation == "Down")

x_max <- max(df$Log2_fold_change, na.rm = TRUE)
x_min <- min(df$Log2_fold_change, na.rm = TRUE)
y_max <- max(df$negLog10P, na.rm = TRUE)



########### Supplementary Figure S3c ############

ggplot(df, aes(x = Log2_fold_change, y = negLog10P, color = regulation)) +
  geom_point(size = 3) +
  geom_text_repel(
    data = top30,
    aes(label = name),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2
  ) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
  theme_classic(base_size = 14) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Regulation",
    title = "Volcano Plot"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  annotate("text", x = x_max, y = y_max, 
           label = paste0("Up: ", up_n),
           hjust = 1, vjust = 1, size = 5, color = "red") +
  annotate("text", x = x_min, y = y_max, 
           label = paste0("Down: ", down_n),
           hjust = 0, vjust = 1, size = 5, color = "blue")


########### Figure 4C ############

df_deg <- lipid_data %>%
  filter(p_value < 0.01)

df_selected <- df_deg %>%
  dplyr::select(name, starts_with("E22_S")) %>%
  as.data.frame()

rownames(df_selected) <- df_selected$name
df_selected$name <- NULL

expr_mat <- as.matrix(df_selected)

expr_mat_scaled <- t(apply(expr_mat, 1, scale))
colnames(expr_mat_scaled) <- colnames(expr_mat)



rownames(meta) <- meta$sample

sample_group <- meta[colnames(expr_mat_scaled), "group", drop = TRUE]
sample_group <- factor(sample_group, levels = c("IRI+FA", "IRI+PM2.5"))

group_colors <- c("IRI+FA" = "#4DAF4A",
                  "IRI+PM2.5" = "#984EA3")

top_anno <- HeatmapAnnotation(
  Group = sample_group,
  col = list(Group = group_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE
)



my_palette <- colorRamp2(c(-2, 0, 2),
                         c("blue", "white", "red"))

Heatmap(
  expr_mat_scaled,
  col = my_palette,
  column_split = sample_group,
  top_annotation = top_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = F,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  row_dend_side = "left",
  row_dend_width = unit(12, "mm"),
  column_dend_side = "bottom",
  column_dend_height = unit(12, "mm"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Expression (Z-score)",
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2")
  )
)



########### Figure 4d ############

library(mixOmics)
library(ggplot2)

X <- as.matrix(plsda_model$input.X)
Y <- droplevels(as.factor(plsda_model$Y))

X_clean <- X[, colSums(!is.finite(X)) == 0, drop = FALSE]
X_clean <- X_clean[, apply(X_clean, 2, sd) > 0, drop = FALSE]

set.seed(123)

plsda_model_clean <- plsda(
  X = X_clean,
  Y = Y,
  ncomp = 2,
  scale = TRUE,
  near.zero.var = TRUE
)

set.seed(123)

perf_plsda <- perf(
  plsda_model_clean,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = FALSE
)

model <- plsda_model_clean

scores <- as.data.frame(model$variates$X[, 1:2, drop = FALSE])
colnames(scores) <- c("Comp1", "Comp2")
scores$Comp1 <- -scores$Comp1
scores$Group <- Y

explained_var <- model$prop_expl_var$X * 100

x_lab <- paste0("PLS-DA Comp1 (", round(explained_var[1], 1), "%)")
y_lab <- paste0("PLS-DA Comp2 (", round(explained_var[2], 1), "%)")

p_plsda <- ggplot(
  scores,
  aes(x = Comp1, y = Comp2, color = Group)
) +
  stat_ellipse(
    type = "norm",
    level = 0.95,
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    size = 2.8,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "IRI+FA" = "#8DD3D1",
      "IRI+PM2.5" = "#7953A9"
    )
  ) +
  labs(
    title = "PLS-DA Group Separation\nfrom Lipidomics",
    x = x_lab,
    y = y_lab,
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(
      color = "#E8E8E8",
      linewidth = 0.3
    ),
    panel.grid.minor = element_line(
      color = "#F2F2F2",
      linewidth = 0.2
    ),
    panel.border = element_rect(
      color = "#BDBDBD",
      linewidth = 0.5
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 12
    ),
    axis.title = element_text(
      size = 11,
      color = "black"
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.key = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(size = 3)
    )
  )

print(p_plsda)


################# Figure 4f ######################


vip_scores <- mixOmics::vip(plsda_model_clean)


vip_df <- data.frame(
  Lipid = rownames(vip_scores),
  VIP_Comp1 = vip_scores[, 1],
  VIP_Comp2 = vip_scores[, 2]
)


vip_df <- vip_df[order(-vip_df$VIP_Comp2), ]


vip_selected <- subset(vip_df, VIP_Comp2 > 1)

top18_vip <- head(
  vip_df[order(-vip_df$VIP_Comp2), ],
  18
)

top18_vip$Lipid <- factor(
  top18_vip$Lipid,
  levels = top20_vip$Lipid
)

p_vip <- ggplot(
  top18_vip,
  aes(x = Lipid, y = VIP_Comp2)
) +
  geom_col(
    fill = "grey60",
    color = "grey30",
    linewidth = 0.3,
    width = 0.8
  ) +
  scale_y_continuous(
    breaks = pretty(c(0, max(top18_vip$VIP_Comp2)), n = 6),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "VIP score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      color = "black",
      size = 9
    ),
    axis.text.y = element_text(
      color = "black",
      size = 10
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 12
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      color = "black"
    ),
    panel.grid = element_blank()
  )

print(p_vip)




################# Supplementary Figure2c ################

metabolites_data_2 <- metabolmic


metabolites_data_2 <- as.data.frame(metabolites_data_2)
metabolites_data_2$negpvalue <- -log10(metabolites_data_2$P_Value)


metabolites_all_sig<-subset(metabolites_data_2, metabolites_data_2$P_Value <=0.05)


metabolites_all_sig_up<-subset(metabolites_all_sig, metabolites_all_sig$Log2_fold_change >= 0)
metabolites_all_sig_down<-subset(metabolites_all_sig, metabolites_all_sig$Log2_fold_change < 0)



metabolites_all_nonsig <- subset(metabolites_data_2, metabolites_data_2$P_Value > 0.05)

num_up_regulated <- length(metabolites_all_sig_up$HMDB.name)
num_down_regulated <- length(metabolites_all_sig_down$HMDB.name)


metabolites_all_sig_up <- metabolites_all_sig_up[order(metabolites_all_sig_up$P_Value), ]
top_up_genes <- head(metabolites_all_sig_up, 20)  

metabolites_all_sig_down <- metabolites_all_sig_down[order(metabolites_all_sig_down$P_Value), ]
top_down_genes <- head(metabolites_all_sig_down, 20)

 
metabolites_vol <- ggplot(metabolites_data_2,
                          aes(x=Log2_fold_change, y=negpvalue)) + 
  geom_point(data=metabolites_all_sig_up, 
             aes(x=Log2_fold_change, y=negpvalue), 
             color="red", 
             size=2) +
  geom_point(data=metabolites_all_sig_down, 
             aes(x=Log2_fold_change, y=negpvalue), 
             color="blue", 
             size=2) +
  geom_point(data=metabolites_all_nonsig, 
             aes(x=Log2_fold_change, y=negpvalue), 
             color="grey", 
             size=2) +
  geom_text_repel(data=top_up_genes, 
                  aes(x=Log2_fold_change, y=negpvalue, label=HMDB.name), 
                  size=4, 
                  box.padding=0.35, 
                  point.padding=0.5, 
                  segment.color='grey50') +
  geom_text_repel(data=top_down_genes, 
                  aes(x=Log2_fold_change, y=negpvalue, label=HMDB.name), 
                  size=4, 
                  box.padding=0.35, 
                  point.padding=0.5, 
                  segment.color='grey50') +
  xlim(-2, 2) +   ylim(0, 2) + 
  scale_color_manual(values=c("red", "blue", "grey")) +
  labs(x="Expression (log2FC)", 
       y="-log10(pvalue)", 
       title="Metabolites") +
  geom_vline(xintercept = 0, lty=3,color = 'black', lwd=0.5) +   geom_hline(yintercept = -log10(0.05), lty=3,color = 'black', lwd=0.5)+
  theme_minimal() + 
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    plot.title   = element_text(size = 18, hjust = 0.5)
  ) +
  guides(color=FALSE) + 
  annotate("text", 
           x=Inf, 
           y=Inf, 
           label=paste("Up:", num_up_regulated), 
           hjust=1.1, 
           vjust=2, 
           color="red", 
           size=6) +
  annotate("text", 
           x=-Inf, 
           y=Inf, 
           label=paste("Down:", num_down_regulated), 
           hjust=-0.1, 
           vjust=2, 
           color="blue", 
           size=6)

print(metabolites_vol)

################# Supplementary Figure2d ################


metabolites_data <- as.data.frame(metabolmic)

sample_cols <- grep("^(IRIFA_|IRIPM_)", colnames(metabolites_data), value = TRUE)

metab_names <- metabolites_data$Analyte.name

metabolites_mat <- as.matrix(
  metabolites_data[, sample_cols, drop = FALSE]
)

rownames(metabolites_mat) <- make.unique(metab_names)

metabolites_mat_scaled <- t(apply(metabolites_mat, 1, scale))
colnames(metabolites_mat_scaled) <- colnames(metabolites_mat)
rownames(metabolites_mat_scaled) <- rownames(metabolites_mat)

sample_group <- factor(
  ifelse(grepl("^IRIFA_", sample_cols), "IRI+FA", "IRI+PM2.5"),
  levels = c("IRI+FA", "IRI+PM2.5")
)

top_anno <- HeatmapAnnotation(
  Group = sample_group,
  col = list(Group = c(
    "IRI+FA" = "#4DAF4A",
    "IRI+PM2.5" = "#984EA3"
  )),
  show_annotation_name = FALSE
)

my_palette <- colorRamp2(
  c(-2, 0, 2),
  c("blue", "white", "red")
)

pathway_df <- as.data.frame(metabolimc_pathway)

pathway_matrix <- as.data.frame.matrix(
  table(
    factor(pathway_df$Metabolite, levels = unique(metab_names)),
    pathway_df$Pathmay
  )
)

pathway_matrix <- pathway_matrix[metab_names, , drop = FALSE]
rownames(pathway_matrix) <- rownames(metabolites_mat_scaled)
pathway_matrix[pathway_matrix > 1] <- 1

primary_pathway <- apply(pathway_matrix, 1, function(x) {
  cols <- names(which(x == 1))
  if (length(cols) == 0) return("Other")
  cols[1]
})

path_order <- c(
  "Glycolysis", "TCA", "Carnitine",
  "Lipid biosynthesis", "Amino Acids",
  "Creatine", "CoA/Serine", "PentosePhosphate",
  "NAD", "Purines", "Pyrimidines",
  "SAMSAH", "Tryptophan", "Urea", "Other"
)

row_order <- order(
  factor(primary_pathway, levels = c(
    path_order,
    setdiff(unique(primary_pathway), path_order)
  )),
  metab_names
)

metabolites_mat_scaled2 <- metabolites_mat_scaled[
  row_order, , drop = FALSE
]

pathway_matrix2 <- pathway_matrix[
  row_order, , drop = FALSE
]

pathway_anno_df <- as.data.frame(
  lapply(pathway_matrix2, function(x) factor(x, levels = c(0, 1)))
)

pathway_col_list <- lapply(pathway_anno_df, function(x) {
  c("0" = "white", "1" = "grey50")
})

pathway_annotation <- rowAnnotation(
  df = pathway_anno_df,
  col = pathway_col_list,
  annotation_name_gp = gpar(col = "black", fontsize = 8),
  show_legend = FALSE,
  border = TRUE
)

heatmap_metabolites <- Heatmap(
  metabolites_mat_scaled2,
  col = my_palette,
  column_split = sample_group,
  top_annotation = top_anno,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 7),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Scaled abundance",
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2")
  ),
  left_annotation = pathway_annotation
)

print(heatmap_metabolites)

########### Figure 3e ###########

library(mixOmics)

metab_plsda <- as.data.frame(metabolmic)


sample_cols <- c(
  paste0("IRIFA_0", 1:6),
  paste0("IRIPM_0", 1:7)
)

X <- t(as.matrix(metab_plsda[, sample_cols]))

colnames(X) <- metab_plsda$Analyte.name

meta_plsda <- as.data.frame(meta)

stopifnot(
  identical(
    as.character(meta_plsda$group),
    c(rep("IRI+FA", 6), rep("IRI+PM2.5", 7))
  )
)

Y <- factor(
  ifelse(meta_plsda$group == "IRI+PM2.5", "IRI+PM", "IRI+FA"),
  levels = c("IRI+FA", "IRI+PM")
)

is_constant <- apply(X, 2, function(x) length(unique(x)) == 1)

X <- X[, !is_constant, drop = FALSE]

X <- scale(X)

dim(X)
table(Y)

colnames(X) <- make.unique(colnames(X))

anyDuplicated(colnames(X))

plsda_model <- plsda(X, Y, ncomp = 2)

set.seed(123)
perf_plsda <- perf(plsda_model, validation = "Mfold", folds = 5,
                   progressBar = FALSE, nrepeat = 10)



score_df <- data.frame(comp1 = plsda_model$variates$X[, 1],
                       comp2 = plsda_model$variates$X[, 2],
                       Group = Y)


explained_var <- plsda_model$prop_expl_var$X * 100



vip_scores <- vip(plsda_model)[, 1]
vip_df <- data.frame(Variable = names(vip_scores), VIP = vip_scores)


top_vars <- vip_df[order(-vip_df$VIP), ][1:15, ]
loadings <- plsda_model$loadings$X
scale_factor <- 60   
arrow_df <- data.frame(
  Variable = rownames(loadings),
  LD1 = loadings[, 1] * scale_factor,
  LD2 = loadings[, 2] * scale_factor
)


arrow_df <- arrow_df[arrow_df$Variable %in% top_vars$Variable, ]


group_colors <- c("IRI+FA" = "#6fc2d7",
                  "IRI+PM" = "#6b519a")

ggplot() +
  geom_point(data = score_df,
             aes(x = comp1, y = comp2, color = Group),
             size = 4) +
  stat_ellipse(data = score_df,
               aes(x = comp1, y = comp2, color = Group),
               level = 0.95, linetype = 2, size = 1) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = LD1, yend = LD2),
               arrow = arrow(length = unit(0.35, "cm")),  
               color = "lightcoral", size = 1.2) +
  geom_text_repel(data = arrow_df,
                  aes(x = LD1, y = LD2, label = Variable),
                  color = "red", size = 4) +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  labs(
    title = "PLS-DA Group Separation with Top Contributors",
    x = paste0("PLS-DA Comp1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PLS-DA Comp2 (", round(explained_var[2], 1), "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
