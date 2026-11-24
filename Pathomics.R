# Environmental Pollution 2026
# Pathomics Analysis
# PMID: 43057449
#
# Description:
# This script  generates figures of pathomics presented in the
# Environmental Pollution manuscript.
#
# Input:
# data/pathomics


library(readxl)       
library(mixOmics)     
library(ggplot2)     
library(ggrepel)     
library(grid)         


########## Figure 1f#############
pathomics_data <- read_excel("data/Pathomics/pathomics_data.xlsx")

data_clean <- pathomics_data[, !colnames(pathomics_data) %in% c("group", "file_name")]


X <- data_clean  
Y <- factor(pathomics_data$group) 

X <- scale(X)


plsda_model <- plsda(X, Y, ncomp = 2)

explained_var <- plsda_model$prop_expl_var$X * 100


score_df <- data.frame(
  comp1 = plsda_model$variates$X[, 1],
  comp2 = plsda_model$variates$X[, 2],
  Group = Y
)

vip_scores <- vip(plsda_model)[, 1]
vip_df <- data.frame(Variable = names(vip_scores), VIP = vip_scores)
top_vars <- vip_df[order(-vip_df$VIP), ][1:6, ]  

scale_factor <- max(abs(score_df$comp1), abs(score_df$comp2)) /
  max(abs(plsda_model$loadings$X[, 1]),
      abs(plsda_model$loadings$X[, 2])) * 1.2

arrow_df <- data.frame(
  Variable = rownames(plsda_model$loadings$X),
  LD1 = plsda_model$loadings$X[, 1] * scale_factor,
  LD2 = plsda_model$loadings$X[, 2] * scale_factor
)
arrow_df <- arrow_df[arrow_df$Variable %in% top_vars$Variable, ]


group_colors <- c(
  "FAsham" = "#f2af3e",
  "Pmsham" = "#93c043",
  "FAIRI"  = "#6fc2d7",
  "PMIRI"  = "#6b519a"
)


ggplot() +

  geom_point(data = score_df,
             aes(x = comp1, y = comp2, color = Group),
             size = 6) +
  stat_ellipse(data = score_df,
               aes(x = comp1, y = comp2, color = Group),
               level = 0.95, linetype = 2, size = 1) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = LD1, yend = LD2),
               arrow = arrow(length = unit(0.35, "cm")),
               color = "lightcoral", size = 0.8) +

  geom_text_repel(data = arrow_df,
                  aes(x = LD1, y = LD2, label = Variable),
                  color = "lightcoral", size = 4) +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  labs(
    title = "PLS-DA Group Separation with Top Contributors",
    x = paste0("PLS-DA Comp1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PLS-DA Comp2 (", round(explained_var[2], 1), "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

