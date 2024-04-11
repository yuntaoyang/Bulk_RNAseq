#---- Load libraries -----------------------------------------------------------
library(ComplexHeatmap)
#---- Load data ----------------------------------------------------------------
exp <- read.csv("exp.csv", row.names = 1)
metadata <- read.csv("metadata.csv")
#---- Parameters ---------------------------------------------------------------
# Color of gene expression level
color_exp <- colorRampPalette(colors = c("blue","white","red"))(100)
# Color of different groups
colAnn <- HeatmapAnnotation(df = data.frame(Group = metadata$group),
                            col = list('Group' = c('group_a' = 'blue', 
                                                   'group_b' = 'red')))
#---- Heatmap ------------------------------------------------------------------
png('heatmap.png', width = 8, height = 5, units = 'in', res = 300)
Heatmap(as.matrix(exp), 
        show_row_names = TRUE,
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        col = color_exp, 
        top_annotation=colAnn,
        heatmap_legend_param = list(title= "Z-Score(logCPM)", 
                                    title_position = "topcenter", 
                                    legend_direction="horizontal"))
dev.off()
