#load library
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel)
library(ggplot2)
library(dplyr)
library(readxl)


data<-read_xlsx("DEG.xlsx",sheet = 1)%>%as.data.frame()

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
data$diffexpressed <- "Not significant"

# if log2Foldchange > 0.15 and pvalue < 0.05, set as "UP"
data$diffexpressed[data$avg_log2FC >= 0.15 & data$p_val < 0.05] <- "Upregulated"

# if log2Foldchange < -0.15 and pvalue < 0.05, set as "DOWN"
data$diffexpressed[data$avg_log2FC <= -0.15 & data$p_val < 0.05] <- "Downregulated"

data$diffexpressed[data$p_val < 0.05 & !(data$diffexpressed %in% c("Upregulated", "Downregulated"))] <- "Significant"

## For volcano plot with top5 up and downgenes only
top_downregulated <- data %>%
  filter(diffexpressed == "Downregulated") %>%
  arrange(avg_log2FC)%>% head(5)
  
         

top_upregulated <- data %>%
  filter(diffexpressed == "Upregulated") %>%
  arrange( avg_log2FC) %>% tail(5)
  


top_genes <- bind_rows(
  top_downregulated,
  top_upregulated
)
data$delabel <- ifelse(data$Genes %in% top_genes$Genes, data$Genes, NA)

# category for legend
data$diffexpressed <- factor(data$diffexpressed, 
                             levels = c("Upregulated","Downregulated","Significant", "Not significant"))


# Define color palette
color_palette <- c(
  "Upregulated" = "red",
  "Downregulated" = "green",
  "Significant" = "lightblue",
  "Not significant" = "grey"
)

# Create a new column for label colors based on diffexpressed
data <- data %>%
  mutate(label_color = case_when(
    diffexpressed == "Upregulated" ~ "red",
    diffexpressed == "Downregulated" ~ "green",
    diffexpressed == "Significant" ~ "lightblue",
    TRUE ~ "black"  # Default color for labels not matching the above categories
  ))


p<-ggplot(data = data, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_vline(xintercept = c(-0.15, 0.15), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + 
  geom_point(aes(col = diffexpressed), size = 2) +  # Plot all genes with color mapping
  geom_label_repel(
    aes(label = delabel, fill = diffexpressed),  # Fill the box with color for categories
    color = data$label_color,        # Set label text color to match the points
    box.padding = 0.75,             
    point.padding = 0.75,           
    segment.color = 'black',         
    size = 4,                       
    min.segment.length = 0,         
    force = 2,                      
    max.overlaps = Inf,             
    force_pull = 0.5,               
    arrow = arrow(length = unit(0.01, "npc")), 
    fill = "white",                # Keep the box fill white
    label.size = 0.3,               # Border thickness of the label box
    fontface = "bold"               # Make the text bold
  ) +  
  scale_fill_manual(values = color_palette) +  # Ensure the fill uses the same color palette
  scale_color_manual(values = color_palette,
                     labels = c("Upregulated", "Downregulated", "Significant (p-value < 0.05)", "Not significant")) +
  labs(x = expression("log"[2]*"(FoldChange)"), y = expression("-log"[10]*"p-value"), color = "Category") +
  ggtitle('A vs. B') +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "right",
    legend.background = element_blank(), # Remove background
    legend.box.background = element_blank(), # Remove box background
    legend.key = element_blank(), # Remove the key (color swatch) background
    legend.title = element_text(face = "bold"), # Bold legend title
    legend.text = element_text(face = "bold") # Bold legend text
  )

ggsave("volcano.pdf", plot = p, width = 15, height = 15) 
