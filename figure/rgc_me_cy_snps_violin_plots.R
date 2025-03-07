#
# Plot non-synonomous, non-splice FC variant frequences as a Violin plot
#
figureTitle = "FC Cytoplasmic, Non-Synonomous SNP \nAllele Frequency by Gene"
# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(scales)   # For number formatting
library(ggrepel)  # For non-overlapping text labels

# Define INPUT file name
variantsFile <- "rgc_me_variant_frequencies_chr1,19_20231004-roifgenes2-dbsnp156-GRCh38.111.2c_onlyTr-onePer_filt_norm.xls"
variantsFile <- "../../cd89paper/rgc_me_variant_frequencies_chr1,19_20231004-roifgenes2-dbsnp156-GRCh38.111.2c_onlyTr-onePer_filt_norm.xls"

# Define OUTPUT filenames
plotFilenameAll <- "rgc_me_cy_nonsyn_vars.all.png"
plotFilenameAllLog <- "rgc_me_cy_nonsyn_vars.all.log.png"

plotFilenamePerGene    <- "rgc_me_cy_nonsyn_vars.per_gene.png"
plotFilenamePerGenePDF <- "rgc_me_cy_nonsyn_vars.per_gene.pdf"

plotFilenamePerGeneNoTitles    <- "rgc_me_cy_nonsyn_vars.per_gene.no_titles.png"
plotFilenamePerGeneNoTitlesPDF <- "rgc_me_cy_nonsyn_vars.per_gene.no_titles.pdf"

plotFilenamePerGeneNoText    <- "rgc_me_cy_nonsyn_vars.per_gene.no_text.png"
plotFilenamePerGeneNoTextPDF <- "rgc_me_cy_nonsyn_vars.per_gene.no_text.pdf"

plotFilenamePerGeneLog <- "rgc_me_cy_nonsyn_vars.per_gene.log.png"

plotFilenameHistLinear <- "rgc_me_cy_nonsyn_vars.hist.linear.png"
plotFilenameHistLog <- "rgc_me_cy_nonsyn_vars.hist.log.png"

# Read in the TSV file
cat(paste0("READING: ", variantsFile, "\n"))
variantsDf <- read_tsv(variantsFile, col_types = cols())
cat(paste0("LOADED: ", paste0(dim(variantsDf),collapse="x"), "\n"))

# Filter data
filteredVariants <- variantsDf %>%
  filter(grepl("CY", DOMAIN, ignore.case = TRUE)) %>%  # Keep rows with "CY" in DOMAIN
  filter(!grepl("splice", `ANN[0].EFFECT`, ignore.case = TRUE)) %>%  # Remove rows with "splice" in ANN[0].EFFECT
  filter(!grepl("CD247", `ANN[0].GENE`, ignore.case = TRUE)) %>%  # Remove rows with "splice" in ANN[0].EFFECT
  filter(!grepl("FCER1G", `ANN[0].GENE`, ignore.case = TRUE))  # Remove rows with "splice" in ANN[0].EFFECT
cat(paste0("FILTERED: ", paste0(dim(filteredVariants),collapse="x"), "\n"))

# Convert AF to numeric if necessary
filteredVariants$AF <- as.numeric(filteredVariants$AF)

# Identify the two rows with the highest AF values
top1Variants <- filteredVariants %>%
  arrange(desc(AF)) %>%
  slice_head(n = 1)  # Select top 2 highest AF rows

top3Variants <- filteredVariants %>%
  arrange(desc(AF)) %>%
  slice_head(n = 3)  # Select top 2 highest AF rows

otherVariants <- filteredVariants[!filteredVariants$ID %in% top3Variants$ID, ]

# Create violin plot of AF for all remaining variants
# plot1Lin <- ggplot(filteredVariants, aes(x = "All Variants", y = AF)) +
#   geom_violin(fill = "blue", alpha = 0.5) +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   labs(title = "Allele Frequency Distribution (All Variants)",
#        x = "Variants",
#        y = "Allele Frequency (%)") +
#   theme_minimal()
# ggsave(plotFilenameAll, plot = plot1Lin, width = 6, height = 4, dpi = 300)
# cat(paste0("PLOT: ", plotFilenameAll, "\n"))

# create log-scale version
# plot1Log <- ggplot(filteredVariants, aes(x = "All Variants", y = AF)) +
#   geom_violin(fill = "blue", alpha = 0.5) +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   scale_y_log10(labels = label_number()) +  # Apply log scale and format labels as decimal
#   labs(title = "Allele Frequency Distribution (All Variants) - Log Scale",
#        x = "Variants",
#        y = "Allele Frequency (%) (log10)") +
#   theme_minimal()
# ggsave(plotFilenameAllLog, plot = plot1Log, width = 6, height = 4, dpi = 300)
# cat(paste0("PLOT: ", plotFilenameAllLog, "\n"))


#
# Create violin plot of AF grouped by ANN[0].GENE
#
plot2Lin <- ggplot(filteredVariants, aes(x = `ANN[0].GENE`, y = AF)) +
  geom_violin(fill = "lightgray",color="darkgray", alpha = 1) +
  geom_jitter(data= otherVariants,width = 0.2, alpha = 1, size = 0.3, color="black", fill="black") +
  geom_point(data = top3Variants, aes(x = `ANN[0].GENE`, y = AF), color = "black", size = 1) +  # Centered top 2 points
  geom_point(data = top1Variants, aes(x = `ANN[0].GENE`, y = AF), color = "black", size = 3) +  # Centered top 2 points
  labs(title = figureTitle,
       cex= 2.3, x = "Gene",
       y = "Allele Frequency (%)") +
  theme_minimal() +
  theme(
    # plot title - big and centered
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),   
    # Axis titles
    # Rotate gene labels for readability, size
    axis.title.y = element_text(size = 16), 
    axis.title.x = element_text(size = 16), 
    # Axis values
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size = 14)   # Increase y-axis tick label size
  )
plot2LinRs <- plot2Lin  +
  geom_text_repel(data = top1Variants, aes(label = ID), size = 4, color = "black", nudge_x=1, nudge_y=-0.01)   # Add labels for top 2 AF values

ggsave(plotFilenamePerGene, plot = plot2LinRs, width = 8, height = 5, dpi = 300)
cat(paste0("PLOT: ", plotFilenamePerGene, "\n"))
ggsave(plotFilenamePerGenePDF, plot = plot2LinRs, width = 8, height = 5, dpi = 300)
cat(paste0("PLOT: ", plotFilenamePerGenePDF, "\n"))
print(plot2Lin)

#
# create a version w/o titles
#
plot2LinNoTitles = plot2LinRs + theme(
  plot.title   = element_blank(),
  axis.title.y = element_blank(),
  axis.title.x = element_blank()
)
ggsave(plotFilenamePerGeneNoTitles, plot = plot2LinNoTitles, width = 8, height = 5, dpi = 300)
cat(paste0("PLOT: ", plotFilenamePerGeneNoTitles, "\n"))
ggsave(plotFilenamePerGeneNoTitlesPDF, plot = plot2LinNoTitles, width = 8, height = 5, dpi = 300)
cat(paste0("PLOT: ", plotFilenamePerGeneNoTitlesPDF, "\n"))

# 
# version with no text
#
plot2LinNoText = plot2Lin + theme(
  plot.title   = element_blank(),
  axis.title.y = element_blank(),
  axis.title.x = element_blank(),
  axis.text.y = element_blank(),
  axis.text.x = element_blank()
)
ggsave(plotFilenamePerGeneNoText, plot = plot2LinNoText, width = 8, height = 5, dpi = 300)
cat(paste0("PLOT: ", plotFilenamePerGeneNoText, "\n"))
ggsave(plotFilenamePerGeneNoTextPDF, plot = plot2LinNoText, width = 8, height = 5, dpi = 300)
cat(paste0("PLOT: ", plotFilenamePerGeneNoTextPDF, "\n"))

print(plot2LinNoText)

# create log-scaled version
# plot2Log = plot2 <- ggplot(filteredVariants, aes(x = `ANN[0].GENE`, y = AF)) +
#   geom_violin(fill = "lightgray", alpha = 0.3) +
#   geom_jitter(width = 0.2, alpha = 0.5, size=3) +
#   scale_y_log10(labels = label_number()) +  # Apply log scale and format labels as decimal
#   labs(title = "Allele Frequency Distribution by Gene - Log Scale",
#        x = "Gene",
#        y = "Allele Frequency (%) (log10)") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate gene labels for readability
# ggsave(plotFilenamePerGeneLog, plot = plot2Log, width = 8, height = 5, dpi = 300)
# cat(paste0("PLOT: ", plotFilenamePerGeneLog, "\n"))
# plot2Log
#
# Histogram of variant frequencies (linear scale)
#
# plot3Lin <- ggplot(filteredVariants, aes(x = AF)) +
#   geom_histogram(bins=100, fill = "purple", alpha = 0.7, color = "black") +
#   labs(title = "Histogram of Allele Frequencies (Linear Scale)",
#        x = "Allele Frequency (%)",
#        y = "Count") +
#   theme_minimal()

# Save third plot (histogram with linear scale)
# ggsave(plotFilenameHistLinear, plot = plot3Lin, width = 6, height = 4, dpi = 300)

# Histogram of variant frequencies (log scale)
# plot3Log <- ggplot(filteredVariants, aes(x = AF)) +
#   geom_histogram(bins=100, fill = "red", alpha = 0.7, color = "black") +
#   scale_x_log10(labels = label_number()) +  # Log scale for x-axis
#   labs(title = "Histogram of Allele Frequencies (Log Scale)",
#        x = "Allele Frequency (%) (log10)",
#        y = "Count") +
#   theme_minimal()

# Save fourth plot (histogram with log scale)
# ggsave(plotFilenameHistLog, plot = plot3Log, width = 6, height = 4, dpi = 300)

# Print plots
# print(plot3Lin)
# print(plot3Log)
# print(plot1Lin)
# print(plot2Lin)
# print(plot1Log)
# print(plot2Log)
# 
# 
