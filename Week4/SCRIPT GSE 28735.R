#INSTALL TOOLS BELOW:
#BiocManager
if (!require("BiocManager", quietly = TRUE))  {
  install.packages("BiocManager") 
}

#Tools from BiocManager
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 

#Install annotation package
BiocManager::install("hugene10sttranscriptcluster.db", ask = FALSE, update = FALSE)

# Install tools with the command below
install.packages(c("pheatmap", "ggplot2", "dplyr"))

#Install umap
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

#CALL LIBRARIES AFTER INSTALLATIONS COMPLETED
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)
library(umap)

#FETCH DATA FROM GEO (Including its matrices and annotations)
gset <- getGEO("GSE28735", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]


# group membership for all samples
gsms <- paste0("01010101010101010101010101010101010101010101010101",
               "0101010101010101010101010101010101010101")
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Tumor","Nontumor"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

# make group names
group_name <- levels(gset$group)
print(group_name)

#put data into groups as per group name
Tumor <- group_name[1]
Non_Tumor <- group_name[2]

#contrast data
contrast_formula <- paste(Tumor, "-", Non_Tumor)
print(paste("analysed contrast:", contrast_formula))

fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#making top table results
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults)

#fetch probe ID
probe_ids <- rownames(topTableResults)

#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#MAKE BOXPLOT
#set colour based on group
group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot of Expression Value per Sample Distribution",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

# MAKE DENSITY PLOT
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Gene Expression Value Distribution",
    x = "Expression Value (log2)",
    y = "Density"
  )


#MAKE UMAP
umap_input <- t(ex)

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Sample Plot Based on Gene Expression",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#MAKE VOLCANO PLOT

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

#classify gene status
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

#Visualise
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot of Tumor vs Nontumor")

#MAKE HEATMAP

#Choose top 50 genes
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

#Pick matrices of the top 50 genes
mat_heatmap <- ex[top50$PROBEID, ]

#use gene symbols
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # if SyMBOL is unavailable, use probe ID
  top50$SYMBOL        # if available, use SYMBOL
)

rownames(mat_heatmap) <- gene_label

#DATA CLEANSING
#removeNA data
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#remove data with the variance walue of 0
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#annotate collums
annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

#Heatmap visualisation
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gene
  annotation_col = annotation_col,
  show_colnames = FALSE,         # gene name deactivated
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)


# SAVE RESULTS
# saving results to CSV
write.csv(topTableResults, "DEG_GSE28735_Result.CSV")

message("file saved

