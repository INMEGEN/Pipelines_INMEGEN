#!/usr/bin/env Rscript

# Cuantificacion y analisis de expresion diferencial utilizando DESeq2
# Inmegen 2024

# Librerias a utilizar
invisible( lapply(c(
"optparse",
"readr",
"dplyr",
"BUSpaRse",
"DESeq2",
"ggplot2",
"ggrepel",
"PCAtools",
"ComplexHeatmap",
"grid",
"AnnotationDbi",
"enrichR",
"tidyr",
"tibble",
"purrr",
"magrittr",
"SummarizedExperiment"), library, character.only=T))

# Lista de opciones opciones(w,i,q,g,k,c,t,l,p,o,m,v,a,f,b,d)
option_list <- list(
    make_option(c("-i","--meta_data"  ), type="character", default=NULL            ,metavar="character" ,help="tsv file that contains the metadata info (-i)"                 ),
    make_option(c("-g","--gtf_file"   ), type="character", default=NULL            ,metavar="character" ,help="path to gtf file (-g)"                                         ),
    make_option(c("-c","--condition1" ), type="character", default="treated"       ,metavar="character" ,help="condition 1 name (default control, -c)"                        ),
    make_option(c("-t","--condition2" ), type="character", default="control"       ,metavar="character" ,help="condition 2 name (default treated, -t)"                        ),
    make_option(c("-l","--Log2FC_th"  ), type="numeric"  , default=1               ,metavar="numeric"   ,help="Log2Fc theshold (default = 1 , -1)"                            ),
    make_option(c("-p","--p_adj_th"   ), type="numeric"  , default=0.5             ,metavar="numeric"   ,help="padj theshold (default = 0.1, -p)"                             ),
    make_option(c("-o","--outdir_pca" ), type="character", default="pca.pdf"       ,metavar="character" ,help="name of PCA plot (default pca.pdf, -o)"                        ),
    make_option(c("-m","--out_p_hm"   ), type="character", default="heatmap.pdf"   ,metavar="character" ,help="name of heatmap plot (default heatmap.pdf, -m)"                ),
    make_option(c("-v","--outdir_vp"  ), type="character", default="volcano.pdf"   ,metavar="character" ,help="name of volcano plot (default volcano.pdf, -v)"                ),
    make_option(c("-x","--out_res"    ), type="character", default="results.csv"   ,metavar="character" ,help="name of tsv results file (default results.csv, -x)"            ),
    make_option(c("-f","--out_deg"    ), type="character", default="fresults.csv"  ,metavar="character" ,help="name of tsv filtered results file (default results.csv, -f)"   ),
    make_option(c("-d","--nsamples"   ), type="numeric"  , default=3               ,metavar="character" ,help=" (min samples number, -d)"                                     ))

# Convertir la lista de opciones a argumentos 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
 
# Archivo con la informacion de las muestras 
samples <- read.table(opt$meta_data,sep="\t",header=T)

# Obtener los ids de genes y transcritos a partir de un archivo GTF
geneIds <- tr2g_gtf(opt$gtf_file, get_transcriptome = F)

# Eliminar las versiones del identificador de genes y transcritos en geneIds
geneIds$gene <- gsub("\\..*","", geneIds$gene)
geneIds$transcript <- gsub("\\..*","", geneIds$transcript)

# Tabla con los IDs de los Genes y los nombres comunes
geneNames <- distinct(geneIds[,c("gene","gene_name")])

# Reemplazar NA en gene_name con el valor de gene
geneNames <- geneNames %>% mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", gene, gene_name))

# Matriz de cuentas de featureCounts
fmcuentas <- read.table("fmcounts.tsv",sep="\t",header=T,stringsAsFactors = FALSE)
rownames(fmcuentas) <- fmcuentas$Geneid
fmcuentas$Geneid <- NULL
mcounts <- as.matrix(fmcuentas)

# Condiciones de las muestras
colData <- merge(data.frame(Sample = colnames(mcounts)), samples ,sort = FALSE) 
sampleTable <- data.frame(Sample = colData$Sample,condition = factor(colData$condition))
rownames(sampleTable) <- sampleTable$Sample
sampleTable$Sample <- NULL


#########################################################################################################
##### Correlación entre muestras
##############################################################################
# Objeto de DESeq2 para correlación entre muestras

dds_all <- DESeqDataSetFromMatrix(countData = mcounts, colData = sampleTable, design = ~condition)

keep_all <- rowSums(counts(dds_all) >= 10) >= opt$nsamples
dds_all <- dds_all[keep_all, ]

dds_all <- DESeq(dds_all)

##############################################################
##########  Expresion diferencial
#############################################################

## Elegir las condiciones de las condiciones a comparar 
filtered_samples <- rownames(sampleTable)[sampleTable$condition %in% c(opt$condition1, opt$condition2)]
head(filtered_samples)

# Seleccionar las columnas correspondientes de la matriz de cuentas
filtered_mcounts <- mcounts[, filtered_samples, drop = FALSE]

# Filtrar las filas del data frame sampleTable
filtered_sampleTable <- sampleTable[filtered_samples, , drop = FALSE]
head(filtered_sampleTable)

# Objeto de DESeq2
dds <- DESeqDataSetFromMatrix(countData = filtered_mcounts, colData = filtered_sampleTable, design = ~condition)

keep <- rowSums(counts(dds) >= 10) >= opt$nsamples
dds <- dds[keep, ]

# Funcion que normaliza los datos y realiza el analisis de expresion diferencial para las muestras elegidas.
dds <- DESeq(dds)

# Obtener los resultados del analisis, nota: Es importante el orden de la comparacion.
res <- results(dds,contrast=c("condition",opt$condition1,opt$condition2))
resOrdered <- res[order(res$pvalue),]
resOrd <- data.frame(gene = rownames(resOrdered), resOrdered)
df_genes <- merge(geneNames,resOrd, by="gene",all.y = TRUE)

# Filtrar los datos con |log2FC| > 1 y un FDR < 0.05
res_subset <- subset(res, abs(res$log2FoldChange) >= opt$Log2FC_th)
resDEG <- subset(res_subset, res_subset$padj < opt$p_adj_th)
resOrdered_f <- resDEG[order(resDEG$pvalue),]
resOrd_f <- data.frame(gene = rownames(resOrdered_f), resOrdered_f)
DEG <- merge(geneNames,resOrd_f, by="gene", all.y = TRUE)

# Transformacion logaritmica de las muestras.
rld <- rlog(dds_all, blind = F)
rlog_matrix <- assay(rld)

# Graficar los datos para el PCA.
ppca <- plotPCA(rld, intgroup = "condition")
ggsave(opt$outdir_pca,ppca,width = 8, height = 6, dpi = 300)

# Info de las condicion de las muestras
sample_names <- colnames(rlog_matrix)
annotation_col <- data.frame(condition = sampleTable$condition)
rownames(annotation_col) <- sample_names

# Clustering de las muestras por correlacion de spearman, las muestras salen ordenadas alfabeticamente 
ordered_sample_names <- rownames(annotation_col)[order(annotation_col$condition)]
rlog_matrix_1 <- rlog_matrix[, ordered_sample_names]

spearman_dist <- cor(rlog_matrix_1, method = "spearman")

# Convertir a NA las entradas superiores a la diagonal principal
spearman_mat <- as.matrix(spearman_dist)
spearman_mat[upper.tri(spearman_mat)] <- NA
clus_annotation <- annotation_col[ordered_sample_names, , drop = FALSE]

png(opt$out_p_hm,width = 2400, height = 1800, res = 300)
set.seed(1)
pheatmap(
  spearman_mat,
  name = "Spearman",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),
  fontsize = 4, fontsize_row = 7, fontsize_col = 7, fontsize_number = 7,
  annotation_col = clus_annotation,
  labels_col = ordered_sample_names,
  labels_row = ordered_sample_names,
  na_col = "white",
  gaps_col = cumsum(table(annotation_col$condition)))
dev.off()

# Heatmap a pariir de la regularización por logaritmos de la matriz de cuentas (regularized logarithm)
# Obtener los nombres de las muestras seleccionadas
sub_muestras <- rownames(filtered_sampleTable)

# Filtrar la matriz de expresión y las anotaciones para el subconjunto de muestras
rld_f <- rlog(dds, blind = F)
rlog_matrix_f <- assay(rld_f)
zscore_matrix_f <- t(apply(rlog_matrix_f, 1, function(x) (x - mean(x)) / sd(x)))
annotation_col_sub <- annotation_col[sub_muestras, , drop = FALSE]
 
# Crear los heatmaps a partir de Log2 y el z-score, si se cambia la opcion show_rownames = F a T se muestran los nombres comunes de los genes.
png("heatmap_log2.png", width = 2400, height = 1800, res = 300)
set.seed(1)
pheatmap(rlog_matrix_f, 
         name = "log2", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         fontsize = 7,
         angle_col = "45", 
         show_rownames = FALSE,
         annotation_col = annotation_col_sub,
         annotation_names_col = FALSE)
dev.off()

png("heatmap_zscore.png", width = 2400, height = 1800, res = 300)
set.seed(1)
pheatmap(zscore_matrix_f,
         name = "z-score",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         fontsize = 7,
         angle_col = "45",
         show_rownames = FALSE,
         annotation_col = annotation_col_sub,
         annotation_names_col = FALSE)
dev.off()

# Hacer grafica de volcan
table_genes <- df_genes
  
df_genes$Type <- ifelse(df_genes$log2FoldChange > 1 & df_genes$padj < 0.05, "Upregulated",
                 ifelse(df_genes$log2FoldChange < -1 & df_genes$padj < 0.05, "Downregulated", "Not Significant"))

volcano_plot <- ggplot(df_genes, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  labs(title = "Volcano Plot",
       x = "Log2FoldChange (Log2FC)",
       y = "-log10(p-adj)") +  
  theme_minimal(base_size = 15) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")

# Esportar la grafica de volcan                 
ggsave(opt$outdir_vp,volcano_plot,width = 8, height = 6, dpi = 300)

# Exportar las tabla con los estadisticos de los genes (hipotesis predeterminadas)
write.table(as.data.frame(table_genes),file = opt$out_res, sep="\t", row.names = FALSE, quote=FALSE)

# Exportar las tabla con los genes diferencialmente expresados
write.table(as.data.frame(DEG),file = opt$out_deg, sep="\t", row.names = FALSE, quote=FALSE)

# Exporta el objeto de resultados de DESeq2
saveRDS(dds, "DESeq2_dds.rds")

### Enriquecimiento de genes 
#setEnrichrSite("Enrichr")
#websiteLive <- TRUE
#dbsAll <- listEnrichrDbs()
#if (is.null(dbsAll)) websiteLive <- FALSE
#if (websiteLive) head(dbsAll)
#Edbs <- c("GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023","Reactome_2022","WikiPathway_2023_Human")

#if (websiteLive) {
#    enriched <- enrichr(DEG$gene, Edbs)
#    printEnrich(enriched)
#}

#enrinch1 <- if (websiteLive) plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
#enrinch2 <- if (websiteLive) plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
#enrinch3 <- if (websiteLive) plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
#enrinch4 <- if (websiteLive) plotEnrich(enriched[[4]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")
#enrinch5 <- if (websiteLive) plotEnrich(enriched[[5]], showTerms = 15, numChar = 40, y = "Count", orderBy = "P.value")

#ggsave("enrichr_GO_Biological_Process_2023_amp.pdf",enrinch1)
#ggsave("enrichr_GO_Cellular_Component_2023_amp.pdf",enrinch2)
#ggsave("enrichr_GO_Molecular_Function_2023_amp.pdf",enrinch3)
#ggsave("enrichr_Reactome_2022.pdf",enrinch4)
#ggsave("enrichr_WikiPathways_2021_Human_amp.pdf",enrinch5)

# R session info
RLogFile <- "R_sessionInfo.log"

sink(RLogFile)
b1  <- res@elementMetadata$description
b2  <- sessionInfo()
print("Tabla utilizada para generar el objeto dds de DESeq2")
print(filtered_sampleTable)
print("Descripcion de los metadatos del objeto results de DESeq2")
print(b1)
print("Resumen del objeto results de DESeq2")
summary(res)
print("Informacion de la sesion de R")
print(b2)
sink()
