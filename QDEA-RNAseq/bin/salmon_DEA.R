#!/usr/bin/env Rscript

# Cuantificación y analisis de expresión diferencial utilizando Tximport & DESeq2
# Inmegen 2024

# librerias a utilizar
invisible( lapply(c(
"optparse",
"readr",
"dplyr",
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
"rtracklayer",
"purrr",
"magrittr",
"SummarizedExperiment"), library, character.only=T))

# Lista de opciones opciones(w,i,q,g,k,c,t,l,p,o,m,v,a,f,b,d)
option_list <- list(
    make_option(c("-i","--meta_data"  ), type="character", default=NULL            ,metavar="character" ,help="tsv file that contains the metadata info (-i)"                 ),
    make_option(c("-g","--gtf_file"   ), type="character", default=NULL            ,metavar="character" ,help="path to gtf file (-g)"                                         ),
    make_option(c("-c","--condition1" ), type="character", default="treated"       ,metavar="character" ,help="condition 1 name (default control, -c)"                        ),
    make_option(c("-a","--condition2" ), type="character", default="control"       ,metavar="character" ,help="condition 2 name (default treated, -a)"                        ),
    make_option(c("-l","--Log2FC_th"  ), type="numeric"  , default=1               ,metavar="numeric"   ,help="Log2Fc theshold (default = 1 , -1)"                            ),
    make_option(c("-p","--p_adj_th"   ), type="numeric"  , default=0.5             ,metavar="numeric"   ,help="padj theshold (default = 0.1, -p)"                             ),
    make_option(c("-o","--outdir_pca" ), type="character", default="pca.pdf"       ,metavar="character" ,help="name of PCA plot (default pca.pdf, -o)"                        ),
    make_option(c("-m","--out_p_hm"   ), type="character", default="heatmap.pdf"   ,metavar="character" ,help="name of heatmap plot (default heatmap.pdf, -m)"                ),
    make_option(c("-v","--outdir_vp"  ), type="character", default="volcano.pdf"   ,metavar="character" ,help="name of volcano plot (default volcano.pdf, -v)"                ),
    make_option(c("-x","--out_res"    ), type="character", default="results.csv"   ,metavar="character" ,help="name of tsv results file (default results.csv, -x)"            ),
    make_option(c("-f","--out_deg"    ), type="character", default="fresults.csv"  ,metavar="character" ,help="name of tsv filtered results file (default results.csv, -f)"   ),
    make_option(c("-d","--nsamples"   ), type="numeric"  , default=3               ,metavar="character" ,help=" (min samples number, -d)"                                     ))

# convertir la lista de opciones a argumentos 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# ruta del directorio de trabajo
dir <- opt$working_dir

#######################################
##### Metadata
#######################################
 
# tsv con la información de las muestras 
samples <- read.table(opt$sample_info,sep="\t",header=T)

# ruta y nombre de los archivos para importar con tximport 
files <- file.path(dir,opt$dir_quants, samples$Sample, "quant.sf")
names(files) <- paste0(samples$Sample)

#####################################
###### Anotacion de genes 
#####################################

gtf_file <- import(opt$gtf_file, format = "gtf")

# tx2gene = data frame con 2 columnas; 1) transcript ID and 2) gene ID
# conservar el orden de las columnas es importante  

tx_entries <- gtf[gtf$type == "transcript"]

tx_names <- unique(data.frame(transcript_id = paste0(tx_entries$transcript_id, ".", tx_entries$transcript_version),
                              gene_id       = paste0(tx_entries$gene_id, ".", tx_entries$gene_version),
                              gene_name     = tx_entries$gene_name,
                              gene_biotype  = tx_entries$gene_biotype, stringsAsFactors = FALSE))

gene_names <- tx_names
gene_names$transcript_id <- NULL

########################################
###### Matriz de cuentas
########################################

txi.salmon <- tximport(files, type = "salmon",tx2gene = tx_names[, c("transcript_id", "gene_id")])

#########################################
########### Analisis de expresión
#########################################

colData <- merge(data.frame(Sample = colnames(txi.salmon$counts)), samples ,sort = FALSE) 
sampleTable <- data.frame(Sample = colData$Sample,condition = factor(colData$condition))
rownames(sampleTable) <- sampleTable$Sample
sampleTable$Sample <- NULL

## Elegir las condiciones de las condiciones a comparar 
filtered_samples <- rownames(sampleTable)[sampleTable$condition %in% c(opt$condition1, opt$condition2)]
head(filtered_samples)

# Filtrar las filas del data frame sampleTable
filtered_mcounts <- txi.salmon$counts[, filtered_samples, drop = FALSE]
filtered_abundance <- txi.salmon$abundance[, filtered_samples, drop = FALSE]
filtered_length <- txi.salmon$length[, filtered_samples, drop = FALSE]

txi.filtered <- list(counts = filtered_mcounts,
                     abundance = filtered_abundance,
                     length = filtered_length)

# Filtrar las filas del data frame sampleTable
filtered_sampleTable <- sampleTable[filtered_samples, , drop = FALSE]
head(filtered_sampleTable)

##### Generar objeto de DESeq2
dds <-  DESeqDataSetFromTximport(txi.filtered, filtered_sampleTable, ~condition)

keep <- rowSums(counts(dds) >= 10) >= opt$nsamples
dds <- dds[keep, ]

dds <- DESeq(dds)

#############################################
####### Resultados de Deseq2
#############################################

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
rld <- rlog(dds, blind = F)
rlog_matrix <- assay(rld)

ppca <- plotPCA(rld, intgroup = "condition")
ggsave("pca_lncRNA.png",ppca,width = 8, height = 6, dpi = 300)

# Info de las condicion de las muestras
sample_names <- colnames(rlog_matrix)
annotation_col <- data.frame(condition = sampleTable$condition)
rownames(annotation_col) <- sample_names

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

# Exportar la matriz de cuentas 
m_counts <- as.data.frame(txi.salmon$counts)
m_counts$gene <- row.names(m_counts)
m_counts1 <- merge(geneNames, m_counts, by="gene", all.y = TRUE)
write.table(m_counts1,file="mcounts_salmon.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# Exportar las tabla con los estadisticos de los genes (hipotesis predeterminadas)
write.table(as.data.frame(table_genes),file = opt$out_res, sep="\t", row.names = FALSE, quote=FALSE)

# Exportar las tabla con los genes diferencialmente expresados
write.table(as.data.frame(DEG),file = opt$out_deg, sep="\t", row.names = FALSE, quote=FALSE)

# Exporta el objeto de resultados de DESeq2
saveRDS(dds, "DESeq2_dds.rds")

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
