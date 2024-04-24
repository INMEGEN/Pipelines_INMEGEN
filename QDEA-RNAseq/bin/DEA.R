#!/usr/bin/env Rscript

# Cuantificación y analisis de expresión diferencial utilizando Tximport & DESeq2
# Visualización de resultados utilizando GeneTonic
# INMEGEN 2022

# librerias a utilizar
invisible( lapply(c(
"optparse",
"tximportData",
"tximport",
"readr",
"dplyr",
"BUSpaRse",
"DESeq2",
"ggplot2",
"ggrepel",
"PCAtools",
"ComplexHeatmap",
"EnhancedVolcano",
"grid",
"AnnotationDbi",
"topGO",
"GeneTonic",
"SummarizedExperiment"), library, character.only=T))

# lista de opciones opciones(w,i,q,g,c,t,l,p,o,m,v,a,f,b,d,s,e,n)
option_list <- list(
    make_option(c("-w","--working_dir"), type="character", default=NULL            ,metavar="character" ,help="path to DEA working directory (-w)"                            ),
    make_option(c("-i","--sample_info"), type="character", default=NULL            ,metavar="character" ,help="csv file that contains the samples info (-i)"                  ),
    make_option(c("-q","--dir_quants" ), type="character", default="kallisto_quant",metavar="character" ,help="directory name with kallisto quants (-q)"                      ),
    make_option(c("-g","--gtf_file"   ), type="character", default=NULL            ,metavar="character" ,help="path to gtf file (-g)"                                         ),
    make_option(c("-c","--condition1" ), type="character", default="treated"       ,metavar="character" ,help="condition 1 name (default control, -c)"                        ),
    make_option(c("-t","--condition2" ), type="character", default="control"       ,metavar="character" ,help="condition 2 name (default treated, -t)"                        ),
    make_option(c("-l","--Log2FC_th"  ), type="numeric"  , default=1               ,metavar="numeric"   ,help="Log2Fc theshold (default = 1 , -1)"                            ),
    make_option(c("-p","--p_adj_th"   ), type="numeric"  , default=0.1             ,metavar="numeric"   ,help="padj theshold (default = 0.1, -p)"                             ),
    make_option(c("-o","--outdir_pca" ), type="character", default="pca.pdf"       ,metavar="character" ,help="name of PCA plot (default pca.pdf, -o)"                        ),
    make_option(c("-m","--out_p_hm"   ), type="character", default="heatmap.pdf"   ,metavar="character" ,help="name of heatmap plot (default heatmap.pdf, -m)"                ),
    make_option(c("-v","--outdir_vp"  ), type="character", default="volcano.pdf"   ,metavar="character" ,help="name of volcano plot (default volcano.pdf, -v)"                ),
    make_option(c("-x","--out_res"    ), type="character", default="results.csv"   ,metavar="character" ,help="name of csv results file (default results.csv, -x)"            ),
    make_option(c("-f","--out_deg"    ), type="character", default="fresults.csv"  ,metavar="character" ,help="name of csv filtered results file (default results.csv, -f)"   ),
    make_option(c("-b","--countmat"   ), type="character", default="countmat.csv"  ,metavar="character" ,help="name of csv with the matrix counts (countmat.csv. -b)"         ),
    make_option(c("-d","--countpm"    ), type="character", default="countpm.csv"   ,metavar="character" ,help="name of csv with the TPM matrix counts (countmaTPM.csv, -d)"   ))

# convertir la lista de opciones a argumentos 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# ruta del directorio de trabajo
dir <- opt$working_dir
 
# csv con la información de las muestras 
samples <- read.table(opt$sample_info,sep="\t",header=T)

# ruta y nombre de los archivos para importar con tximport 
files <- file.path(dir,opt$dir_quants, samples$Sample_name, "abundance.h5")
names(files) <- paste0(samples$Sample_id)

# tx2gene = data frame con al menos 2 columnas; 1) transcript ID and 2) gene ID
# conservar el orden de las columnas es importante  
geneIds <- tr2g_gtf(opt$gtf_file, get_transcriptome = F)
tx2gene <- geneIds[,c("transcript","gene")]

# Tabla con los Genes IDs y los nombres comunes
geneNames <- distinct(geneIds[,c("gene","gene_name")])

# Importar datos con tximport
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene)

# Importar los datos escalados usando la longitud promedio del transcripto y el tamaño de la biblioteca
txi.countpm  <- tximport(files, type = "kallisto", countsFromAbundance = "lengthScaledTPM",tx2gene = tx2gene)


# Las condiciones a comparar, DEPENDEN TOTALMENTE DEL DISEÑO EXPERIMENTAL y se toman del csv con la información de las muestras
sampleTable <- data.frame(condition = factor(samples$condition))
rownames(sampleTable) <- colnames(txi.kallisto$counts)

# Objeto de DESeq2
dds <- DESeqDataSetFromTximport(txi.kallisto, 
	                        sampleTable, 
	                        ~condition)

# Prefiltrado para quitar genes con cuentas bajas (definir bien los thresholds)
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]

# Función que normaliza los datos y realiza el análisis de expresión diferencial
dds <- DESeq(dds)

# Obtener los resultados del análisis de DE, nota: Es importante el orden de la comparación
# Definir simpre las condiciones a tratar porf default es treated vs control
# Resultados asuminendo las hipotesis predeterminadas (|Log2FC| >= 0 y padj = 0.1)
res <- results(dds,contrast=c("condition",opt$condition1,opt$condition2))
resOrdered <- res[order(res$pvalue),]
resOrd <- data.frame(gene = rownames(resOrdered), resOrdered)
statistics_genes <- merge(geneNames, resOrd, by="gene")

# Filtrar los datos con un |log2FC| > 1 y un FDR < 0.05
res_subset <- subset(res, res$log2FoldChange >= 1 | res$log2FoldChange <= -1 )
resDEG <- subset(res_subset, res_subset$padj < 0.05)
resOrdered_f <- resDEG[order(resDEG$pvalue),]
resOrd_f <- data.frame(gene = rownames(resOrdered_f), resOrdered_f)
DEG <- merge(geneNames, resOrd_f, by="gene")

# Transformación recomendada para sets de con menos de 30 muestras
rld <- rlog(dds, blind = F)

# Graficar los datos para el PCA
ppca<-plotPCA(rld, intgroup = "condition")
ggsave(opt$outdir_pca,ppca)

# Heatmap a pariir de la regularización por logaritmos de la matriz de cuentas (regularized logarithm) 
Z <- t(scale(t(assay(rld))))
mat <- merge(geneNames,Z,by.x='gene',by.y=0)
mat<- na.omit(mat)
mat_gn <- as.matrix(mat[,-(1:2)])
rownames(mat_gn) <- mat[,2]

# Crear el heatmap, si se cambia la opción show_rownames = F a T se muestran los nombres comunes de los genes
pdf(opt$out_p_hm)
pheatmap(mat_gn, 
         name = "Row Z-Score", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D2", 
         fontsize = 11, 
         angle_col = "45", 
         show_rownames = F)
dev.off()

# Hacer gráfica de volcán  
volc_plot <- EnhancedVolcano(statistics_genes,
                             lab = statistics_genes$gene_name,
                             x = 'log2FoldChange',
                             y = 'pvalue',
                             title = 'Differential expression genes',
                             pCutoff = opt$p_adj_th,
                             FCcutoff = opt$Log2FC_th,
                             cutoffLineType = 'twodash',
                             cutoffLineWidth = 0.8,
                             pointSize = 3.0,
                             labSize = 6.0,
                             col = c('gray', 'gray', 'gray', 'red3'),
                             legendLabels=c('n.s.','Log2FC','p-value','Log2FC & p-value (DEG)'),
                             legendPosition = 'right',
                             subtitle = bquote(italic("Treated vs Control")))

# Esportar la gráfica de volcán                 
ggsave(opt$outdir_vp,volc_plot)

# Exportar las tabla con los estadisticos de los genes (hipotesis predeterminadas)
write.table(as.data.frame(statistics_genes),file = opt$out_res, sep="\t", row.names = FALSE, quote=FALSE)

# Exportar las tabla con los genes diferencialmente expresados
write.table(as.data.frame(DEG),file = opt$out_deg, sep="\t", row.names = FALSE, quote=FALSE)

# Exportar la matriz de cuentas 
m_counts1 <- as.data.frame(txi.kallisto$counts)
m_counts1$gene <- row.names(m_counts1)
m_counts2 <- merge(geneNames, m_counts1, by="gene")
write.table(m_counts2,file=opt$countmat, sep="\t", row.names = FALSE, quote=FALSE)

# Exportar la matriz de cuentas TPM 
m_countst1 <- as.data.frame(txi.countpm$counts)
m_countst1$gene <- row.names(m_countst1)
m_countst2 <- merge(geneNames, m_countst1, by="gene")
write.table(m_countst2,file=opt$countpm, sep="\t", row.names = FALSE, quote=FALSE)

# R session info
RLogFile <- "R_sessionInfo.log"

sink(RLogFile)
b1  <- res@elementMetadata$description
b2  <- sessionInfo()
print("Tabla utilizada para generar el objeto dds de DESeq2")
print(sampleTable)
print("Descripción de los metadatos del objeto results de DESeq2")
print(b1)
print("Resumen del objeto results de DESeq2")
summary(res)
print("Información de la sesión de R")
print(b2)
sink()
