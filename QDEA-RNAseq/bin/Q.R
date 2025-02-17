#!/usr/bin/env Rscript

# Cuantificación y analisis de expresión diferencial utilizando Tximport & DESeq2
# Inmegen 2024

# librerias a utilizar
invisible( lapply(c(
"optparse",
"tximportData",
"tximport",
"readr",
"dplyr",
"BUSpaRse",
"SummarizedExperiment"), library, character.only=T))

# lista de opciones opciones(w,i,q,g,b,d)
option_list <- list(
    make_option(c("-w","--working_dir"), type="character", default=NULL            ,metavar="character" ,help="path to DEA working directory (-w)"                            ),
    make_option(c("-i","--sample_info"), type="character", default=NULL            ,metavar="character" ,help="csv file that contains the samples info (-i)"                  ),
    make_option(c("-q","--dir_quants" ), type="character", default="kallisto_quant",metavar="character" ,help="directory name with kallisto quants (-q)"                      ),
    make_option(c("-g","--gtf_file"   ), type="character", default=NULL            ,metavar="character" ,help="path to gtf file (-g)"                                         ),
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
files <- file.path(dir,opt$dir_quants, samples$Sample, "quant.sf")
names(files) <- paste0(samples$Sample)

# tx2gene = data frame con al menos 2 columnas; 1) transcript ID and 2) gene ID
# conservar el orden de las columnas es importante  
geneIds <- tr2g_gtf(opt$gtf_file, get_transcriptome = F)
tx2gene <- geneIds[,c("transcript","gene")]

# Tabla con los Genes IDs y los nombres comunes
geneNames <- distinct(geneIds[,c("gene","gene_name")])

geneNamesT <- distinct(geneIds[,c("transcript","gene","gene_name")])

# Reemplazar NA en gene_name con el valor de gene
geneNames <- geneNames %>% mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", gene, gene_name))

#### Matrices de cuentas a nivel de gen Salmon
txi.salmon <- tximport(files, type = "salmon",tx2gene = tx2gene)

#### Matrices de cuentas a nivel de transcrito Salmon
txi.salmonT <- tximport(files, type = "salmon", txOut = TRUE)

# Importar los datos escalados usando la longitud promedio del transcripto y el tamaño de la biblioteca
txi.salmontpm  <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

# Importar los datos escalados usando la longitud promedio del transcripto y el tamaño de la biblioteca
txi.salmontpmT  <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "lengthScaledTPM")

# Exportar la matriz de cuentas 
m_counts <- as.data.frame(txi.salmon$counts)
m_counts$gene <- row.names(m_counts)
m_counts1 <- merge(geneNames, m_counts, by="gene", all.y = TRUE)
write.table(m_counts1,file=opt$countmat, sep="\t", row.names = FALSE, quote=FALSE)

# Exportar la matriz de cuentas TPM 
m_countstpm <- as.data.frame(txi.salmontpm$counts)
m_countstpm$gene <- row.names(m_countstpm)
m_countstpm1 <- merge(geneNames, m_countstpm, by="gene", all.y = TRUE)
write.table(m_countstpm1,file=opt$countpm, sep="\t", row.names = FALSE, quote=FALSE)

# Exportar la matriz de cuentas
m_countsT <- as.data.frame(txi.salmonT$counts)
m_countsT$transcript <- row.names(m_countsT)
m_counts1T <- merge(geneNamesT, m_countsT, by="transcript", all.y = TRUE)
write.table(m_counts1T,"mcuentastx.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# Exportar la matriz de cuentas TPM
m_countstpmT <- as.data.frame(txi.salmontpmT$counts)
m_countstpmT$transcript <- row.names(m_countstpmT)
m_countstpm1T <- merge(geneNamesT, m_countstpmT, by="transcript", all.y = TRUE)
write.table(m_countstpm1T,file="mcuentastx_tpm.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# R session info
RLogFile <- "R_sessionInfo.log"

sink(RLogFile)
b2  <- sessionInfo()
print("Información de la sesión de R")
print(b2)
sink()
