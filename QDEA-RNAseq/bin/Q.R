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
files <- file.path(dir,opt$dir_quants, samples$Sample_name, "abundance.h5")
names(files) <- paste0(samples$Sample_id)

# tx2gene = data frame con al menos 2 columnas; 1) transcript ID and 2) gene ID
# conservar el orden de las columnas es importante  
geneIds <- tr2g_gtf(opt$gtf_file, get_transcriptome = F)
tx2gene <- geneIds[,c("transcript","gene")]

# Importar datos con tximport
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene)

# Importar los datos escalados usando la longitud promedio del transcripto y el tamaño de la biblioteca
txi.countpm  <- tximport(files, type = "kallisto", countsFromAbundance = "lengthScaledTPM",tx2gene = tx2gene)

# Nombres comunes de los genes_ids
geneNames <- distinct(geneIds[,c("gene","gene_name")])

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
b2  <- sessionInfo()
print("Información de la sesión de R")
print(b2)
sink()
