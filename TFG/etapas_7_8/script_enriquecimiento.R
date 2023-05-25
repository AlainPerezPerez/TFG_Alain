##Paquetes
if (!require("Biostrings")) {
  install.packages("Biostrings", dependencies = TRUE)
  library(Biostrings)
}

if (!require("universalmotif")) {
  install.packages("universalmotif", dependencies = TRUE)
  library(universalmotif)
}

if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if (!require("seqinr")) {
  install.packages("seqinr", dependencies = TRUE)
  library(seqinr)
}

if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}

##Variables
ruta_motif_caracteristico <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_6\\human\\aff4\\starved\\motif_caracteristico\\motif_starved_caracteristico.motif"
ruta_salida <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_7\\human\\aff4\\starved"
ruta_fasta_raw <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_5\\entrada\\human\\aff4\\secuencias_fasta\\starved\\secuencias_AFF4_STARVED.fasta"
fasta_filtrado <- "secuencias_starved_filtradas.fasta"
ruta_bed_entrada <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_7\\human\\aff4\\starved\\secuencias_bed_entrada\\GSE30267.AFF4.HCT-116_STARVED.bed"
ruta_bed_salida <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_7\\human\\aff4\\starved\\secuencias_bed_salida"
nombre_bed_salida <- "peaks_filtrados_aff4_starved"
ruta_fasta_salida <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_7\\human\\aff4\\starved\\secuencias_fasta_filtradas"
condicion <- "starved"
tf <- "aff4"
##Script

set_target <- readDNAStringSet(ruta_fasta_raw)

motif_enriquecimiento <- read_motifs(ruta_motif_caracteristico)
resultado <- scan_sequences(motif_enriquecimiento, RC=F, sequences = set_target, threshold = 0.25, threshold.type = "logodds", calc.qvals = FALSE)


##Cribado en las secuencias raw

total_secuencias <- as.character(set_target)

separado <- strsplit(total_secuencias, ":", fixed = TRUE)
nombres <- names(separado)
df <- do.call(data.frame, separado)
df <- t(df)
df <- as.data.frame(df)
secuencias <- df$V1
motif <- resultado$match
inicio_motif <- resultado$start
final_motif <- resultado$stop
tabla_secuencias <- data.frame(nombres, secuencias)


##Exportar el archivo .bed filtrado

tabla_bed <- read.table(ruta_bed_entrada)
tabla_bed$V2 <- tabla_bed$V2+1


resultado$inicio <- as.numeric(sub(".*-(\\d+).*", "\\1", resultado$sequence))
coordenadas_seleccionadas <- resultado$inicio
tabla_bed_filtrada <- tabla_bed[tabla_bed[, 2] %in% coordenadas_seleccionadas, ]
tabla_bed_filtrada$V2 <- tabla_bed_filtrada$V2 - 1 
ruta_archivo <- file.path(ruta_bed_salida, nombre_bed_salida)
write.table(tabla_bed_filtrada, ruta_archivo, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##DATAFRAME DE LAS COORDENADAS DE LOS MOTIFS

expresiones <- resultado$sequence
extraer_numero <- function(expresion) {
  parte <- sub(".*-(\\d+):.*", "\\1", expresion)
  return(parte)
}

partes <- sapply(expresiones, extraer_numero)
partes <- str_extract(partes, "\\d+")
partes <- as.numeric(partes)
inicio_motif_bed <- partes + resultado$start
final_motif_bed <- partes + resultado$stop
summit <- match(partes, tabla_bed$V2)
summit_2 <- tabla_bed$V8[summit]

tabla_motifs_bed <- data.frame(resultado$sequence, inicio_motif_bed, final_motif_bed, summit_2)

##Referencia del summit respecto a la coordenada del motif


# Crear un vector vacío para almacenar los resultados
resultado_summit <- vector(length = nrow(tabla_motifs_bed))

# Realizar las comparaciones y asignar los resultados al vector
for (i in 1:nrow(tabla_motifs_bed)) {
  if (tabla_motifs_bed$summit_2[i] >= tabla_motifs_bed$inicio_motif_bed[i] &&
      tabla_motifs_bed$summit_2[i] <= tabla_motifs_bed$final_motif_bed[i]) {
    resultado_summit[i] <- "yes"
  } else if (tabla_motifs_bed$summit_2[i] < tabla_motifs_bed$inicio_motif_bed[i]) {
    resultado_summit[i] <- tabla_motifs_bed$summit_2[i] - tabla_motifs_bed$inicio_motif_bed[i]
  } else {
    resultado_summit[i] <- tabla_motifs_bed$summit_2[i] - tabla_motifs_bed$final_motif_bed[i]
  }
}

# Calcular los porcentajes de "yes"
porcentaje_yes <- sum(resultado_summit == "yes") / length(resultado_summit) * 100

##Exportar secuencias formato fasta

filtrado <- filter(tabla_secuencias, nombres %in% resultado$sequence)
cabecera <- paste0(">", filtrado$nombres, "_", tf, "_", condicion, "_", "coordenadas_motif_secuencia_fasta", "-", tabla_motifs_bed$inicio_motif_bed, ":", tabla_motifs_bed$final_motif_bed, "-orientación:", resultado$strand, "-summit:", resultado_summit, "_", porcentaje_yes, "%")
secuencia_filtrado <- paste(cabecera, filtrado$secuencias, sep = "\n")
ruta_fasta_filtrado <- file.path(ruta_fasta_salida, fasta_filtrado)
writeLines(secuencia_filtrado, ruta_fasta_filtrado)
