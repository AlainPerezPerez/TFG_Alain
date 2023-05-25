##PAQUETES

if (!require("Biostrings")) {
  install.packages("Biostrings", dependencies = TRUE)
  library(Biostrings)
}

if (!require("universalmotif")) {
  install.packages("universalmotif", dependencies = TRUE)
  library(universalmotif)
}

if (!require("ggseqlogo")) {
  install.packages("ggseqlogo", dependencies = TRUE)
  library(ggseqlogo)
}

if (!require("heatmaply")) {
  install.packages("heatmaply", dependencies = TRUE)
  library(heatmaply)
}

#Variables
nombre_motif <- "motif_starved_caracteristico"

#RUTAS

ruta_motifs_entrada <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_5\\resultados\\human\\aff4\\starved\\motif"
ruta_logos_salida <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_5\\resultados\\human\\aff4\\starved\\logo"
ruta_Motifs_salida <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_6\\human\\aff4\\starved\\motif_caracteristico"
ruta_logo_caracteristico <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_6\\human\\aff4\\starved\\logo_caracteristico"

##Generar lista de motifs refinados

archivos_motifs <- list.files(ruta_motifs_entrada, pattern = ".motif", full.names = TRUE)
lista_motifs_refinados <- lapply(archivos_motifs, read_motifs)

##Cambiar el bkg

fasta_file <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\Background_seq\\seq_genome_background_regions_both_filter_subset.fasta"
background_seqs <- readDNAStringSet(fasta_file)
background_freq <- get_bkg(background_seqs, k=1:2)
frecuencias_uno <- background_freq$probability/sum(background_freq$probability)

## Generar el nuevo bkg con 20 k-mers
new_bkg <- setNames(background_freq$probability, background_freq$klet)

## Asignar el nuevo bkg a cada uno de los motifs
for (i in seq_along(lista_motifs_refinados)) {
  motif <- lista_motifs_refinados[[i]]
  motif@bkg <- new_bkg
  lista_motifs_refinados[[i]] <- motif
}

##Generar secuencias logo de los motifs refinados

get_motif_matrix <- function(x) {
  return(x@motif)
}
matrix_list <- lapply(lista_motifs_refinados, get_motif_matrix)  
matrix_list <- setNames(lapply(lista_motifs_refinados, function(x) x@motif), sapply(lista_motifs_refinados, function(x) x@name))
get_motif_logo <- function(matrix_list) {
  logo_list <- list()
  for (i in 1:length(matrix_list)) {
    motif_name <- names(matrix_list)[i]
    motif_matrix <- matrix_list[[i]]
    motif_logo <- ggseqlogo(motif_matrix, alphabet = "DNA", 
                            stack.width = "fixed", padding.width = 0.1) +
      ggtitle(motif_name) + theme(plot.title = element_text(hjust = 0.5))
    logo_list[[i]] <- motif_logo
    
  }
  return(logo_list)
}

logo_list <- get_motif_logo(matrix_list)
for (i in seq_along(logo_list)) {
  filename <- paste0(ruta_logos_salida, "/", names(matrix_list)[i], ".pdf") # especificar la ruta completa del archivo
  ggsave(filename, logo_list[[i]], width = 5, height = 4, dpi = 300) # guardar el gráfico
}

##Generar el motif característico

merged_motif <- merge_motifs(lista_motifs_refinados, method = "PCC", use.type = "PPM", tryRC = F)
trimmed_motif <- trim_motifs(merged_motif, min.ic = 0.25)
trimmed_motif@name <- nombre_motif
ruta_motif <- file.path(ruta_Motifs_salida, paste0(nombre_motif, ".motif"))
write_motifs(trimmed_motif, file = ruta_motif)

##Generar secuencia logo del motif característico

logo_motif_trimmed <- view_motifs(trimmed_motif)
filename_2 <- paste0(ruta_logo_caracteristico, "/", trimmed_motif@name, ".pdf")
ggsave(filename_2, plot = logo_motif_trimmed + labs(title = trimmed_motif@name))

