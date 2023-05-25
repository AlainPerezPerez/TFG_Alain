##Paquetes
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

##variables 

carpeta_logo <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\Etapa_4\\motifs_logo\\human\\aff4\\starved"

##decimales y notacion cientifica

options(scipen = 999)
options(digits = 22)

##Cambiar dirección de carpeta
dirección_carpeta <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\Etapa_4\\motifs_entrada\\human\\aff4\\starved"
ruta_carpeta <- gsub("\\\\", "/", dirección_carpeta)


##Generar lista de motifs en formato unviersalmotif

motif_files <- list.files(ruta_carpeta, pattern = ".motif", full.names = TRUE)
motif_list <- lapply(motif_files, read_homer)

# Asignar nombres a cada motif en la lista
motif_names <- tools::file_path_sans_ext(basename(motif_files))
for (i in seq_along(motif_list)) {
  motif <- motif_list[[i]]
  motif@name <- motif_names[i]
  motif_list[[i]] <- motif
}

##Cambiar el pseudocount

for (i in seq_along(motif_list)) {
  motif <- motif_list[[i]]
  motif@pseudocount <- 1/3
  motif_list[[i]] <- motif
}

##Cambiar el bkg

fasta_file <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\Background_seq\\seq_genome_background_regions_both_filter_subset.fasta"
background_seqs <- readDNAStringSet(fasta_file)
background_freq <- get_bkg(background_seqs, k=1:2)
frecuencias_uno <- background_freq$probability/sum(background_freq$probability)

## Generar el nuevo bkg con 20 k-mers
new_bkg <- setNames(background_freq$probability, background_freq$klet)

## Asignar el nuevo bkg a cada uno de los motifs
for (i in seq_along(motif_list)) {
  motif <- motif_list[[i]]
  motif@bkg <- new_bkg
  motif_list[[i]] <- motif
}

##Hacer el mapa de calor del conjunto de motifs
comparisons <- compare_motifs(motif_list, use.type = "ICM", method = "PCC", tryRC = FALSE, min.mean.ic = 0)
##heatmap_data <- as.dist(1 - comparisons)
##heatmap_clust <- hclust(comparisons, method = "ward.D2")
heatmap <- ggheatmap(as.matrix(comparisons))

##Exportar los motifs en formato universalmotif

for (i in seq_along(motif_list)) {
  motif <- motif_list[[i]]
  motif_name <- motif@name
  file_path <- paste0("C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\Etapa_4\\motifs_universalmotif\\human\\aff4\\starved", motif_name, ".motif")
  write_motifs(motif, file_path)
}

##Obtención de las secuencias logo
get_motif_matrix <- function(x) {
  return(x@motif)
}
matrix_list <- lapply(motif_list, get_motif_matrix)  
matrix_list <- setNames(lapply(motif_list, function(x) x@motif), sapply(motif_list, function(x) x@name))
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
  filename <- paste0(carpeta_logo, "/", names(matrix_list)[i], ".pdf") # especificar la ruta completa del archivo
  ggsave(filename, logo_list[[i]], width = 5, height = 4, dpi = 300) # guardar el gráfico
}





