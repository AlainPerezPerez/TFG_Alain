##VARIABLES

ruta_salida <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\Etapa_4\\motifs_entrada\\human\\GATA1\\BIO"
dirección_carpeta <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\resultados_Homer\\Human\\GATA1\\BIO\\motifs\\motifs_BIO_100"
ruta_carpeta <- gsub("\\\\", "/", dirección_carpeta)
condicion <- "bio"
longitud <- "100"


##Lista para guardar los resultados 
motif_folder <- ruta_carpeta
motif_files <- list.files(ruta_carpeta, pattern = ".motif", full.names = TRUE)

motif_data <- list()

resultados <- list() 

for (file_name in motif_files) {
  file_lines <- readLines(file_name)
  
  # Extraer el logaritmo del p_value y las secuencias target de la cabecera
  header_line <- file_lines[1]
  header_info <- strsplit(header_line, "\t")[[1]]
  col4_value <- as.numeric(header_info[4])
  col6_value <- as.numeric(gsub(".*T:(\\d+\\.\\d+).*", "\\1", header_info[6]))
  
  # Guardar la información en la lista
  motif_data[[file_name]] <- list(col4 = col4_value, col6 = col6_value)
  
  # Calcular el valor de la función y agregarlo a la lista de resultados
  resultado <- -col4_value*log(col6_value)
  resultados[[file_name]] <- resultado
}


# Ordenar los resultados en orden descendente y seleccionar los tres primeros
resultados_ordenados <- sort(unlist(resultados), decreasing = TRUE)
top_resultados <- resultados_ordenados[1:5]

##Exportar los motifs
nombres_motifs_seleccionados <- names(top_resultados)

# Copiar los tres motifs seleccionados y cambiarles el nombre
for (i in 1:5) {
  motif_name <- names(resultados)[which(resultados == top_resultados[i])]
  new_name <- paste0(motif_name, "_", condicion, "_", longitud, ".motif")
  file.copy(motif_name, file.path(ruta_salida, basename(new_name)))
}
