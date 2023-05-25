##Paquetes

if (!require("universalmotif")) {
  install.packages("universalmotif", dependencies = TRUE)
  library(universalmotif)
}

##Script

ruta_entrada <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_6\\arabidopsis\\zat10\\ja\\motif_caracteristico\\motif_ja_caracteristico.motif"
ruta_salida <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_7\\arabidopsis\\zat10\\ja\\motif_para_validar\\motif_ja_validar.motif"
motif <- read_motifs(ruta_entrada)

resultado <- write_meme(motif, file = ruta_salida)