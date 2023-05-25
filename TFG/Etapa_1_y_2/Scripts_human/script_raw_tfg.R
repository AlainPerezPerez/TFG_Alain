##PAQUETES NECESARIOS EN EL SCRIPT

if (!require("httr")) {
  install.packages("httr", dependencies = TRUE)
  library(httr)
}

if (!require("jsonlite")) {
  install.packages("jsonlite", dependencies = TRUE)
  library(jsonlite)
}

if (!require("xml2")) {
  install.packages("xml2", dependencies = TRUE)
  library(xml2)
}

##SCRIPT PARA GENERAR LAS SECUENCIAS DE ADN

peaks <- read.table(file="/Users/abskiel/Desktop/peaks_tfg/H.sapiens/GATA1/GSE29194/GSE29194_GATA1_ERYTH_BMP.bed", header = FALSE, sep = '')
cromosomas_limpio <- substring(peaks[,1], 4)
peaks$V1 <- cromosomas_limpio
coordenada_inicio <- peaks$V2
coordenada_final <- peaks$V3


df <- data.frame(matrix(ncol = 5, nrow = length(cromosomas_limpio)))
nombres_columnas <- c("Cromosoma", "Coordenada_inicial", "Coordenada_final","Longitud","peak")
colnames(df) <- nombres_columnas

df$Cromosoma <- peaks$V1 
df$Coordenada_inicial <- peaks$V2 + 1
df$Coordenada_final <- peaks$V3
df$Longitud <- df$Coordenada_final - df$Coordenada_inicial
texto <- paste0("/sequence/region/human/", df$Cromosoma, ":", df$Coordenada_inicial, "..", df$Coordenada_final, ":1?")

vector_secuencias <- character()
for (i in 1:length(texto)) {
  server <- "https://rest.ensembl.org"
  r <- GET(paste(server, texto[i], sep = ""), content_type("text/plain")) 
  stop_for_status(r)
  print(content(r))
  vector_secuencias[i] <- content(r)
}

df$peak <- vector_secuencias

##GENERAR LOS RESULTADOS EN FORMATO FASTA

encabezados <- paste(">Secuencia_", 1:length(vector_secuencias),"-", df$Cromosoma, "-", df$Coordenada_inicial, ":", df$Coordenada_final, sep = "")
fasta <- paste(encabezados, vector_secuencias, sep = "\n" )
writeLines(fasta, "secuencias_GATA1_BMP.fasta")






