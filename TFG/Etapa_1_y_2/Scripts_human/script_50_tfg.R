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

##SCRIPT

peaks <- read.table(file="/Users/abskiel/Desktop/peaks_tfg/H.sapiens/GATA1/GSE29194/GSE29194_GATA1_ERYTH_BMP.bed", header = FALSE, sep = '')
cromosomas_limpio <- substring(peaks[,1], 4)
peaks$V1 <- cromosomas_limpio
coordenada_inicio <- peaks$V2
coordenada_final <- peaks$V3

df50 <- data.frame(matrix(ncol=4, nrow=length(cromosomas_limpio)))
nombres_columnas50 <- c("Cromosoma50","Coordenada_inicial50","Coordenada_final50", "peak50")
colnames(df50) <- nombres_columnas50

df50$Cromosoma50 <- peaks$V1
df50$Coordenada_inicial50 <- peaks$V8 - 50
df50$Coordenada_final50 <- peaks$V8 + 50

texto50 <- paste0("/sequence/region/human/", df50$Cromosoma50, ":", df50$Coordenada_inicial50, "..", df50$Coordenada_final50, ":1?")

vector_secuencias50 <- character()

for (i in 1:length(texto50)) {
  server <- "https://rest.ensembl.org"
  r <- GET(paste(server, texto50[i], sep = ""), content_type("text/plain"))
  stop_for_status(r)
  print(content(r))
  vector_secuencias50[i] <- content(r)
}

df50$peak50 <- vector_secuencias50

##GENERAR LOS RESULTADOS EN FORMATO FASTA

encabezados <- paste(">Secuencia_", 1:length(vector_secuencias50),"-", df50$Cromosoma50, "-", df50$Coordenada_inicial50, ":", df50$Coordenada_final50, sep = "")
fasta <- paste(encabezados, vector_secuencias50, sep = "\n" )
writeLines(fasta, "secuencias50_GATA1_BMP.fasta")