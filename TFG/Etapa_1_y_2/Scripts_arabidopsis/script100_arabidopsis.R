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

peaks <- read.table(file="/Users/abskiel/Desktop/peaks_tfg/A.thaliana/ZAT10/GSE133408/GSE133408_ZAT10_JA.bed", header = FALSE, sep = '')
coordenada_inicio <- peaks$V2
coordenada_final <- peaks$V3

df100 <- data.frame(matrix(ncol=4, nrow=length(peaks$V1)))
nombres_columnas100 <- c("Cromosoma100","Coordenada_inicial100","Coordenada_final100","peak100")
colnames(df100) <- nombres_columnas100

df100$Cromosoma100 <- peaks$V1
df100$Coordenada_inicial100 <- peaks$V8 - 100
df100$Coordenada_final100 <- peaks$V8 + 100

options(scipen=999)
texto100 <- paste0("/sequence/region/arabidopsis_thaliana/", df100$Cromosoma100, ":", df100$Coordenada_inicial100, "..", df100$Coordenada_final100, ":1?")

vector_secuencias100 <- character()

for (i in 1:length(texto100)) {
  server <- "https://rest.ensembl.org"
  r <- GET(paste(server, texto100[i], sep = ""), content_type("text/plain"))
  stop_for_status(r)
  print(content(r))
  vector_secuencias100[i] <- content(r)
}

df100$peak100 <- vector_secuencias100

##GENERAR LOS RESULTADOS EN FORMATO FASTA

encabezados <- paste(">Secuencia_", 1:length(vector_secuencias100),"-", df100$Cromosoma100, "-", df100$Coordenada_inicial100, ":", df100$Coordenada_final100, sep = "")
fasta <- paste(encabezados, vector_secuencias100, sep = "\n" )
writeLines(fasta, "secuencias100_ZAT10_JA.fasta")