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

df25 <- data.frame(matrix(ncol=4, nrow=length(peaks$V1)))
nombres_columnas25 <- c("Cromosoma25","Coordenada_inicial25","Coordenada_final25","peak25")
colnames(df25) <- nombres_columnas25

df25$Cromosoma25 <- peaks$V1
df25$Coordenada_inicial25 <- peaks$V8 - 25
df25$Coordenada_final25 <- peaks$V8 + 25


texto25 <- paste0("/sequence/region/arabidopsis_thaliana/", df25$Cromosoma25, ":", df25$Coordenada_inicial25, "..", df25$Coordenada_final25, ":1?")

vector_secuencias25 <- character()
for (i in 1:length(texto25)) {
  server <- "https://rest.ensembl.org"
  r <- GET(paste(server, texto25[i], sep = ""), content_type("text/plain"))
  stop_for_status(r)
  print(content(r))
  vector_secuencias25[i] <- content(r)
}

df25$peak25 <- vector_secuencias25

##GENERAR LOS RESULTADOS EN FORMATO FASTA

encabezados <- paste(">Secuencia_", 1:length(vector_secuencias25),"-", df25$Cromosoma25, "-", df25$Coordenada_inicial25, ":", df25$Coordenada_final25, sep = "")
fasta <- paste(encabezados, vector_secuencias25, sep = "\n" )
writeLines(fasta, "secuencias25_ZAT10_JA.fasta")
