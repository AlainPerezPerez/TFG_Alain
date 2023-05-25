##Paquetes
if (!require("Biostrings")) {
  install.packages("Biostrings", dependencies = TRUE)
  library(Biostrings)
}

if (!require("universalmotif")) {
  install.packages("universalmotif", dependencies = TRUE)
  library(universalmotif)
}


##variables

ruta_motif <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_5\\entrada\\human\\aff4\\starved\\starvedmotif1.motif_starved_100.motif"
ruta_secuencias_fasta <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_5\\entrada\\human\\aff4\\secuencias_fasta\\starved\\secuencias100_AFF4_STARVED.fasta"
nombre <- "motif_refined_starved_100"
salida_motif_refinado <- "C:\\Users\\Usuario\\OneDrive\\Escritorio\\TFG\\etapa_5\\resultados\\human\\aff4\\starved\\"
  
##Cambiar el bkg
fasta_file <- "C:/Users/Usuario/OneDrive/Escritorio/TFG/Background_seq/seq_genome_background_regions_both_filter_subset.fasta"
background_seqs <- readDNAStringSet(fasta_file)
background_freq <- get_bkg(background_seqs, k=1:2)
new_bkg <- setNames(background_freq$probability, background_freq$klet)

##Script
set_target <- readDNAStringSet(ruta_secuencias_fasta)

motif_etapa5 <- read_motifs(ruta_motif)
resultado <- scan_sequences(motif_etapa5, RC=T, sequences = set_target, threshold = 0, threshold.type = "logodds", calc.qvals = FALSE)

motif_refinado <- create_motif(DNAStringSet(resultado$match), name = nombre, pseudocount = 1/3, bkg = new_bkg)


##Exportar el motif
ruta_salida_motif_refinado <- paste0(salida_motif_refinado, nombre, ".motif")
write_motifs(motif_refinado, ruta_salida_motif_refinado)

