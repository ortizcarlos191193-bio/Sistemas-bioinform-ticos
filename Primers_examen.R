############################################################
# EJERCICIO: DISEÑO SIMPLE DE PRIMERS DESDE UNA SECUENCIA NCBI
############################################################

# 1. Función para leer la secuencia FASTA ------------------

leer_fasta <- function(archivo) {
  lineas <- readLines(archivo, warn = FALSE)
  sec <- lineas[!grepl("^>", lineas)]        # quitar encabezado
  sec <- paste(sec, collapse = "")          # unir todo
  sec <- gsub("\\s+", "", sec)              # limpiar espacios
  sec <- toupper(sec)                       # mayúsculas
  return(sec)
}

# 2. Funciones para evaluar primers ------------------------

calc_gc <- function(primer) {
  bases <- strsplit(primer, "")[[1]]
  n_gc <- sum(bases == "G" | bases == "C")
  round(100 * n_gc / length(bases), 2)
}

calc_tm <- function(primer) {
  bases <- strsplit(primer, "")[[1]]
  A <- sum(bases == "A"); T <- sum(bases == "T")
  G <- sum(bases == "G"); C <- sum(bases == "C")
  2 * (A + T) + 4 * (G + C)          # Regla de Wallace
}

reverse_complement <- function(seq) {
  bases <- strsplit(seq, "")[[1]]
  comp <- bases
  comp[bases == "A"] <- "T"
  comp[bases == "T"] <- "A"
  comp[bases == "G"] <- "C"
  comp[bases == "C"] <- "G"
  rev <- rev(comp)
  paste(rev, collapse = "")
}

# 3. Leer secuencia descargada de NCBI ---------------------

secuencia <- leer_fasta("Btub1.fasta")   # <-- el alumno descarga este archivo
cat("Longitud total:", nchar(secuencia), "bases\n\n")

cat("Primeras 80 bases de la secuencia:\n")
cat(substr(secuencia, 1, 80), "\n\n")

# 4. Tomar dos regiones de la secuencia --------------------

primer_forward  <- substr(secuencia, 101, 120)   # 20 nt desde posición 101
region_reverse  <- substr(secuencia, 501, 520)   # 20 nt desde posición 501

primer_reverse <- reverse_complement(region_reverse)

# 5. Evaluación de primers --------------------------------

cat("Primer Forward:", primer_forward, "\n")
cat("GC% Forward:", calc_gc(primer_forward), "\n")
cat("Tm Forward:", calc_tm(primer_forward), "\n\n")

cat("Región Reverse (sentido genómico):", region_reverse, "\n")
cat("Primer Reverse (reverse complement):", primer_reverse, "\n")
cat("GC% Reverse:", calc_gc(primer_reverse), "\n")
cat("Tm Reverse:", calc_tm(primer_reverse), "\n\n")
