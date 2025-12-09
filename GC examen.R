############################################################
# EJERCICIO DE GC%
# Leer un archivo FASTA descargado del NCBI
# Calcular longitud, frecuencias y GC%
############################################################

# 1. Función simple para leer el FASTA ---------------------
leer_fasta <- function(archivo) {
  lineas <- readLines(archivo, warn = FALSE)
  
  # eliminar encabezados (líneas que comienzan con ">")
  sec <- lineas[!grepl("^>", lineas)]
  
  # unir líneas en una sola secuencia
  sec <- paste(sec, collapse = "")
  
  # limpiar espacios y pasar a mayúsculas
  sec <- gsub("\\s+", "", sec)
  sec <- toupper(sec)
  
  return(sec)
}

# 2. Leer la secuencia descargada --------------------------

secuencia <- leer_fasta("GFP.fasta")

cat("Secuencia cargada:\n")
print(substr(secuencia, 1, 60))  # imprimir solo los primeros 60 caracteres
cat("...\n\n")

# 3. Calcular longitud -------------------------------------

longitud <- nchar(secuencia)
cat("Longitud total:", longitud, "bases\n")

# 4. Calcular frecuencias de bases -------------------------

bases <- strsplit(secuencia, "")[[1]]

frecuencias <- table(bases)
cat("\nFrecuencias de bases:\n")
print(frecuencias)

# 5. Calcular porcentaje de GC ------------------------------

G <- sum(bases == "G")
C <- sum(bases == "C")

GC_porcentaje <- 100 * (G + C) / longitud

cat("\nContenido GC:", round(GC_porcentaje, 2), "%\n")

# 6. Visualización rápida de las bases ---------------------

barplot(
  frecuencias,
  main = "Frecuencia de nucleótidos",
  ylab = "Conteo",
  col = c("skyblue", "tomato", "gold", "forestgreen")
)
