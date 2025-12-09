############################################################
# EXAMEN – EJERCICIO PCA APLICADO A QSAR
############################################################

# 1. Cargar datos QSAR ------------------------------------

data <- read.csv("qsar_data.csv")

cat("Dimensiones del dataset:", dim(data), "\n")
cat("Columnas disponibles:\n")
print(colnames(data))
cat("\n")

# 2. Seleccionar solo descriptores numéricos --------------

# Eliminamos la columna Activity si está presente
descriptores <- data[, sapply(data, is.numeric)]

cat("Número de descriptores numéricos:", ncol(descriptores), "\n\n")

# 3. Escalamiento de datos (muy importante en PCA) --------

descriptores_scaled <- scale(descriptores)

# 4. Aplicar PCA usando prcomp() --------------------------

pca_res <- prcomp(descriptores_scaled)

# 5. Mostrar varianza explicada ---------------------------

cat("Varianza explicada por los primeros componentes:\n")
print(summary(pca_res)$importance[2, 1:3])   # PC1-PC3

# 6. Mostrar cargas (loadings) ----------------------------

cat("\nCargas (loadings) de los descriptores en PC1 y PC2:\n")
print(pca_res$rotation[, 1:2])

# 7. Coordenadas de las moléculas en PC1 y PC2 ------------

scores <- pca_res$x

cat("\nPrimeras 5 moléculas proyectadas en PC1 y PC2:\n")
print(scores[1:5, 1:2])
############################################################
# 8. Graficar PC1 vs PC2 para interpretar el PCA
############################################################

pc1 <- scores[, 1]
pc2 <- scores[, 2]

plot(
  pc1, pc2,
  pch = 19,
  col = "blue",
  xlab = paste0("PC1 (", round(100 * summary(pca_res)$importance[2,1], 1), "%)"),
  ylab = paste0("PC2 (", round(100 * summary(pca_res)$importance[2,2], 1), "%)"),
  main = "PCA de QSAR: Moléculas proyectadas en PC1 vs PC2"
)

# Opcional: mostrar etiquetas con nombre de moléculas
text(pc1, pc2, labels = data$Mol, pos = 4, cex = 0.7)


############################################################
