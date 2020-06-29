set.seed(999)
mat = matrix(runif(8), 8, 1) 
mat = matrix(c(0.1302376,
             0.11563469999999998,
             0.040762999999999994,
             0.07287690000000002,
             0.11317050000000001,
             0.14295190000000002,
             0.12426769999999998,
             0.2600977), 8, 1) 
rownames(mat) = paste0("db", 1:8)
colnames(mat) = paste0("peptides")

chordDiagram(mat)
title("Default")
chordDiagram(mat, grid.col = 1:9, scale = TRUE)
title("scale = TRUE")
