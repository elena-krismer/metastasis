library(preprocessCore)
data <- read.csv("~/Documents/GitHub/metastasis/data/gene_expression/raw/lung_metastasis_raw.csv", sep = ";", dec = ",")
data_mat <- data.matrix(data[,-1])
data_norm <- normalize.quantiles(data_mat, copy = TRUE)
samplenames <- colnames(data)[2:6]
genenames <- data[,1]
colnames(data_norm) <- samplenames
rownames(data_norm) <- genenames

write.csv(data_norm, "~/Documents/GitHub/metastasis/data/gene_expression/lung_metastasis_quantile_normalized.csv", row.names = TRUE
          )