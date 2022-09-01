#Clean dataframe
df_clean <- subset(omim_analysis, select = -c(Symbol, HGNCID, name, avg, std, median))

#omim_analysis is the raw dataset with extra columns than wnated

row.names(df_clean) = omim_analysis$Symbol


#Forloop to match two datasets
#batch is the dataset containing samples with corresponding batch numbers

rev_df <- t(df_clean)


batch <- read.csv(file.choose())

for(i in 1:nrow(batch)) {
  row <- batch[i,2]
  print(row)
  B2 <- paste(c("X", batch[i,2]), collapse = "")
  batch[i,2] <- gsub("-" , ".", B2)
  print(batch[i,2])
}


match(data.frame(t(x)), data.frame(t(y)))


#merge two datasets
merged <- merge(batch, rev_df)

my_cols <- c("Batch", "Sample")

M2 <- do.call(paste, c(merged[my_cols], sep = "-"))

new_merged <- merge(M2, rev_df)

summary(new_merged)



#alphabetical order

batch[,order(colnames(batch))]

a=unlist(batch$Sample)

data_new1 <- df_clean[, match(a, colnames(df_clean))] 

#heatmap 

heatmap_test <- heatmap(as.matrix(data_new1), Rowv=NA, Colv=NA, col = cm.colors(256), scale="column")


heatmap_test <- heatmap(as.matrix(data_new1), ColSideColors=colors)


colors=c()


heatmap(as.matrix(data_new1),ColSideColors=colors)

for(i in 1:nrow(batch)) {
  row <- batch[i, ]
  if (row$Batch == 200128){
    colors=c(colors,"#7fc97f")
  }
  else if ( row$Batch == 324337) {
    colors=c(colors,"#beaed4")
  } 
  else if ( row$Batch == 378455) {
    colors=c(colors,"#fdc086")
  } 
  else if ( row$Batch == 445390) {
    colors=c(colors,"#ffff99")
  } 
  else if ( row$Batch == 502514) {
    colors=c(colors,"#386cb0")
  } 
}

heatmap(as.matrix(data_new1),ColSideColors=colors)
heatmap(as.matrix(merged[, -1]))


#HEATMAP RIN by using the same code used for batch number

RIN <- read.csv("/Users/amira/RINH.csv", sep="\t")


RIN[,order(colnames(RIN))]


data_new1 <- df_clean[,order(colnames(df_clean))]


heatmap_test <- heatmap(as.matrix(data_new1), Rowv=NA, Colv=NA, col = cm.colors(256), scale="column")

heatmap(as.matrix(data_new1), ColSideColors=colors)

colors=c()
for(i in 1:nrow(RIN)) {
  row <- RIN[i, ]
  
  if (row$CG.RIN < 6) {
    colors=c(colors,"#e5f5f9")
  }
  else if ( row$CG.RIN > 8) {
    colors=c(colors,"#2ca25f")
  } 
  else{
    colors=c(colors,"#99d8c9")
  } 

  
}




heatmap(as.matrix(data_new1),ColSideColors=colors)
heatmap(as.matrix(merged[, -1]))


#PCA plot

pc <- prcomp(df_clean,
             center = TRUE,
             scale. = TRUE)
attributes(pc)

print(pc)

str(df_clean)

head(df_clean)

df_clean_pca <- prcomp(df_clean[,c(1:104)],
                       center = TRUE,
                       scale. = TRUE)

summary(df_clean_pca)


str(df_clean_pca)

install.packages("ggfortify")

library(ggfortify)

df_pca_plot <- autoplot(df_clean_pca,
                        data = df_clean,x=1,y=2, label = TRUE, label.size = 3, shape = FALSE)

biplot(df_clean_pca)


#PCA with excluded HBA1,HBA2, HBB, HBD genes to get a more precise vizualization

df_PCA2 <- df_clean[-c(1643, 1644, 1645, 1646),]

pc <- prcomp(df_PCA2,
             center = TRUE,
             scale. = TRUE)
attributes(pc)

print(pc)

str(df_PCA2)

head(df_PCA2)

df_PCA2_plot <- prcomp(df_PCA2[,c(1:104)],
                       center = TRUE,
                       scale. = TRUE)

summary(df_PCA2_plot)


str(df_PCA2_plot)

install.packages("ggfortify")

library(ggplot2)

library(ggfortify)

PCA2_plot <- autoplot(df_PCA2_plot,
                      data = df_PCA2,x=1,y=2, label = TRUE, label.size = 3, 
                      shape = FALSE, ylim=c(-0.025, 0.03), xlim=c(0.01, 0.02))

biplot(df_PCA2_plot, expand=10, xlim=c(-0.3, 0.8), ylim=c(-0.45, 0.27))


#the xlim and ylim were adjuste to zoom in and zoom out in the plot


par(mfrow=c(2,2))









