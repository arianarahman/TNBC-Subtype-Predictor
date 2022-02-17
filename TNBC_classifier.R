

### Install Package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("CancerSubtypes")
BiocManager::install("ComplexHeatmap")
install.packages("tidyverse", dependencies = TRUE)
install.packages("gplots", dependencies = TRUE)
install.packages("tidytext", dependencies = TRUE)
install.packages("ggfortify", dependencies = TRUE)
install.packages("preprocessCore", dependencies = TRUE)
install.packages("caret", dependencies = TRUE)
install.packages("randomForest", dependencies = TRUE)
install.packages("BBmisc", dependencies = TRUE)
install.packages("pROC", dependencies = TRUE)
install.packages("hrbrthemes", dependencies = TRUE)

## ----global_palette, results = 'asis'------------------------------------
rm(list=ls())
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
start_time <- Sys.time()


### Load Necessary Libraries
library(CancerSubtypes)
library(tidyverse)
library(gplots)
library(tidytext)
library(ggfortify)
library(preprocessCore)
library(caret)
library(randomForest)
library(BBmisc)
library(pROC)
library(hrbrthemes)
library(ComplexHeatmap)


### Define the input and output folders
data_folder <- "D:/TNBC-Project/data/"
output_folder <- "D:/TNBC-Project/output/"

### Load the mRNA expression table into a dataframe
input_file_name = paste0(data_folder,"mrna_processed_GSE76274.txt")
mrna_processed_GSE76274 <- read.table(input_file_name, header = TRUE, sep = "\t")

### Load the sample characteristics table into a dataframe
input_file_name = paste0(data_folder,"GSE76124_series_matrix_mod.txt")
mrna_series_tmat_GSE76274 <- read.table(input_file_name, header = TRUE, sep = "\t")

### Load the matching gene ids for probe ids into a dataframe
input_file_name = paste0(data_folder,"mrna_processed_GSE76274_gene_id_GEO_sel.txt")
mrna_processed_GSE76274_gene_id_GEO_sel <- read.table(input_file_name, header = TRUE, sep = "\t")

### Define phenotype data, expression data and feature data
pdata=mrna_series_tmat_GSE76274
edata=mrna_processed_GSE76274[,-1]
fdata=data.frame(probe_id = mrna_processed_GSE76274[,1])

# define row names of edata
row.names(edata) = mrna_processed_GSE76274[,1]

### Check the dimensions to see if the dimensions match
dim(edata)
dim(pdata)
dim(fdata)

### Define the TNBC subgroups of the phenotype data pdata
pdata$group <- factor(ifelse(grepl("LAR", paste0(pdata[,16],pdata[,17],pdata[,18],pdata[,19])), "LAR", 
                             (ifelse(grepl("MES", paste0(pdata[,16],pdata[,17],pdata[,18],pdata[,19])), "MES",
                                     (ifelse(grepl("BLIA", paste0(pdata[,16],pdata[,17],pdata[,18],pdata[,19])), "BLIA",
                                             (ifelse(grepl("BLIS", paste0(pdata[,16],pdata[,17],pdata[,18],pdata[,19])), "BLIS",
                                                     "NONE"))))))))


### expand one probe ids matching with multiple gene symbols
temp_data <- mrna_processed_GSE76274_gene_id_GEO_sel[which(mrna_processed_GSE76274_gene_id_GEO_sel$ID %in% row.names(edata)),] %>% 
              unnest_tokens(Gene.Symbol, Gene.Symbol, token = "regex", pattern = " /// ") %>% as_tibble()

### add probe id column to edata
edata$probe_id <- tolower(row.names(edata))

### match probe ids in edata with gene symbols in temp_data
edata_gene <- left_join(temp_data, edata, by = c("ID" = "probe_id")) %>% dplyr::select(-ID)

### Average the multiple measurements of expression values the same gene into one.
edata_agg <- aggregate(. ~ Gene.Symbol, data = edata_gene, mean)

### House keeping and make the Gene symbols to row names for edata_agg
edata_agg$Gene.Symbol <- factor(toupper(edata_agg$Gene.Symbol))
row.names(edata_agg) <- edata_agg$Gene.Symbol
edata_agg <- edata_agg[,-1]

### filter the data by 5000 most variable genes
edata_var <- FSbyVar(edata_agg, cut.type = "topk",value = 5000)
dim(edata_var)

### Transpose the data for PCA and further downstream analysis
edata_filt_t <- as.data.frame(t(edata_var))

pc1 = prcomp(edata_filt_t, center = TRUE)
edata_filt_t$group = pdata$group
edata_filt_t <- dplyr::select(edata_filt_t, group, everything())

autoplot(pc1, data = edata_filt_t, colour = "group")

### Plot the PCA results
plot_file_name = paste0(output_folder,"TNBC_7000_var_genes_pca_plot.pdf")

# Step 1: Call the pdf command to start the plot
pdf(file = plot_file_name,   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
autoplot(pc1, data = edata_filt_t, colour = "group", size = 5, main ="PCA") + theme_bw() +
         theme(text = element_text(size=18)) + theme(plot.title=element_text(hjust=0.5))

# Step 3: Run dev.off() to create the file!
dev.off()

### Define a functino RF_classifier for multiple calls
RF_classifier <- function(input_df){
  ### make a 70/30 split of train and test data
  set.seed(3456)
  train.index <- createDataPartition(input_df$group, p = .6, list = FALSE)
  train <- input_df[train.index,]
  test  <- input_df[-train.index,]

  train.data = train[,-1] #Set the training set
  test.data = test[,-1] #Set the testing set
  
  train.output = train[,1]  #store the labels of train data
  test.output = test[,1]  # store the true labels of test data

  ### Fitting Random Forest to the train dataset
  set.seed(234)
  classifier_RF = randomForest(x = train.data,
                               y = train.output,
                               ntree = 500,
                               importance = TRUE)
  
  ### Print out details of the classiier
  #print(classifier_RF)

  ### Draw plots of accuracy
  plot1 <- plot (classifier_RF, main='Random forest plot', lty=1, cex=1.5)
  #print(plot1)
  
  ### Prediction on test.data as a whole
  test.predict = predict(classifier_RF, test.data)
  
  result <- confusionMatrix(as.factor(test.predict), test.output,
                            positive = "pos")
  
  feat_imp_df <- importance(classifier_RF) %>% 
    data.frame() %>% 
    mutate(feature = row.names(.)) %>% arrange(desc(MeanDecreaseAccuracy)) 
  
  return(list(classifier_RF, result, feat_imp_df))    
}

RF_accuracy_metric = c()
No_top_sel_feat = 5000
edata_filt_t <- as.data.frame(t(edata_var))
edata_filt_t$group = pdata$group
edata_filt_t <- dplyr::select(edata_filt_t, group, everything())
dim(edata_filt_t)
i = 1
feat <- list()

### Repeat RF_classifer calls until all the features are used up or until number of important genes stays same 
repeat { 
  prev_No_top_sel_feat = No_top_sel_feat
  ### call the RF classifier
  RF_full <- RF_classifier(edata_filt_t)
  
  ### plot the classifier accuracy per class
  plot (RF_full[[1]], main='Random forest plot', lty=1, cex=1.5)
  
  ### confusion matrix results
  result <- RF_full[[2]]
  
  ### Feature importance results
  feat_imp_df <- RF_full[[3]]
  
  ### Filter genes/features with zero importance
  feat_imp_df[feat_imp_df$MeanDecreaseAccuracy>0, ]
  top_sel_feat <- factor(feat_imp_df[feat_imp_df$MeanDecreaseAccuracy>0, ]$feature)
  No_top_sel_feat <- length(top_sel_feat)
  
  ### Store the features in a list for retrieval
  feat[[i]] <-top_sel_feat
  i <- i + 1
  
  ### subset edata for those genes with greater than zero importance value estimated
  edata_var_sel_feat <- as.data.frame(edata_filt_t[,top_sel_feat])
  edata_var_sel_feat$group <- edata_filt_t$group
  edata_filt_t <- dplyr::select(edata_var_sel_feat, group, everything())
  
  ### collect the accuracy metrics
  temp_acc_metric <- cbind(No_sel_feat = No_top_sel_feat,
                              Accuracy_BLIA = result$byClass[41],
                              Accuracy_BLIS = result$byClass[42],
                              Accuracy_LAR = result$byClass[43],
                              Accuracy_MES = result$byClass[44], 
                              Accuracy_over_all = result$overall[[1]], 
                              Kappa = result$overall[[2]])
  
  RF_accuracy_metric <- rbind(temp_acc_metric, RF_accuracy_metric)
  ### print the accuracy metrics
  print (temp_acc_metric)
  
  ### Exit condition
  if(No_top_sel_feat == prev_No_top_sel_feat){ break }
  if(dim(edata_filt_t)[2]<3){break}
  }

### Make the accuracy metrics as dataframe
RF_accuracy_metric <- as.data.frame(RF_accuracy_metric)

### Find the best accuracy metrics and identify the classifier 
max_accuracy_index <- RF_accuracy_metric[RF_accuracy_metric$Accuracy_over_all+RF_accuracy_metric$Kappa == max(RF_accuracy_metric$Accuracy_over_all+RF_accuracy_metric$Kappa),] %>% 
  row.names() %>% as.numeric()

### Define the maximum accuracy metrics for plotting
x_max <- RF_accuracy_metric[RF_accuracy_metric$Accuracy_over_all+RF_accuracy_metric$Kappa == max(RF_accuracy_metric$Accuracy_over_all+RF_accuracy_metric$Kappa),]$No_sel_feat

### Plot the accuracy metric against number of important features
ggplot(RF_accuracy_metric, aes(x=No_sel_feat, y=Accuracy_over_all)) +
  geom_line( color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=3) +
  scale_x_reverse() +
  geom_vline(xintercept=x_max, linetype="dashed", color = "red") +
  theme_ipsum() +
  ggtitle("RF Classifier Accuracy")

### Random forest classification for the best important features.
top_important_features <- data.frame(genes = feat[dim(RF_accuracy_metric)[1] - max_accuracy_index][[1]])
top_sel_feat <- top_important_features$genes

#write.table(top_sel_feat, file = "/Users/Raghav/Pine Biotech/Research Fellows/Ariana/ISEF project/script/Github/data/Gene_signature_265.txt", sep = "\t",
            #quote = FALSE,  row.names = FALSE, col.names = TRUE)

edata_filt_t <- as.data.frame(t(edata_var))
#edata_var_sel_feat <- as.data.frame(edata_filt_t[,top_sel_feat])
edata_var_sel_feat <- dplyr::select(edata_filt_t, top_sel_feat)

pc1 = prcomp(edata_var_sel_feat, center = TRUE)
edata_var_sel_feat$group = pdata$group
edata_var_sel_feat <- dplyr::select(edata_var_sel_feat, group, everything())

autoplot(pc1, data = edata_var_sel_feat, colour = "group")

RF_full <- RF_classifier(edata_var_sel_feat)
result <- RF_full[[2]]

### Plot the heatmap comparison
plot_file_name = paste0(output_folder,"Random_forest_significant_genes_accuracy.pdf")

# Step 1: Call the pdf command to start the plot
pdf(file = plot_file_name,   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
# Specific color for each bar? Use a well known palette
library(RColorBrewer)
coul <- brewer.pal(4, "Set2") 
barplot(height=result$byClass[41:44], 
        names=row.names(result$table), 
        col=coul,
        xlab="TNBC subtypes", 
        ylab="Accuracy", 
        main="RF Accuracy",
        ylim=c(0,1),
        las = 2,
        cex.axis=1, cex.names=1,
        cex.lab=1.5
)
# Step 3: Run dev.off() to create the file!
dev.off()

end_time <- Sys.time()
