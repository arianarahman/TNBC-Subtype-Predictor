#### Test in data for prediction combined data with GSE76275_id, GSE48390_id, GSE21653_id GSE20711_id data

### Define the input and output folders
data_folder <- "D:/TNBC-Project/data/"
output_folder <- "D:/TNBC-Project/output/"

### Load the processed training dataset for training the RF classifier
input_file_name = paste0(data_folder,"mrna_train_trans_265.txt")
mrna_processed_train <- read.table(input_file_name, header = TRUE, row.names = 1, sep = "\t")

### Load the mrna_expresion table for prediction
input_file_name = paste0(data_folder,"mrna_data_for_pred_265gene.txt")
mrna_processed_test <- read.table(input_file_name, header = TRUE, row.names = 1, sep = "\t")

### Define test data for testing
tdata=mrna_processed_test[,-1]

pc1 = prcomp(tdata, center = TRUE)
autoplot(pc1, data = tdata_var, colour = "red")

### Plot the PCA results
plot_file_name = paste0(output_folder,"mrna_prediction_data_pca_plot.pdf")

# Step 1: Call the pdf command to start the plot
pdf(file = plot_file_name,   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
autoplot(pc1, data = tdata, colour = "red", size = 5, main ="PCA") + theme_bw() +
  theme(text = element_text(size=18)) + theme(plot.title=element_text(hjust=0.5))

# Step 3: Run dev.off() to create the file!
dev.off()

### Train the model with training data set
train.data = mrna_processed_train[,-1] #Set the training set
test.data = tdata #Set the testing set

dim(train.data)
dim(test.data)

train.output = factor(mrna_processed_train[,1])  #store the labels of train data

### Fitting Random Forest to the train dataset
set.seed(234)
classifier_RF = randomForest(x = train.data,
                             y = train.output,
                             ntree = 500,
                             importance = TRUE)

### Print out details of the classiier
print(classifier_RF)

### Draw plots of accuracy
plot (classifier_RF, main='Random forest plot', lty=1, cex=1.5)

### Prediction on test.data as a whole
test.predict = predict(classifier_RF, test.data)

### Plot heat map of transcriptomic profiles for visual verification
# scale the training data for the heatmap
train.data.scaled <- as.matrix(train.data) %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

#Set annotation
ann = data.frame(group=mrna_processed_train$group)
colours <- list('group' = c('MES' = 'red', 'LAR' = 'royalblue', 'BLIS' = 'limegreen', 'BLIA' = 'gold'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

split <- factor(mrna_processed_train$group, levels=c("MES","LAR","BLIA","BLIS"))
hmap1 <- Heatmap(
  t(train.data.scaled),
  name = "TNBC mRNA top 250 genes", 
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  #column_title = "Train Data",
  column_split=split, cluster_row_slices = FALSE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  width = unit(100, "mm"),
  top_annotation=colAnn)


test.data$group <- test.predict
test.data <- dplyr::select(test.data, group, everything())

ann = data.frame(group=test.data$group)
colours <- list('group' = c('MES' = 'red', 'LAR' = 'royalblue', 'BLIS' = 'limegreen', 'BLIA' = 'gold'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

split <- factor(test.data$group, levels=c("MES","LAR","BLIA","BLIS"))
hmap2 <- Heatmap(
  t(test.data[,-1]),
  name = "TNBC mRNA top 250 genes", 
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_row_dend = TRUE,
  row_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  #column_title = "Test Data",
  column_split=split, cluster_row_slices = FALSE,
  width = unit(100, "mm"),
  top_annotation=colAnn)

### Plot the heatmap comparison
plot_file_name = paste0(output_folder,"mrna_prediction_data_vs_train_data_heatmap.pdf")

# Step 1: Call the pdf command to start the plot
pdf(file = plot_file_name,   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plot with R code
draw(hmap1 + hmap2, heatmap_legend_side="left")

# Step 3: Run dev.off() to create the file!
dev.off()

end_time <- Sys.time()
