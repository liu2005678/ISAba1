BiocManager::install("DEP")
setwd("D:\\37662\\manuscript")
library("DEP")
#run_app("TMT")
data=read.csv("D:\\37662\\manuscript/P20191202271_proteins.csv",header = T,na.strings=c(" ","",0,"NA"))
colnames(data)

data$ID=data$locus_tag
data$name=data$gene.names
data$name %>% duplicated() %>% any()
data %>% group_by(gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "gene.names", "locus_tag", delim = ";")
data$name %>% duplicated() %>% any()

# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("Abundances_", colnames(data_unique)) # get LFQ column numbers
label=c("WT.1","WT.2","RM1.1","RM1.2","RM2.1","RM2.3")
condition=c("WT","WT","RM1","RM1","RM2","RM2")
replicate=c("1","2","1","2","1","2")
experimental_design =data.frame(label=label,condition=condition,replicate=replicate)
data_se <- make_se(data_unique, LFQ_columns, experimental_design)
data_se
plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 1)
plot_numbers(data_filt)
plot_coverage(data_filt)
data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)
plot_missval(data_filt)
plot_detect(data_filt)
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)
plot_imputation(data_norm, data_imp)
data_diff <- test_diff(data_imp, type = "control", control = "WT")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep, significant = TRUE, lower = 0.5, upper = 1, pal = "Spectral")
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))
plot_volcano(dep, contrast = "RM1_vs_WT", label_size = 4, add_names = TRUE)
plot_volcano(dep, contrast = "RM2_vs_WT", label_size = 4, add_names = TRUE)
plot_single(dep, proteins = c("queG","Nnr1","GNAT","bioB"))
plot_cond(dep)
df_wide <- get_df_wide(dep)
write.csv(df_wide,file="deg.csv")
