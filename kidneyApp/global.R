################################## IMPORTING LIBRARIES ##################################
library(dplyr)
library(tidyr)
library(devtools)
library(ggplot2)
library(stringr)
library(stringdist)
library(GEOquery) 
library(R.utils)
library(Biobase)
library(reshape2)
library(ggplot2)
library(limma)
library(stats)
library(blorr)
library(viridis)
library(plotly)
library(ggrepel)
library(CPOP)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(pathview)
library(enrichplot)
library(msigdbr)
library(stringr)
library(DT)
library(rpart)
library(randomForest)

################################## HANDLING GSE36059 DATASET ##################################
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5) 
readr::local_edition(1)
# Reading in Data
GEO_GSE36059 = getGEO("GSE36059")
GSE36059 = GEO_GSE36059$GSE36059_series_matrix.txt.gz

# Load relevant matrices
eMat_GSE36059 = exprs(GSE36059)
p_GSE36059 = dplyr::rename(pData(GSE36059), diagnosis = `diagnosis:ch1`)
p_GSE36059$diagnosis = gsub("non-rejecting", "NR", p_GSE36059$diagnosis)

# Rows & patient id with control diagnosis of nephrectomy
control_idx = which(p_GSE36059$diagnosis == "nephrectomy")
remove_patient_id = rownames(p_GSE36059)[control_idx]

# Remove controls from patient and expression dataset
p_GSE36059 = p_GSE36059[!(row.names(p_GSE36059) %in% remove_patient_id),]
eMat_GSE36059 = eMat_GSE36059[, !(colnames(eMat_GSE36059) %in% remove_patient_id)]

# List of the main genes
all_gene_symbols = unlist(lapply(strsplit(fData(GSE36059)$`Gene Symbol`, ' /// ', 1), `[`, 1))

# Get index of probe_ids to keep which are not duplicates or have NA values
#idx = which(!duplicated(fData(GSE36059)$`Gene Symbol`) & !is.na(fData(GSE36059)$`Gene Symbol`))
idx = which(!duplicated(all_gene_symbols) & !is.na(fData(GSE36059)$`Gene Symbol`) & !(fData(GSE36059)$`Gene Symbol` == ""))
# Remove duplicate gene probes & probes w/o gene names
eMat_GSE36059 = eMat_GSE36059[idx,]

# Print gene symbols
gene_symbols = fData(GSE36059) %>%
  dplyr::select(ID, `Gene Symbol`) 
kept_gene_symbols = gene_symbols[idx,] %>% dplyr::pull(`Gene Symbol`)

# Map Affymetrix gene probe to gene symbol
AFFX_gene_symbols = gene_symbols[idx,] %>% mutate(`Gene Symbol` = all_gene_symbols[idx] )

# Replace probe ids with gene names
rownames(eMat_GSE36059) = lapply(strsplit(kept_gene_symbols, ' /// ', 1), `[`, 1)


################################## HANDLING GSE48581 DATSET ##################################

GEO_GSE48581 = getGEO("GSE48581")
GSE48581 = GEO_GSE48581$GSE48581_series_matrix.txt.gz

# Load relevant matrices
eMat_GSE48581 = exprs(GSE48581)
p_GSE48581 = dplyr::rename(pData(GSE48581), diagnosis = `diagnosis (tcmr, abmr, mixed, non-rejecting, nephrectomy):ch1` )

# Rename non-rejecting level
p_GSE48581$diagnosis = gsub("non-rejecting", "NR", p_GSE48581$diagnosis)


# Remove rows & patient id with control diagnosis of nephrectomy or mixed
control_idx = which(p_GSE48581$diagnosis == "nephrectomy" | p_GSE48581$diagnosis == "mixed")
remove_patient_id = rownames(p_GSE48581)[control_idx]

# Remove from patient and expression dataset
p_GSE48581 = p_GSE48581[!(row.names(p_GSE48581) %in% remove_patient_id),]
eMat_GSE48581 = eMat_GSE48581[, !(colnames(eMat_GSE48581) %in% remove_patient_id)]

# List of the main genes
all_gene_symbols = unlist(lapply(strsplit(fData(GSE48581)$`Gene Symbol`, ' /// ', 1), `[`, 1))

# Get index of probe_ids to keep which are not duplicates or have NA values
#idx = which(!duplicated(fData(GSE48581)$`Gene Symbol`) & !is.na(fData(GSE48581)$`Gene Symbol`))
idx = which(!duplicated(all_gene_symbols) & !is.na(fData(GSE48581)$`Gene Symbol`) & !(fData(GSE48581)$`Gene Symbol` == ""))
# Remove duplicate gene probes & probes w/o gene names
eMat_GSE48581 = eMat_GSE48581[idx,]

# Print gene symbols
gene_symbols = fData(GSE48581) %>%
  dplyr::select(ID, `Gene Symbol`) 
kept_gene_symbols = gene_symbols[idx,] %>% dplyr::pull(`Gene Symbol`)

# Replace probe ids with gene names
rownames(eMat_GSE48581) = lapply(strsplit(kept_gene_symbols, ' /// ', 1), `[`, 1)

################################## SELECTING GENES ##################################

# Takes lmFit model, contrast matrix, condition (ABMR or TCMR), returns top table dataframe & genes selected
ORA_contrasts = function(fit, contrast_matrix, condition){
  # For linear model fit of microarray data, compute estimated coefficients and standard errors for given contrasts
  constrast_fit = contrasts.fit(fit, contrast = contrast_matrix[, condition])
  # Adjust standard errors towards empirical derived global value
  efit = eBayes(constrast_fit, robust = TRUE)
  # Show how many genes were significant
  #summary(decideTests(efit, adjust.method = "BH", p.value = 0.05))
  tT = topTable(efit, n = 200, p.value = 0.05)
  # Rounding because significant figures
  DT::datatable(round(tT[1:100,], 2))
  return(tT)
}

# Takes pData with diagnosis column at levels ("ABMR", "TCMR", "NR"), expression matrix, returns list of topTables
get_ORA_genes = function(patient_data, expression_matrix) {
  design = model.matrix(~0 + diagnosis, data = patient_data)
  colnames(design) = gsub("diagnosis", "", colnames(design))
  fit = lmFit(expression_matrix, design)
  CM = makeContrasts(ABMR = ABMR - NR,
                     TCMR = TCMR - NR,
                     levels = design)
  tT_ABMR = ORA_contrasts(fit, CM, "ABMR")
  tT_TCMR = ORA_contrasts(fit, CM, "TCMR")
  return(list("tT_ABMR" = tT_ABMR, "tT_TCMR" = tT_TCMR))
}


# Finds only the names of the ABMR vs Healthy genes significant to both datasets
tT_ABMR_48581 = get_ORA_genes(p_GSE48581, eMat_GSE48581)$tT_ABMR
tT_ABMR_36059 = get_ORA_genes(p_GSE36059, eMat_GSE36059)$tT_ABMR
genes_selected_ABMR = intersect(rownames(tT_ABMR_36059), rownames(tT_ABMR_48581))

# Finds only the names of the TCMR vs Healthy genes significant to both datasets
tT_TCMR_48581 = get_ORA_genes(p_GSE48581, eMat_GSE48581)$tT_TCMR
tT_TCMR_36059 = get_ORA_genes(p_GSE36059, eMat_GSE36059)$tT_TCMR
genes_selected_TCMR = intersect(rownames(tT_TCMR_36059), rownames(tT_TCMR_48581))

################################## INPUT TOPTABLE OUTPUT LOG FOLD TABLE ##################################

# Input two toptables to be combined & significant genes
# Output single average log fold ratio table with ENTREZID
get_average_FC = function(topTable1, topTable2, genes_selected) {
  
  FC_1 = topTable1 %>% 
    dplyr::select(logFC) %>% 
    filter(rownames(topTable1) %in% genes_selected)
  
  FC_2 = topTable2 %>% 
    dplyr::select(logFC) %>% 
    filter(rownames(topTable2) %in% genes_selected)
  
  FC_avg = (FC_1+FC_2)/2
  
  translated = bitr(rownames(FC_avg),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db"
  )
  FC_ENTREZID = FC_avg %>%
    filter(rownames(FC_avg) %in% translated$SYMBOL)
  
  row.names(FC_ENTREZID) = translated$ENTREZID
  
  FC_ENTREZID = setNames(FC_ENTREZID$logFC, rownames(FC_ENTREZID))
  return(FC_ENTREZID)
}

# Input single toptable & significant genes
# Output single log fold ratio table with ENTREZID
get_single_FC = function(topTable, genes_selected) {
  
  FC = topTable %>% 
    dplyr::select(logFC) %>% 
    filter(rownames(topTable) %in% genes_selected)
  
  translated = bitr(rownames(FC),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db"
  )
  
  FC_ENTREZID = FC %>%
    filter(rownames(FC) %in% translated$SYMBOL)
  
  row.names(FC_ENTREZID) = translated$ENTREZID
  FC_ENTREZID = setNames(FC_ENTREZID$logFC, rownames(FC_ENTREZID))
  return(FC_ENTREZID)
}

# Get average log fold changes of significant genes between datasets
FC_ABMR_ENTREZID = get_average_FC(tT_ABMR_48581, tT_ABMR_36059, genes_selected_ABMR)
FC_TCMR_ENTREZID = get_average_FC(tT_TCMR_48581, tT_TCMR_36059, genes_selected_TCMR)

# Alternatively, get log fold changes of a single dataset
FC_ABMR_ENTREZID = get_single_FC(tT_ABMR_48581, genes_selected_ABMR)

################################## VISUALISATIONS ##################################

################################## Molecular Signatures Functional Enrichment Network Plot ##################################

# Input vector with names as ENTREZID & log-fold change values, return enrich result
get_enrich = function(geneList) {
  
  # Retrieves human molecular signatures database v7.5.1
  # Category, H -> Hallmark Gene Sets, C7 -> Immunologic Signature Gene Sets
  msigdbr_df = msigdbr(species = "human", category = "H")
  
  # Prepare term to gene input
  msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
  
  # Get over-representation results
  enrich = enricher(gene = names(geneList), TERM2GENE = msigdbr_t2g)
  return(enrich)
  
}

# Input vector with names as ENTREZID & log-fold change values, return network plot
make_network_plot = function(geneList){
  
  enrich = get_enrich(geneList)
  # ENTREZID to gene symbol 
  enrich_readable = setReadable(enrich, 'org.Hs.eg.db', 'ENTREZID')
  enrich_readable@result$Description = enrich_readable@result$Description %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
  
  # Plot a cnetplot
  p1 = cnetplot(enrich_readable, 
                showCategory = 6,
                categorySize="pvalue", 
                foldChange = geneList, 
                cex_label_gene = 0.3,
                colorEdge = TRUE,
                cex_label_category = 0.5
  )
  return(p1)
  
}

################################## DOT PLOT ##################################

# Makes strings returned by msigdb look prettier
label_func = function(string){
  
  return(str_to_sentence(str_replace_all(string, "_", " ")))
  
}

# Input vector with names as ENTREZID & log-fold change values, return dotplot
make_dotplot = function(geneList){
  
  enrich = get_enrich(geneList)
  dotplot(enrich, label_format = label_func)
}

################################## TREE PLOT ##################################

# Input vector with names as ENTREZID & log-fold change values, return tree plot
make_treeplot = function(geneList){
  
  enrich = get_enrich(geneList)
  # ENTREZID to gene symbol 
  enrich_readable = setReadable(enrich, 'org.Hs.eg.db', 'ENTREZID')
  # Calculates pairwsie similarity using Jaccard's similarity index
  similiarity_enrich = pairwise_termsim(enrich_readable)
  # Tree plot where default agglomeration method uses ward.D
  p1 = treeplot(similiarity_enrich)
  return(p1)
  
}

################################## KEGGGGGGG ##################################

# #search_kegg_organism('hsa', by='kegg_code')
# 
# # Searches for relevant pathways based on selected genes
# kk = enrichKEGG(gene = translated$ENTREZID,
#                 organism = 'hsa',
#                 pvalueCutoff = 0.05)
# 
# # Shows upregulated pathways, pay attention to ID column
# head(kk)
# 
# # Displays upregulated pathway of ID = 'hsa05332' in browser
# browseKEGG(kk, 'hsa05332')
# 
# # Saves an image of the pathway in root as {ID}.pathview.png e.g. hsa05332.pathview.png
# hsa05332 = pathview(gene.data  = geneList,
#                     pathway.id = "hsa05332",
#                     species    = "hsa",
#                     limit      = list(gene=max(abs(geneList)), cpd=1))

################################## SUMMARY VISUALISATIONS ##################################

# make_network_plot(FC_ABMR_ENTREZID)
# make_network_plot(FC_TCMR_ENTREZID)
# make_dotplot(FC_ABMR_ENTREZID)
# make_dotplot(FC_TCMR_ENTREZID)
# make_treeplot(FC_ABMR_ENTREZID)
# make_treeplot(FC_TCMR_ENTREZID)

################################## SELECTING STABLE GENES ##################################

# Input two expression matrices, phenotype dataframes, significant genes, positive class (ABMR, TCMR)
# Output cpop result
get_cpop_result = function(exp_matrix1, exp_matrix2, pdata1, pdata2, genes_selected, positive_class){
  
  # Selects only significant genes from each expression matrix
  exp_matrix1_reduced = exp_matrix1[genes_selected,]
  exp_matrix2_reduced = exp_matrix2[genes_selected,]
  
  # Finds patients whom are either non-rejection or ABMR
  indx = which(pdata1$diagnosis == 'NR' | pdata1$diagnosis == positive_class)
  
  # Selects the correct patients from gene expression ratios & outcomes
  x1 = t(exp_matrix1_reduced[, indx])
  y1 = factor(pdata1$diagnosis[indx])
  
  indx = which(pdata2$diagnosis == 'NR' | pdata2$diagnosis == positive_class)
  
  x2 = t(exp_matrix2_reduced[, indx])
  y2 = factor(pdata2$diagnosis[indx])
  
  # Perform CPOP
  cpop_result = cpop_model(
    x1 = x1, x2 = x2,
    y1 = y1, y2 = y2,
    family = "binomial",
    alpha = 1,
    n_features = 10
  )
  
  # Find combined feature matrix
  features = bind_rows(data.frame(cpop_result$z1, check.names = FALSE), data.frame(cpop_result$z2, check.names = FALSE)) %>%
    dplyr::select(cpop_result$feature)
  
  return(list("cpop_result" = cpop_result, "features" =  features, "outcome" = c(y1, y2)))
  
}

################################## CPOP ON ABMR ##################################

abmr_cpop_results = get_cpop_result(eMat_GSE48581, eMat_GSE36059, p_GSE48581, p_GSE36059, genes_selected_ABMR, "ABMR")
cpop_result = abmr_cpop_results$cpop_result
abmr_nonrej_features = abmr_cpop_results$features
abmr_nonrej_outcome = abmr_cpop_results$outcome

plot_cpop(cpop_result = cpop_result, type = "ggraph")

################################## CPOP ON TCMR ##################################

tcmr_cpop_results = get_cpop_result(eMat_GSE48581, eMat_GSE36059, p_GSE48581, p_GSE36059, genes_selected_TCMR, "TCMR")
cpop_result = tcmr_cpop_results$cpop_result
tcmr_nonrej_features = tcmr_cpop_results$features
tcmr_nonrej_outcome = tcmr_cpop_results$outcome

plot_cpop(cpop_result = cpop_result, type = "ggraph")

################################## Model Prediction ##################################

# How to Ingest a Test Case
#test_case = read.csv('example_gene_input.csv')
#colnames(test_case) = c('probe_id', 'expr')

# Input feature names df from CPOP, test_case with probe ID & expression
# Output feature vector of relevant pairwise differences
get_pairwise_differences_probe_id = function(features, test_case_probe_id){
  colnames(test_case_probe_id) = c('probe_id', 'expr')
  
  feature_names = data.frame(names(features)) %>%
    tidyr::separate(`names.features.`,c("from", "to"), "--")
  
  
  feature_differences = feature_names %>%
    dplyr::select(from, to) %>%
    apply(., 1, function(x){
      from = AFFX_gene_symbols$ID[AFFX_gene_symbols$`Gene Symbol` == x[1]]
      to = AFFX_gene_symbols$ID[AFFX_gene_symbols$`Gene Symbol` == x[2]]
      return(test_case_probe_id$expr[test_case_probe_id$probe_id == from] - test_case_probe_id$expr[test_case_probe_id$probe_id == to])
    })
  
  names(feature_differences) = names(features)
  
  return(feature_differences)
}

# Input feature names df from CPOP, test_case with gene symbol & expression
# Output feature vector of relevant pairwise differences
get_pairwise_differences_gene_symbol = function(features, test_case_gene_symbol){
  colnames(test_case_gene_symbol) = c('gene_symbol', 'expr')
  
  feature_names = data.frame(names(features)) %>%
    tidyr::separate(`names.features.`,c("from", "to"), "--")
  
  feature_differences = feature_names %>%
    dplyr::select(from, to) %>%
    apply(., 1, function(x){
      from = test_case_gene_symbol$expr[test_case_gene_symbol$gene_symbol == x[1]]
      to = test_case_gene_symbol$expr[test_case_gene_symbol$gene_symbol == x[2]]
      return(as.numeric(from)-as.numeric(to))
    })
  names(feature_differences) = colnames(features)
  
  return(feature_differences)
}

# Single Prediction Classifiers

# Logistic Regression
log_pred = function(train_features, train_outcomes, test_features) {
  df = data.frame(cbind(train_features, "Outcome" = train_outcomes))
  log_reg = glm(Outcome ~ ., family = binomial(link = "logit"), data = df)
  prediction = predict(log_reg, data.frame(t(test_features)), type = "response")
  if (prediction <= 0.5) {
    classifier_prediction = levels(df$Outcome)[1]
  } else {
    classifier_prediction = levels(df$Outcome)[2]
  }
  return(classifier_prediction)
}

# kNN
knn_pred = function(train_features, train_outcomes, test_features){
  fit5 = class::knn(
    train = train_features,
    test = test_features,
    cl = train_outcomes,
    k=5
  )
  return(fit5)
}

# SVM
svm_pred = function(train_features, train_outcomes, test_features){
  svm_res = e1071::svm(x = train_features, y = as.factor(train_outcomes))
  classifier_prediction = predict(svm_res, data.frame(t(test_features)))
  return(classifier_prediction)
}

# Trees
tree_pred = function(train_features, train_outcomes, test_features){
  df = data.frame(cbind(train_features, "Outcome" = train_outcomes))
  rpart_res = rpart(Outcome ~ ., data = df)
  prediction = predict(rpart_res, data.frame(t(test_features)))
  classifier_prediction = colnames(prediction)[which(prediction == max(prediction))]
  return(classifier_prediction)
}


# Ensembles of Trees
rf_pred = function(train_features, train_outcomes, test_features){
  rf_res = randomForest::randomForest(x = train_features, y = as.factor(train_outcomes))
  classifier_prediction = predict(rf_res, test_features)
  return(classifier_prediction)
}



# Input CPOP pairwise feature matrix, outcome vector, test features, classifier prediction, class ("ABMR" or "TCMR")
# Output ggplotly object
# class_model = {"log", "svm", "tree", "rf", "knn"}

get_PCA_plot = function(train_features, train_outcomes, test_features, positive_class, class_model) {
  
  # New x projections on rotated set of axes
  pca_df = prcomp(train_features, retx = TRUE)
  # Get variance explained by PC1
  PC1_var = round(summary(pca_df)$importance[2] * 100, 2)
  # Get variance explained by PC2
  PC2_var = round(summary(pca_df)$importance[5] * 100, 2)
  # Get PCA of test features
  PCA_pred = predict(pca_df, t(test_features))
  
  # Get classifier prediction
  # Calling function using class_model variable name
  classifier_prediction = do.call(paste(class_model, "_pred", sep=''), list(train_features, train_outcomes, test_features))
  
  
  # Proper names for classifiers
  classifier_names = as.data.frame(matrix(
    c("log", "svm", "tree", "rf", "knn", "by Logistic Regression",
      "by Simple Vector Machine", "by Tree", "by Random Forest", "by KNN"),
    ncol = 2, byrow=FALSE))
  
  classifier_name = classifier_names$V2[classifier_names$V1 == class_model]
  
  
  # Create plot dataframe with test patient
  Outcome = c(paste("Predicted as", as.character(classifier_prediction), classifier_name))
  PCA_pred_classified = cbind(PCA_pred, Outcome)
  pca_plot_df = data.frame(pca_df$x) %>%
    cbind("Outcome" = train_outcomes) %>%
    rbind(., "Test_Patient" = PCA_pred_classified)
  
  
  # Custom color palette
  rej_healthy_col = c("#4daf4a", "#E41A1C", "#377eb8")
  names(rej_healthy_col) = c("NR", positive_class, paste("Predicted as", as.character(classifier_prediction), classifier_name))
  
  # Create PCA plot with annotated prediction from classifier
  p = ggplot(pca_plot_df,
             aes(x=as.numeric(PC1),
                 y=as.numeric(PC2),
                 color=Outcome
             )) +
    geom_point(alpha=0.7) +
    scale_color_manual(values = rej_healthy_col) +
    theme_bw() +
    labs(title = paste("PCA of Pairwise Gene Expression by Kidney Graft Outcome (", positive_class, "vs No Rejection )"), 
         color = "Outcome",
         x = paste("PC1", " (", PC1_var,"%)", sep=""),
         y = paste("PC2", " (", PC2_var,"%)", sep=""))
  
  return(ggplotly(p))
}


# n: number of samples
# cvK: fold number
# n_sim: number of repeats

cross_validation = function(n, cvK, n_sim, X, y){
  
  cv_rep_acc_log = cv_rep_acc_knn = cv_rep_acc_svm = cv_rep_acc_tree = cv_rep_acc_rf = c()
  
  for (i in 1:n_sim) {
    cvSets =  cvTools::cvFolds(n, cvK)  # permute all the data, into 5 folds
    cv_acc_log = cv_acc_knn = cv_acc_svm = cv_acc_tree = cv_acc_rf = c()
    
    for (j in 1:cvK) {
      test_id = cvSets$subsets[cvSets$which == j]
      X_test = X[test_id, ]
      X_train = X[-test_id, ]
      y_test = y[test_id]
      y_train = y[-test_id]
      
      ## Binary Logistic
      df = as.data.frame(cbind(X_train, "Outcome" = y_train))
      log_reg = glm(Outcome ~ ., family = binomial(link = "logit"), data = df)
      responses = predict(log_reg, X_test, type = "response")
      preds = sapply(responses, function(x){
        if (x <= 0.5) {
          levels(df$Outcome)[1]
        } else {
          levels(df$Outcome)[2]
        }
      })
      cv_acc_log[j] = mean(y_test == preds)
      
      ## KNN - 5 Nearest Neighbours
      fit5 = class::knn(
        train = X_train,
        test = X_test,
        cl = y_train,
        k=5
      )
      cv_acc_knn[j] = mean(y_test == fit5)
      
      ## SVM
      svm_res = e1071::svm(x = X_train, y = as.factor(y_train))
      fit = predict(svm_res, X_test)
      cv_acc_svm[j] = mean(y_test == fit)
      
      ## Trees
      
      #df = as.data.frame(cbind(X_train, "Outcome" = y_train))
      rpart_res = rpart(Outcome ~ ., data = df)
      prediction = predict(rpart_res, X_test)
      levels = colnames(prediction)
      preds = apply(prediction, 1, function(x){
        levels[which(x == max(x))]
      })
      cv_acc_tree[j] = mean(y_test == preds)
      
      
      #classifier_prediction = colnames(prediction)[which(prediction == max(prediction))]
      #return(classifier_prediction)
      
      ## Random Forest
      rf_res = randomForest::randomForest(x = X_train, y = as.factor(y_train))
      fit = predict(rf_res, X_test)
      cv_acc_rf[j] = mean(y_test == fit)
      
    }
    cv_rep_acc_log = append(cv_rep_acc_log, mean(cv_acc_log))
    cv_rep_acc_knn = append(cv_rep_acc_knn, mean(cv_acc_knn))
    cv_rep_acc_svm = append(cv_rep_acc_svm, mean(cv_acc_svm))
    cv_rep_acc_tree = append(cv_rep_acc_tree, mean(cv_acc_tree))
    cv_rep_acc_rf = append(cv_rep_acc_rf, mean(cv_acc_rf))
  }
  return(cbind(cv_rep_acc_log, cv_rep_acc_knn, cv_rep_acc_svm, cv_rep_acc_tree, cv_rep_acc_rf))
}


# Calls cross-validation function and then plots results
# cvK = 5
# n_sim = 10
# X = tcmr_nonrej_features
# y = tcmr_nonrej_outcome
# n = nrow(X)

get_cross_val_plot = function(n, cvK, n_sim, X, y,ab) {
  
  # Get cross validation accuracy
  acc = cross_validation(n, cvK, n_sim, X, y)
  acc = as.data.frame(acc)
  
  # Positive class either TCMR or ABMR
  positive_class = levels(y)[2]
  
  # For pretty renaming
  names(acc) = c("Logistic Regression", "k-Nearest-Neighbours", "Simple Vector Machine", "Tree", "Random Forest")
  
  # Producing boxplot
  p = acc %>%
    pivot_longer(cols = 1:5,
                 names_to = "Classifier",
                 values_to = "Accuracy") %>%
    mutate( type=ifelse(Classifier==ab,
                        "Highlighted","Normal")) %>%
    ggplot() +
    aes(x = Classifier,
        y = Accuracy,fill=type, alpha=type) +
    geom_boxplot() +
    scale_fill_manual(values=c("#377eb8", "grey")) +
    scale_alpha_manual(values=c(1,0.1)) +
    theme_bw() +
    theme(legend.position = "none")+
    labs(
      title = paste("Accuracy Comparison of Classifiers on Distinguishing", positive_class, "vs Non-Rejection"),
      subtitle = "Repeated 5-fold Cross Validation (N=10)",
      xlab = "Classifiers",
      ylab = "Accuracy"
    )
  return(ggplotly(p))
}

get_genes_for_sliders = function(features){
  names(features) %>%
  sapply(., function(x){str_split(x,'--')}) %>%
  unlist() %>%
  unique()
}
