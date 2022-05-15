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

################################## HANDLING GSE36059 DATASET ##################################

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
