###############################################################################
###############################################################################
###############################################################################
rm(list = ls(all.names = TRUE))
set.seed(123)
###############################################################################
###############################################################################
###############################################################################
source("funzione_filtraggio.R")
library(stringr)
library(dplyr)
library(tibble)

# 1) librerie
load_libraries(get_default_packages())

# 2) lettura file
counts     <- read.table("GSE68086_TEP_data_matrix.txt", header = TRUE)
colData    <- read.table("patients_GSE68086.csv", header = TRUE, sep = ";", row.names = 1)
genes_info <- read.delim("Human.GRCh38.p13.annot.tsv")


# 3) preprocessing + DESeq2
prep <- preprocess_data(
  counts = counts,
  colData = colData,
  genes_info = genes_info,
  condition = "GBM", # can put also Lung or Breast
  min_row_sum = 100,
  padj_cutoff = 0.01,
  lfc_cutoff = 0
)

res_sig <- prep$res_sig



###############################################################################
###############################################################################
##############################################################################


# 7) network + top genes
net <- run_network_analysis(
  norm_mat = prep$norm_mat,
  group = prep$group,
  degs_ensembl = rownames(prep$res_sig), 
  annot = prep$annot,
  thr = 0.7,
  signed = T,     
  top_n = 100,    
  select = "intersection"      # you can use the following: "intersection", "union", "degree", "betweenness"
)


top_all_df <- net$top_df

top_symbols <- na.omit(top_all_df$SYMBOL)
#write.csv(file="D&B_top_100_union_symbol_LUNG.csv", top_symbols)

top_all_df_ensemble <- unique(top_all_df$ENSEMBL)

################################################################################
################################################################################
################################################################################

# PERMANOVA E PERMDISP

ensembl_interest <- top_all_df_ensemble 


out <- run_permanova_permdisp_sets(
  dds = prep$dds,
  group = prep$group,
  degs_use = rownames(prep$res_sig),
  genes_network = ensembl_interest,
  vst_blind = TRUE,
  n_iter = 100,
  n_perm_disp = 999
)


out$permanova_table
out$permdisp_table


################################################################################
################################################################################
################################################################################

# TRAIN AND TEST.  

patientsGSE156902 <- read.csv("patientGSE156902s.csv")
dim(patientsGSE156902)

load("GSE156902_TEP_Count_Matrix.RData")
tep.exprGSE156902.raw <- as.data.frame(dgeIncludedSamples$raw.counts) 
tep.exprGSE156902 <- as.data.frame(dgeIncludedSamples$counts)
dim(tep.exprGSE156902)
rm(dgeIncludedSamples)

patients <- patientsGSE156902
counts2 <- tep.exprGSE156902.raw
rm(patientsGSE156902,tep.exprGSE156902,tep.exprGSE156902.raw)

patients$name <- sub(".*\\[(.*?)\\].*", "\\1", patients$title)
colnames(counts2) <- sub("^[^-]*-(.*)", "\\1", colnames(counts2))

patientsGBM <- patients[patients$group.ch1 == "GBM",]
patientsGBM <- patientsGBM %>% filter(!str_detect(title, "t1|t2|t3|t4|t5|t6|t7|t8|t9|t10"))

patientsHC <- patients[patients$group.ch1 == "asymptomaticControls",]

patients.names.GBM <- patientsGBM$name
patients.names.HC <-  patientsHC$name

countsGBM <- counts2[,colnames(counts2) %in% patients.names.GBM]
countsHC <- counts2[,colnames(counts2) %in% patients.names.HC]

colnames(counts2[colnames(counts2) == "Maas-GBM-NICT-035G-TR2170"])
counts2 <-  counts2[,colnames(counts2) %in% c(patients.names.HC,patients.names.GBM)]

dim(counts2)
colnames(counts2[colnames(counts2) == "Maas-GBM-NICT-035G-TR2170"])
counts2 <-  counts2[,colnames(counts2) %in% c(patients.names.HC,patients.names.GBM)]
dim(counts2)

colData2 <- rbind(patientsGBM, patientsHC)

# Ricrea un vettore ordinato di cancer.type.ch1 in base ai rownames
cancer_type <- colData2[match(colnames(counts2), colData2$name), "group.ch1"]

#NORMALIZZAZIONE CPM SAMPLES INDIPENDENTI TEST DATA

oggettoDGEList2 <- DGEList(counts = counts2[,1,drop = FALSE]) 
normalized_counts.test <- cpm(oggettoDGEList2,normalized.lib.sizes = F,prior.count= 1 ,log = T)
normalized_counts.test <- as.data.frame(normalized_counts.test)

for (i in 2:ncol(counts2)){
  # Estrarre la iesima riga di TEST e farne un oggetto edgeR
  oggettoDGEList2 <- DGEList(counts = counts2[,i,drop = FALSE]) 
  ## drop = FALSE permette di continuare a considerare l'oggetto come un dataframe e quindi mantenere i nomi di riga
  
  # Normalizzare la iesima riga di test
  normalized_counts.test_row <- cpm(oggettoDGEList2,normalized.lib.sizes = F,prior.count= 1 ,log = T)
  
  # Rendere la iesima riga normalizzata un dataframe
  normalized_counts.test_row <- as.data.frame(normalized_counts.test_row)
  
  # Salvare la iesima riga
  normalized_counts.test <- cbind(normalized_counts.test, normalized_counts.test_row) 
  
}

counts_test <- normalized_counts.test

test_data <- counts_test[rownames(counts_test) %in% ensembl_interest,] 
test_data<- t(test_data)

# Ricrea un vettore ordinato di cancer.type.ch1 in base ai rownames
cancer_type <- colData2[match(rownames(test_data), colData2$name), "group.ch1"]

test_data <- as.data.frame(test_data)

test_data$cancer.type.ch1 <- cancer_type

test_data$cancer_binary <- ifelse(test_data$cancer.type.ch1 == "GBM", 0, 1)

test_data <- test_data[, !(names(test_data) %in% "cancer.type.ch1")]
dim(test_data)

write.csv(file="test.csv", test_data, row.names = F)



# TRAIN
# note: for breast cancer, you can change here directly to Breast. also to the accordingly test set


colData    <- read.table("patients_GSE68086.csv", header = TRUE, sep = ";", row.names = 1)
colData$sample_name <- sapply( colData$source_name_ch1, # Vettore al quale applicare la funzione
                               # Funzione da applicare ad ogni valore nel vettore
                               function(x) gsub("-", ".", x)
)

colnames(counts)<- sapply(colnames(counts), function(x) sub("X", "", x))

# Specificare la condizione da analizzare
condizione.da.analizzare = "GBM"

# Estrarre MALATI con la condizione specificata
condition.samples <- colData[colData$cancer.type.ch1==condizione.da.analizzare,]$sample_name

# Dimensioni del samples (controllo per la correttezza dei risultati)
n.condition.samples <- length(condition.samples)

# Estrarre il nome/codice dei sample SANI
healhty.samples <- colData[colData$cancer.type.ch1=="HC",]$sample_name
## Dimensioni del samples (controllo per la correttezza dei risultati)
length(healhty.samples)

# Unire SANI e MALATI
tutti.samples <- c(condition.samples,healhty.samples)
length(tutti.samples)
## [1] 95

# Filtrare la matrice delle conte di nuovo perchè sono state selezionate solo alcune patologie
counts <- counts[, tutti.samples]
dim(counts)

colData <- colData[colData$sample_name %in% colnames(counts),c("sample_name","cancer.type.ch1")]
colData[1:5,]

colData <- colData[match(colnames(counts), colData$sample_name), ]
all(colData$sample_name == colnames(counts))  # deve dare TRUE

group <- ifelse(colData$cancer.type.ch1 == condizione.da.analizzare,
                "Cancer", "Healthy")
group <- factor(group)
table(group)  

degs_ens <- rownames(prep$res)
counts_degs <- counts[rownames(counts) %in% degs_ens, ]

oggettoDGEList <- DGEList(counts = counts_degs[,1,drop = FALSE]) 
normalized_counts_train <- cpm(oggettoDGEList,normalized.lib.sizes = F,prior.count= 1 ,log = T)
normalized_counts_train <- as.data.frame(normalized_counts_train)

for (i in 2:ncol(counts_degs)){
  # Estrarre la iesima riga di TEST e farne un oggetto edgeR
  oggettoDGEList <- DGEList(counts = counts_degs[,i,drop = FALSE]) 
  ## drop = FALSE permette di continuare a considerare l'oggetto come un dataframe e quindi mantenere i nomi di riga
  
  # Normalizzare la iesima riga di test
  normalized_counts_train_row <- cpm(oggettoDGEList,normalized.lib.sizes = F,prior.count= 1 ,log = T)
  
  # Rendere la iesima riga normalizzata un dataframe
  normalized_counts_train_row <- as.data.frame(normalized_counts_train_row)
  
  # Salvare la iesima riga
  normalized_counts_train <- cbind(normalized_counts_train, normalized_counts_train_row) 
  
}

normalized_counts <- normalized_counts_train

filtro <- ensembl_interest[ensembl_interest %in% rownames(counts2)] #controllo che ci siano negli altri dati

# Dati per export
train_data <- normalized_counts[rownames(normalized_counts) %in% filtro,]
train_data<- t(train_data)

# Ricrea un vettore ordinato di cancer.type.ch1 in base ai rownames
cancer_type <- colData[match(rownames(train_data), colData$sample_name), "cancer.type.ch1"]

train_data <- as.data.frame(train_data)

train_data$cancer.type.ch1 <- cancer_type

train_data$cancer_binary <- ifelse(train_data$cancer.type.ch1 == "GBM", 0, 1)

train_data <- train_data[, !(names(train_data) %in% "cancer.type.ch1")]
dim(train_data)

write.csv(train_data ,file="train_data_new.csv", row.names = F)


###########################################################################
#############################################################################
##############################################################################

# TRAIN AND TEST. FOR TEST, RUN LINES 175 -260


## ---- INPUT ----
load("data/GSE183635_TEP_Count_Matrix.RData")
counts2  <- as.data.frame(TEP_Count_Matrix)
colData2 <- readxl::read_excel("data/GSE183635_patients.xlsx") %>% as.data.frame()

table(colData2$Group)

## ---- CONTROLLI ----
stopifnot(!is.null(rownames(counts2)))
stopifnot(nrow(counts2) > 0, ncol(counts2) > 0)
stopifnot(all(c("Group", "Sample ID") %in% colnames(colData2)))

## ---- 1) tieni SOLO: (Asymptomatic controls) + (Glioma) ----
is_control <- colData2$Group == "Asymptomatic controls"
is_glioma  <- str_detect(tolower(colData2$Group), "\\bglioma\\b")

colData2 <- colData2[is_control | is_glioma, , drop = FALSE]
if(nrow(colData2) == 0) stop("Dopo filtro (Asymptomatic controls + Glioma) non rimane nulla.")

table(colData2$Group)

## ---- 2) match Sample ID -> colnames(counts2) togliendo prefissi ----
cn <- colnames(counts2)

counts_key <- cn %>%
  str_replace("^\\d+-", "") %>%                 # rimuove "1-"
  str_replace("^countMatrix\\.\\d+-", "")       # rimuove "countMatrix.3488-"

map_counts <- setNames(cn, counts_key)          # key -> vero colname counts

colData2$counts_col <- unname(map_counts[as.character(colData2$`Sample ID`)])

colData2 <- colData2[!is.na(colData2$counts_col), , drop = FALSE]
if(nrow(colData2) == 0) stop("Nessun match tra colData2$`Sample ID` e colonne di counts2 (dopo rimozione prefissi).")

colData2 <- colData2 %>% distinct(counts_col, .keep_all = TRUE)

## ---- 3) subset counts e allinea ordine ----
counts2_sub <- counts2[, colData2$counts_col, drop = FALSE]
stopifnot(identical(colnames(counts2_sub), colData2$counts_col))

## ---- 4) counts numerici ----
counts2_sub[] <- lapply(counts2_sub, function(x) as.numeric(as.character(x)))

## ---- 5) CPM log ----
dge <- DGEList(counts = counts2_sub)
counts_logcpm <- cpm(dge, normalized.lib.sizes = FALSE, prior.count = 1, log = TRUE)
counts_logcpm <- as.data.frame(counts_logcpm)

## ---- 6) FILTRO GENI: ensembl_interest ----
ensembl_interest <- as.character(ensembl_interest)

gene_sub <- counts_logcpm[rownames(counts_logcpm) %in% ensembl_interest, , drop = FALSE]
if(nrow(gene_sub) == 0) stop("Nessun gene di ensembl_interest trovato nelle rownames della matrice.")

test_data <- as.data.frame(t(gene_sub))  # righe=campioni, colonne=geni
stopifnot(identical(rownames(test_data), colData2$counts_col))

## ---- 7) LABEL: cancro=0, healthy control=1 ----
test_data$cancer_binary <- ifelse(colData2$Group == "Asymptomatic controls", 1, 0)
dim(test_data)

## ---- 8) export ----
write.csv(test_data, file = "test.csv", row.names = FALSE)


################################################################################
################################################################################
################################################################################

# TRAIN AND TEST. FOR TEST, RUN LINES 175 -260

counts2 <- read.table("data/GSE107868_TEP_Count_Matrix.txt")
sample_cols <- colnames(counts2)

# --- 2) Leggi i NOMI dei file e parsali (formato: ID;Name;...;Group;Age;Gender) ---
folder <- "data/patients_GSE107868"
files  <- list.files(folder, full.names = FALSE)

# rimuove eventuale estensione (se presente) e prende solo il basename
files_base <- tools::file_path_sans_ext(basename(files))

# IMPORTANT: dai tuoi output servono 6 campi (non 5)
# ID ; Name ; (campo extra) ; Group ; Age ; Gender
parts <- str_split_fixed(files_base, pattern = ";", n = 6)

patients <- tibble(
  ID     = parts[, 1],
  Name   = parts[, 2],
  Extra  = parts[, 3],   
  Group  = parts[, 4],   
  Age    = parts[, 5],
  Gender = parts[, 6],
  file   = files
) %>%
  mutate(
    ID     = suppressWarnings(as.integer(ID)),
    Age    = suppressWarnings(as.numeric(Age)),
    Group  = str_trim(Group),
    Gender = str_trim(Gender),
    Name   = str_trim(Name)
  )

# --- 3) Normalizzazione per matchare colonne counts (con '.') e Name file (con '-') ---
norm_key <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\.", "-") %>%   # . -> -
    str_replace_all("--+", "-") %>%   # collassa doppi trattini
    str_trim()
}

map_counts <- tibble(sample_col = sample_cols) %>%
  mutate(sample_key = norm_key(sample_col))

patients2 <- patients %>%
  mutate(sample_key = norm_key(Name))

# --- 4) Join: colonna counts -> info paziente/condizione ---
sample_annotation <- map_counts %>%
  left_join(patients2, by = "sample_key") %>%
  select(sample_col, ID, Name, Group, Age, Gender, file)

# --- 5) Controlli ---
unmatched <- sample_annotation %>% filter(is.na(Group))

duplicates <- patients2 %>%
  group_by(sample_key) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

cat("Tot colonne counts:", nrow(map_counts), "\n")
cat("Match trovati:", sum(!is.na(sample_annotation$Group)), "\n")
cat("Non matchati:", nrow(unmatched), "\n")

if (nrow(duplicates) > 0) {
  cat("ATTENZIONE: ci sono sample_key duplicati nei file!\n")
  print(duplicates)
}

# guarda rapidamente distribuzione gruppi
print(table(sample_annotation$Group, useNA = "ifany"))

# --- Output utili ---
# tabella finale (colonna counts -> gruppo ecc.)
head(sample_annotation, 20)


# 0) tieni solo campioni con annotazione valida
sa <- sample_annotation %>%
  filter(!is.na(Group)) %>%
  distinct(sample_col, .keep_all = TRUE)

# 1) tieni in counts2 SOLO le colonne presenti in sample_annotation$sample_col
common_samples <- intersect(colnames(counts2), sa$sample_col)

counts2_sub <- counts2[, common_samples, drop = FALSE]

# 2) riordina sample_annotation nello stesso identico ordine delle colonne di counts2_sub
sa_sub <- sa %>%
  filter(sample_col %in% common_samples) %>%
  mutate(sample_col = factor(sample_col, levels = colnames(counts2_sub))) %>%
  arrange(sample_col) %>%
  mutate(sample_col = as.character(sample_col))

# 3) riordina anche counts2_sub per sicurezza (stesso ordine di sa_sub)
counts2_sub <- counts2_sub[, sa_sub$sample_col, drop = FALSE]

# 4) controllo "corrispondono alla perfezione"
stopifnot(identical(colnames(counts2_sub), sa_sub$sample_col))



## -----------------------------
## 1) tieni SOLO: Controls (HD) + Glioma (LGG/GBM)
##    (equivalente al tuo filtro Asymptomatic + Glioma)
## -----------------------------
is_control <- sa_sub$Group == "HD"
is_glioma  <- sa_sub$Group %in% c("LGG", "GBM")

sa_sub2 <- sa_sub[is_control | is_glioma, , drop = FALSE]
if (nrow(sa_sub2) == 0) stop("Dopo filtro (HD + LGG/GBM) non rimane nulla.")

counts2_sub2 <- counts2_sub[, sa_sub2$sample_col, drop = FALSE]
stopifnot(identical(colnames(counts2_sub2), sa_sub2$sample_col))

table(sa_sub2$Group)

## -----------------------------
## 2) (opzionale) crea un colData 'stile GEO' con una colonna 'Sample ID'
##    Qui puoi decidere cosa vuoi come Sample ID:
##    - sample_col (es: Vumc.HD.130.TR926)  -> sempre univoco e già matchato
##    - ID numerico (es: 22)               -> ok se non duplicato
## -----------------------------
colData2 <- sa_sub2 %>%
  transmute(
    `Sample ID` = sample_col,  # consigliato: usa questo come chiave univoca
    Group = Group,
    Age = Age,
    Gender = Gender,
    patient_ID = ID
  ) %>% as.data.frame()

stopifnot(identical(colnames(counts2_sub2), colData2$`Sample ID`))

## -----------------------------
## 3) counts numerici (se serve)
## -----------------------------
counts2_sub2 <- as.data.frame(counts2_sub2)
counts2_sub2[] <- lapply(counts2_sub2, function(x) as.numeric(as.character(x)))

## -----------------------------
## 4) CPM log (come nel tuo esempio)
## -----------------------------
dge <- DGEList(counts = counts2_sub2)
counts_logcpm <- cpm(dge, normalized.lib.sizes = FALSE, prior.count = 1, log = TRUE)
counts_logcpm <- as.data.frame(counts_logcpm)

## -----------------------------
## 5) FILTRO GENI: ensembl_interest (come nel tuo esempio)
## -----------------------------
ensembl_interest <- as.character(ensembl_interest)

gene_sub <- counts_logcpm[rownames(counts_logcpm) %in% ensembl_interest, , drop = FALSE]
if (nrow(gene_sub) == 0) stop("Nessun gene di ensembl_interest trovato nelle rownames della matrice.")

test_data <- as.data.frame(t(gene_sub))  # righe=campioni, colonne=geni
stopifnot(identical(rownames(test_data), colData2$`Sample ID`))

## -----------------------------
## 6) LABEL binaria (come nel tuo esempio):
##    nel tuo esempio: control=1, cancer=0
##    qui: HD=1, glioma(LGG/GBM)=0
## -----------------------------
test_data$cancer_binary <- ifelse(colData2$Group == "HD", 1, 0)
dim(test_data)

## -----------------------------
## 7) export
## -----------------------------
write.csv(test_data, file = "test_GSE107868.csv", row.names = FALSE)



################################################################################
################################################################################
################################################################################


# LUNG TEST INDIPENDENTE 
# TRAIN AND TEST. FOR TEST, RUN LINES 175 -260

counts_LUNG     <- read.table("data/GSE89843_TEP_Count_Matrix.txt", header = TRUE)
colData_LUNG    <- readxl::read_excel("data/GSE89843_patients.xlsx")

## ---- 1) Oggetti base ----
counts2  <- as.data.frame(counts_LUNG)
patients <- as.data.frame(colData_LUNG)

stopifnot("Gene.ID" %in% colnames(counts2))
stopifnot("Sample name" %in% colnames(patients))
stopifnot("Patient group" %in% colnames(patients))

rownames(counts2) <- counts2$Gene.ID
counts2$Gene.ID <- NULL

## ---- 2) Match: colData -> colonne count matrix ----
# Colonne in counts2 sono tipo "Vumc.HD.161.TR1522"
# Sample name in colData è tipo "Vumc-HD-161"
norm_id <- function(x){
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("-", "\\.")
}

cn <- colnames(counts2)
patients$name_base <- norm_id(patients$`Sample name`)

patients$name <- vapply(patients$name_base, function(b){
  m <- cn[str_detect(cn, fixed(b))]
  if(length(m) == 0) NA_character_ else m[1]
}, character(1))

patients_matched <- patients %>% filter(!is.na(name))
if(nrow(patients_matched) == 0){
  stop("Nessun match tra colData_LUNG$`Sample name` e colnames(counts_LUNG). Controlla i formati.")
}

## ---- 3) Filtra SOLO: Healthy Control e NSCLC ----
group_control <- "Healthy Control"
group_case    <- "NSCLC"

patientsHC <- patients_matched %>% filter(`Patient group` == group_control)
patientsNSCLC <- patients_matched %>% filter(`Patient group` == group_case)

if(nrow(patientsHC) == 0) stop("Nessun campione 'Healthy Control' dopo matching.")
if(nrow(patientsNSCLC) == 0) stop("Nessun campione 'NSCLC' dopo matching.")

patients.names.HC    <- patientsHC$name
patients.names.NSCLC <- patientsNSCLC$name

# Sotto-matrice counts con SOLO HC + NSCLC
counts2_sub <- counts2[, colnames(counts2) %in% c(patients.names.HC, patients.names.NSCLC), drop = FALSE]
if(ncol(counts2_sub) == 0) stop("counts2_sub vuoto: nessuna colonna selezionata.")

# ColData allineato alle colonne
colData2 <- rbind(patientsNSCLC, patientsHC)
colData2 <- colData2[match(colnames(counts2_sub), colData2$name), , drop = FALSE]

## ---- 4) Normalizzazione CPM log (prior.count=1) ----
dge <- DGEList(counts = counts2_sub)
counts_test <- cpm(dge, normalized.lib.sizes = FALSE, prior.count = 1, log = TRUE)
counts_test <- as.data.frame(counts_test)

## ---- 5) Selezione geni -> trasposizione -> label binaria -> export ----
ensembl_interest <- as.character(ensembl_interest)

test_data <- counts_test[rownames(counts_test) %in% ensembl_interest, , drop = FALSE]
if(nrow(test_data) == 0) stop("Nessun gene di ensembl_interest trovato nelle rownames della matrice.")

test_data <- t(test_data)
test_data <- as.data.frame(test_data)

# Etichetta gruppo (testuale)
test_data$patient_group <- colData2[match(rownames(test_data), colData2$name), "Patient group"]

# Binario: NSCLC = 0, Healthy Control = 1
test_data$cancer_binary <- ifelse(test_data$patient_group == group_case, 0, 1)

# Rimuovi label testuale se vuoi solo geni + binario
test_data <- test_data[, !(names(test_data) %in% "patient_group"), drop = FALSE]
dim(test_data)
write.csv(test_data, file = "test.csv", row.names = FALSE)



################################################################################

# TRAIN AND TEST. FOR TEST, RUN LINES 175 -260

counts2  <- as.data.frame(counts_LUNG)
patients <- as.data.frame(colData_LUNG)

stopifnot("Gene.ID" %in% colnames(counts2))
stopifnot("Sample name" %in% colnames(patients))
stopifnot("Patient group" %in% colnames(patients))

## ---- 1) Gene.ID -> rownames ----
rownames(counts2) <- counts2$Gene.ID
counts2$Gene.ID <- NULL

## ---- 2) Match colData -> colonne matrice count ----
# Colonne in counts2: "Vumc.HD.161.TR1522"
# Sample name: "Vumc-HD-161"
norm_id <- function(x){
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("-", "\\.")
}

cn <- colnames(counts2)
patients$name_base <- norm_id(patients$`Sample name`)

patients$name <- vapply(patients$name_base, function(b){
  m <- cn[str_detect(cn, fixed(b))]
  if(length(m) == 0) NA_character_ else m[1]
}, character(1))

patients_matched <- patients %>% filter(!is.na(name))
if(nrow(patients_matched) == 0){
  stop("Nessun match tra colData_LUNG$`Sample name` e colnames(counts_LUNG). Controlla i formati.")
}

## ---- 3) Usa TUTTI i campioni matchati ----
counts2_sub <- counts2[, colnames(counts2) %in% patients_matched$name, drop = FALSE]
if(ncol(counts2_sub) == 0) stop("counts2_sub vuoto: nessuna colonna selezionata.")

# Allinea colData2 all’ordine delle colonne
colData2 <- patients_matched[match(colnames(counts2_sub), patients_matched$name), , drop = FALSE]

## ---- 4) Normalizzazione CPM log ----
dge <- DGEList(counts = counts2_sub)
counts_test <- cpm(dge, normalized.lib.sizes = FALSE, prior.count = 1, log = TRUE)
counts_test <- as.data.frame(counts_test)

## ---- 5) Selezione geni -> trasposizione ----
ensembl_interest <- as.character(ensembl_interest)

test_data <- counts_test[rownames(counts_test) %in% ensembl_interest, , drop = FALSE]
if(nrow(test_data) == 0) stop("Nessun gene di ensembl_interest trovato nelle rownames della matrice.")

test_data <- t(test_data)
test_data <- as.data.frame(test_data)

## ---- 6) Label binaria: NSCLC=0, RESTO=1 ----
# (ci teniamo patient_group solo per costruire la label; poi lo togliamo)
test_data$patient_group <- colData2[match(rownames(test_data), colData2$name), "Patient group"]

test_data$cancer_binary <- ifelse(test_data$patient_group == "NSCLC", 0, 1)

# rimuovi patient_group (se vuoi SOLO geni + binario)
test_data <- test_data[, !(names(test_data) %in% "patient_group"), drop = FALSE]
dim(test_data)

## ---- 7) Export ----
write.csv(test_data, file = "test.csv", row.names = FALSE)


####################################################################################
####################################################################################

# TRAIN AND TEST. FOR TEST, RUN LINES 175 -260

# 2) lettura file 
counts_LUNG     <- read.table("data/GSE207586_NSCLC_CountMatrix.txt", header = TRUE)
colData_LUNG    <- readxl::read_excel("data/GSE207586_patients.xlsx")

counts2  <- as.data.frame(counts_LUNG)
patients <- as.data.frame(colData_LUNG)

stopifnot(!is.null(rownames(counts2)))
stopifnot("Sample Name" %in% colnames(patients))
stopifnot("Patient Group" %in% colnames(patients))

# counts numerici
counts2[] <- lapply(counts2, function(x) as.numeric(as.character(x)))

cn <- colnames(counts2)

norm_id <- function(x){
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("-", "\\.")   # Vumc-HD-161 -> Vumc.HD.161
}

# match robusto: preferisci la versione "base." senza seq2 se esiste, altrimenti seq2
pick_col <- function(base, cn){
  base_esc <- str_replace_all(base, "\\.", "\\\\.")   # escape dei punti per regex
  # candidate: colonne che iniziano con base.
  cand <- cn[str_detect(cn, regex(paste0("^", base_esc, "\\.")))]
  if(length(cand) == 0) return(NA_character_)
  
  # preferenza: (1) base.<qualcosa> ma NON .seq2. (2) base.*seq2*
  non_seq2 <- cand[!str_detect(cand, regex("\\.seq2\\.", ignore_case = TRUE))]
  if(length(non_seq2) > 0) return(non_seq2[1])
  
  seq2 <- cand[str_detect(cand, regex("\\.seq2\\.", ignore_case = TRUE))]
  if(length(seq2) > 0) return(seq2[1])
  
  return(cand[1])
}

patients_matched <- patients %>%
  mutate(
    name_base = norm_id(`Sample Name`),
    name = vapply(name_base, pick_col, character(1), cn = cn)
  ) %>%
  filter(!is.na(name))

if(nrow(patients_matched) == 0){
  stop("Nessun match tra colData_LUNG$`Sample Name` e colnames(counts_LUNG). Controlla i formati.")
}

patients_matched <- patients_matched %>% distinct(name, .keep_all = TRUE)

# subset counts e allineamento ordine
counts2_sub <- counts2[, patients_matched$name, drop = FALSE]
idx <- match(colnames(counts2_sub), patients_matched$name)
colData2 <- patients_matched[idx, , drop = FALSE]
rownames(colData2) <- colnames(counts2_sub)
stopifnot(identical(rownames(colData2), colnames(counts2_sub)))

table(colData2$`Patient Group`)

# CPM log
dge <- DGEList(counts = counts2_sub)
counts_test <- cpm(dge, normalized.lib.sizes = FALSE, prior.count = 1, log = TRUE)
counts_test <- as.data.frame(counts_test)

# geni di interesse
ensembl_interest <- as.character(ensembl_interest)
test_data <- counts_test[rownames(counts_test) %in% ensembl_interest, , drop = FALSE]
if(nrow(test_data) == 0) stop("Nessun gene di ensembl_interest trovato nelle rownames della matrice.")

test_data <- as.data.frame(t(test_data))
stopifnot(identical(rownames(test_data), rownames(colData2)))

# label
test_data$cancer_binary <- ifelse(colData2$`Patient Group` == "NSCLC", 0, 1)
dim(test_data)

# export
write.csv(test_data, file = "test.csv", row.names = FALSE)
































