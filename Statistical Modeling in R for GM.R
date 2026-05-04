
# Install CRAN packages
install.packages(c("tidyverse", "ggplot2", "plotly", "sva", "ggrepel", "writexl"))
install.packages(c('dplyr'))
install.packages(c('pbapply', 'lmerTest', 'dplyr', 'devtools', 
                   'magrittr', 'gtools', 'nlme', 'numDeriv', 'lavaSearch2'))
devtools::install_github("stop-pre16/lmerSeq", build_vignettes = TRUE)
# You can skip this first two line if you have BiocManager already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install these
BiocManager::install("lmerSeq")
BiocManager::install("sva")


library(tidyverse)
library(dplyr)
library(ggplot2)
library(plotly)

#File path
setwd("C:/Users/tanji/Desktop/Proteomic R Analysis/GM 40Prec Threshold")
getwd()
list.files()

metadata.df = read.csv("metadata_MATCHED.csv", sep=",")
matrix_raw_df = read.csv("FINAL_GM_QC_matrix.csv", sep=",", check.names = FALSE)


# OPTIONAL Column Name check
#setdiff(
#  c("ColumnName","TissueRegionID","ApoE4","TMA","Sex","Race",
#    "Braak","AMY","Thal"),
#  names(metadata.df)
#)


# Transform the factor variables 
metadata.df[,c("LC_stage", "TMA", "ApoE4", "Sex", "Race", "AMY",
               "Braak", "TissueDiagnosisID", "Thal", "SupraTau")] <- 
  lapply(metadata.df[,c("LC_stage","TMA","ApoE4","Sex", "Race",
                        "AMY","Braak","TissueDiagnosisID","Thal", "SupraTau")], factor)

metadata.df$Braak <- factor(metadata.df$Braak, ordered=TRUE,
                            levels = c("0","0a","0b","1a","1b","I","II","III","IV","V"))
metadata.df$Thal <- factor(metadata.df$Thal, ordered=TRUE,
                           levels = c("0","1","2","3","4","5"))
metadata.df$Plaque <- factor(metadata.df$Plaque, ordered=TRUE,
                             levels = c("0","Diffuse","Neuritic"))
metadata.df$ApoE4 <- relevel(metadata.df$ApoE4, "-")
metadata.df$LC_stage <- factor(metadata.df$LC_stage, ordered=TRUE,
                               levels = c("0","0a","0b"))
metadata.df$AMY <- factor(metadata.df$AMY, levels = c("NEG","POS"))


# ============================================================
# Build RAW intensity matrix (no log yet)
# ============================================================

peptide_ids <- matrix_raw_df$Peptide_ID                    # pull peptide IDs for row names
matrix_raw_all <- as.matrix(matrix_raw_df[, -1])           # drop Peptide_ID column, keep intensity values only
rownames(matrix_raw_all) <- peptide_ids                    # set peptide IDs as rownames of matrix


# ============================================================
# Base alignment only (NO phenotype subsetting here)
# ============================================================

metadata.df$ColumnName <- as.character(metadata.df$ColumnName)     # ensure sample IDs are plain strings
rownames(metadata.df) <- metadata.df$ColumnName                    # set rownames so we can match by sample ID

keep_ids <- intersect(rownames(metadata.df), colnames(matrix_raw_all))  # keep only sample IDs that exist in the matrix columns

matrix_raw_base <- matrix_raw_all[, keep_ids, drop = FALSE]              # subset matrix columns to matched IDs
metadata_base <- metadata.df[colnames(matrix_raw_base), , drop = FALSE]  # reorder metadata rows to match matrix column order exactly

stopifnot(identical(rownames(metadata_base), colnames(matrix_raw_base))) # hard check
metadata_base <- droplevels(metadata_base)                               # drop unused factor levels after filtering


# Identify all factor columns
factor_cols <- names(metadata_base)[sapply(metadata_base, is.factor)]

for (col in factor_cols) {
  cat("\n===== ", col, " =====\n")
  print(table(metadata_base[[col]], useNA = "ifany"))
}

# ============================================================
# Metadata and sample alignment checks (BASE matrix)
# ============================================================

anyDuplicated(colnames(matrix_raw_base))                         # check for duplicated sample IDs in matrix columns (expect 0)
summary(metadata_base)                                           # inspect metadata structure and factor levels after filtering
idx <- match(colnames(matrix_raw_base), metadata_base$ColumnName) # map matrix columns to metadata rows by sample ID
all(colnames(matrix_raw_base) == metadata_base$ColumnName)        # confirm matrix columns and metadata rows are perfectly aligned




#-------------------------------------------------------------------ONLY for PCA-------------------------------------------
# ============================================================
# LOG2 TRANSFORMATION (ONLY for PCA)
# ============================================================

matrix_log_base <- log2(matrix_raw_base + 1)   # <-- BASE MATRIX FOR PCA/ComBat VIS


# ============================================================
# Metadata and sample alignment checks (log-scale matrix)
# ============================================================

anyDuplicated(colnames(matrix_log_base))                         # This is different file name from above. This is checking for the Log file (basically the same result will be generated, but just to make the code logical when you want to run the PCA)
summary(metadata_base)                                           # 
idx <- match(colnames(matrix_log_base), metadata_base$ColumnName) # 
all(colnames(matrix_log_base) == metadata_base$ColumnName)        # 

# ============================================================
# PCA before batch correction (log scale)
# ============================================================

pca.before <- prcomp(t(matrix_log_base), scale.=FALSE)
pc.var <- pca.before$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pca.before.df <- as_tibble(pca.before$x)

pca.plot <- ggplot(pca.before.df) +
  aes(PC1, PC2, color = metadata_base$TMA) +
  geom_point(size=4) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)


# ============================================================
# ComBat batch correction (log scale)
# ============================================================

library(sva)
Batch_matrix <- ComBat(dat = matrix_log_base,
                       batch = metadata_base$TMA,
                       par.prior = TRUE)


# ============================================================
# PCA after batch correction
# ============================================================

pca.after <- prcomp(t(Batch_matrix), scale.=FALSE)
pc.var <- pca.after$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100, 1)
pca.after.df <- as_tibble(pca.after$x)

pca.plot <- ggplot(pca.after.df) +
  aes(PC1, PC2, color = metadata_base$TMA) +
  geom_point(size=4) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

#--------------------------------------------------------------------------------------------------------------------------
# ===========================================================================================================================
#                                                     START OF ANALYSIS
# ===========================================================================================================================

# Explaination of the following code format: You dont have to run this section
if (FALSE) {
  "
  You start with a “master” aligned dataset:
    
    metadata_base = all samples that match the matrix
    matrix_raw_base = all columns (samples) that match the metadata
    
    For each analysis variable (Braak, Thal, AMY, etc.) you:
      
    - filter metadata_base to keep only rows where that variable is not NA
    - save that filtered metadata into metadata.(variables)
    - subset matrix_raw_base to those same samples into matrix_raw
    - log-transform into matrix_log
    - fit the model using metadata.(variables) and matrix_log
    
    IMPORTANT: metadata.(variables) is different for each analysis. Pay attention to the file names. 
"
}


# ============================================================
# lmerSeq ANALYSES (EACH USES ITS OWN SUBSET)
# ============================================================

library(lmerSeq)
library(writexl)

out_dir <- "Analysis Summary"

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}
# ---- Braak 
keep_rows <- (!is.na(metadata_base$Braak) & metadata_base$Braak != "") #TRUE for samples that have a defined Braak value
metadata.Braak <- metadata_base[keep_rows, , drop = FALSE] #Subsets metadata to only Braak-valid samples
matrix_raw <- matrix_raw_base[, rownames(metadata.Braak), drop = FALSE] #Subsets the expression matrix by columns (samples)
metadata.Braak <- metadata.Braak[colnames(matrix_raw), , drop = FALSE] #Reorders metadata rows to match the exact column order of the matrix
stopifnot(identical(rownames(metadata.Braak), colnames(matrix_raw))) #Hard failure if alignment is broken
metadata.Braak <- droplevels(metadata.Braak)

# Recode Braak: collapse 0/0a/0b -> 0
metadata.Braak$Braak <- as.character(metadata.Braak$Braak)
metadata.Braak$Braak[metadata.Braak$Braak %in% c("0","0a","0b")] <- "0" #0, 0a, and 0b should be pre-braak
metadata.Braak$Braak[metadata.Braak$Braak %in% c("1a","1b")] <- "1a/b" # combine 1a and 1b
metadata.Braak$Braak[metadata.Braak$Braak %in% c("IV","V")] <- "IV+" # combine IV and V
metadata.Braak$Braak <- factor(metadata.Braak$Braak, ordered=TRUE,
                               levels = c("0","1a/b","I","II","III","IV+")) 
table(metadata.Braak$Braak, useNA = "ifany") #Sanity check for all the groups

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ Braak + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.Braak,
                            REML = TRUE)
colnames(model.matrix(~ Braak, data = metadata.Braak))    # inspect design matrix column names to identify the coefficient used in the model
# add in a sanity check, graph the peptide that is significant to see if it is driven a limited number of data
model_Braak <- lmerSeq.summary(fit.analysis,
                               coefficient = "Braak.L",
                               p_adj_method = "BH",
                               ddf = "Satterthwaite",
                               sort_results = TRUE)

head(model_Braak$summary_table)
write_xlsx(
  model_Braak$summary_table,
  file.path(out_dir, "GM_Braak_Linear.xlsx")
)


# ---- SupraTau
keep_rows <- (!is.na(metadata_base$Braak) & metadata_base$Braak != "")
metadata.SupraTau <- metadata_base[keep_rows, , drop = FALSE]
matrix_raw <- matrix_raw_base[, rownames(metadata.SupraTau), drop = FALSE]
metadata.SupraTau <- metadata.SupraTau[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.SupraTau), colnames(matrix_raw)))
metadata.SupraTau <- droplevels(metadata.SupraTau)
matrix_log <- log2(matrix_raw + 1)

metadata.SupraTau$BraakSupra <- ifelse(as.character(metadata.SupraTau$Braak) %in% c("0","0a","0b"), # define non-supratentorial set
                                       "NEG","POS")                                          # NEG is 0/0a/0b, POS is everything else
metadata.SupraTau$BraakSupra <- factor(metadata.SupraTau$BraakSupra, levels=c("NEG","POS"))         # set NEG as reference
table(metadata.SupraTau$Braak, metadata.SupraTau$BraakSupra, useNA = "ifany")  #Check if the level is correct
table(metadata.SupraTau$SupraTau, useNA = "ifany") 

fit.analysis <- lmerSeq.fit(form = ~ BraakSupra + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.SupraTau,
                            REML = TRUE)

colnames(model.matrix(~ BraakSupra, data = metadata.SupraTau))

model_SupraTau <- lmerSeq.summary(fit.analysis,
                                  coefficient = "BraakSupraPOS",
                                  p_adj_method = "BH",
                                  ddf = "Satterthwaite",
                                  sort_results = TRUE)

head(model_SupraTau$summary_table)

write_xlsx(
  model_SupraTau$summary_table,
  file.path(out_dir, "GM_SupraTau.xlsx")
)



# ---- PercTau (continuous)

keep_rows <- (!is.na(metadata_base$PercTau))
metadata.PercTau <- metadata_base[keep_rows, , drop = FALSE]
matrix_raw <- matrix_raw_base[, rownames(metadata.PercTau), drop = FALSE]
metadata.PercTau <- metadata.PercTau[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.PercTau), colnames(matrix_raw)))
metadata.PercTau <- droplevels(metadata.PercTau)
matrix_log <- log2(matrix_raw + 1)
c(
  Value    = sum(!is.na(metadata.PercTau$PercTau)),
  NA_Counts = sum(is.na(metadata.PercTau$PercTau))
)#sanity check how many samples remained

hist(metadata.PercTau$PercTau,
     breaks = 30,
     main = "Raw PercTau distribution",
     xlab = "PercTau") #Check for data distribution

hist(log10(metadata.PercTau$PercTau + 1e-6),
     breaks = 30,
     main = "Log10-transformed PercTau",
     xlab = "log10(PercTau)") #Check data distribution if we put a log scale

metadata.PercTau$logPercTau <- log10(metadata.PercTau$PercTau + 1e-6)



fit.analysis <- lmerSeq.fit(form = ~ logPercTau + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.PercTau,
                            REML = FALSE)

colnames(model.matrix(~ logPercTau, data = metadata.PercTau))

model_PercTau <- lmerSeq.summary(fit.analysis,
                                 coefficient = "logPercTau",
                                 p_adj_method = "BH",
                                 ddf = "Satterthwaite",
                                 sort_results = TRUE)

head(model_PercTau$summary_table)
write_xlsx(
  model_PercTau$summary_table,
  file.path(out_dir, "GM_PercTau.xlsx")
)


# ---- ActivatedMicroglia (continuous)

keep_rows <- (!is.na(metadata_base$ActivatedMicroglia))
metadata.ActivatedMicroglia <- metadata_base[keep_rows, , drop = FALSE]
matrix_raw <- matrix_raw_base[, rownames(metadata.ActivatedMicroglia), drop = FALSE]
metadata.ActivatedMicroglia <- metadata.ActivatedMicroglia[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.ActivatedMicroglia), colnames(matrix_raw)))
metadata.ActivatedMicroglia <- droplevels(metadata.ActivatedMicroglia)
matrix_log <- log2(matrix_raw + 1)
nrow(metadata.ActivatedMicroglia) #sanity check how many samples remained

hist(metadata.ActivatedMicroglia$ActivatedMicroglia,
     breaks = 30,
     main = "Raw ActivatedMicroglia distribution",
     xlab = "ActivatedMicroglia") #Check for data distribution

hist(log10(metadata.ActivatedMicroglia$ActivatedMicroglia + 1),
     breaks = 30,
     main = "Log10-transformed ActivatedMicroglia",
     xlab = "log10(ActivatedMicroglia)") #Check data distribution if we put a log scale

metadata.ActivatedMicroglia$logActivatedMicroglia <- log10(metadata.ActivatedMicroglia$ActivatedMicroglia + 1)


fit.analysis <- lmerSeq.fit(form = ~ logActivatedMicroglia + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.ActivatedMicroglia,
                            REML = TRUE)

colnames(model.matrix(~ logActivatedMicroglia, data = metadata.ActivatedMicroglia))

model_ActivatedMicroglia <- lmerSeq.summary(fit.analysis,
                                            coefficient = "logActivatedMicroglia",
                                            p_adj_method = "BH",
                                            ddf = "Satterthwaite",
                                            sort_results = TRUE)

head(model_ActivatedMicroglia$summary_table)
write_xlsx(
  model_ActivatedMicroglia$summary_table,
  file.path(out_dir, "GM_ActivatedMicroglia.xlsx")
)

# ---- PercMicro (continuous)

keep_rows <- (!is.na(metadata_base$PercMicro))
metadata.PercMicro <- metadata_base[keep_rows, , drop = FALSE]
matrix_raw <- matrix_raw_base[, rownames(metadata.PercMicro), drop = FALSE]
metadata.PercMicro <- metadata.PercMicro[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.PercMicro), colnames(matrix_raw)))
metadata.PercMicro <- droplevels(metadata.PercMicro)
matrix_log <- log2(matrix_raw + 1)
nrow(metadata.PercMicro) #sanity check how many samples remained

hist(metadata.PercMicro$PercMicro,
     breaks = 30,
     main = "Raw PercMicro distribution",
     xlab = "PercMicro") #Check for data distribution


fit.analysis <- lmerSeq.fit(form = ~ PercMicro + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.PercMicro,
                            REML = TRUE)

colnames(model.matrix(~ PercMicro, data = metadata.PercMicro))

model_PercMicro <- lmerSeq.summary(fit.analysis,
                                   coefficient = "PercMicro",
                                   p_adj_method = "BH",
                                   ddf = "Satterthwaite",
                                   sort_results = TRUE)

head(model_PercMicro$summary_table)
write_xlsx(
  model_PercMicro$summary_table,
  file.path(out_dir, "GM_PercMicro.xlsx")
)


# ---- Microglia (continuous)
keep_rows <- (!is.na(metadata_base$Microglia))
metadata.Microglia <- metadata_base[keep_rows, , drop = FALSE]
matrix_raw <- matrix_raw_base[, rownames(metadata.Microglia), drop = FALSE]
metadata.Microglia <- metadata.Microglia[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.Microglia), colnames(matrix_raw)))
metadata.Microglia <- droplevels(metadata.Microglia)
matrix_log <- log2(matrix_raw + 1)


hist(metadata.Microglia$Microglia,
     breaks = 30,
     main = "Raw Microglia distribution",
     xlab = "Microglia") # Check raw Microglia distribution


hist(log10(metadata.Microglia$Microglia + 1),
     breaks = 30,
     main = "Log10-transformed Microglia",
     xlab = "log10(Microglia)")# Check log-transformed Microglia distribution

metadata.Microglia$logMicroglia <- log10(metadata.Microglia$Microglia + 1) # Log-transform Microglia


fit.analysis <- lmerSeq.fit( # Fit mixed model
  form = ~ logMicroglia + (1 | TMA),
  expr_mat = matrix_log,
  sample_data = metadata.Microglia,
  REML = TRUE
) #Model failed to converge with max|grad| = 0.00221625 (tol = 0.002, component 1)

colnames(model.matrix(~ logMicroglia, data = metadata.Microglia)) # Check model matrix

model_Microglia <- lmerSeq.summary(
  fit.analysis,
  coefficient = "logMicroglia",
  p_adj_method = "BH",
  ddf = "Satterthwaite",
  sort_results = TRUE
) # Summarize results

head(model_Microglia$summary_table)
write_xlsx(
  model_Microglia$summary_table,
  file.path(out_dir, "GM_Microglia.xlsx")
)


# ---- LC_stage (binary: 0 + 0a = NEG, 0b = POS)

keep_rows <- (!is.na(metadata_base$LC_stage) & metadata_base$LC_stage != "")  # Filtering metadata_base to keep samples where LC_stage is not NA and not empty
metadata.LC_stage <- metadata_base[keep_rows, , drop = FALSE]  # Creating metadata.LC_stage by keeping only filtered rows from metadata_base
matrix_raw <- matrix_raw_base[, rownames(metadata.LC_stage), drop = FALSE] # Subsetting matrix_raw_base to keep only columns matching samples in metadata.LC_stage
metadata.LC_stage <- metadata.LC_stage[colnames(matrix_raw), , drop = FALSE]  # Reordering metadata.LC_stage rows to match matrix_raw column order exactly
stopifnot(identical(rownames(metadata.LC_stage), colnames(matrix_raw))) # Verifying metadata.LC_stage and matrix_raw are perfectly aligned by sample ID
metadata.LC_stage <- droplevels(metadata.LC_stage) # Removing unused factor levels in metadata.LC_stage after filtering (metadata_base remains unchanged)
matrix_log <- log2(matrix_raw + 1) # Creating matrix_log by log2-transforming matrix_raw (which came from matrix_raw_base subset)

# Define binary LC_stage
metadata.LC_stage$LC_stage_bin <- ifelse(
  as.character(metadata.LC_stage$LC_stage) %in% c("0", "0a"),
  "NEG",
  "POS"
) # Creating binary LC_stage_bin in metadata.LC_stage based on LC_stage values

metadata.LC_stage$LC_stage_bin <- factor(metadata.LC_stage$LC_stage_bin, # Converting LC_stage_bin into a factor inside metadata.LC_stage
                                         levels = c("NEG", "POS"))  # NEG as reference

# Sanity check
table(metadata.LC_stage$LC_stage, metadata.LC_stage$LC_stage_bin, useNA = "ifany") # Checking that LC_stage values were grouped correctly in metadata.LC_stage

fit.analysis <- lmerSeq.fit(form = ~ LC_stage_bin + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.LC_stage,
                            REML = TRUE)

colnames(model.matrix(~ LC_stage_bin, data = metadata.LC_stage))

model_LC_stage <- lmerSeq.summary(fit.analysis,
                                  coefficient = "LC_stage_binPOS",
                                  p_adj_method = "BH",
                                  ddf = "Satterthwaite",
                                  sort_results = TRUE)

head(model_LC_stage$summary_table)
write_xlsx(
  model_LC_stage$summary_table,
  file.path(out_dir, "GM_LC_stage_Binary.xlsx")
)


# ---- AMY
keep_rows <- (!is.na(metadata_base$AMY) & metadata_base$AMY != "") # TRUE for samples that have a defined AMY value
metadata.AMY <- metadata_base[keep_rows, , drop = FALSE] # Subsets metadata to only AMY-valid samples
matrix_raw <- matrix_raw_base[, rownames(metadata.AMY), drop = FALSE] # Subsets the expression matrix by columns (samples)
metadata.AMY <- metadata.AMY[colnames(matrix_raw), , drop = FALSE] # Reorders metadata rows to match matrix column order
stopifnot(identical(rownames(metadata.AMY), colnames(matrix_raw))) # Hard failure if alignment is broken
metadata.AMY <- droplevels(metadata.AMY)
table(metadata.AMY$AMY, useNA = "ifany")  # Sanity check: sample count per AMY group after NA removal

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ AMY + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.AMY,
                            REML = TRUE)
colnames(model.matrix(~ AMY, data = metadata.AMY))

model_Amy <- lmerSeq.summary(fit.analysis,
                             coefficient = "AMYPOS",
                             p_adj_method = "BH",
                             ddf = "Satterthwaite",
                             sort_results = TRUE)

head(model_Amy$summary_table)
write_xlsx(
  model_Amy$summary_table,
  file.path(out_dir, "GM_AMY.xlsx")
)


# ---- Thal
# Group Thal group 3, ,4 ,5 into one group as a thal 3+ 
keep_rows <- !is.na(metadata_base$Thal) & metadata_base$Thal != ""
metadata.Thal <- metadata_base[keep_rows, , drop = FALSE]

matrix_raw <- matrix_raw_base[, rownames(metadata.Thal), drop = FALSE]
metadata.Thal <- metadata.Thal[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.Thal), colnames(matrix_raw)))
metadata.Thal <- droplevels(metadata.Thal)

# Collapse Thal: 3/4/5 -> "3+"
metadata.Thal$Thal3p <- as.character(metadata.Thal$Thal)
metadata.Thal$Thal3p[metadata.Thal$Thal3p %in% c("3","4","5")] <- "3+"
metadata.Thal$Thal3p <- factor(metadata.Thal$Thal3p, ordered = TRUE,
                               levels = c("0","1","2","3+"))

# Sanity check: sample counts per grouped stage
table(metadata.Thal$Thal3p, useNA = "ifany")

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ Thal3p + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.Thal,
                            REML = TRUE)

colnames(model.matrix(~ Thal3p, data = metadata.Thal))

model_Thal <- lmerSeq.summary(fit.analysis,
                              coefficient = "Thal3p.L",
                              p_adj_method = "BH",
                              ddf = "Satterthwaite",
                              sort_results = TRUE)

head(model_Thal$summary_table)
write_xlsx(
  model_Thal$summary_table,
  file.path(out_dir, "GM_Thal_Linear.xlsx.xlsx")
)

# ---- Race
keep_rows <- (!is.na(metadata_base$Race) & metadata_base$Race != "") # TRUE for samples that have a defined Race value
metadata.Race <- metadata_base[keep_rows, , drop = FALSE]
matrix_raw <- matrix_raw_base[, rownames(metadata.Race), drop = FALSE]
metadata.Race <- metadata.Race[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.Race), colnames(matrix_raw)))
metadata.Race <- droplevels(metadata.Race)
table(metadata.Race$Race, useNA = "ifany")  # Sanity check: sample count per Race group after NA removal

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ Race + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.Race,
                            REML = TRUE)
colnames(model.matrix(~ Race, data = metadata.Race))

model_Race <- lmerSeq.summary(fit.analysis,
                              coefficient = "RaceCaucasian",
                              p_adj_method = "BH",
                              ddf = "Satterthwaite",
                              sort_results = TRUE)

head(model_Race$summary_table)
write_xlsx(
  model_Race$summary_table,
  file.path(out_dir, "GM_Race.xlsx")
)


# ---- Sex
keep_rows <- (!is.na(metadata_base$Sex) & metadata_base$Sex != "") # TRUE for samples that have a defined Sex value
metadata.Sex <- metadata_base[keep_rows, , drop = FALSE]
matrix_raw <- matrix_raw_base[, rownames(metadata.Sex), drop = FALSE]
metadata.Sex <- metadata.Sex[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.Sex), colnames(matrix_raw)))
metadata.Sex <- droplevels(metadata.Sex)
table(metadata.Sex$Sex, useNA = "ifany")  # Sanity check: sample count per Sex group after NA removal

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ Sex + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.Sex,
                            REML = TRUE)
colnames(model.matrix(~ Sex, data = metadata.Sex))

model_Sex <- lmerSeq.summary(fit.analysis,
                             coefficient = "SexM",
                             p_adj_method = "BH",
                             ddf = "Satterthwaite",
                             sort_results = TRUE)

head(model_Sex$summary_table)
write_xlsx(
  model_Sex$summary_table,
  file.path(out_dir, "GM_Sex.xlsx")
)

# ---- ApoE4
keep_rows <- (!is.na(metadata_base$ApoE4) & metadata_base$ApoE4 != "") # TRUE for samples that have a defined ApoE4 value
metadata.ApoE4 <- metadata_base[keep_rows, , drop = FALSE]
matrix_raw <- matrix_raw_base[, rownames(metadata.ApoE4), drop = FALSE]
metadata.ApoE4 <- metadata.ApoE4[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.ApoE4), colnames(matrix_raw)))
metadata.ApoE4 <- droplevels(metadata.ApoE4)
table(metadata.ApoE4$ApoE4, useNA = "ifany")  # Sanity check: sample count per ApoE4 group after NA removal

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ ApoE4 + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.ApoE4,
                            REML = TRUE)
colnames(model.matrix(~ ApoE4, data = metadata.ApoE4))

model_ApoE4 <- lmerSeq.summary(fit.analysis,
                               coefficient = "ApoE4+",
                               p_adj_method = "BH",
                               ddf = "Satterthwaite",
                               sort_results = TRUE)

head(model_ApoE4$summary_table)
write_xlsx(
  model_ApoE4$summary_table,
  file.path(out_dir, "GM_ApoE4.xlsx")
)

# AGE!!!!! As continues and grouped, such as by decades, 10 years, or young vs old (65-16)/3 = 3 groups

# ---- AGE (continuous)
keep_rows <- !is.na(metadata_base$Age) & metadata_base$Age != ""
metadata.Age <- metadata_base[keep_rows, , drop = FALSE]

matrix_raw <- matrix_raw_base[, rownames(metadata.Age), drop = FALSE]
metadata.Age <- metadata.Age[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.Age), colnames(matrix_raw)))
metadata.Age <- droplevels(metadata.Age)

# Make sure Age is numeric (important if CSV read it as character)
metadata.Age$Age <- as.numeric(as.character(metadata.Age$Age))

# Sanity checks: counts + distribution
c(Value = sum(!is.na(metadata.Age$Age)),
  NA_count = sum(is.na(metadata.Age$Age)))
range(metadata.Age$Age, na.rm = TRUE)
matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ Age + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.Age,
                            REML = TRUE)

colnames(model.matrix(~ Age, data = metadata.Age))

model_Age_cont <- lmerSeq.summary(fit.analysis,
                                  coefficient = "Age",
                                  p_adj_method = "BH",
                                  ddf = "Satterthwaite",
                                  sort_results = TRUE)

head(model_Age_cont$summary_table)

write_xlsx(
  model_Age_cont$summary_table,
  file.path(out_dir, "GM_Age_continuous.xlsx")
)


# ---- AGE (grouped into 3 bins from 16 to 65)
keep_rows <- !is.na(metadata_base$Age) & metadata_base$Age != ""
metadata.Age3 <- metadata_base[keep_rows, , drop = FALSE]

matrix_raw <- matrix_raw_base[, rownames(metadata.Age3), drop = FALSE]
metadata.Age3 <- metadata.Age3[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.Age3), colnames(matrix_raw)))
metadata.Age3 <- droplevels(metadata.Age3)

# Ensure numeric
metadata.Age3$Age <- as.numeric(as.character(metadata.Age3$Age))
range(metadata.Age$Age, na.rm = TRUE)
# Define 3 groups based on (65-16)/3
breaks <- c(16, 16 + (65 - 16)/3, 16 + 2*(65 - 16)/3, 65)  # 16, 32.333..., 48.666..., 65
metadata.Age3$AgeGroup3 <- cut(metadata.Age3$Age,
                               breaks = breaks,
                               include.lowest = TRUE,
                               right = TRUE,
                               labels = c("Young", "Mid", "Old"))

# Make ordered factor so you can test trend (L/Q)
metadata.Age3$AgeGroup3 <- factor(metadata.Age3$AgeGroup3,
                                  ordered = TRUE,
                                  levels = c("Young", "Mid", "Old"))

# Sanity check: counts per group + NA
table(metadata.Age3$AgeGroup3, useNA = "ifany")

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ AgeGroup3 + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.Age3,
                            REML = TRUE)

colnames(model.matrix(~ AgeGroup3, data = metadata.Age3))

# With 3 ordered levels, you typically get AgeGroup3.L and AgeGroup3.Q
model_Age_group3 <- lmerSeq.summary(fit.analysis,
                                    coefficient = "AgeGroup3.Q",
                                    p_adj_method = "BH",
                                    ddf = "Satterthwaite",
                                    sort_results = TRUE)

head(model_Age_group3$summary_table)

write_xlsx(
  model_Age_group3$summary_table,
  file.path(out_dir, "GM_Age_group3_Quadratic.xlsx")
)


# ---- AGE (Young vs Old extremes: 16–26 vs 55–65)

# Start from base metadata
metadata.AgeYO <- metadata_base
metadata.AgeYO$Age <- as.numeric(as.character(metadata.AgeYO$Age))

# Define Young/Old labels (middle ages remain NA)
metadata.AgeYO$YoungOld <- NA
metadata.AgeYO$YoungOld[metadata.AgeYO$Age >= 16 & metadata.AgeYO$Age <= 30] <- "Young"
metadata.AgeYO$YoungOld[metadata.AgeYO$Age >= 60 & metadata.AgeYO$Age <= 65] <- "Old"

# Keep only extremes (this automatically drops NA Age and middle ages)
metadata.AgeYO <- metadata.AgeYO[!is.na(metadata.AgeYO$YoungOld), , drop = FALSE]

# Align ONCE after final filtering
matrix_raw <- matrix_raw_base[, rownames(metadata.AgeYO), drop = FALSE]
metadata.AgeYO <- metadata.AgeYO[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.AgeYO), colnames(matrix_raw)))
metadata.AgeYO <- droplevels(metadata.AgeYO)

# Factor with Young as reference
metadata.AgeYO$YoungOld <- factor(metadata.AgeYO$YoungOld, levels = c("Young", "Old"))

# Sanity check
table(metadata.AgeYO$YoungOld, useNA = "ifany")


matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ YoungOld + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.AgeYO,
                            REML = TRUE)

colnames(model.matrix(~ YoungOld, data = metadata.AgeYO))

model_Age_young_vs_old <- lmerSeq.summary(fit.analysis,
                                          coefficient = "YoungOldOld",
                                          p_adj_method = "BH",
                                          ddf = "Satterthwaite",
                                          sort_results = TRUE)

head(model_Age_young_vs_old$summary_table)

write_xlsx(
  model_Age_young_vs_old$summary_table,
  file.path(out_dir, "GM_Age_YOUNGvsOLD.xlsx")
)



# ---- Amy_Extent

keep_rows <- !is.na(metadata_base$Amy_Extent) & metadata_base$Amy_Extent != ""
metadata.AmyExtent <- metadata_base[keep_rows, , drop = FALSE]

matrix_raw <- matrix_raw_base[, rownames(metadata.AmyExtent), drop = FALSE]
metadata.AmyExtent <- metadata.AmyExtent[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.AmyExtent), colnames(matrix_raw)))
metadata.AmyExtent <- droplevels(metadata.AmyExtent)

# Ensure numeric so we can format consistently, then convert to STRING labels
metadata.AmyExtent$Amy_Extent_num <- as.numeric(as.character(metadata.AmyExtent$Amy_Extent))

# Convert to ordinal integer score by *3:
# 0, 0.333333, 0.666667, 1  ->  0, 1, 2, 3
metadata.AmyExtent$Amy_Extent_int <- as.integer(round(metadata.AmyExtent$Amy_Extent_num * 3))

# Sanity check
table(metadata.AmyExtent$Amy_Extent_int, useNA = "ifany")

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ Amy_Extent_int + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.AmyExtent,
                            REML = TRUE)

# Inspect coefficient names (important)
colnames(model.matrix(~ Amy_Extent_int, data = metadata.AmyExtent))

# Coefficient will be the slope per +1 step (i.e., per 1/3 increase)
model_Amy_slope <- lmerSeq.summary(fit.analysis,
                                   coefficient = "Amy_Extent_int",
                                   p_adj_method = "BH",
                                   ddf = "Satterthwaite",
                                   sort_results = TRUE)
head(model_Amy_slope$summary_table)

write_xlsx(
  model_Amy_slope$summary_table,
  file.path(out_dir, "GM_AmyExtent.xlsx")
)




# ---- Amy_Extent (binary: 0 vs >0)

keep_rows <- !is.na(metadata_base$Amy_Extent) & metadata_base$Amy_Extent != ""
metadata.AmyExtent_bin <- metadata_base[keep_rows, , drop = FALSE]

matrix_raw <- matrix_raw_base[, rownames(metadata.AmyExtent_bin), drop = FALSE]
metadata.AmyExtent_bin <- metadata.AmyExtent_bin[colnames(matrix_raw), , drop = FALSE]
stopifnot(identical(rownames(metadata.AmyExtent_bin), colnames(matrix_raw)))
metadata.AmyExtent_bin <- droplevels(metadata.AmyExtent_bin)

# Numeric for robust thresholding (does not change interpretation)
x <- as.numeric(as.character(metadata.AmyExtent_bin$Amy_Extent))

# Create binary group: 0 = NoAmy, >0 = CortexPos
metadata.AmyExtent_bin$Amy_Extent_bin <- ifelse(x == 0, "NoAmy", "CortexPos")
metadata.AmyExtent_bin$Amy_Extent_bin <- factor(metadata.AmyExtent_bin$Amy_Extent_bin,
                                                levels = c("NoAmy", "CortexPos"))

# Sanity check
table(metadata.AmyExtent_bin$Amy_Extent_bin, useNA = "ifany")

matrix_log <- log2(matrix_raw + 1)

fit.analysis <- lmerSeq.fit(form = ~ Amy_Extent_bin + (1|TMA),
                            expr_mat = matrix_log,
                            sample_data = metadata.AmyExtent_bin,
                            REML = TRUE)

colnames(model.matrix(~ Amy_Extent_bin, data = metadata.AmyExtent_bin))

# Coefficient should be CortexPos vs NoAmy (reference)
model_AmyExtent_bin <- lmerSeq.summary(fit.analysis,
                                       coefficient = "Amy_Extent_binCortexPos",
                                       p_adj_method = "BH",
                                       ddf = "Satterthwaite",
                                       sort_results = TRUE)

head(model_AmyExtent_bin$summary_table)

write_xlsx(
  model_AmyExtent_bin$summary_table,
  file.path(out_dir, "GM_AmyExtent_CortexPos_vs_NoAmy.xlsx")
)



#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)

df <- model_Microglia$summary_table
colnames(df)[colnames(df) == "gene"] <- "peptide"
# If p has zeros, -log10 will create Inf, so fix that
df$p_val_raw[df$p_val_raw <= 0] <- .Machine$double.xmin

# Significance color labels (same logic as your original)
df$Color <- "NS"
df$Color[df$p_val_raw < 0.05] <- "P < 0.05"
df$Color[df$p_val_adj < 0.05] <- "FDR < 0.05"
df$Color[df$p_val_adj < 0.001] <- "FDR < 0.001"
df$Color <- factor(df$Color, levels = c("NS", "P < 0.05", "FDR < 0.05", "FDR < 0.001"))

ggplot(df, aes(x = Estimate, y = -log10(p_val_raw), color = Color)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_point(size = 2.5) +
  geom_text_repel(
    data = subset(df, p_val_adj < 0.05),
    aes(label = peptide),
    size = 2.8,
    box.padding = 0.2,
    point.padding = 0.6,
    segment.size = 0.4,
    min.segment.length = 0,
    force = 8,
    max.iter = 20000,
    max.time = 5,
    max.overlaps = Inf,
    direction = "both",
    color = "black"
  )+
  scale_color_manual(values = c("FDR < 0.001" = "dodgerblue",
                                "FDR < 0.05"  = "red",
                                "P < 0.05"    = "orange2",
                                "NS"          = "gray")) +
  labs(x = "Total Microglial Density (mix model with TMA random intercept)",
       y = "-log10(raw p)",
       color = "Significance") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")





# ==============================================================================================================================
# ==============================================================================================================================

#---------------Peptide for Braak Graph----------------
library(ggplot2)
library(dplyr)

peptide_id <- "18707_2555.2622_AVEILGNTEAAHPPSPIRCCWLR_M+Na"

# safety: peptide exists and matrix/metadata aligned
stopifnot(peptide_id %in% rownames(matrix_log))
stopifnot(identical(colnames(matrix_log), rownames(metadata.Braak)))

# plotting df (use the SAME Braak variable used in the model)
plot_df <- data.frame(
  Sample = colnames(matrix_log),
  Intensity = as.numeric(matrix_log[peptide_id, ]),
  Braak = droplevels(metadata.Braak$Braak),
  TMA = droplevels(metadata.Braak$TMA)
)

ggplot(plot_df, aes(x = Braak, y = Intensity)) +
  geom_boxplot(aes(fill = Braak),
               width = 0.4,
               outlier.shape = NA,
               alpha = 0.35) +
  geom_jitter(aes(color = TMA),
              width = 0.12,
              size = 2,
              alpha = 0.85) +
  labs(
    title = peptide_id,
    x = "Braak stage (GM)",
    y = "Peptide intensity (log2(raw + 1))",
    color = "TMA"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")


#----activated Microglia--------------

library(ggplot2)
library(dplyr)

peptide_id <- "20394_2534.2724_VQDPPGTSTDCYLLPVLKPGHFK_M+Na"

# safety
stopifnot(peptide_id %in% rownames(matrix_log))
stopifnot(identical(colnames(matrix_log), rownames(metadata.ActivatedMicroglia)))

# plotting df
plot_df <- data.frame(
  Sample = colnames(matrix_log),
  Intensity = as.numeric(matrix_log[peptide_id, ]),
  ActivatedMicroglia = as.numeric(metadata.ActivatedMicroglia$ActivatedMicroglia),
  logActivatedMicroglia = as.numeric(metadata.ActivatedMicroglia$logActivatedMicroglia),
  TMA = droplevels(metadata.ActivatedMicroglia$TMA)
)

# scatter: logActivatedMicroglia vs peptide intensity
ggplot(plot_df, aes(x = logActivatedMicroglia, y = Intensity)) +
  geom_point(aes(color = TMA), size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = peptide_id,
    x = "log10(ActivatedMicroglia + 1)",
    y = "Peptide intensity (log2(raw + 1))",
    color = "TMA"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")


#----------Bar graph for Activated Microglia-----------

mean(matrix_raw[" 11176_1222.6588_GGIVGMTLPIAR_Hydroxylation_M+Na", ] > 0)

peptide_id <- " 11176_1222.6588_GGIVGMTLPIAR_Hydroxylation_M+Na"

plot_df <- data.frame(
  Sample = colnames(matrix_raw),
  ActivatedMicroglia = as.numeric(metadata.ActivatedMicroglia$ActivatedMicroglia),
  Intensity = log2(as.numeric(matrix_raw[peptide_id, ]) + 1)
)

plot_df$ActivatedMicroglia_bin <- cut(
  plot_df$ActivatedMicroglia,
  breaks = c(0,20,40,60,Inf),
  include.lowest = TRUE
)

detection_summary <- plot_df %>%
  group_by(ActivatedMicroglia_bin) %>%
  summarise(
    n_samples = n(),
    mean_intensity = mean(Intensity, na.rm = TRUE)
  )

ggplot(detection_summary,
       aes(x = ActivatedMicroglia_bin, y = mean_intensity)) +
  geom_bar(stat = "identity", fill = "brown") +
  geom_text(aes(label = paste0("n=", n_samples)),
            vjust = -0.4,
            size = 4) +
  labs(
    x = "Activated Microglia range",
    y = "Log2 Mean peptide intensity",
    title = peptide_id
  ) +
  theme_bw(base_size = 12)

#-----------------------Total Microglial Density----------------------
peptide_id <- "11176_1222.6588_GGIVGMTLPIAR_Hydroxylation_M+Na"

# safety
stopifnot(peptide_id %in% rownames(matrix_log))
stopifnot(identical(colnames(matrix_log), rownames(metadata.Microglia)))

# plotting df
plot_df <- data.frame(
  Sample = colnames(matrix_log),
  Intensity = as.numeric(matrix_log[peptide_id, ]),
  Microglia = as.numeric(metadata.Microglia$Microglia),
  logMicroglia = as.numeric(metadata.Microglia$logMicroglia),
  TMA = droplevels(metadata.Microglia$TMA)
)

# scatter: logMicroglia vs peptide intensity
ggplot(plot_df, aes(x = logMicroglia, y = Intensity)) +
  geom_point(aes(color = TMA), size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = peptide_id,
    x = "log10(Microglia + 1)",
    y = "Peptide intensity (log2(raw + 1))",
    color = "TMA"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
#MicroglialDensity_11176_1222.6588
#--------Total Microglial Density Bar Graph----------
library(dplyr)
library(ggplot2)

# Check what fraction of samples have nonzero intensity for this peptide
mean(matrix_raw["20394_2534.2724_VQDPPGTSTDCYLLPVLKPGHFK_M+Na", ] > 0)

# Peptide to plot
peptide_id <- "20394_2534.2724_VQDPPGTSTDCYLLPVLKPGHFK_M+Na"

# Build plotting data frame using log2-transformed peptide intensity
plot_df <- data.frame(
  Sample = colnames(matrix_raw),
  Microglia = as.numeric(metadata.Microglia$Microglia),
  Intensity = log2(as.numeric(matrix_raw[peptide_id, ]) + 1)
)

# Bin Microglia into ranges
plot_df$Microglia_bin <- cut(
  plot_df$Microglia,
  breaks = c(0, 20, 40, 60, 80, 100, 120, 140, Inf),
  include.lowest = TRUE
)

# Summarize mean peptide intensity within each Microglia bin
intensity_summary <- plot_df %>%
  group_by(Microglia_bin) %>%
  summarise(
    n_samples = n(),
    mean_intensity = mean(Intensity, na.rm = TRUE)
  )

print(intensity_summary)

# Plot mean log2 peptide intensity by Microglia range
ggplot(intensity_summary, aes(x = Microglia_bin, y = mean_intensity)) +
  geom_bar(stat = "identity", fill = "brown") +
  geom_text(
    aes(label = paste0("n=", n_samples)),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    x = "Microglia range",
    y = " log2 Mean peptide intensity",
    title = peptide_id
  ) +
  theme_bw(base_size = 12)

