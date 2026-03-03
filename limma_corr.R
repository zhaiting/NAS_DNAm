##############################################################################
##  EWAS with limma + duplicateCorrelation  (random intercept on ID)
##############################################################################
message("Started at: ", Sys.time())

suppressPackageStartupMessages({
  library(limma)      # fast linear models + empirical-Bayes
  library(bacon)      # bias / genomic-control correction
})

## --- paths and data --------------------------------------------------------
path_data <- "../Data/"
path_out  <- "../Output/"

message("Reading data...")

dnam   <- readRDS(file.path(path_data, "dnam.RDS"))        # CpG × samples (M-values)
pheno  <- readRDS(file.path(path_data, "pheno_methy.RDS")) # sample metadata
stopifnot(all(colnames(dnam) == rownames(pheno)))          # sanity check
message("Data loaded: ", nrow(dnam), " CpGs × ", ncol(dnam), " samples")

## --- (1) design matrix: identical fixed effects ---------------------------
design <- model.matrix(
  ~ CRII.mov30_IQR + AGE + SMK + packyrs + bmi + educmax + Statin_flag + med + diabete +
    temp_c.mov30 + VISIB.mov30 + year + season +
    CD8T_est450k + CD4T_est450k + NK_est450k + Bcell_est450k + Mono_est450k + Gran_est450k +
    sv1 + sv2 + sv3 + sv4 + sv5 + sv6 + sv7 + sv8 ,
    data = pheno
)

## --- (2) estimate subject-level correlation on a RANDOM subset ------------
message("Sampling random subset for duplicateCorrelation...")

set.seed(1996) 
n_sub  <- min(10000, nrow(dnam))                     # ≤ 10 000 rows or all
idx    <- sample(nrow(dnam), n_sub)
dt <- system.time({
  corfit <- duplicateCorrelation(dnam[idx, ], design, block = pheno$ID)
})
message("duplicateCorrelation elapsed: ", round(dt["elapsed"], 1), " sec. Memory used: ",
        format(object.size(dnam), units="GB"))
message("Estimated consensus correlation ρ = ", signif(corfit$consensus, 3))

## save the full corfit object for future reference
saveRDS(corfit, file.path(path_out, "corfit_dupCor_crii_mov30.rds"))

## --- (3) fit model with consensus correlation -----------------------------
message("Fitting limma model with consensus correlation...")
fit <- lmFit(dnam, design,
             block       = pheno$ID,
             correlation = corfit$consensus)

## --- (4) empirical-Bayes moderation ---------------------------------------
message("Running empirical Bayes moderation (eBayes)...")
fit <- eBayes(fit)

## --- (5) pull results for exposure of interest ----------------------------
message("Extracting topTable results for coefficient: CRII.mov30_IQR")
tt <- limma::topTable(fit,
               coef   = "CRII.mov30_IQR",
               number = nrow(dnam),
               sort.by = "none")

ewas <- data.frame(
  CpG  = rownames(tt),
  beta = tt$logFC,
  t    = tt$t,
  p    = tt$P.Value,
  row.names = NULL
)

## --- (6) bias-correction + FDR --------------------------------------------
message("Applying bacon correction and BH FDR adjustment")
bc            <- bacon(ewas$t)
ewas$p_bacon  <- pval(bc)[, 1]
ewas$q_bh     <- p.adjust(ewas$p_bacon, method = "BH")

message("Saving EWAS results (limma) ...")
saveRDS(ewas, file.path(path_out, "EWAS_CRII_mov30_results_limma.rds"))
message("Finished at: ", Sys.time())
