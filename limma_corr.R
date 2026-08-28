##############################################################################
##  EWAS with limma and duplicateCorrelation (random intercept for ID)
##
##  All five space weather indicators. Both inputs are written by
##  NAS_DNAm_Analysis.Rmd: dnam_mval.RDS holds the M-values, and
##  pheno_methy_Aug2025.RDS holds the *_IQR exposure columns and the est450k
##  cell fractions this design uses.
##############################################################################
message("Started at: ", Sys.time())

suppressPackageStartupMessages({
  library(limma)      # linear models and empirical Bayes moderation
  library(bacon)      # bias and inflation correction
})

## --- paths and data --------------------------------------------------------
path_data <- "../Data/"
path_out  <- "../Output/"

message("Reading data...")

dnam   <- readRDS(file.path(path_data, "dnam_mval.RDS"))             # CpG × samples (M-values)
pheno  <- readRDS(file.path(path_data, "pheno_methy_Aug2025.RDS"))   # sample metadata

## Align samples using samplename_450k
pheno  <- as.data.frame(pheno)
common <- intersect(pheno$samplename_450k, colnames(dnam))
stopifnot(length(common) > 0)
pheno  <- pheno[match(common, pheno$samplename_450k), , drop = FALSE]
dnam   <- dnam[, match(common, colnames(dnam)), drop = FALSE]
rownames(pheno) <- pheno$samplename_450k

stopifnot(identical(colnames(dnam), pheno$samplename_450k))           # verify sample alignment
message("Data loaded: ", nrow(dnam), " CpGs × ", ncol(dnam), " samples")

## --- indicators and output filename tags ----------------------------------
indicators <- c(CRII = "CRII.mov30_IQR",
                ntr  = "ntr.mov30_IQR",
                SSN  = "SSN.mov30_IQR",
                IMF  = "IMF.mov30_IQR",
                kp   = "Kp_index.mov30_IQR")

stopifnot(all(indicators %in% names(pheno)))

covariates <- paste(
  "AGE + SMK + packyrs + bmi + educmax + Statin_flag + med + diabete +",
  "temp_c.mov30 + VISIB.mov30 + year + season +",
  "CD8T_est450k + CD4T_est450k + NK_est450k + Bcell_est450k + Mono_est450k + Gran_est450k +",
  "sv1 + sv2 + sv3 + sv4 + sv5 + sv6 + sv7 + sv8"
)

run_ewas <- function(tag, coef_name) {

  message("\n=== ", tag, " (", coef_name, ") ===")

  ## --- (1) design matrix ---------------------------------------------------
  design <- model.matrix(as.formula(paste("~", coef_name, "+", covariates)),
                         data = pheno)

  ## --- (2) estimate subject-level correlation on a random subset ----------
  message("Sampling random subset for duplicateCorrelation...")

  set.seed(1996)
  n_sub  <- min(10000, nrow(dnam))                   # ≤ 10 000 rows or all
  idx    <- sample(nrow(dnam), n_sub)
  dt <- system.time({
    corfit <- duplicateCorrelation(dnam[idx, ], design, block = pheno$ID)
  })
  message("duplicateCorrelation elapsed: ", round(dt["elapsed"], 1), " sec. Memory used: ",
          format(object.size(dnam), units = "GB"))
  message("Estimated consensus correlation ρ = ", signif(corfit$consensus, 3))

  ## Save the estimated correlation object for reproducibility
  saveRDS(corfit, file.path(path_out, sprintf("corfit_dupCor_%s_mov30.rds", tag)))

  ## --- (3) fit model with consensus correlation ---------------------------
  message("Fitting limma model with consensus correlation...")
  fit <- lmFit(dnam, design,
               block       = pheno$ID,
               correlation = corfit$consensus)

  ## --- (4) empirical-Bayes moderation ------------------------------------
  message("Running empirical Bayes moderation (eBayes)...")
  fit <- eBayes(fit)

  ## --- (5) extract results for the exposure of interest -------------------
  message("Extracting topTable results for coefficient: ", coef_name)
  tt <- limma::topTable(fit,
                        coef    = coef_name,
                        number  = nrow(dnam),
                        sort.by = "none")

  ewas <- data.frame(
    CpG  = rownames(tt),
    beta = tt$logFC,
    t    = tt$t,
    p    = tt$P.Value,
    row.names = NULL
  )

  ## --- (6) bias correction and FDR adjustment -----------------------------
  message("Applying bacon correction and BH FDR adjustment")
  bc            <- bacon(ewas$t)
  ewas$p_bacon  <- pval(bc)[, 1]
  ewas$q_bh     <- p.adjust(ewas$p_bacon, method = "BH")

  message(tag, " inflation: ", signif(inflation(bc), 4),
          " | FDR < 0.05: ", sum(ewas$q_bh < 0.05, na.rm = TRUE), " CpGs")

  message("Saving EWAS results (limma) ...")
  saveRDS(ewas, file.path(path_out, sprintf("EWAS_%s_mov30_results_limma.rds", tag)))
}

for (tag in names(indicators)) run_ewas(tag, indicators[[tag]])

message("Finished at: ", Sys.time())
