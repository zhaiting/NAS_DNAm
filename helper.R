## Shared helpers for the analysis documents.

library(ggplot2)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(stringr)
library(viridis)

## Fits y ~ exposure/IQR + covariates + (1|ID) for every outcome-exposure pair.
## The IQR is calculated over the complete column before lmer applies listwise deletion.

run_lmer_loop_iqr <- function(data, outcomes, exposures, covariates, random_effect = "ID") {
  results <- list()
  
  for (x in exposures) {
    iqr_val <- IQR(data[[x]], na.rm = TRUE)
    if (iqr_val == 0 || is.na(iqr_val)) next
    
    scaled_var <- paste0(x, "_scaled")
    data[[scaled_var]] <- data[[x]] / iqr_val
    
    for (y in outcomes) {
      rhs <- paste(c(scaled_var, covariates), collapse = " + ")
      formula_str <- paste0(y, " ~ ", rhs, " + (1 | ", random_effect, ")")
      formula <- as.formula(formula_str)
      
      mod <- tryCatch(
        lmer(formula, data = data, REML = FALSE),
        error = function(e) return(NULL)
      )
      
      if (!is.null(mod)) {
        tidy_mod <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE) %>%
          filter(term == scaled_var) %>%
          mutate(
            outcome = y,
            exposure = x,
            iqr = iqr_val,
            star = case_when(
              p.value < 0.001 ~ "***",
              p.value < 0.01  ~ "**",
              p.value < 0.05  ~ "*",
              p.value < 0.1   ~ ".",
              TRUE            ~ ""
            )
          ) %>%
          dplyr::select(outcome, exposure, iqr, estimate, std.error, conf.low, conf.high, p.value, star)
        
        results[[paste0(y, "_", x)]] <- tidy_mod
      }
    }
  }
  
  final_df <- bind_rows(results)
  return(final_df)
}


## circular Hallmark plot.
circ_hallmark <- function(df,
                          id_col   = "ID_clean",
                          nes_col  = "NES",
                          padj_col = "padj",
                          # radial and fill encodings
                          y_var    = c("NES", "neglogFDR"),
                          fill_var = c("neglogFDR", "NES"),
                          title    = "Enriched Hallmark pathways (GSEA)",
                          # display settings
                          sort_by  = c("padj", "NES", "ID"),   # pathway order
                          wrap     = 28,                        # label width; 0 disables wrapping
                          inner_offset = 0                      # radial baseline offset
) {
  
  y_var    <- match.arg(y_var)
  fill_var <- match.arg(fill_var)
  sort_by  <- match.arg(sort_by)
  
  need <- c(id_col, nes_col, padj_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  dat <- df %>%
    mutate(
      term      = as.character(.data[[id_col]]),
      NES       = as.numeric(.data[[nes_col]]),
      padj_num  = as.numeric(.data[[padj_col]]),
      neglogFDR = -log10(padj_num)
    )
  
  # Replace non-finite values resulting from zero or missing FDR values
  if (any(!is.finite(dat$neglogFDR))) {
    max_finite <- max(dat$neglogFDR[is.finite(dat$neglogFDR)], na.rm = TRUE)
    dat$neglogFDR[!is.finite(dat$neglogFDR)] <- max_finite + 0.5
  }
  
  # Wrap long pathway labels when requested
  if (wrap > 0) dat$term <- stringr::str_wrap(dat$term, width = wrap)
  
  # Order pathways around the circle
  dat <- switch(
    sort_by,
    padj = dat %>% arrange(padj_num, term),
    NES  = dat %>% arrange(desc(NES), term),
    ID   = dat %>% arrange(term)
  )
  dat$term <- factor(dat$term, levels = rev(unique(dat$term)))
  
  # Select the requested radial and fill variables
  dat$yval   <- dat[[y_var]]
  dat$fillv  <- dat[[fill_var]]
  if (!is.numeric(dat$yval))  stop("Chosen y_var is not numeric.")
  if (!is.numeric(dat$fillv)) stop("Chosen fill_var is not numeric.")
  
  # Define the fill-color scale
  if (fill_var == "NES") {
    fill_scale <- scale_fill_gradient2(
      name = "NES",
      low = "#2166AC", mid = "#F7F7F7", high = "#EDA35A",
      midpoint = 1
    )
  } else {
    # −log10(FDR)
    fill_scale <- scale_fill_viridis_c(
      option = "mako",
      name = expression(-log[10](italic(FDR)))
    )
  }
  
  ggplot(dat, aes(x = term, y = yval + inner_offset, fill = fillv)) +
    geom_col(width = 1, color = "black", linewidth = 0.25) +
    coord_polar(start = 0, clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    fill_scale +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank(),
      panel.grid  = element_blank(),
      plot.margin = margin(12, 12, 12, 12),
      axis.text.x = element_text(size = 9, face = "bold"),
      legend.position = "right"
    )
}
