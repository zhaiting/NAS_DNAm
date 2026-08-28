
# NAS_DNAm

Analysis code for:

> Zhai T, Zilli Vieira CL, Vokonas P, Baccarelli AA, Nagel ZD, Schwartz J, Koutrakis P. **Space weather fluctuations and epigenetic aging in an elderly male cohort from Massachusetts, USA.** *npj Aging* (2026). <doi:10.1038/s41514-026-00498-z>

## Files

|   | File | Contents |
|------------------------|------------------------|------------------------|
| 1 | `NAS_DNAm_Analysis.Rmd` | BMIQ normalization, M-values, SmartSVA surrogate variables, 450K annotation, epigenetic clock age acceleration, and the phenotype-level mixed models |
| 2 | `limma_corr.R` | Epigenome-wide association analysis for all five indicators: `limma` with `duplicateCorrelation` blocked on participant, then `bacon` correction, Benjamini–Hochberg FDR, and inflation factors |
| 3 | `NAS_DNAm_Regions_Pathways.Rmd` | CpG overlap and direction of change, DMR-highlighted Manhattan plots, and Hallmark pathway enrichment |
|  | `helper.R` | Plotting theme, the mixed-model fitting loop, and the circular Hallmark plot |

