#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GWASTools)
})

"%&%" <- function(a, b) paste0(a, b)

## ---------------- args ----------------
args <- commandArgs(trailingOnly = TRUE)
p <- as.numeric(args[1])

## ---------------- directories ----------------
dir_PGS <- "/projects/rpci/wchs/pnfioric/pgs_wchs/WCHS_BrCa_pgs/pgs_weights/bmi_weights"
PGS.files <- list.files(dir_PGS, full.names = TRUE)
PGS.scoreID <- sub("\\.txt$", "", basename(PGS.files))

setwd(dir_PGS)

dir.create(paste0("./../../", PGS.scoreID[p]), showWarnings = FALSE)
dir.create(paste0("./../../", PGS.scoreID[p], "/TotalScores"), showWarnings = FALSE)

## ---------------- load weights ----------------
weight <- fread(PGS.files[p])
setnames(weight, old = "other_allele", new = "reference_allele", skip_absent = TRUE)

weight[, `:=`(
  ea = toupper(effect_allele),
  oa = toupper(reference_allele),
  cpos = paste0(chr_name, ":", chr_position)
)]
weight1 <- weight
weight <- weight[!is.na(effect_weight) & effect_weight != 0]


removed_variants = nrow(weight1)- nrow(weight)
print(removed_variants %&% " variants were removed because of a ZERO EFFECT SIZE")
print("This corresponds to " %&% (100*removed_variants/ nrow(weight1)) %&% "% of variants!!")

## ---------------- helpers ----------------
strand_flip <- function(x) chartr("ATCG", "TAGC", x)

is_ambiguous <- function(a1, a2) {
  (a1 %in% c("A","T") & a2 %in% c("A","T")) |
    (a1 %in% c("C","G") & a2 %in% c("C","G"))
}

## ---------------- main loop ----------------
for (pop in "main") {
  
  outdir <- paste0("./../../", PGS.scoreID[p], "/", pop)
  dir.create(outdir, showWarnings = FALSE)
  
  for (chr in intersect(1:22, unique(weight$chr_name))) {
    
    message("Processing chr ", chr)
    
    outfile <- paste0(outdir, "/chr.", chr, ".scores.txt")
    outfile1 <- paste0(outfile, ".nvar")
    outfile2 <- paste0(outfile, ".valid.nvar")
    outfile3 <- paste0(outfile, ".n.unmatched.alleles")
    
    gdsDir <- "/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/split_files/" %&%
      pop %&% "/wchs_chr" %&% chr
    
    files <- list.files(gdsDir, pattern = "\\.gds$", full.names = TRUE)
    
    scores <- NULL
    nvar <- valid.nvar <- umchd.Alle <- 0
    
    for (gdsfile in files) {
      
      gds <- GdsGenotypeReader(gdsfile)
      genoData <- GenotypeData(gds)
      
      anno <- data.table(
        cpos = paste0(getChromosome(gds), ":", getPosition(gds)),
        A1 = toupper(getAlleleA(gds)),
        A2 = toupper(getAlleleB(gds))
      )
      
      geno <- getGenotype(genoData, use.names = TRUE)
      close(gds)
      
      m <- merge(
        anno,
        weight,
        by = "cpos",
        allow.cartesian = FALSE
      )
      
      if (nrow(m) == 0) next
      
      m[, allele_match :=
          (ea == A1 & oa == A2) |
          (ea == A2 & oa == A1)
      ]
      
      m[, `:=`(
        ea_flip = strand_flip(ea),
        oa_flip = strand_flip(oa),
        ambiguous = is_ambiguous(ea, oa)
      )]
      
      m[, strand_match :=
          (ea_flip == A1 & oa_flip == A2) |
          (ea_flip == A2 & oa_flip == A1)
      ]
      
      m[, usable := allele_match | (strand_match & !ambiguous)]
      
      m[, effect_weight_aligned :=
          fifelse(
            allele_match,
            effect_weight,
            fifelse(strand_match & !ambiguous, -effect_weight, NA_real_)
          )
      ]
      
      m <- m[usable & !is.na(effect_weight_aligned)]
      
      if (nrow(m) == 0) next
      
      idx <- match(m$cpos, anno$cpos)
      g <- geno[idx, , drop = FALSE]
      
      nvar <- nvar + nrow(m)
      valid.nvar <- valid.nvar + nrow(m)
      umchd.Alle <- umchd.Alle + sum(!m$allele_match)
      
      s <- t(g) %*% m$effect_weight_aligned
      scores <- if (is.null(scores)) s else cbind(scores, s)
    }
    
    if (is.null(scores)) {
      write.table(0, outfile, col.names = FALSE, row.names = FALSE, quote = FALSE)
    } else {
      scores <- rowSums(scores)
      write.table(scores, outfile, col.names = FALSE, quote = FALSE)
    }
    
    write.table(nvar, outfile1, col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(valid.nvar, outfile2, col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(umchd.Alle, outfile3, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
}
