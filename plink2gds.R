## Convert PLINK to GDS


library(GWASTools)
library(SNPRelate)
library(dplyr)
library(data.table)
#library(argparse)

"%&%"=function(a,b) paste(a,b,sep="")

#parser <- ArgumentParser()
#parser$add_argument("--bed", help = "PLINK genotype bed file", required = T)
#parser$add_argument("--fam", help = "PLINK genotype fam file", required = T)
#parser$add_argument("--bim", help = "PLINK genotype bim file", required = T)
#parser$add_argument("--out", help = "outfile for GDS", required = T)

#args <- parser$parse_args()

#bed_file<-args$bed
#bim_file<-args$bim
#fam_file<-args$fam
#out_file<-args$out

file_path<- "/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/split_files/subfiles/wchs_3048_chr"
#bim_file<-"/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/wchs_3048_chr"
#fam_file<-"/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/wchs_3048_chr"
out_file<-"/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/split_files/main/wchs_chr"
#wd <-""

max_subgroups_per_chromosome <- 10

for (i in c(9)) {
  num_subgroups <- 0  # Initialize the number of subgroups found for this chromosome
  
  # Test for how many files there are
  for (j in 1:max_subgroups_per_chromosome) {
    bed_file <- file_path  %&% i %&% "/chr" %&% i %&% "."%&% j%&% ".bed"
    fam_file <- file_path  %&% i %&%  "/chr" %&% i %&%"."%&% j%&% ".fam"
    bim_file <- file_path  %&% i %&%  "/chr" %&% i %&%"."%&% j%&% ".bim"
    print(bim_file)
    
    # Create the directory if it doesn't exist
    if (!file.exists(out_file %&% i)) {
      dir.create(out_file %&% i , recursive = TRUE)  # Create the directory and any missing parent directories
    }
    
    # Check if the files exist
    if (file.exists(bed_file) && file.exists(fam_file) && file.exists(bim_file)) {
      num_subgroups <- num_subgroups + 1
      
      # Here, you can read in the subfile and perform any necessary operations
      snpgdsBED2GDS(
        bed.fn = bed_file,
        fam.fn = fam_file,
        bim.fn = bim_file,
        out.gdsfn = out_file %&% i %&% "/chr" %&% i %&% "." %&% j %&% ".gds",
        family = FALSE,
        compress.annotation = "ZIP.max",
        option = NULL,
        verbose = TRUE
      )
    } else {
      break  # Stop searching for subgroups if a file is missing
    }
  }
  
  cat(sprintf("Chromosome %d has %d subgroups.\n", i, num_subgroups))
}
