
library(data.table)
library(dplyr)
library(GWASTools)

"%&%"=function(a,b) paste(a,b,sep="")

pgs.dir <- "/projects/rpci/wchs/pnfioric/pgs_wchs/WCHS_BrCa_pgs/pgs_weights/"
out.dir <- "/projects/rpci/wchs/pnfioric/pgs_wchs/WCHS_BrCa_pgs/troubleshooting/"
PGS_list<-c("PGS000004", "PGSAABCGS_black", "PGSGao2022_BrCa", "PGSGao2022_ERNEG", "PGSGao2022_ERPOS", "PGSShieh2023")

# Read in an combined all PGS into one data.table
# PGS SNPs will be labeled according to PGS

pgs_total<-data.table()
for (PGS in PGS_list){
  a<-fread(pgs.dir %&% PGS %&% ".txt", header = T)
  a$pgs <- PGS
  pgs_total <- bind_rows(pgs_total, a)
}


bed_dir<- "/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/split_files/"

for (i in 1:22){
  bim <- fread(bed_dir %&% "wchs_3048_chr" %&% i %&% ".bim", header = F)
  gdsDir = "/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/split_files/main/wchs_chr" %&% i 
  files = list.files(gdsDir, pattern = paste("chr", i, ".*.gds$", sep = ""), full.names = T)
  
  print("There are " %&% length(files) %&%" gds subfiles for chromosome " %&% i %&% ".")
  
  chr_files<-data.table()
  for (j in 1:length(files)){
    input_file <- gdsDir %&%"/chr" %&% i %&%  "."%&% j %&%".gds"
    file<-GdsGenotypeReader(input_file)
    
    b <- data.table()
    b$chromosome <- getChromosome(file)
    b$position <- getPosition(file)
    b$alleleA <- getAlleleA(file)
    b$alleleB <- getAlleleB(file)
    close(file)
    
    chr_files <- rbind(chr_files, b)
  }
  if(nrow(bim)==nrow(chr_files)){
    print("chromosome "%&%  i %&% " bim files contains the same number of SNPs at the .gds file.")
  }else{
    print("chromosome "%&%  i %&% " .gds files contain " %&% nrow(chr_files) %&%" SNPs.")
    print("chromosome "%&%  i %&% " .BIM files contain " %&% nrow(bim) %&%" SNPs.")
    print( (nrow(bim)-nrow(chr_files)) %&% " SNPs may have been lost in the conversion.")
  }
  pgs_chr <- pgs_total %>% filter(chr_name == i)
  for (PGS in PGS_list){
    pgs_chr_name <- pgs_chr %>% filter(pgs == PGS)
    missing <- pgs_chr_name %>% filter(!chr_pos %in% bim$V4)
    print(PGS %&% ": " %&% nrow(missing) %&% " of " %&% nrow(pgs_chr_name) %&% " SNPs (" %&% (nrow(missing)/nrow(pgs_chr_name)) %&% ") on chromosome " %&% i %&% " are missing from the BIM genotype file.")
    if(nrow(missing)==0){}
    else{
      fwrite(missing, out.dir %&% PGS %&% "_snps_missing_from_BIM_chr" %&% i %&% ".txt", col.names= T, sep = "\t" , quote = F, row.names=F)
    }
    missing1 <- pgs_chr_name %>% filter(!chr_pos %in% chr_files$position)
    print(PGS %&% ": " %&% nrow(missing) %&% " of " %&% nrow(pgs_chr_name) %&% " SNPs (" %&% (nrow(missing)/nrow(pgs_chr_name)) %&% ") on chromosome " %&% i %&% " are missing from the .gds genotype file.")
    if(nrow(missing)==0){}
    else{
      fwrite(missing, out.dir %&% PGS %&% "_snps_missing_from_GDS_chr" %&% i %&% ".txt", col.names= T, sep = "\t" , quote = F, row.names=F)
    }
  }
}