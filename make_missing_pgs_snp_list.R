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
  
for (i in 1:22){
  # Read in GDS genotype files one file at a time
  gdsDir = "/projects/rpci/wchs/pnfioric/geno_filter_WCHS_Merged_0.01_0.3_AABC_AMBER_FULL/split_files/main/wchs_chr" %&% i 
  files = list.files(gdsDir, pattern = paste("chr", i, ".*.gds$", sep = ""), full.names = T)
  # Note the number of sub files per chromosome
  print("There are " %&% length(files) %&%" gds subfiles for chromosome " %&% i %&% ".")
  
  
  geno_data<-data.table()
  for (j in 1:length(files)){
    input_file <- gdsDir %&%"/chr" %&% i %&%  "."%&% j %&%".gds"
    file<-GdsGenotypeReader(input_file)
    
    b <- data.table()
    b$chromosome <- getChromosome(file)
    b$position <- getPosition(file)
    b$alleleA <- getAlleleA(file)
    b$alleleB <- getAlleleB(file)
    
    geno_data <- rbind(geno_data, b)
    close(file)
  }
  print(nrow(geno_data))

  for (PGS in PGS_list){
    spec<-pgs_total %>% filter(pgs == PGS) %>% filter(chr_name == i)
    missing <- spec %>% filter(!chr_pos %in% geno_data$position)
    print(100*nrow(missing)/nrow(spec))
    if(nrow(missing)>0){
    	fwrite(missing, out.dir %&% PGS %&%"_chr"%&%i%&% "_missing_from_genotypes.txt", quote = F, sep = "\t", row.names = F, col.names = T)
	}else{}
  }
}


