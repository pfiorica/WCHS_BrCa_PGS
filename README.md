# WCHS_BrCa_PGS
Calculating breast cancer PGS for the Women's Circle of Health Study

This is a repository for calculating PGS for the Women's Circle of Health Study (WCHS). Here, we calculate 6 different PGS for breast cancer and its various subtypes. The main document that can be followed for notes on how I calculated PGS is `WCHS_PGS_NOTES.Rmd`. This document provides background on the calculation, PGS Weight cleaning, genotype filtering, and PGS calculation. 

The [troubleshooting](https://github.com/pfiorica/WCHS_BrCa_PGS/tree/main/troubleshooting) folder contains a list of SNPS for each PGS and chromosome that did not match to the genotype. Given that this was checked on both the PLINK.bim files and .gds files, the files in this folder are redundant.

## BMI PGS

In January 2026, 5 PGS for BMI (PGS003844-PGS003848) were calculated from the WCHS genotypes. A few additional scripts were made to perform these calculations. The documentation of this process can be found in `BMI_PGS_Calculation_Notes.Rmd`. It follows the same pipeline as above with a slight change to how the `% Matched` was calculated when intersecting the PGS variants and the genotype variants..

## Contact
For further questions about these PGS and their calculation, please contact me, Peter Fiorica, at peter.fiorica@roswellpark.org.
