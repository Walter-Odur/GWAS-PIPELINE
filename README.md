<h1 align = "center">
    GWAS PIPELINE - Rmarkdown
</h1>

<h3 align = 'center'>
    Author: Walter Odur <br>
    Affiliation: ACE-Uganda
</h3>

This pipeline is divided in 10 steps
(1) Reading data into R to create an R object; 
(2) SNP-level filtering (part 1); 
(3) Sample-level filtering; 
(4) SNP-level filtering (part 2); 
(5) Principal component analysis (PCA); 
(6) Imputation of non-typed genotypes; 
(7) Association analysis of typed SNPs; 
(8) Association analysis of imputed data; 
(9) Integration of imputed and typed SNP results; and
(10) Visualization and quality control of association findings.
### Loading Libraries
```{r}
library(snpStats)
library(LDheatmap)
library(devtools)
library(plyr)
library(doParallel)
library(coin)
library(igraph)
library(downloader)
library(SNPRelate)
library(GenABEL.data)
library(tibble)
library(GenABEL)
```
## Data pre-processing
.ped: The .ped file contains information on each study participant including family ID, participant ID, father ID, mother ID, sex, phenotype, and the full typed genotype.
.map: The .map file contains a row for each SNP with rsNumber (SNP) and corresponding chromosome (chr) and coordinate (BPPos) based on the current genome build.
.bim: The .bim file contains a row for each SNP and six columns, containing information for the chromosome number, rsNumber, genetic distance,position identifier, allele 1, and allele 2.
.bed: The .bed file contains a binary version of the genotype data. This is the largest of the three files(.bim, .bed & .fam) because it contains every SNP in the study, as well as the genotype at this SNP for each individual.
.fam: The .fam file contains the participant identification information, including a row for each individual and six columns, corresponding the same columns described for the .ped file with the exception of the genotype data.
Clinical data file: An additional ascii .txt or .csv file is typically available, which includes clinical data on each study subject. The rows of this file represent each subject, and the columns correspond to available covariates and phenotypes.
### 1(a). Loading data into R
```{r}
#-------step1---------
#setwd("C:/Users/DELL/Downloads/GWAS_ASSIGNMENT_2024_GROUP1")
data.dir <- "C:/Users/DELL/Downloads/GWAS_ASSIGNMENT_2024_GROUP1/"
out.dir <- "C:/Users/DELL/Downloads/GWAS_ASSIGNMENT_2024_GROUP1/"
gwas.fn <- lapply(c(bed='bed', bim='bim',fam='fam',gds='gds'), function(x) sprintf("%s/GWAStutorial.%s", data.dir, x))
clinical.fn <- sprintf("%s/GWAStutorial_clinical.csv", data.dir)
onethou.fn <- lapply(c(info='info',ped='ped'), function(x) sprintf("%s/chr16_1000g_CEU.%s", data.dir, x))
protein.coding.coords.fname <- sprintf("%s/ProCodgene_coords.csv", out.dir)

# Output files
gwaa.fname <- sprintf("%s/GWAStutorialout.txt", out.dir)
gwas.unadj.fname <- sprintf("%s/GWAStutorialoutUnadj.txt", out.dir)
impute.out.fname <- sprintf("%s/GWAStutorial_imputationOut.csv", out.dir)
CEPT.fname <- sprintf("%s/CEPT_GWASout.csv", out.dir)
```


### (b) Reading data into r
```{r}
#-------step1-------
# Creating a list with thge data files
geno <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam)

# Obtaining a snpmatrix for the data
gwas_genotype <- geno$genotypes
print(gwas_genotype)

# Obtaining the snp information table
genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
head(genoBim)

```
![image](https://hackmd.io/_uploads/SyCE62gGA.png)


### Or the non replicative way
```{r}
# Setting working directories
setwd("C:/Users/DELL/Downloads/GWAS_ASSIGNMENT_2024_GROUP1/")
getwd()
list.files()

# Reading data into R
geno <- read.plink("GWAStutorial.bed", "GWAStutorial.bim", "GWAStutorial.fam", na.strings = ("-9"))
head(geno)

# Creating a genotype object from the data list
gwas_genotype <- geno$genotypes
print(gwas_genotype)

# SNP information table
gwas_bim <- geno$map
colnames(gwas_bim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
names(gwas_bim)
head(gwas_bim)
```

### c). Clinical data
```{r}
#------Step1------
# Clinical data
setwd("C:/Users/DELL/Downloads/GWAS_ASSIGNMENT_2024_GROUP1/")
clinical <- read.csv("GWAStutorial_clinical.csv", colClasses = c("character", "factor", "factor", "numeric", "numeric", "numeric", "numeric"))
#class(clinical)
#head(clinical)
rownames(clinical) <- clinical$FamID
head(clinical)

# Subsetting data with only individuals who have both genotype and phenotype data
gwas_genotype_filtered_1 <- gwas_genotype[clinical$FamID, ]
print(gwas_genotype_filtered_1)

# Freeing up space
rm(geno)

```
![image](https://hackmd.io/_uploads/BJZ_pheMC.png)
![image](https://hackmd.io/_uploads/rkOs6nlGC.png)


### (2). SNP-level filtering (part 1)
SNP-level filtering based on a large amount of missing data and lower variability is performed first.
Call rate: The call rate for a given SNP is defined as the proportion of individuals in the study for which the corresponding SNP information is not missing.
Minor Allele Frequency (MAF): A large degree of homogeneity at a given SNP across study participants generally results in inadequate power to infer a statistically significant relationship between the SNP and the trait under study.
Considered MAF of 1% and calling rate of 95% (meaning we retain SNPs for which there is less than 5% missing data)
```{r}

#---------step2------
gwas_snpsum_col <- col.summary(gwas_genotype)
print(head(gwas_snpsum_col))

gwas_genotype_filtered_1 <- with(gwas_snpsum_col, (!is.na(MAF) & MAF > "0.01") & Call.rate >= "0.95")
head(gwas_genotype_filtered_1)
gwas_genotype_filtered_1[is.na(gwas_genotype_filtered_1)] <- FALSE
head(gwas_genotype_filtered_1)

cat(ncol(gwas_genotype)-sum(gwas_genotype_filtered_1), "SNPs will be removed due to low MAF or call rate. \n")

# Subsetting genotype and SNP summarry that pass the filters
gwas_genotype <- gwas_genotype[, gwas_genotype_filtered_1]
gwas_snpsum_col <- gwas_snpsum_col[gwas_genotype_filtered_1,]

print(gwas_genotype)

```
 ![image](https://hackmd.io/_uploads/rJ5W0hgzA.png)
   
    ***Report the number of SNPs removed due to low MAF or call rate and those that pass call rate and MAF criteria.***
    203287 SNPs will be removed due to low MAF or call rate 
    658186 SNPs pass call rate and MAF criteria.
    


### (3)(a). Sample-level filtering

Criteria for sample-level filtering are generally based on missing data, sample contamination, correlation (for population-based investigations), and racial, ethnic, or gender ambiguity or discordance.
Heterozygosity: Heterozygosity refers to the presence of each of the two alleles at a given SNP within an individual. This is expected under HWE to occur with probability 2∗p∗(1 − p), where p is the dominant allele frequency at that SNP (assuming a bi-allelic SNP). Excess heterozygosity across typed SNPs within an individual may be an indication of poor sample quality, while deficient heterozygosity can indicate inbreeding or other substructure in that person.
Individuals who are missing genotype data for more than 5% of the typed SNPs are removed.

```{r}
#-------step3------

# Creating sample statistics (call rate, heterozygosity)
gwas_snpsum_row <- row.summary(gwas_genotype)
head(gwas_snpsum_row)

# Adding the f-stat (inbreeding coefficient) to gwas_snpsum_row
MAF <- gwas_snpsum_col$MAF
memory.limit(size = 9999999999)
callmatrix <- !is.na(gwas_genotype)
heterozygosity_expected <- callmatrix %*% (2*MAF*(1-MAF))
heterozygosity_observed <- with(gwas_snpsum_row, Heterozygosity*(ncol(gwas_genotype))*Call.rate)
gwas_snpsum_row$heterozygosity_final <- 1-(heterozygosity_observed/heterozygosity_expected)
head(gwas_snpsum_row)

# Given the sample call rate cut-off is 0.95 and inbreeding cutoff being 0.1
sampleuse <- with(gwas_snpsum_row, !is.na(Call.rate) & Call.rate > 0.95 & abs(heterozygosity_final) <= 0.1)
sampleuse[is.na(sampleuse)] <- FALSE
cat(nrow(gwas_genotype) - sum(sampleuse), "subjects will be removed due to low sample call rate or inbreeding coefficient.\n")
# Zero were removed


# Subsetting genotypes and clinical data for subjects who pass call rate and heterozygosity criteria
gwas_genotype <- gwas_genotype[sampleuse,]
clinical <- clinical[rownames(gwas_genotype),]
```
![image](https://hackmd.io/_uploads/BJHtA3efA.png)

    How many subjects were removed due to low sample call rate or inbreeding coefficient?
    No samples were removed at this step 

### (b). Checking for relatedness through Identity By Descent (IBD) and Principal Component Analysis (PCA)
A common measure of relatedness (or duplication) between pairs of samples is based on identity by descent (IBD). An IBD kinship coefficient of greater than 0.10 may suggest relatedness, duplicates, or sample mixture.
Typically, the individual of a related pair with lower genotype call rate is removed. We note that gender identity can also be checked at this stage to confirm that self-reported gender is consistent with the observed X and Y chromosomes
```{r}
# Given LD cutoff of 0.2 and kinship cutoff of 0.1
# Creating gds file, required for snprelate functions
snpgdsBED2GDS(gwas.fn$bed, gwas.fn$fam, gwas.fn$bim, gwas.fn$gds)
genofile <- snpgdsOpen(gwas.fn$gds, readonly = FALSE)
#genofile
#snpgdsClose(genofile)

# Automatically added "-1" to be removed
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
#View(gds.id)
gds.id <- sub("-1","",gds.ids)
# View(gds.id)
add.gdsn(genofile, "sample.id", gds.id, replace = TRUE)
```
![image](https://hackmd.io/_uploads/HyyYk6xGC.png)
 
    

### 4). SNP-level filtering (part 2) 
### (a). Cleaning for LD
Applying linkage disequilibrium (LD) pruning using a threshold value of 0.2, eliminates a large degree of redundancy in the data and reduces the influence of chromosomal artifacts.
An alternative to this approach is the ‘HapMap rooted’ analysis, which involves first performing PCA in a reference panel, for example, HapMap or 1000 Genomes, and then projecting the study sample onto the resulting space
```{r}

# Prunning snps for IBD analysis
set.seed(1000)
geno.sample.ids <- rownames(gwas_genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = 0.2, sample.id = geno.sample.ids, snp.id = colnames(gwas_genotype))

snpset.ibd <- unlist(snpSUB, use.names = FALSE)
cat(length(snpset.ibd), "will be used in IBD analysis\n") # Expect 72812 snps

```
![image](https://hackmd.io/_uploads/By0w4pefR.png)

How many SNPs will be used in the IBD analysis?
    
    72812 SNPs  
    
### (b). Cleaning for IBD
```{r}
IBD <- snpgdsIBDMoM(genofile, kinship = TRUE, sample.id = geno.sample.ids, snp.id = snpset.ibd, num.thread = 4)
IBDcoeff <- snpgdsIBDSelection(IBD)
head(IBDcoeff)

# Checking for any candidates that are related
IBDcoeff <- IBDcoeff[IBDcoeff$kinship >= 0.1, ]

# Using a while loop to remove sample with a high pairing
related.samples <- NULL
while (nrow(IBDcoeff) > 0) {
  # Counting the number of  occurrences of each and taking the top one
  sample.counts <- arrange(count(c(IBDcoeff$ID1, IBDcoeff$ID2)), -freq)
  rm.sample <- sample.counts[1, 'x']
  cat("Removing sample", as.character(rm.sample), 'too closely related to', sample.counts[1, 'freq'], 'other samples.\n')
  
  # Removing samples from ibdcoeff and to list
  IBDcoeff <- IBDcoeff[IBDcoeff$ID1 != rm.sample & IBDcoeff$ID2 != rm.sample, ]
  related.samples <- c(as.character(rm.sample), related.samples)
}

# Filter genotype and clinical to include only unrelated samples
gwas_genotype <- gwas_genotype[!(rownames(gwas_genotype) %in% related.samples), ]
clinical <- clinical[!(clinical$FamID %in% related.samples), ]

geno.sample.ids <- rownames(gwas_genotype)

cat(length(related.samples), "similar samples removed due to correlation coeffecient >= 0.1")
print(gwas_genotype) # Expect all 1401 subjects to remain
```
![image](https://hackmd.io/_uploads/ryJcNalM0.png)

    Report the number of similar samples removed (if any) removed due to correlation coefficient >= 0.1. How many subjects remain? Do we have to remove any samples? If yes, how many and why and if not, why not?
    None were removed because all the samples were from different individuals given the correlation coeffecient.
    

### (5)(a). Principal component analysis (PCA)
Ancestry PCA is one approach to visualizing and classifying individuals into ancestry groups based on their observed genetic makeup. 
We do this for two reasons: 
First, self-reported race and ethnicity can differ from clusters of individuals that are based solely on genetic
information.
Second, the presence of an individual not appearing to fall within a racial/ethnic cluster may be suggestive of a sample-level error.
```{r}
# Find PCA matrix
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids, snp.id = snpset.ibd, num.thread = 4)

# Create a dataframe of first two principle components
pctab <- data.frame(sample.id = pca$sample.id, PC1 = pca$eigenvect[,1], PC2 = pca$eigenvect[,2], stringsAsFactors = FALSE)

# Plot the first two principle components
plot(pctab$PC2, pctab$PC1, xlab = "Principle component 2", ylab = "Principle component 1", main = "Ancestry Plot")
```
![image](https://hackmd.io/_uploads/SJP3BpeGR.png)
![image](https://hackmd.io/_uploads/HyIaHTgGR.png)



### (b). Checking for HWE
Violations of HWE can be an indication of the presence of population substructure or the occurrence of a genotyping error.
While they are not always distinguishable, it is a common practice to assume a genotyping error and remove SNPs for which HWE is violated. If case-control status is available, we limit this filtering to analysis of controls as a violation in cases may be an indication of association. Departures from HWE are generally measured at a given SNP using a CHI-Square (𝜒2) goodness-of-fit test between the observed and expected genotypes. 
We remove SNPs for which the HWE test statistic has a corresponding p-value of less than 1 × 10−6 in controls.

```{r}
hardy <- 10^-6

CADcontrols <- clinical[clinical$CAD == 0, 'FamID']
snpsum.colCont <- col.summary(gwas_genotype[CADcontrols, ])
HWE_use <- with(snpsum.colCont, !is.na(z.HWE) & (abs(z.HWE) < abs(qnorm(hardy/2))))
rm(snpsum.colCont)

HWE_use[is.na(HWE_use)] <- FALSE
cat(ncol(gwas_genotype) - sum(HWE_use), "SNPs will be removed due to high HWE. \n") # 1296 SNPS removed

# Subset genotypes and SNP summary data for SNPs that pass HWE criteria
gwas_genotype <- gwas_genotype[, HWE_use]
print(gwas_genotype) # 656890 SNPs remain
```
![image](https://hackmd.io/_uploads/H1dBFpezA.png)

Why is it important to perform a check for HWE?
    
    Detecting Genotyping Errors: HWE deviations can often indicate genotyping errors within a dataset. These errors can be introduced during various stages of the genotyping process, leading to inaccurate results. This helps detect potential errors and improve the data quality.

Reasons for using a lenient cut-off (e.g., 1x10^-6):

    Multiple Testing Correction: Performing a statistical test on a large number of SNPs necessitates adjusting the significance threshold to account for multiple testing. A stricter cut-off would likely result in discarding a substantial number of valid SNPs due to chance alone.
    Focus on Larger Deviations: While some deviations from HWE might occur due to random fluctuations, a lenient cut-off allows us to focus on SNPs with substantial HWE deviations, which are more likely to indicate true genotyping errors or other underlying population genetic factors.

Why HWE is tested only on CAD controls:

    Association Studies with Disease: Case groups (individuals with CAD) might exhibit deviations from HWE due to the disease itself or selection bias. Therefore, focusing on control groups, which are not expected to be influenced by the disease, provides a better assessment of genotyping quality and potential technical issues.

How many SNPs were removed due to high HWE?

    1296 SNPS
Report the number of SNPs that remain?

    656890 SNPs

### (6). Imputation of non-typed genotypes (New data generation)
The first are PCs that are intended to capture information of latent population substructure that is typically not available in self-reported race and ethnicity variables. Substructure, also referred to as population admixture and population stratification, refers to the presence of genetic diversity (e.g., different allele frequencies) within an apparently homogenous population that is due to population genetic history (e.g., migration, selection, and/or ethnic integration)
The second are genotypes of untyped SNPs that may have a functional relationship to the outcome and therefore provide additional power for identifying association.
### (a). Creating pcs for capturing population substructure
```{r}
# Setting LD
geno.sample.ids <- rownames(gwas_genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = 0.2, sample.id = geno.sample.ids, snp.id = colnames(gwas_genotype))
snpset.pca <- unlist(snpSUB, use.names = FALSE)
cat(length(snpset.pca))

pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids, snp.id = snpset.pca, num.thread = 4)
# Find and record first 10 principle components
# PCs will be a N:10 matrix each column with a pc
pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1:10], stringsAsFactors = FALSE)
colnames(pcs)[2:11] <- paste("pc", 1:10, sep = "")
head(pcs)

```
![image](https://hackmd.io/_uploads/HydFYalMA.png)
![image](https://hackmd.io/_uploads/Hy_ct6gzC.png)

What is the importance of performing PCA at this step?

    To reduce the dimensionality by identifying a smaller number of new, uncorrelated variables (principal components) that capture the majority of the variance in the original data.
    
How many SNPs will be used in the IBD analysis?

    72812 SNPs
### (b). Imputing non-typed SNPs using 1000 genome data
Typed data is one generated from the sequencing machines (chip array technology)
Analyzing association of genotypes of non-typed SNPs with disease outcomes because functional (causal) variants may not be measured.
Some of the stand alone packages/software that can be used for this analysis include - BEAGLE, IMPUTE2, MACH

```{r}
# Read in 1000genome data for chromosome 16
thougeno <- read.pedfile(onethou.fn$ped, snps = onethou.fn$info, which=1)

# Obtain genotype data for given chromosome
genoMatrix <- thougeno$genotypes

# Obtain the chromosome position for each SNP
gwas_support <- thougeno$map
colnames(gwas_support) <- c("SNP", "Position", "A1", "A2")
head(gwas_support)

# Imputation of non-typed snps
presSNPs <- colnames(gwas_genotype)

# Subset for SNPs on given chromosome
presDatChr <- genoBim[genoBim$SNP %in% presSNPs & genoBim$chr == 16, ]
targetSNPs <- presDatChr$SNP
#targetSNPs
# Subset 1000 genome data for our snps
# 'missing' and 'present' are snpMatrix objects needed for imputation rules 
is.present <- colnames(genoMatrix) %in% targetSNPs

missing <- genoMatrix[, !is.present]
#print(missing)

present <- genoMatrix[, is.present]
#print(present)

# Obtain position of snps to be used for imputatio rules
pos.present <- gwas_support$Position[is.present]
#pos.present
pos.miss <- gwas_support$Position[!is.present]
#pos.miss
# Calculate and store imputation rules
rules <- snp.imputation(present, missing, pos.present, pos.miss)
#head(rules)
```


```{r}
# Remove failed imputations
rules <- rules[can.impute(rules)]
cat("Imputation rules for", length(rules), "SNPs were estimated\n") # 197888 snps 

# Quality control for imputation certainity and MAF
rules <- rules[imputation.r2(rules) >= 0.7]
cat(length(rules), "Imputation rules remain after uncertain impuation were removed \n") # 162565 

rules <- rules[imputation.maf(rules) >= 0.01]
cat(length(rules), "Imputation rules remain after MAF filtering \n") # 162565

# Obtain posterior expectation of genotypes of imputed snps
target <- gwas_genotype[, targetSNPs]
imputed <- impute.snps(rules,target, as.numeric = FALSE)
print(imputed) # 162565
```
![image](https://hackmd.io/_uploads/ryNd6alfR.png)

Briefly describe these two projects (HapMap and 1000 Genomes).
HapMap Project:

    Goal: To create a catalog of common human genetic variations (Single Nucleotide Polymorphisms - SNPs) across diverse populations.
    Methodology: Sequenced the genomes of 270 individuals from four populations (African, European, Japanese, and Chinese).
    Impact: Revolutionized human genetic studies by providing a reference database for identifying associations between SNPs and disease risk.
            

1000 Genomes Project:

    Goal: To create a comprehensive map of human genetic variation, including both common and rare SNPs, structural variations, and other forms of genetic diversity.
    Methodology: Sequenced the genomes of over 2,500 individuals from 26 populations worldwide, utilizing advanced sequencing technologies.
    Impact: Provided a deeper understanding of human genetic diversity and enabled the discovery of new variants associated with complex traits and diseases.

How many SNPs are missing? 

    400000 SNPs
How many SNPs were typed? 

    656890 SNPs   
Report the number of SNPs tagged by a single SNP and those tagged by multiple tag haplotypes. 
    
    SNPs tagged by a single SNP: 82119
    SNPs tagged by multiple tag haplotypes (saturated model): 115769
Imputation rules for how many SNPs were estimated? 

    197888 SNPs
Use 0.7 for the R 2 threshold and 0.01 for the MAF. How many imputation rules remain after imputations with low certainty were removed? 

    162565 SNPs
Imputation rules remain after MAF filtering? 

    162565 SNPs
How many SNPs were imputed?

    162565 SNPs
## Genome-wide association analysis
### (7). Association analysis of typed snps
Association analysis typically involves regressing each SNP separately on a given trait, adjusted for patient-level clinical, demographic, and environmental factors.
Each SNP is represented as the corresponding number of minor alleles (0, 1, or 2).
A Bonferonni-corrected genome-wide significance threshold of 5 × 10−8 is used for control of the family-wise error rate. This cutoff is based on research, suggesting approximately one-million independent SNPs across the genome, so tends be applied regardless of the actual number of typed or imputed SNPs under investigation.

```{r}
# Merge clinical data and principal components to create phenotype table
phenoSub <- merge(clinical, pcs)

# Performing a rank-based inverse normal transformation of hdl
phenoSub$phenotype <- rntransform(phenoSub$hdl, family="gaussian")

# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(phenoSub$hdl, main="Histogram of HDL", xlab="HDL")
hist(phenoSub$phenotype, main="Histogram of Transformed HDL", xlab="Tranformed HDL")

# Removed unnecessary columns from table
phenoSub$hdl <- NULL
phenoSub$ldl <- NULL
phenoSub$tg <- NULL
phenoSub$CAD <- NULL

# Rename columns to match names necessary for GWAS() function
phenoSub <- rename(phenoSub, replace = c(FamID="id"))

# Include only subjects with hdl data
phenoSub <- phenoSub[!is.na(phenoSub$phenotype), ] # 1309
head(phenoSub)

```
![image](https://hackmd.io/_uploads/S1EIIRlzC.png)
![image](https://hackmd.io/_uploads/HJIvICefC.png)

```{r}
source("GWAA.R")

start <- Sys.time()

GWAA(genodata = gwas_genotype, phenodata = phenoSub, filename = gwaa.fname)
end <- Sys.time()
print(end - start)

```
![Screenshot 2024-04-30 194215](https://hackmd.io/_uploads/HJiY1AxMC.png)

How many subjects were included with phenotype data? 

    1309 subjects
Using this phenotype data, perform model fitting on each of the typed SNPs. How much time did it take to perform the parallel model fitting? 

    2 hours and 15 minutes (2.226162 hours)
State the platform used (e.g., windows), number of cores and RAM on your computer.

    Windows, 4 cores and 8GB RAM

### (8). Association analysis of imputed data
```{r}
# Carryout association testing for imputed SNPs 
rownames(phenoSub) <- phenoSub$id

imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, family = "Gaussian", data = phenoSub, snp.data = target, rules = rules)

# Obtaining p - values for imputed SNPS by calling methods on the returned GlmTests object
results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results <- results[!is.na(results$p.value), ]

# Write a file containing the results
write.csv(results, impute.out.fname, row.names = FALSE)

# Merge imputation testing results with support to obtain coordinates
imputeOUt <- merge(results, gwas_support[, c("SNP", "Position")])
imputeOUt$chr <- 16

imputeOUt$type <- "imputed"

# Find the -log_10 of the p_values
imputeOUt$Neg_logP <- -log10(imputeOUt$p.value)

# Order by p-value
imputeOUt <- arrange(imputeOUt, p.value)
head(imputeOUt)
```
![7](https://hackmd.io/_uploads/BkW_j0lG0.png)


```{r}
map2gene <- function(gene, coords, SNPs, extend.boundary = 5000){
  coordsSub <- coords[coords$gene == gene, ] #Subset coordinate file for specified gene
  
  coordsSub$start <- coordsSub$start - extend.boundary # Extend gene boundaries 
  coordsSub$stop <- coordsSub$stop + extend.boundary
  
  SNPsub <- SNPs[SNPs$Position >= coordsSub$start & SNPs$Position <= coordsSub$stop & SNPs$chr == coordsSub$chr,] # Subset for SNPs in gene
  return(data.frame(SNPsub, gene = gene, stringsAsFactors = FALSE))
}
 
```

```{r}
#source("map2gene.R")
protein.coding.coords.fname <- "ProCodgene_coords.csv"
genes <- read.csv(protein.coding.coords.fname, stringsAsFactors = FALSE)

# Subset for CETP SNPs
impCETP <- map2gene("CETP", coords = genes, SNPs = imputeOUt)

# Filter only the imputed CETP SNP genotypes
impCETPgeno <- imputed[ , impCETP$SNP]
head(impCETPgeno)
```
![2](https://hackmd.io/_uploads/rkQQCAlMA.png)

State the SNP with the lowest p-value in both the typed and imputed SNP analysis. 

    Typed - rs1532625
    Imputed - rs1532624
    
In the boundaries of which gene does the SNP with the lowest p-value for imputed SNP analysis lie? 

    Position 57005479
What about the typed? 

    Position 57005301
    
State the gene(s).

    The cholesteryl ester transfer protein (CETP) gene

## Post-analytic visualisation and genomic interrogation
### (9). Integration of imputed and typed SNP results
The reference used in this study was GRCh37 (hg19)
```{r}
GWASout <- read.table(gwaa.fname, header = T, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))

# Find the -log_10 of the p-value
GWASout$Neg_logP <- -log10(GWASout$p.value)

# Merge output with genoBim by SNP name to add position and chromosome number
GWASout <- merge(GWASout, genoBim[, c("SNP", "chr", "position")])

# Order SNPs by significance
GWASout <- arrange(GWASout, -Neg_logP)
head(GWASout)
```
![4](https://hackmd.io/_uploads/rJEu-y-f0.png)


```{r}
map2gene <- function(gene, coords, SNPs, extend.boundary = 5000){
  coordsSub <- coords[coords$gene == gene, ] #Subset coordinate file for specified gene
  
  coordsSub$start <- coordsSub$start - extend.boundary # Extend gene boundaries 
  coordsSub$stop <- coordsSub$stop + extend.boundary
  
  SNPsub <- SNPs[SNPs$position >= coordsSub$start & SNPs$position <= coordsSub$stop & SNPs$chr == coordsSub$chr,] # Subset for SNPs in gene
  return(data.frame(SNPsub, gene = gene, stringsAsFactors = FALSE))
}
```


```{r}
GWASout$type <- "typed"

GWAScomb <- rbind.fill(GWASout, imputeOUt)
head(GWAScomb)
tail(GWAScomb)

# Subset for CETP snps
typCETP <- map2gene("CETP", coords = genes, SNPs = GWASout)

# Combine CETP snps from imputed and typed analysis
CETP <- rbind.fill(typCETP, impCETP)[, c("SNP", "p.value", "Neg_logP", "chr", "position", "type", "gene")]
print(CETP)
```
![5](https://hackmd.io/_uploads/H164W1ZzR.png)

![6](https://hackmd.io/_uploads/HJsoeyWGR.png)

![7](https://hackmd.io/_uploads/HyM5Wk-fC.png)

## Manhattan plots
Manhattan plots are used to visualize GWA significance level by chromosome location. Here, each dot corresponds to a single SNP. 
The x-axis represents gene coordinates, and the numbers shown correspond to chromosome numbers. 
The y-axis is the negative of the log p-value, so that large values correspond to small p-values.
The solid horizontal line indicates the Bonferonni corrected significance threshold (− log(5 × 10−8)). 
The dotted horizontal line is a less stringent suggestive association threshold (− log(5 × 10−6)) that we use as an indicator of a suggestive association and requiring further validation.
```{r}
#    ---- manhattan ----
# Receives a data.frame of SNPs with Neg_logP, chr, position, and type.
# Plots Manhattan plot with significant SNPs highlighted.
GWAS_Manhattan <- function(GWAS, col.snps=c("black","gray"),
                           col.detected=c ("blue"), col.imputed=c("red"), col.text="black", 
                           title="GWAS Tutorial Manhattan Plot", display.text=TRUE, 
                           bonferroni.alpha=0.05, bonferroni.adjustment=1000000,
                           Lstringent.adjustment=10000){
  bonferroni.thresh <- -log10(bonferroni.alpha / bonferroni.adjustment)
  Lstringent.thresh <- -log10(bonferroni.alpha / Lstringent.adjustment) 
  xscale <- 10000000
  
  manhat <- GWAS[!grepl("[A-z]",GWAS$chr),]
  
  #sort the data by chromosome and then location
  manhat.ord <- manhat[order(as.numeric(manhat$chr),manhat$position),] 
  manhat.ord <- manhat.ord[!is.na(manhat.ord$position),]
  
  ##Finding the maximum position for each chromosome
  max.pos <- sapply(1:21, function(i) { max(manhat.ord$position[manhat.ord$chr==i],0) }) 
  max.pos2 <- c(0, cumsum(max.pos))
  
  #Add spacing between chromosomes
  max.pos2 <- max.pos2 + c(0:21) * xscale * 10
  
  #defining the positions of each snp in the plot
  manhat.ord$pos <- manhat.ord$position + max.pos2[as.numeric(manhat.ord$chr)]
  
  # alternate coloring of chromosomes
  manhat.ord$col <- col.snps[1 + as.numeric(manhat.ord$chr) %% 2]
  
  # draw the chromosome label roughly in the middle of each chromosome band 
  text.pos <- sapply(c(1:22), function(i) { mean(manhat.ord$pos[manhat.ord$chr==i]) })
  
  
  # PLOT THE DATA
  plot(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=='typed'], 
       pch=20, cex=.3, col= manhat.ord$col[manhat.ord$type=="typed"], xlab="Chromosome", 
       ylab="Negative Log P-Value", axes=F, ylim=c(0,max(manhat$Neg_logP)+1))
  
  points(manhat.ord$pos[manhat.ord$type=="imputed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="imputed"], 
         pch=20, cex=.4, col = col.imputed)
  
  points(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"], 
         pch=20, cex=.3, col = manhat.ord$col[manhat.ord$type=="typed"])
  
  axis(2)
  abline(h=0)
  
  SigNifSNPs <- as.character(GWAS[GWAS$Neg_logP > Lstringent.thresh & GWAS$type=="typed", "SNP"]) 
  
  #Add legend
  legend("topright",c("Bonferroni Corrected Threshold*", "Less Stringent Threshold**"), 
         border="black", col=c("gray60", "gray60"), pch=c(0, 0), lwd=c(1,1), 
         lty=c(1,2), pt.cex=c(0,0), bty="o", cex=0.7)
  
  #Add chromosome number
  text(text.pos/xscale, -.3, seq(1,22,by=1), xpd=TRUE, cex=1)
  
  #Add bonferroni line
  abline(h=bonferroni.thresh, untf = FALSE, col = "gray60")
  
  #Add "less stringent" line
  abline(h=Lstringent.thresh, untf = FALSE, col = "gray60", lty = 2 )
  
  #Plotting detected genes 
  #Were any genes detected? 
  if (length(SigNifSNPs)>0){
    sig.snps <- manhat.ord[,'SNP'] %in% SigNifSNPs
    
    points(manhat.ord$pos[sig.snps]/xscale, 
           manhat.ord$Neg_logP[sig.snps], 
           pch=15,col=col.detected, bg=col.detected, cex=0.5)
    
    text(manhat.ord$pos[sig.snps]/xscale, 
         manhat.ord$Neg_logP[sig.snps],
         as.character(manhat.ord[sig.snps,1]), col=col.text, offset=1, adj=-.1, cex=.7)
  }
}


```

```{r}
#source("GWAS_ManhattanFunction.R")

par(mfrow = c(1,2))

GWAS_Manhattan(GWAScomb)

```
![8](https://hackmd.io/_uploads/ByGLfyZf0.png)

Create a Manhattan plot to visualize GWA significant results by chromosome location. Use a “Bonferroni” adjusted significance cut-off of −log10(5×10 −8 ) and a “Candidate” cut-off of −log10(5×10 −6 ). What do you observe?

    None of the SNPs (both typed and imputed) reached the adjusted signficant cut-off.
    
## Quantile - Quantile (Q-Q) plots and the lambda statistic.
Q-Q plots are used to visualize the relationship between the expected and observed distributions of SNP level test statistics. Here we compare these statistics for the unadjusted model (left) compared with the model adjusted for confounders by incorporating the first ten principal components along with clinical covariates. Create QQ plots for adjusted and unadjusted model outputs.
```{r}
# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[, c("id", "phenotype")]

GWAA(genodata = gwas_genotype, phenodata = phenoSub2, filename = gwas.unadj.fname)
GWASoutUndaj <- read.table(gwas.unadj.fname, header = T, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))

# Create QQ plots for adjusted and unadjusted model outputs
par(mfrow= c(1,2))
lambdaAdj <- estlambda(GWASout$t.value^2, plot = TRUE, method = "median")
lambdaUnadj <- estlambda(GWASoutUndaj$t.value^2, plot = TRUE, method = "median")
cat(sprintf("Unadjusted lambda: %s\nAdjusted lambda: %s\n", lambdaUnadj$estimate, lambdaAdj$estimate))

```
![59](https://hackmd.io/_uploads/HkeahgbMC.png)

```{r}
# ------Step10-c--------
# Calculate standardized lambda
lambdaAdj_1000 <- 1+(lambdaAdj$estimate - 1)/nrow(phenoSub)*1000
lambdaUnadj_1000 <- 1+(lambdaUnadj$estimate -1)/nrow(phenoSub)*1000
cat(sprintf("Standardized unadjusted lambda: %s\nStandardized adjusted lambda: %s\n", lambdaUnadj_1000, lambdaAdj_1000))

```
![57](https://hackmd.io/_uploads/rJIlhx-zC.png)


```{r}
#--------Step10-d--------
library(LDheatmap)
library(rtracklayer)
  
# Add nra247617" to CETP
CETP <- rbind.fill(GWASout[GWASout$SNP = "rs247617",], CETP)

# Combine genotypes and imputed genotypes for CETP region
subgen <- cbind(gwas_genotype[.colnames(gwas_genotype) %in% CETP$SNP], impCETPgeno)    # CETP subsets from typed and imputed SNPs

# Subset SNPs for only certain genotypes
certain <- apply(as(subgen, 'numeric'), 2, function(x) { all(x %in% c (0,1,2,NA)) }) 
subgen <- subgen[,certain]

# Subset and order CETP SNPB by position 
CETP <- CETP[CETP$SNP %in% colnamas(subgan),]
CETP <- arrange(CETP, position)
subgen <- subgen[, order(match(colnames(subgen), CETP$SNP)) ]

# Create LDheatmap
ld <- ld(subgen, subgen, stats-"R.squared") # Find LD map of CETP SNPs
ll <- LDheatmap(ld, CETP$position, flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination
llplusgenes <- LDheatmap.addGenes(ll, chr = "chr16", genome = "hg19", genesLocation = 0.01)

# Add plot of -log(p)
library(ggplot2)
plot.new()
llplot2<-LDheatmap.addGrob(llplusgenes, rectGrob(gp = gpar(col = "white")), height = .34) 
pushViewport(viewport(x = 0.483, y= 0.76, width = .91, height = .4))

grid.draw(ggplotGrob({
  qplot(position, Neg_logP, data = CETP, xlab = "", ylab = "Negative Log P-value", xlim = range(CETP$position), 
        asp = 1/10, color = factor(type), colour=c("#000000", "#D55E00")) + 
    theme(axis.text.x = element_blank(),
          axis.title.y = element.text(size = rel(0.75)), legend.position = "none", 
          panel.background = element.blank(), 
          axis.line = element.line(colour = "black") + 
    scale.color.manual(values = c("red", "black"))
)}))


# ---- steplO-e ----
# Create regional association plot

library (postgwas)

# Create data.frame of most significant SNP only 
snps<-data.frame(SNP=c("rs1532625"))
                                    
# Change column names necessary to run regionalplot function
GWAScomb <- rename(GWAScomb, replace=c(p.value="P", chr="CHR", position="BP"))

# Edit biomartConfigs so regionalplot function
# pulls from human genome build 37/hgl9

myconfig <- biomartConfigs$hsapiens 
myconfig$hsapiens$gene$host <- "grch37.ensembl.org" 
myconfig$hsapiens$gene$mart <- "ENSEMBL_MART_ENSEMBL" 
myconfig$hsapiens$snp$host <- "grch37.ensembl.org" 
myconfig$hsapiens$snp$mart <- "ENSEMBL_MART_SNP"

# Run regionalplot using HAPMAP data (pop = CEU)
regionalplot(snps, GWAScomb, biomart.config = myconfig, window.size = 400000, draw.snpname = data.frame(
  snps = c("rs1532625", "rs247617"), 
  text = c("rs1532625", "rs247617"), 
  angle = c(20, 160), 
  length = c(1, 1), 
  cex = c(0.8)
),
ld.options = list(
  gts.source = 2,
  max.snps.per.window = 2000,
  rsquare.min = 0.8,
  show.rsquare.text = FALSE
),
out.format = list(file = "png", panels.per.page = 4))
```
![61](https://hackmd.io/_uploads/rysOhgWfA.png)




















































































































































































