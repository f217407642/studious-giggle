# studious-giggle
I am working on my R scripting skills ~
Included in this repository is my R code for a group project I worked on in a bioinformatics class.
It is broken down into 3 parts.
I share it with you all so you may offer suggestions for making my code more efficient.
Be constructive, please and thank you !

# PART 1

Aim: Clean up a patient sample annotation csv file by:
1. Reading in a csv file
2. Removing characters within each field (column) that precede and include a semicolon
3. Creating a new file that includes the revised data frame
4. Removing rows within the data frame with NA's in specified columns


# load the package
library(readr)
# read the file
GSEa <- read.csv("GSE4115_clinical_annotation.csv",head=TRUE)
# look at the beginning contents
head(GSEa)
# get the dimensions of the data frame
dim(GSEa) # 192 rows, 11 cols

# parse through each column in GSEa data frame 
# use regular expression 
df_age <- gsea_age_parsed <- gsub(".*:","",GSEa$AGE)
df_gender <- gsea_age_parsed <- gsub(".*:","",GSEa$GENDER)
df_smoking_status <- gsea_age_parsed <- gsub(".*:","",GSEa$SMOKING_STATUS)
df_pack_years <- gsea_age_parsed <- gsub(".*:","",GSEa$PACK_YEARS)
df_diag_bronch <- gsea_age_parsed <- gsub(".*:","",GSEa$DIAG_BRONCH)
df_mass_size <- gsea_age_parsed <- gsub(".*:","",GSEa$MASS_SIZE)
df_biomarker_score <- gsea_age_parsed <- gsub(".*:","",GSEa$BIOMARKER_SCORE)
df_pretest_probability <- gsea_age_parsed <- gsub(".*:","",GSEa$PRETEST_PROBABILITY)
df_cancer_status <- gsub(".*:","",GSEa$CANCER_STATUS)
df_geo_accession <- gsub(".*:","",GSEa$GEO_ACCESSION)

# combine all columns into one data frame
df_combo <-cbind(df_age,df_gender,df_smoking_status,df_pack_years,df_diag_bronch,df_mass_size,df_biomarker_score,df_pretest_probability,df_cancer_status,df_geo_accession)

# write the combined columns to a file and edit columns names manually!
write.csv(df_combo,"gsea4115_gsea_split.csv",row.names=F,quote=FALSE)

# read the file
gsea_split <- read.csv("gsea4115_gsea_split.csv",head=TRUE) 
# remove the rows with NA in any of the fields:
gsea_split_cleaned <- na.omit(gsea_split,cols=c("CANCER_STATUS","PACK_YEARS","SMOKING_STATUS"))

# write it to a file!
write.csv(gsea_split_cleaned,"gsea4115_gsea_split_cleaned.csv",row.names=F,quote=FALSE)

# read it!
gsplitcleaned <- read.csv("gsea4115_gsea_split_cleaned.csv",head=TRUE)

# PART 2
Aim: Pre-process microarray expression data (22,215 probes) for 192 samples by:
1. Calculating the z-score of each probe per sample
2. Finding the absolute value of each probe's z-score per sample
3. Calculating the mean probe z-score of each sample
4. Writing the data frame to a file

# Read the file
zf <- read.csv("GSE4115_probeset_expression.csv",header=TRUE) # read the file
# Convert the data frame to a matrix!
zm <- data.matrix(zf) # converting a data frame to a numeric matrix
# Calculate the z-score of each value per field 
zm_zscore <- scale(zm,center=TRUE,scale=TRUE)# find z score of each point per vector plus deets on std and mean of column 
# Convert all z-score values to + values
zm_zscore_posit <- abs(zm_zscore)

# Conver the matrix to a data frame
zf_zscore_posit <- data.frame(zm_zscore_posit)
# Write the df to a file 
write.csv(zf_zscore_posit,"zf_zscore_posit.csv",col.names=NA,row.names=TRUE) 
# read the file
zf_zscore_posit <- read.csv("zf_zscore_posit.csv",header=TRUE) 
# get dimensions of file
dim(zf_zscore_posit) 
# Create the df to a matrix 
zm_zscore_posit <- data.matrix(zf_zscore_posit) 
# Find the mean of each column
zn_zscore_posit_mean <- colMeans(zm_zscore_posit,na.rm=TRUE) 

#  Convert the numeric to a df
zf_zscore_posit_mean <- data.frame(zn_zscore_posit_mean) # BEST 
# Write the df to a file 
write.csv(zf_zscore_posit_mean,"zf_zscore_posit_mean.csv",col.names=NA,row.names=TRUE)

# Read the file
c <- read.csv("zf_zscore_posit_mean.csv",header=TRUE)
# Load the tidyverse package
library(tidyverse)
# Arrange the row values of a field by values ascending
c_ordered <- c %>% arrange(zn_zscore_posit_mean)
# Write the df to a file
write.csv(c_ordered,"zf_zscore_posit_mean_ordered.csv",col.names=NA,row.names=TRUE)
# Read the file
c <- read.csv("zf_zscore_posit_mean_ordered.csv",header=TRUE)

# PART 3
Aim: Create a histogram of mean probeset z-score values of all samples

# Read the file
zf <- read.csv("GSE4115_probeset_expression.csv",header=TRUE) 
# Check the class of the var
class(zf)
# Convert the df to a matrix
zm <- data.matrix(zf, rownames.force = NA) 
# Find the z-score of each probe per sample
zm_zscore <- scale(zm,center=TRUE,scale=TRUE)
# convert to dataframe
zf_zscore<-data.frame(zm_zscore) 
# convert all zscores to + values
zf_zscore_posit <- abs(zf_zscore)
# Find the mean z-score per sample 
zn_zscore_posit_mean <- colMeans(zf_zscore_posit)

# histogram of the mean z-score per sample 
low<- 0.792655197257325 #mean(zm_zscore[,c("GSM94118")])

hist(zn_zscore_posit_mean,
     breaks=100,
     xlim=c(.778,.805),
     ylim = c(0,25),
     xlab="Z-score",
     ylab="Frequency",
     main="")
abline(v = low,
       col = "red",
       lwd = 2)
