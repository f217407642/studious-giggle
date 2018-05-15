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

