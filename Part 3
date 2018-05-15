
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
