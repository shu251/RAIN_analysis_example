# RAIN analysis
# last updated 07-04-2018

# Import/install R libraries
# install.packages("compositions")
# install.packages("pracma")
# install.packages("rain")
library(compositions)
library(pracma)
library(rain)

# Import test data (derived from OTU table)
load("test_RAIN_df.RData",verbose=T)
head(in_df) # test OTU results. Data frame columns are equal to time points and rows equate to OTUs.
row.names(in_df)<-in_df$OTU.ID; in_df[1]<-NULL # Set OTU.IDs as row names

# Centered log-ratio transformation. This is necessary due to the compositional nature of OTU results (see Gloor et al. 2017).
# ?clr()
df_clr<-as.data.frame(clr(in_df))

# Detrend
# ?detrend()
df_detr<-detrend(t(df_clr))

# RAIN analysis
# ?rain()
# Set up RAIN parameters. This analysis will treat each time point independently (rather than each 'day' as a replicate)
t<-(1:19)
ft <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
# Run RAIN
output_RAIN<-rain(as.matrix(df_detr), period=24, measure.sequence=ft, deltat=4, method="independent", na.rm = TRUE)

# Report results
head(output_RAIN[1:3,]) # Output data frame
hist(output_RAIN$pVal) # distribution of p-value results
sig<-subset(output_RAIN, pVal < 0.05); dim(sig)[1] # Total number of OTUs with significant rhymicity (p<0.05)
sig$OTU.ID<-row.names(sig); list_of_sig_OTUs<-unique(sig$OTU.ID) # list of OTUs with significant rhythmicity

# Assessing significance
# Bonferroni correction
# This is the LEAST forgiving form of multiple testing correction (least power and least type I error)
bonf_sig<-subset(output_RAIN, pVal <= (0.05/nrow(output_RAIN)))

# Benjamini-Hochberg Correction
#?p.adjust()
fdr_ps<-p.adjust(output_RAIN$pVal,method='BH')
bh_sig<-output_RAIN[which(fdr_ps<=0.05),]

# Adaptive Benjamini-Hochberg Correction (see Benjamini and Hochberg 2001)
# Step 1: Sort p-values in ascending order
abh_ps<-output_RAIN$pVal[order(output_RAIN$pVal,decreasing=FALSE)]

# Step 2: Find p-values smaller than those expected at 5% FDR.
# Rejecting all null hypotheses which have p-values less than their corresponding abh_metric will
# yield the same results as p.adjust(method='BH') at the 0.05 significance level.
q<-0.05 #Significance level
m<-length(abh_ps) #Number of hypotheses tested
abh_metric<-q*(1:m)/m

#For a visualization of what this procedure does -- we consider all red points that fall below the black
# line significant
plot(1:m,abh_metric,type='l',
  xlab='p-value rank',
  ylab='p-value',lwd=2,
  main='Benjamini-Hochberg FDR Control')
points(1:m,abh_ps,col='red')
legend(650,0.05,fill=c('red','black'),legend=c('p-value','FDR controlled p threshold'))

# Step 3: Estimate the number of false null hypotheses according to the procedure outlined
# in Benjamini and Hochberg 2001
s=0
s_i=0
i=1
while(s>=s_i){
  s_i=s
  s=(1-abh_ps[i])/(1+m-i)
  i=i+1
}
m_0<-ceiling(1/s)+1

# Step 4: Find p-values smaller than those expected at 5% FDR for adjusted number of false hypotheses
abh_updated<-(1:length(abh_ps))*0.05/m_0
# For a similar visualization
plot(1:m,abh_updated,type='l',
     xlab='p-value rank',
     ylab='p-value',lwd=2,
     main='Adaptive Benjamini-Hochberg FDR Control')
points(1:m,abh_ps,col='red')
legend(650,0.05,fill=c('red','black'),legend=c('p-value','FDR controlled p threshold'))

# Step 5: Determine which tests resulted in a p-value significant at the 5% level considering FDR
abh_sigs<-output_RAIN[which(output_RAIN$pVal %in% abh_ps[which(abh_ps<=abh_updated)]),]
