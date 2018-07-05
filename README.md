# RAIN analysis

Below are explanations, test data, and example code for how to treat data ahead of RAIN analysis. Rhythmicity Analysis Incorporating Nonparametric method (RAIN) analysis [Thaben and Westermark 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4266694/) is a robust non-parametric method to detect rhythms in time series data. 

### Contributors:
Daniel Muratore (Georgia Tech) & Sarah K. Hu (Univ Southern California)

## Origin of test data

Data was derived from clustering OTUs at 97% (see Hu et al. 2016 and Hu et al. submitted) - link to github for OTU analysis. Test data provided for below example has been subset to that each OTU has 20 or more sequences. ([Origin of data](https://github.com/shu251/18Sdiversity_diel)). Steps below are applicable to 18S or 16S tag-sequencing results (from HTS methods) or metagenome/metatranscriptome results.

## Explanation of steps in R script (RAINanalysis_wOTUdata.R)

*Centered log-ratio transformation*

Count data from HTS (tag-sequence or meta-genome/transcriptome) is a randomly generated subsample of the relative abundances of the actual molecules targeted for sequencing. This means that absolute abundances cannot be described from this data - this is way we always perform careful conversions to relative abundances, rarefied, or employ other library normalization tools (i.e. edgeR or DEseq). This is the definition of 'compositional data' ([see more](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full)).
While there will continue to be debate on how best to normalize HTS results, it is important that one thinks critically about their data and interprets ecological findings with caution.


One way to deal with compositional data is to perform a centered-log ratio normalization ([see more](https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.4710300705)). This step essentially converts count or relative abundance data into ratios that more clearly represent the relationship between samples in your data. The log-transformation serves to make the data symmetrical and linear. Note that this is still not obtaining total counts, but this normalized data is suitable for most of the statistical tools used downstream in HTS analyses.

```
library(compositions)
df_clr<-as.data.frame(clr(in_df))
```

*Detrend*

Since our goal is to look for significant rhythmicity in our data, we need to detrend ahead of running RAIN analysis. RAIN is rank-based method that uses consecutive Mann-Whitney tests to determine if one sampling time of the specified period tends to have the highest rank over your time series (ie, is there a repeating hump?). If there is a global trend across your entire time series, then later periods will systematically have higher or lower ranks than their counterpart from a period earlier in the series, and that signal obfuscates the signal of the periodicity and substantially reduces RAIN's power.

```
library(pracma)
df_detr<-detrend(t(df_clr))
```

*RAIN analysis*

First set up RAIN parameters. This particular test treats each time point (columns) independently (rather than assigning a day as a 'replicate').
```
# Set up RAIN parameters. This analysis will treat each time point independently (rather than each 'day' as a replicate)
t<-(1:19)
ft <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
```

Use detrended data as input for RAIN.
* period - best match of number of measurements (over a 24 hour period)
* measure.sequence - Number of replicates for each time point (in this case we have none, so all are set to 1)
* deltat - sampling interval (4 hours in this case)
* method - recommended to run with "independent"
* na.rm - treatment of NAs

```
output_RAIN<-rain(as.matrix(df_detr), period=24, measure.sequence=ft, deltat=4, method="independent", na.rm = TRUE)
``` 

*Interpreting results from RAIN analysis*

Output from RAIN will report pvalues, phase, period, and peak shape. Most likely, you will want to determine OTUs (or genes) with significant rhythmicity based on the p-value.
```
sig<-subset(output_RAIN, pVal < 0.05); dim(sig)[1] # Total number of OTUs with significant rhymicity (p<0.05)
```

However, while RAIN's output are corrected p-values, those values are corrected for the multiple testing of determining the significance of different possible phases/shapes, not the multiple testing of running RAIN on many time series. R script highlights a few possible options for multiple testing correction for analyzing large numbers of OTUs. 

*including:*
* **Bonferroni correction**: this is the LEAST forgiving form of multiple testing correction (least power and least type I error)
* **Benjamini-Hochberg Correction**: this approach is still conservative but less than any familywise-error rate correction (more power, but more type I error) - (see Benjamini and Hochberg 1995)
* **Adaptive Benjamini-Hochberg Correction** (see Benjamini and Hochberg 2001): This method is a more powerful version of the above correction (R script demonstrates how to calculate this and provides a visualization

##### Last updated 07-06-2018
