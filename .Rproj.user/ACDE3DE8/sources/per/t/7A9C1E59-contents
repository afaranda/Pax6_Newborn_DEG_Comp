################################################################################
# File: GSE29402_Process.R                                                     #
# Created: May 16, 2019                                                        #
# Author: Adam Faranda                                                         #
# Purpose: Calculate differential expression statistics for GSE13244,          #
#      	                                                                       #
#                                                                              #
################################################################################


# Setup Environment
wd<-getwd()
library(limma)
library(GEOquery)
library(dplyr)
library(reshape2)
library(biomaRt)

# Import Normalized Expression values from GSE13244, clean up Grouping Column
gse<-getGEO(filename='GSE13244_series_matrix.txt.gz', destdir=wd)
exprs(gse)<-log2(exprs(gse))
pData(gse)$characteristics_ch1<-gsub(" ", "_", pData(gse)$characteristics_ch1)
pData(gse)$characteristics_ch1<-gsub("-", "", pData(gse)$characteristics_ch1)

# Process Expression Data, get DEGList
print(apply(exprs(gse), 2, summary))
dm<-model.matrix(~0+characteristics_ch1, pData(gse))
colnames(dm)<-gsub('characteristics_ch1', '', colnames(dm))
cmat<-makeContrasts(`Pax6_heterozygous_RNA` - `Pax6_wildtype_RNA` , levels=dm)
fit<-lmFit(gse, dm)
fit<-contrasts.fit(fit, cmat)
fit<-eBayes(fit)

cols<-c(
	'ID', 'Gene.Symbol', 'ENTREZ_GENE_ID', 'AveExpr',
	'logFC', 'P.Value', 'adj.P.Val'
)
deg<-topTable(fit, n=Inf)[, cols]
print(head(deg))

# Join Group Averages on the degList
Gr1<-pData(gse)[pData(gse)$characteristics_ch1 == 'Pax6_wildtype_RNA', 'geo_accession']
Gr2<-pData(gse)[pData(gse)$characteristics_ch1 == 'Pax6_heterozygous_RNA', 'geo_accession']

ex<-as.data.frame(exprs(gse))
ex$ID<-row.names(ex)
ex<-melt(ex, id.vars='ID')

gr.avg<-as.data.frame(
	inner_join(
		ex %>%
		   filter(variable %in% Gr1) %>%
		   group_by(ID) %>%
		   summarize(Pax6_wildtype_Mean=mean(value)),
		ex %>%
		   filter(variable %in% Gr2) %>%
		   group_by(ID) %>%
		   summarize(Pax6_Heterozygous_Mean=mean(value)),
		by='ID'
	)
)
deg<-as.data.frame(inner_join(deg, gr.avg,by='ID'))

# Choose a representative Probe for each gene
x<-c(1,1,2,0.5, -0.3, 1)
n_up<-function(x){
  tally()
}


