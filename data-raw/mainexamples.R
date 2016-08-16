library(OncoPhase)
library(limSolve)

mymode="PhasedSNP"

# Example 1: Loading a simple example data set with two somatic mutations, 5 germlines SNP, and 3 tumor samples
data(simpleExample2)
se=simpleExample2
prevalence_df=getPrevalenceMultiSamples(se$snp_allelecount_df, se$ref_allelecount_df,  se$major_copynumber_df,se$minor_copynumber_df,phasing_association_df=se$phasing_association_df, mode=mymode)
print(prevalence_df)



cat("\n\n\n")

#Chrom     End IsGermline  Tumour1        Tumour2        Tumour3
#mutation2  chr2 3003000          0 C2:0|0|1 C2:0.15|0|0.85 C2:0.12|0|0.88
#mutation6  chr2 4008000          0 C1:1|0|0       C1:1|0|0 C2:0|0.24|0.76

#stop(1)



# Example 2: Running a case study as illustrated in the accompanying paper. Available case studies: A, B, C, 1, 2, . . ., 9
data(CaseStudy_6)
cs=CaseStudy_6
prevalence_CaseStudy6=getPrevalenceMultiSamples(cs$snp_allelecount_df, cs$ref_allelecount_df,  cs$major_copynumber_df,cs$minor_copynumber_df,phasing_association_df=cs$phasing_association_df,mode=mymode)
print(prevalence_CaseStudy6)
#Chrom  End IsGermline          Tumour1
#somaticM  chr3 1000          0 C2:0.25|0.25|0.5

cat("\n\n\n")

data(CaseStudy_A)
cs=CaseStudy_A
prevalence_CaseStudy_A=getPrevalenceMultiSamples(cs$snp_allelecount_df, cs$ref_allelecount_df, phasing_association_df=cs$phasing_association_df, cs$major_copynumber_df,cs$minor_copynumber_df,detail=FALSE,mode=mymode)
print(prevalence_CaseStudy_A)
# 	Chrom  End IsGermline   Tumour1
# somaticM  chr3 1000          0 0.66
cat("\n\n\n")
LongExamples=F
if(LongExamples){
  
  


#Example 3 : Computing somatic mutation cellular prevalence on chromosome 15 of  patient 11152 (data retrieved from a parallel study)

data("chr15_OP1019")
ds=chr15_OP1019
masterprevalence_df=getPrevalenceMultiSamples(ds$snp_allelecount_df, ds$ref_allelecount_df, phasing_association_df=ds$phasing_association_df, ds$major_copynumber_df,ds$minor_copynumber_df,cnv_fraction=ds$cnv_fraction,nbFirstColumns=6,detail=FALSE,mode=mymode)
print(head(masterprevalence_df))
cat("\n\n\n")
data("chr10_OP1019")
df=chr10_OP1019
masterprevalence_df=getPrevalenceMultiSamples(df$snp_allelecount_df, df$ref_allelecount_df, phasing_association_df=df$phasing_association_df, df$major_copynumber_df,df$minor_copynumber_df,cnv_fraction=df$cnv_fraction,nbFirstColumns=6, region="chr10:50000000-180000000",mode=mymode)
print(head(masterprevalence_df))

}


cat("\n\n\n")


# Example 4 : Creating a simple example with one somatic mutation and one germline mutation on a single tumor sample

#Empty dataframe
snp_allelecount_df=as.data.frame(matrix(ncol=4,nrow=2))
names(snp_allelecount_df) = c("Chrom","End","IsGermline","Tumour1")
rownames(snp_allelecount_df) = c("mutation1","mutation2")
ref_allelecount_df = snp_allelecount_df
major_copynumber_df= as.data.frame(matrix(ncol=1,nrow=2))
names(major_copynumber_df) = "Tumour1"
rownames(major_copynumber_df) = c("mutation1","mutation2")
minor_copynumber_df = major_copynumber_df
cnv_fraction = major_copynumber_df

#Filling the dataframes
snp_allelecount_df["mutation1",] = c("chr1", 200100,0,40)
snp_allelecount_df["mutation2",] = c("chr1",  200900,1,60)
ref_allelecount_df["mutation1",] = c("chr1",  200100,0,20)
ref_allelecount_df["mutation2",] = c("chr1",  200900,1,40)
major_copynumber_df["Tumour1"] = c(1,1)
minor_copynumber_df["Tumour1"] = c(1,1)
cnv_fraction["Tumour1"] = c(0.2,0.2)

#Phasing association
phasing_association_df = as.data.frame(matrix(ncol=1,nrow=1))
colnames(phasing_association_df) = c("PhasedGermline")
rownames(phasing_association_df) = c("mutation1")
phasing_association_df["mutation1","PhasedMutations"] = "mutation2"

#Computing the prevalence
prevalence_df=getPrevalenceMultiSamples(snp_allelecount_df, ref_allelecount_df, phasing_association_df=phasing_association_df, major_copynumber_df, minor_copynumber_df,cnv_fraction=cnv_fraction,detail=TRUE,mode=mymode)

print(prevalence_df)

#Chrom    End IsGermline      Tumour1
#mutation1  chr1 200100          0 C2:0|0.5|0.5
