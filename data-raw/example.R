#Simple case I : 1 somatic , 1 germline, 1 tumour

nbTumour = 1
chrom="chr1"

snp_allelecount_df=as.data.frame(matrix(ncol=3+nbTumour,nrow=2))
names(snp_allelecount_df) = c("Chrom","End","IsGermline",paste("Tumour",1:nbTumour,sep=""))
rownames(snp_allelecount_df) = c("mutation1","mutation2")
ref_allelecount_df = snp_allelecount_df
major_copynumber_df= snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
minor_copynumber_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
normalfraction_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]

snp_allelecount_df["mutation1",] = c(chrom, 200100,0,40)
snp_allelecount_df["mutation2",] = c(chrom, 200900,1,60)
ref_allelecount_df["mutation1",] = c(chrom, 200100,0,20)
ref_allelecount_df["mutation2",] = c(chrom, 200900,1,40)


major_copynumber_df["Tumour1"] = c(1,1)
minor_copynumber_df["Tumour1"] = c(1,1)
normalfraction_df["Tumour1"] = c(0.2,0.2)

phasing_association_df = as.data.frame(matrix(ncol=1,nrow=1))
colnames(phasing_association_df) = c("PhasedGermline")
rownames(phasing_association_df) = c("mutation1")
phasing_association_df["mutation1","PhasedGermlines"] = "mutation2"


nbFirstColumns = 3
tumoursamples = "Tumour1"
region = chrom

simpleExample1=list(snp_allelecount_df=snp_allelecount_df, ref_allelecount_df=ref_allelecount_df,
                     phasing_association_df=phasing_association_df, major_copynumber_df=major_copynumber_df,
                     minor_copynumber_df=minor_copynumber_df,normalfraction_df=normalfraction_df )

devtools::use_data(simpleExample1)
prevalence_I=getPrevalence(snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df,nbFirstColumns,tumoursamples, region)







#Simple case II : 2 somatic , 5 germline, 3 tumour

nbTumour = 3
nbmutations=7
chrom="chr2"

snp_allelecount_df=as.data.frame(matrix(ncol=3+nbTumour,nrow=nbmutations))
names(snp_allelecount_df) = c("Chrom","End","IsGermline",paste("Tumour",1:nbTumour,sep=""))
rownames(snp_allelecount_df) = paste("mutation",1:nbmutations,sep="")

ref_allelecount_df = snp_allelecount_df
major_copynumber_df= snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
minor_copynumber_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
normalfraction_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]

snp_allelecount_df["mutation1",] = c(chrom, 3001000,1,40,60,80)
snp_allelecount_df["mutation2",] = c(chrom, 3003000,0,30,20,30)
snp_allelecount_df["mutation3",] = c(chrom, 3009000,1,42,55,72)
snp_allelecount_df["mutation4",] = c(chrom, 4001000,1,30,20,25)
snp_allelecount_df["mutation5",] = c(chrom, 4003000,1,28,18,28)
snp_allelecount_df["mutation6",] = c(chrom, 4008000,0,10,16,20)
snp_allelecount_df["mutation7",] = c(chrom, 4010000,1,32,21,30)


ref_allelecount_df["mutation1",] = c(chrom, 3001000,1,20,40,30)
ref_allelecount_df["mutation2",] = c(chrom, 3001000,0,10,15,20)
ref_allelecount_df["mutation3",] = c(chrom, 3001000,1,25,20,50)
ref_allelecount_df["mutation4",] = c(chrom, 4001000,1,40,25,15)
ref_allelecount_df["mutation5",] = c(chrom, 4003000,1,18,22,10)
ref_allelecount_df["mutation6",] = c(chrom, 4008000,0,8,4,10)
ref_allelecount_df["mutation7",] = c(chrom, 4010000,1,15,20,25)




major_copynumber_df["Tumour1"] = c(2,2,2,3,3,3,3)
minor_copynumber_df["Tumour1"] = c(1,1,1,1,1,1,1)
normalfraction_df["Tumour1"] = c(0.2,0.2,0.2,0.3,0.3,0.3,0.3)

major_copynumber_df["Tumour2"] = c(2,2,2,3,3,3,3)
minor_copynumber_df["Tumour2"] = c(1,1,1,1,1,1,1)
normalfraction_df["Tumour2"] = c(0.2,0.2,0.2,0.3,0.3,0.3,0.3)


major_copynumber_df["Tumour3"] = c(2,2,2,2,2,2,2)
minor_copynumber_df["Tumour3"] = c(1,1,1,1,1,1,1)
normalfraction_df["Tumour3"] = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2)




phasing_association_df = as.data.frame(matrix(ncol=1,nrow=2))
colnames(phasing_association_df) = c("PhasedGermlines")
rownames(phasing_association_df) = c("mutation2","mutation6")
phasing_association_df["mutation2","PhasedGermlines"] = "mutation1:mutation3"
phasing_association_df["mutation6","PhasedGermlines"] = "mutation4:mutation5:mutation7"


nbFirstColumns = 3
tumoursamples = paste("Tumour",1:nbTumour,sep="")
region = chrom

simpleExample2=list(snp_allelecount_df=snp_allelecount_df, ref_allelecount_df=ref_allelecount_df,
                    phasing_association_df=phasing_association_df, major_copynumber_df=major_copynumber_df,
                    minor_copynumber_df=minor_copynumber_df,normalfraction_df=normalfraction_df )

devtools::use_data(simpleExample2)

prevalence_II=getPrevalence(snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df,nbFirstColumns,tumoursamples, region)





