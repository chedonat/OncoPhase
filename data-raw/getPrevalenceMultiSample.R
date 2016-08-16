library(limSolve)


#data("chr15_OP1019")
#ds=chr15_OP1019
#data=ds


fromdata=F


if(fromdata){
  snp_allelecount_df = data$snp_allelecount_df
  ref_allelecount_df = data$ref_allelecount_df
  major_copynumber_df = data$major_copynumber_df
  minor_copynumber_df = data$minor_copynumber_df
  if(!is.null(data$phasing_association_df))
    phasing_association_df = data$phasing_association_df
  
  if(!is.null(data$cnv_fraction))
    cnv_fraction= data$cnv_fraction
}


mode="PhasedSNP"
#mode="FlankingSNP"
#mode="OptimalSNP"
#mode="SNVOnly"
cnv_fraction=NULL
#phasing_association_df=NULL
NormalcellContamination_df=NULL
tumoursamples=NULL
nbFirstColumns=3
region=NULL
min_cells=2
min_alleles=4
detail=F
LocusRadius = 100000
NoPrevalence.action="Skip"





{
  
  # Extract the somatic mutations 
  
  
  compulsory_columns=c("Chrom","End","IsGermline")
  
  if (length(setdiff(compulsory_columns,colnames(snp_allelecount_df)))>0){
    stop(" The allele count master matrices should have at least the following headers
         columns : ")
    print(compulsory_columns)
  }
  
  
  if (length(setdiff(compulsory_columns,colnames(ref_allelecount_df)))>0){
    stop(" The allele count master matrices should have at least the following 
         headers columns : ")
    print(compulsory_columns)
  }
  
  if (is.null(tumoursamples)){
    tumoursamples = colnames(snp_allelecount_df[(nbFirstColumns+1):ncol(snp_allelecount_df)])
  }
  
  #print(colnames(snp_allelecount_df))
  
  
  tumoursamples = Reduce(intersect,list(tumoursamples,colnames(snp_allelecount_df),
                                        colnames(ref_allelecount_df),
                                        colnames(major_copynumber_df), 
                                        colnames(minor_copynumber_df)
  ))
  
  if(!is.null(cnv_fraction))
    tumoursamples =intersect(tumoursamples,colnames(cnv_fraction) )
  
  
  if(length(tumoursamples) ==0)
  {
    stop(" None of the tumour samples provided is present in the five  matrices :
         snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df, cnv_fraction")
  }
  
  snp_allelecount_df=numeric_column(snp_allelecount_df,tumoursamples)
  ref_allelecount_df=numeric_column(ref_allelecount_df,tumoursamples)  
  major_copynumber_df=numeric_column(major_copynumber_df,tumoursamples)
  minor_copynumber_df=numeric_column(minor_copynumber_df,tumoursamples)
  if(!is.null(cnv_fraction)) cnv_fraction=numeric_column(cnv_fraction,tumoursamples)
  
 
  masterprevalence = getPrevalence_Matrice(snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,mode,cnv_fraction, phasing_association_df,NormalcellContamination_df,tumoursamples,  nbFirstColumns, region,detail,  LocusRadius,NoPrevalence.action)
  
  masterprevalence
  
  }



print(head(masterprevalence))
