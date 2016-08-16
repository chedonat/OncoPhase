#Simple case I : 1 somatic , 1 germline, 1 tumour



SavingData=F

#Case Study 0
#============
cat(" \n \n \t Case STudy 1 \n")
#Build the input data 
cs = build_casestudy(lambda_G=c(8,12,10), mu_G=c(5,4,8),lambda_S=c(3,6,10),mu_S=c(10,8,12),phi_G=c(3/5,4/5,3/5),major_cn=c(2,2,3),minor_cn=c(2,1,2)
                     , depthOfCoverage = c(60,100,200)
                     )
#Save the case
CaseStudy_10=cs
if(SavingData) devtools::use_data(CaseStudy_10, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)

stop(9)


#Case Study 1
#============
cat(" \n \n \t Case STudy 1 \n")
#Build the input data 
cs = build_casestudy(lambda_G=8, mu_G=5,lambda_S=3,mu_S=10,phi_G=3/5,major_cn=2,minor_cn=2  )
#Save the case
CaseStudy_1=cs
if(SavingData) devtools::use_data(CaseStudy_1, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)


#Case Study 2
#============
cat(" \n \n \t Case STudy 2 \n")
#Build the input data 
cs = build_casestudy(lambda_G=20, mu_G=8,lambda_S=3,mu_S=25,phi_G=6/8,major_cn=3,minor_cn=1  )
#Save the case
CaseStudy_2=cs
if(SavingData) devtools::use_data(CaseStudy_2, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)


#Case Study 3
#============
cat(" \n \n \t Case STudy 3 \n")
#Build the input data 
cs = build_casestudy(lambda_G=5, mu_G=3,lambda_S=3,mu_S=5,phi_G=2/5,major_cn=1,minor_cn=0  )
#Save the case
CaseStudy_3=cs
if(SavingData) devtools::use_data(CaseStudy_3, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)



#Case Study 4
#============
cat(" \n \n \t Case STudy 4 \n")
#Build the input data 
cs = build_casestudy(lambda_G=9, mu_G=6,lambda_S=7,mu_S=8,phi_G=3/6,major_cn=2,minor_cn=1  )
#Save the case
CaseStudy_4=cs
if(SavingData) devtools::use_data(CaseStudy_4, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)





#Case Study 5
#============
cat(" \n \n \t Case STudy 5 \n")
#Build the input data 
cs = build_casestudy(lambda_G=4, mu_G=4,lambda_S=2,mu_S=6,phi_G=0,major_cn=1,minor_cn=1  )
#Save the case
CaseStudy_5=cs
if(SavingData) devtools::use_data(CaseStudy_5, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)



#Case Study 6
#============
cat(" \n \n \t Case STudy 6 \n")
#Build the input data 
cs = build_casestudy(lambda_G=16, mu_G=8,lambda_S=14,mu_S=10,phi_G=4/8,major_cn=3,minor_cn=1  )
#Save the case
CaseStudy_6=cs
if(SavingData) devtools::use_data(CaseStudy_6, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)




#Case Study 7
#============
cat(" \n \n \t Case STudy 7 \n")
#Build the input data 
cs = build_casestudy(lambda_G=8, mu_G=0,lambda_S=1,mu_S=7,phi_G=1,major_cn=2,minor_cn=0  )
#Save the case
CaseStudy_7=cs
if(SavingData) devtools::use_data(CaseStudy_7, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)



#Case Study 8
#============
cat(" \n \n \t Case STudy 8 \n")
#Build the input data 
cs = build_casestudy(lambda_G=6, mu_G=2,lambda_S=3,mu_S=5,phi_G=4/6,major_cn=1,minor_cn=0  )
#Save the case
CaseStudy_8=cs
if(SavingData) devtools::use_data(CaseStudy_8, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)







#Case Study 9
#============
cat(" \n \n \t Case STudy 9 \n")
#Build the input data 
cs = build_casestudy(lambda_G=7, mu_G=1,lambda_S=6,mu_S=2,phi_G=3/4,major_cn=2,minor_cn=0  )
#Save the case
CaseStudy_9=cs
if(SavingData) devtools::use_data(CaseStudy_9, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)


#Case Study A
#============
cat(" \n \n \t Case STudy A \n")
#Build the input data 
cs = build_casestudy(lambda_G=8, mu_G=6,lambda_S=6,mu_S=8,phi_G=2/6,major_cn=2,minor_cn=1  )
#Save the case
CaseStudy_A=cs
if(SavingData) devtools::use_data(CaseStudy_A, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)


#Case Study B
#============
cat(" \n \n \t Case STudy B \n")
#Build the input data 
cs = build_casestudy(lambda_G=8, mu_G=4,lambda_S=4,mu_S=8,phi_G=2/6,major_cn=2,minor_cn=0  )
#Save the case
CaseStudy_B=cs
if(SavingData) devtools::use_data(CaseStudy_B, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)



#Case Study C
#============
cat(" \n \n \t Case STudy C \n")    
#Build the input data 
cs = build_casestudy(lambda_G=8, mu_G=12,lambda_S=6,mu_S=12,phi_G=4/8,major_cn=2,minor_cn=1  )
#Save the case
CaseStudy_C=cs
if(SavingData) devtools::use_data(CaseStudy_C, overwrite=TRUE)
#Run the case
prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
                         cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#print the result
print(prevalence)































