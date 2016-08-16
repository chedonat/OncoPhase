





#' Build the input data matrices for a case study
#' 
#' This is a generic function to automatically build the five input data frame (snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,CNVfraction_df ) for a case study with one somatic mutation, one germline mutation and one or more tumor sample.
#' 
#' @param lambda_G : A count  or a vector of counts (In the case of  multiple tumor samples) of 
#' alleles supporting the variant sequence of  the Germline SNP
#' @param mu_G : A count or a vector of counts (In the case of  multiple tumor samples)  of 
#' alleles  supporting the reference sequence  of  the Germline SNP
#' @param lambda_S : A count or a vector of counts (In the case of  multiple tumor samples)  of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param mu_S : A count or a vector of counts (In the case of  multiple tumor samples)  of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param cnv_fraction: Estimated fraction (or a vector of fractions if multiple tumor samples) 
#' of cells affected by the CNV (1- normal genotype cell fraction). 
#' @param major_cn: Major copy number (or a vector of copy number if multiple tumor samples)
#' at the locus of the mutation
#' @param minor_cn : Minor copy number (or a vector of copy number if multiple tumor samples)
#'  at the locus of the mutation 
#' @param depthOfCoverage : Coverage depth (or a vector of depth coverage  if multiple tumor samples) 
#' at the locus of the mutation. If not provided the exact value of the counts passed as parameters 
#' are considered. If provided then a binomial sampling with replacement is performed to generate the 
#' counts. For the germline, the sampling is done with the parameters p=lambda_G / (lambda_G + mu_G) 
#' and N= depthOfCoverage and will yield  the count of allele supporting the variant sequence of the 
#' germline and the count of allele supporting the reference.  The same sampling is  apply to the 
#' somatic mutation  with the parameters : p=lambda_S/(lambda_S + mu_S) and N = depthOfCoverage.
#'  
#' @return A list containing the following data frames:
#'  \describe{
#'        \item{snp_allelecount_df}{A data frame containing the  count of alleles at each tumour samples supporting the variant sequence at the somatic and germline  mutations. Chrom is set to chr3, Position of the germline and somatic mutations are respectively set to 100 and 1000}
#'        \item{ref_allelecount_df}{A data frame containing the  count of alleles  at each tumour sample  supporting the reference sequence at the somatic and germline  mutation. Chrom is set to chr3, Position of the germline and somatic mutations are respectively set to 100 and 1000}
#'        \item{phasing_association_df}{A data frame containing the phasing association between the somatic and the germline mutation}
#'         \item{major_copynumber_df}{A data frame containing the major copy number  at each tumour sample  at the mutation locus}
#'         \item{minor_copynumber_df}{A data frame containing the minor copy number  at each tumour sample  at the mutation locus}
#'        \item{normalfraction_df}{A data frame containing the proportion of cells with a normal genotype  at each tumour sample. Present only if the method is "PhasedSNPgeneral" }
#'      }
#'      
#' @examples
#' 
#' #Example 1
#' # We reproduce here the case study No 6 of the paper
#'  #Build the input data 
#'  cs = build_casestudy(lambda_G=16, mu_G=8,lambda_S=14,mu_S=10,cnv_fraction=4/8,major_cn=3,minor_cn=1  )
#'  #Run the case
#' prevalence=getPrevalenceMultiSamples(cs$snp_allelecount_df, cs$ref_allelecount_df,cs$major_copynumber_df,cs$minor_copynumber_df, 
#' phasing_association_df=cs$phasing_association_df,cnv_fraction=cs$cnv_fraction)
#' #print the result
#' print(prevalence)
#' #Chrom  End IsGermline Tumour1
#' #somaticM  chr3 1000          0    0.75
#' 
#' 
#' #Example 2
#' #Multiple tumours and stochastic generation of the counts.
#' CaseStudy_10 = build_casestudy(lambda_G=c(8,12,10), mu_G=c(5,4,8),lambda_S=c(3,6,10),mu_S=c(10,8,12),
#' major_cn=c(2,2,3),minor_cn=c(1,1,2) , depthOfCoverage = c(60,100,200))
#' 
#' cs = CaseStudy_10
#' prevalence=getPrevalenceMultiSamples(cs$snp_allelecount_df, cs$ref_allelecount_df,cs$major_copynumber_df,cs$minor_copynumber_df, 
#' phasing_association_df=cs$phasing_association_df,cnv_fraction=cs$cnv_fraction)
#' #print the result
#' print(prevalence)
#' 
#' @export
build_casestudy<-function(lambda_G, mu_G,lambda_S, mu_S, major_cn, minor_cn,  cnv_fraction= NULL, depthOfCoverage=NULL)
{
  
  chrom="chr3"
  N=length(lambda_G) # Number of samples
  if((length(mu_G)!=N) || (length(lambda_S) !=N) || (length(mu_S) !=N) ||
     (length(major_cn) !=N) || (length(minor_cn)!=N) ||
     (!is.null(cnv_fraction) && length(cnv_fraction) !=N))
    stop("\n\nThe vectors passed as input should have the same size\n\n")
  
  if(!is.null(depthOfCoverage))
    if(length(depthOfCoverage)!=N)
      stop("\n\n Vector of coverage should have same size as input vector")
  
  
  
  nbTumour = N
  snp_allelecount_df=as.data.frame(matrix(ncol=3+nbTumour,nrow=2))
  names(snp_allelecount_df) = c("Chrom","End","IsGermline",paste("Tumour",1:nbTumour,sep=""))
  rownames(snp_allelecount_df) = c("germlineM","somaticM")
  ref_allelecount_df = snp_allelecount_df
  major_copynumber_df= snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
  minor_copynumber_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
  cnv_fraction_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
  
  lambda_G_prime=vector("numeric",nbTumour)
  lambda_S_prime=vector("numeric",nbTumour)
  mu_G_prime =vector("numeric",nbTumour)
  mu_S_prime =vector("numeric",nbTumour)
  
  
  
  if(!is.null(depthOfCoverage))
  {
    
    p_G= lambda_G/(lambda_G + mu_G)
    p_S= lambda_S/(lambda_S + mu_S)
    
    
    for(isample in 1:nbTumour){
      lambda_G_prime[isample] = sum(sample(c(1,0), depthOfCoverage[isample], replace = TRUE,prob=c(p_G[isample],1-p_G[isample])))
      mu_G_prime[isample] = depthOfCoverage[isample] - lambda_G_prime[isample] 
      lambda_S_prime[isample] = sum(sample(c(1,0), depthOfCoverage[isample], replace = TRUE,prob=c(p_S[isample],1-p_S[isample])))
      mu_S_prime[isample] = depthOfCoverage[isample] - lambda_S_prime[isample]
      
    }
    
  }else{
    lambda_G_prime=lambda_G
    lambda_S_prime=lambda_S
    mu_G_prime =mu_G
    mu_S_prime =mu_S
  }
  
  
  snp_allelecount_df["germlineM",] = c(chrom, 100,1,lambda_G_prime)
  snp_allelecount_df["somaticM",] = c(chrom, 1000,0,lambda_S_prime)
  ref_allelecount_df["germlineM",] = c(chrom, 100,1,mu_G_prime)
  ref_allelecount_df["somaticM",] = c(chrom, 1000,0,mu_S_prime)
  
  
  major_copynumber_df["germlineM",paste("Tumour",1:nbTumour,sep="")] = major_cn
  minor_copynumber_df["germlineM",paste("Tumour",1:nbTumour,sep="")] = minor_cn
  if(!is.null(cnv_fraction)) cnv_fraction_df["germlineM",paste("Tumour",1:nbTumour,sep="")] = cnv_fraction
  
  major_copynumber_df["somaticM",paste("Tumour",1:nbTumour,sep="")] = major_cn
  minor_copynumber_df["somaticM",paste("Tumour",1:nbTumour,sep="")] = minor_cn
  if(!is.null(cnv_fraction))cnv_fraction_df["somaticM",paste("Tumour",1:nbTumour,sep="")] = cnv_fraction
  
  
  
  phasing_association_df = as.data.frame(matrix(ncol=1,nrow=1))
  colnames(phasing_association_df) = c("PhasedMutations")
  rownames(phasing_association_df) = c("somaticM")
  phasing_association_df["somaticM","PhasedMutations"] = "germlineM"
  
  #   stop(1)
  
  if(!is.null(cnv_fraction)){
    cs=list(snp_allelecount_df=snp_allelecount_df, ref_allelecount_df=ref_allelecount_df,
            phasing_association_df=phasing_association_df, major_copynumber_df=major_copynumber_df,
            minor_copynumber_df=minor_copynumber_df, cnv_fraction=cnv_fraction_df)
  }else {
    cs=list(snp_allelecount_df=snp_allelecount_df, ref_allelecount_df=ref_allelecount_df,
            phasing_association_df=phasing_association_df, major_copynumber_df=major_copynumber_df,
            minor_copynumber_df=minor_copynumber_df)
  }
  
  
  
  cs
}








add_depthcoverage<-function(lambda_S, mu_S, major_cn, minor_cn,lambda_G, mu_G, depthOfCoverage)
{
  
  chrom="chr3"
  N=length(lambda_G) # Number of samples
  if((length(mu_G)!=N) || (length(lambda_S) !=N) || (length(mu_S) !=N) ||
     (length(major_cn) !=N) || (length(minor_cn)!=N) )
    stop("\n\nThe vectors passed as input should have the same size\n\n")
  
  if(length(depthOfCoverage)!=N)
    stop("\n\n Vector of coverage should have same size as input vector")
  
  
  
  nbTumour = N
  
  lambda_G_prime=vector("numeric",nbTumour)
  lambda_S_prime=vector("numeric",nbTumour)
  mu_G_prime =vector("numeric",nbTumour)
  mu_S_prime =vector("numeric",nbTumour)
  
  
  
  if(!is.null(depthOfCoverage))
  {
    
    p_G= lambda_G/(lambda_G + mu_G)
    p_S= lambda_S/(lambda_S + mu_S)
    
    
    for(isample in 1:nbTumour){
      lambda_G_prime[isample] = sum(sample(c(1,0), depthOfCoverage[isample], replace = TRUE,prob=c(p_G[isample],1-p_G[isample])))
      mu_G_prime[isample] = depthOfCoverage[isample] - lambda_G_prime[isample] 
      lambda_S_prime[isample] = sum(sample(c(1,0), depthOfCoverage[isample], replace = TRUE,prob=c(p_S[isample],1-p_S[isample])))
      mu_S_prime[isample] = depthOfCoverage[isample] - lambda_S_prime[isample]
      
    }
    
  }else{
    lambda_G_prime=lambda_G
    lambda_S_prime=lambda_S
    mu_G_prime =mu_G
    mu_S_prime =mu_S
  }
  
  
  list(lambda_S_prime,mu_S_prime,major_cn, minor_cn,lambda_G_prime,mu_G_prime)
  
}









