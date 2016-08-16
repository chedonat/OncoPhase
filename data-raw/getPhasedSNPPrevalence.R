#getPhasedSNPPrevalence<-function( lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G, detail=0,Trace=FALSE )
{
  #if(length(lambda_G)>1)
  {
    tumoursamples= names(lambda_G)
    if (is.null(tumoursamples)){
      tumoursamples = paste("Sample",c(1:length(lambda_G)),sep="_")
    }
    
    #To avoid some side error, we transform them into vector
    lambda_G=as.vector(lambda_G)
    mu_G=as.vector(mu_G)
    lambda_S=as.vector(lambda_S)
    mu_S = as.vector(mu_S)
    major_cn = as.vector(major_cn)
    minor_cn=as.vector(minor_cn)
    
    names(lambda_G) =  tumoursamples
    names(mu_G) =  tumoursamples
    names(lambda_S) =  tumoursamples
    names(mu_S) =  tumoursamples
    names(major_cn) =  tumoursamples
    names(minor_cn) =  tumoursamples
    
  }
  
  Normalise=T
  if(Normalise){
    Total_S=mu_S +lambda_S
    Total_G=mu_G + lambda_G
    mu_S = (mu_S/Total_S) *Total_G 
    lambda_S = (lambda_S/Total_S) *Total_G 
  }
  
  if(detail==1)
  { prev_S = list()
  }else{
    prev_S =vector("numeric", length=length(lambda_G))
    names(prev_S) = names(lambda_G)
    prev_S[prev_S==0]<-NA  
  }
  
  
  for(sample in tumoursamples)
  {
    if(Trace) cat("\n\n\n Computing the prevalence on sample :", sample)
    args_list=list(
      lambda_S=lambda_S[sample],
      mu_S=mu_S[sample],
      major_cn=major_cn[sample],
      minor_cn=minor_cn[sample],
      lambda_G=lambda_G[sample],
      mu_G=mu_G[sample],#/ omega_G[sample] - lambda_G[sample],
      # NoPrevalence.action=NoPrevalence.action,
      detail=1,
      Trace=Trace)
    
    if(anyNA(args_list))
      next
    if(lambda_G[sample]+mu_G[sample] ==0)
      next
    if(lambda_S[sample]+ mu_S[sample] ==0)
      next
    
    prevalence=do.call(getPhasedSNPPrevalence_on_singlemutation, args_list)
    
    # print(prevalence)
    
    #  if detail, the context and  tree type of prevalence are collapsed else only the prevalence is outputed
    if(detail==2){
      prev_S[sample] = prevalence$CondensedPrevalence
    } else if(detail==1){
      prev_S[sample] =list(prevalence)
    }else{
      prev_S[sample] =prevalence$Prevalence
    }
  }
  
  if(Trace) {cat("\n\n\n Computed prevalences are  :\n")
    print(prev_S)}
  
  
  prev_S
}