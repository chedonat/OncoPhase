

#' @export
#getFlankingSNPPrevalence<-function( lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G, detail=FALSE,Trace=False,SameTumour=T )
{
  
  #For each case, we compute the prevalence twice. By considering the somatic to be phased to the germline SNP or phased with the alternative chromosome harboring the reference of the Germline. The latter is achieved just by 
  
  prevalence_phasedSNP = getPhasedSNPPrevalence(lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G, detail=1,Trace=0)
  prevalence_phasedREF= getPhasedSNPPrevalence(lambda_S,mu_S,major_cn,minor_cn, mu_G,lambda_G, detail=1, Trace=0)
  
  if(Trace){
    cat("\n\n Prevalence case phased with SNP \n ")
    print( prevalence_phasedSNP)
    cat("\n Prevalence case phased with REF \n ")
    print( prevalence_phasedREF)
  }
  
  #Initialisation
  prevalence_flankingSNP=prevalence_phasedSNP
  
  
  # we now have choose the prevalence having the highest residual, if Same tumour we sum the residual across all the sample and we compare, else We go sample by sample.
  #if both have the same residual, we check into sample prevalence if they have the same prevalence, in that case the value of that prevalence will be assigned. Otherwise NA will be assigned.
  
  
  if(SameTumour)
  {
    #We sum the residuals and we choose the prevalence with the lower residuals.
    SNP_residual=sum(as.numeric(unlist(lapply(prevalence_phasedSNP, function(x)  if(!is.na(x[["ResidualNorm"]])) as.numeric(x[["ResidualNorm"]]) else NA ))) ,na.rm=T)
    REF_residual=sum(as.numeric(unlist(lapply( prevalence_phasedREF, function(x)  if(!is.na(prevalence_phasedREF[1])) as.numeric(x[["ResidualNorm"]]) else NA ))),na.rm=T)
    
    if(SNP_residual<1e-8)
      SNP_residual=0
    if(REF_residual<1e-8)
      REF_residual=0
    
   if(SNP_residual<REF_residual){
      prevalence_flankingSNP=prevalence_phasedSNP
    }else if(SNP_residual>REF_residual) {
      prevalence_flankingSNP=prevalence_phasedREF
    } else{
      for(sample in names(prevalence_flankingSNP))
      if (prevalence_phasedSNP[[sample]]$Prevalence != prevalence_phasedREF[[sample]]$Prevalence)
        prevalence_flankingSNP[[sample]] =NA 
    }
    
    
  }else{
    for(sample in names(prevalence_flankingSNP))
    {
      #First any residual less than 1e-8 is considered to be 0.
      
      if(prevalence_phasedSNP[[sample]]$ResidualNorm <1e-8)
        prevalence_phasedSNP[[sample]]$ResidualNorm=0
      if(prevalence_phasedREF[[sample]]$ResidualNorm <1e-8)
        prevalence_phasedREF[[sample]]$ResidualNorm=0  
      
      if(prevalence_phasedSNP[[sample]]$ResidualNorm < prevalence_phasedREF[[sample]]$ResidualNorm){
        prevalence_flankingSNP[[sample]] = prevalence_phasedSNP[[sample]]  
      } else if (prevalence_phasedSNP[[sample]]$ResidualNorm > prevalence_phasedREF[[sample]]$ResidualNorm) {
        prevalence_flankingSNP[[sample]] = prevalence_phasedREF[[sample]]  
      }else if (prevalence_phasedSNP[[sample]]$Prevalence != prevalence_phasedREF[[sample]]$Prevalence){
        prevalence_flankingSNP[[sample]] =NA
        
        warning("The prevalence is not resolved without the knowledge of the Phased Germline at this sample")
      }
    }
  }
  
  
  
  
  
  if(detail==1){
    prev_S= prevalence_flankingSNP
  }else{
    prev_S =vector("numeric", length=length(prevalence_flankingSNP))
    names(prev_S) = names(prevalence_flankingSNP)
    prev_S[prev_S==0]<-NA 
    for(sample in names(prevalence_flankingSNP))
    {
      if(is.na(prevalence_flankingSNP[[sample]]))
        next
      
      if(detail==2){
        prev_S[sample] = prevalence_flankingSNP[[sample]]$CondensedPrevalence             
      }else{
        prev_S[sample] = prevalence_flankingSNP[[sample]]$Prevalence               
      }
    }
  }
  
  if(Trace){
    cat("\n Final prevalence \n ")
    print( prev_S)
  }
  
  prev_S
}




