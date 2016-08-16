#getPrevalence<-function(lambda_S,mu_S,major_cn,minor_cn, lambda_G=NULL, mu_G=NULL,  detail=0, mode="PhasedSNP",Trace=FALSE )

#mode="FlankingSNP"
#mode="PhasedSNP"
#Trace=T
lambda_S=lambda_somatic
mu_S=mu_somatic
lambda_G=lambda_LinkedGermline
mu_G=mu_LinkedGermline
  {
  
  N=length(lambda_S) # Number of samples
  
  if((length(mu_S)!=N) || 
     (!is.null(lambda_G) && ((length(lambda_G) !=N) || (length(mu_G) !=N))) ||
     (length(major_cn) !=N) || (length(minor_cn)!=N) )
    stop("\n\nThe vectors passed as input should have the same size\n\n")
  
  
  
  if(mode=="PhasedSNP"){
    prev_somatic=getPhasedSNPPrevalence( lambda_S,mu_S , major_cn,minor_cn, lambda_G , mu_G,detail,Trace=Trace )
  }else if(mode=="FlankingSNP"){
    prev_somatic=getFlankingSNPPrevalence( lambda_S,mu_S ,  major_cn,minor_cn,lambda_G , mu_G, detail,Trace=Trace)
  }else if(mode=="SNVOnly"){
    prev_somatic=getSNVOnlyPrevalence(lambda_S,mu_S ,major_cn,minor_cn, detail, Trace=Trace )
  }else {
    stop("parameter mode should be either FlankingCNP, PhasedSNP or SNVOnly")
  }
  
  
  
  
  
  prev_somatic
  }

print(prev_somatic)
