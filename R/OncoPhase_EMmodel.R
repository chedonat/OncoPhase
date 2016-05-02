

complete_likelihood<-function(lambda_S, lambda_G,theta, C, alpha,beta)
{
  
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  N= length(lambda_S)
  m=nrow(lambda_G)
  prob=0
  for (i in 1:N)
  {
    
    
    if(!is.na(hatlambda_S[i] ) && !is.na(hatlambda_G[i])&& hatlambda_G[i]>0 )
    {
      # cat (" here PC ")
      #Constrained part
      if(C=="C1")
      {445
        if(is.na(alpha[i])) next
        if(hatlambda_S[i] <= alpha[i]  * hatlambda_G[i])
        {
          if (prob==0)
          {
            prob = dpois2(lambda_S[i], hatlambda_S[i])
          }else
          {
            prob = prob * dpois2(lambda_S[i], hatlambda_S[i])
          }
          
          for (j in 1:m)
          {
            if(is.na(lambda_G[j,i]) ) next
            prob = prob * dpois2(lambda_G[j,i],hatlambda_G[i])
          }
          
          
        }else
        {
          prob = prob * 0
        }
      }else#(C=="C2")
      {
        if (is.na(beta[i] )) next
        if((hatlambda_S[i] >= beta[i]  * hatlambda_G[i]) && (hatlambda_S[i] <=   hatlambda_G[i]))
        {if (prob==0)
        {
          prob = dpois2(lambda_S[i], hatlambda_S[i])
        }else
        {
          prob = prob * dpois2(lambda_S[i], hatlambda_S[i])
        } 
          for (j in 1:m)
          {
            if(is.na(lambda_G[j,i]) ) next
            prob = prob * dpois2(lambda_G[j,i],hatlambda_G[i])
          }
          
          
        }else
        {
          prob = prob * 0
        }
      }
      
    }
    #  cat(prob)
  }
  
  prob
  
  
}



complete_loglikelihood<-function(lambda_S, lambda_G,theta, C, alpha,beta)
{
  
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  N= length(lambda_S)
  #m=nrow(lambda_G)
  prob=-Inf
  for (i in 1:N)
  {
    if(!is.na(hatlambda_S[i] ) && !is.na(hatlambda_G[i] ) && hatlambda_G[i]>0 )
    {
      #Constrained part
      if(C=="C1")
      {
        if(is.na(alpha[i] )) next
        
        if(hatlambda_S[i] <= alpha[i]  * hatlambda_G[i])
        {
          
          if (prob==-Inf)
          {
            prob =  dpois2(lambda_S[i], hatlambda_S[i], log=T)
          }else{
            prob = prob + dpois2(lambda_S[i], hatlambda_S[i], log=T)            
          }
          
          
          #for (j in 1:m)
          #{
          #if(is.na(lambda_G[j,i]) ) next
          #prob = prob + dpois2(lambda_G[j,i],hatlambda_G[i], log=T)
          #}
          
          if(is.na(lambda_G[i]) ) next
          prob = prob + dpois2(lambda_G[i],hatlambda_G[i], log=T)
          
          
          
        }else
        {
          prob = prob  + -Inf 
        }
      }else
      {
        if(is.na(beta[i] )) next
        if((hatlambda_S[i] >= beta[i]  * hatlambda_G[i]) && (hatlambda_S[i] <=   hatlambda_G[i]))
        {
          if (prob==-Inf)
          {
            prob =  dpois2(lambda_S[i], hatlambda_S[i], log=T)
          }else{
            prob = prob + dpois2(lambda_S[i], hatlambda_S[i], log=T)
          }
          #for (j in 1:m)
          # {
          if(is.na(lambda_G[i]) ) next
          prob = prob + dpois2(lambda_G[i],hatlambda_G[i], log=T)
          # }
          
          
        }else
        {
          prob = prob + -Inf
        }
        
      }
      
      
    }
  }
  
  prob
  
  
}

incomplete_likelihood<-function(lambda_S, lambda_G, theta)
{
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  
  N= length(lambda_S)
  m=nrow(lambda_G)
  prob=1
  for (i in 1:N)
  {
    if(is.na(lambda_S[i]) ) next
    #if(sum(!is.na(lambda_G[,i]) )==0) next
    if(is.na(lambda_G[i]) ) next
    
    prob = prob * dpois2(lambda_S[i], hatlambda_S[i])
    #for (j in 1:m)
    # {
    if(is.na(lambda_G[i]) ) next
    prob = prob * dpois2(lambda_G[i],hatlambda_G[i])
    # }
    
  }
  
  prob
}

incomplete_loglikelihood<-function(lambda_S, lambda_G, theta)
{
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  
  N= length(lambda_S)
  m=nrow(lambda_G)
  prob=0
  for (i in 1:N)
  {
    if(is.na(lambda_S[i]) ) next
    if(is.na(lambda_G[i]) ) next
    #if(sum(!is.na(lambda_G[,i]) )==0) next
    
    prob = prob + dpois2(lambda_S[i], hatlambda_S[i])
    #for (j in 1:m)
    # {
    if(is.na(lambda_G[i]) ) next
    prob = prob + dpois2(lambda_G[i],hatlambda_G[i])
    # }
    
  }
  
  prob
  
  
}




Qfunction<-function(lambda_S, lambda_G, theta, theta_m, alpha,beta)
{
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  
  N= length(lambda_S)
  # m=nrow(lambda_G)
  logprob=0
  PC=1/2 # Prior
  PC1= complete_loglikelihood(lambda_S, lambda_G, theta_m,"C1", alpha,beta)/ PC
  PC2= complete_loglikelihood(lambda_S, lambda_G, theta_m,"C2", alpha,beta)/ PC
  gamma=1/(PC1 + PC2)
  
  
  #C1_term
  prob =0
  for (i in 1:N)
  {
    
    if(!is.na(hatlambda_S[i] ) && !is.na(hatlambda_G[i]) && !is.na(alpha[i])  && hatlambda_G[i] > 0 )
    {
      
      if(hatlambda_S[i] <= alpha[i]  * hatlambda_G[i])
      {
        prob = prob + dpois2(lambda_S[i], hatlambda_S[i], log=T)
        #for (j in 1:m)
        #  {
        if(is.na(lambda_G[i]) ) next
        prob = prob + dpois2(lambda_G[i],hatlambda_G[i], log=T)
        #  }
        
        #print(dpois2(lambda_S[i], hatlambda_S[i], log=T))
        if (dpois2(lambda_S[i], hatlambda_S[i], log=T)==-Inf) stop(paste("lambda_S[i] ", lambda_S[i], "hatlambda_S[i]", hatlambda_S[i]), "hatlambda_G[i]", hatlambda_G[i])
        
      }else
      {
        prob = prob + -Inf
      }
    }
  }
  
  C1_logprob=prob
  if (PC1==-Inf)
    C1_logprob=0
  
  #C2_term
  prob =0
  for (i in 1:N)
  {
    
    if(!is.na(hatlambda_S[i] ) && !is.na(hatlambda_G[i]) && !is.na(beta[i])  && hatlambda_G[i] > 0)
    {
      if((hatlambda_S[i] >= beta[i]  * hatlambda_G[i]) && (hatlambda_S[i] <=   hatlambda_G[i]))
      {
        prob = prob + dpois2(lambda_S[i], hatlambda_S[i], log=T)
        
        #for (j in 1:m)
        #{
        if(is.na(lambda_G[i]) ) next
        prob = prob + dpois2(lambda_G[i],hatlambda_G[i], log=T)
        #}
        
      }else
      {
        prob = prob + -Inf 
      }
    }
  }
  
  C2_logprob=prob
  if (PC2==-Inf)
    C2_logprob=0
  
  
  Qvalue =C1_logprob +    C2_logprob
  Qvalue
  
}



bestAllele<-function(lambda_S,lambda_G, alpha,beta){
  
  
  traceoutput=0
  
  N=length(lambda_S)
  #M=nrow(lambda_G)
  
  #if alpha[i] = beta[i], we add a slight increase on beta[i]
  beta[beta==alpha] = beta[beta==alpha]* 1.01
  
  if(N>0)
  {
    
    # Expectation Maximization algorithm
    
    theta<-list()
    theta[[1]] = list(hatlambda_S=lambda_S,hatlambda_G= lambda_G)
    lik0=1
    lik1= incomplete_loglikelihood(lambda_S,lambda_G,theta[[1]])
    
    lik<-list()
    i=0
    while(lik1 -lik0 !=0 )
    {
      i=i+1
      lik[i]=lik1
      lik0=lik1
      #We find theta_C1
      hatlambda_S=c()
      hatlambda_G=c()
      for (isample in 1:N)
      {
        hatlambda_S[isample] = lambda_S[isample]
        #hatlambda_G[isample] = as.numeric(colMeans(lambda_G,na.rm=T)[isample])
        hatlambda_G[isample] = lambda_G[isample]
        if(is.na(hatlambda_S[isample] ) || is.na(hatlambda_G[isample])|| is.na(alpha[isample] )) next
        #Out of boundary cases
        if (hatlambda_S[isample]>alpha[isample] *hatlambda_G[isample])
          hatlambda_S[isample]=alpha[isample] * hatlambda_G[isample]
      }
      theta_C1=list(hatlambda_S=hatlambda_S,hatlambda_G=hatlambda_G )
      
      
      #We find theta_C2
      hatlambda_S=c()
      hatlambda_G=c()
      for (isample in 1:N)
      {
        hatlambda_S[isample] = lambda_S[isample]
        hatlambda_G[isample] = lambda_G[isample]
        #hatlambda_G[isample] = as.numeric(colMeans(lambda_G,na.rm=T)[isample])
        if(is.na(hatlambda_S[isample] ) || is.na(hatlambda_G[isample]) || is.na(beta[isample] )) next
        #Out of boundary cases
        if (hatlambda_S[isample]<beta[isample] *hatlambda_G[isample]) {
          hatlambda_S[isample]=beta[isample] *hatlambda_G[isample]
        }else  if(hatlambda_S[isample]>hatlambda_G[isample])
        {
          hatlambda_S[isample]= hatlambda_G[isample]
        }
        
        
        
      }
      
      theta_C2=list(hatlambda_S=hatlambda_S,hatlambda_G=hatlambda_G )
      
      
      QC2= Qfunction(lambda_S,lambda_G,theta_C2, theta[[i]], alpha,beta )
      QC1= Qfunction(lambda_S,lambda_G,theta_C1, theta[[i]],alpha,beta ) 
      
      if (QC2>=QC1)
        theta[[i+1]] =theta_C2
      else
        theta[[i+1]] =theta_C1
      
      
      lik1= incomplete_loglikelihood(lambda_S,lambda_G,theta[[i+1]])
      
      if(traceoutput)
      {
        cat("\n\n\n\n \t\t One step done \n\t\t############")
        
        cat("\n\t theta[i] :\n"); print(theta[[i]]  )
        cat("\n\t lik[i] : \n"); print(lik[[i]]  )
        cat("\n\t theta_C1 : \n"); print(theta_C1  )
        cat("\n\t theta_C2 :  \n"); print(theta_C2  )
        cat("\n\t QC1 : \n"); print( QC1  )
        cat("\n\t QC2 :  \n"); print(QC2  )
        cat("\n\t theta[i+1] : \n"); print(theta[[i+1]]  )
        cat("\n\t Lik[i+1]:  \n"); print(lik1  )
        
        
        
      }
      
      
      #We find theta_C2 
    }
    
    
    
    LikC2=complete_loglikelihood(lambda_S,lambda_G, theta[[i+1]],"C2",alpha, beta)
    LikC1=complete_loglikelihood(lambda_S,lambda_G, theta[[i+1]],"C1",alpha, beta)
    if(LikC2>=LikC1)
    {
      bestC="C2"
    }else
    {
      bestC="C1"
    }
    
    
    if(traceoutput)
    { 
      cat("\n\t LikC2 :",LikC2  )
      cat("\n\t LikC1 : ", LikC1  )
      
    }
    
    
    list(hatlambda_S=theta[[i+1]]$hatlambda_S, hatlambda_G=theta[[i+1]]$hatlambda_G,bestC=bestC)
    
  }else
  {
    if(traceoutput)
      cat("\n\n No input, check your inputs")
    
    list(hatlambda_S=c(), hatlambda_G=c(),bestC="")
    
  }
  
  
}



dpois2<-function(x, lambda, log = FALSE){
  dpois(round(x, digits=0), lambda, log = log)
}










