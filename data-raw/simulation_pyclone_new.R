
library(OncoPhase)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(reshape2)
jBrewColors <- brewer.pal(n = 8, name = "Set2")

#Load the somatic mutation data
somatic_folder="data-raw/doni2/output/"
germline_folder= "data-raw/doni2/output2/"

depth<-c(30,60,120,1000)

n_replicates=30
n_depth=4


mymode1="SNVOnly"
mymode2="PhasedSNP"





getAllMutationPlot<-function(mymode){
  Pprev<-list()
  Presolv<-list()
 
  
  Compute=T
  if(Compute)
    for(idepth in 1:length(depth)){
      results<-as.data.frame(matrix(nrow=50,ncol=n_depth))
      colnames(results) = c("depth","mutation_id","replicate","prevalence")
      
      cat("\n\n Depth : ", depth[idepth],"")
      for (irep in 1:n_replicates){
        somatic_output<-read.table(paste(somatic_folder,"pyclone_readdepth-",depth[idepth],"_n-1_rep-",irep,".txt",sep=""),header=T)
        germline_output<-read.table(paste(germline_folder,"pyclone_germline_readdepth-",depth[idepth],"_n-1_rep-",irep,".txt",sep=""),header=T)
        germline_output["normalgenotype"] = c(1,1,0.7,1,0.8)
        
        for(imutation in 1:5){
          lambda_G=germline_output[imutation,"var_counts"]
          mu_G=germline_output[imutation,"ref_counts"]
          lambda_S=somatic_output[imutation,"var_counts"]
          mu_S=somatic_output[imutation,"ref_counts"]
          phi_G=1-germline_output[imutation,"normalgenotype"]
          major_cn=germline_output[imutation,"major_cn"]
          minor_cn=germline_output[imutation,"minor_cn"]
          prevalence = getPrevalence(lambda_S,mu_S,major_cn,minor_cn,lambda_G,mu_G,mode=mymode)
          
          myprevalence=as.numeric(prevalence)
          
          cat ( imutation + 5* (irep-1), " " )
          results[imutation + 5* (irep-1),"depth"] = depth[idepth]
          results[imutation + 5* (irep-1),"replicate"] = irep
          results[imutation + 5* (irep-1),"mutation_id"] = somatic_output[imutation,"mutation_id"]
          results[imutation + 5* (irep-1),"prevalence"] = myprevalence
          
        }
      }
      
      write.table(results[2:4],file=paste("data-raw/doni2/prevalence_depth-",depth[idepth],".txt",sep=""),quote=F,sep="\t")
      
      if(idepth==1){
        wholeresults=results
      }else{
        wholeresults=rbind( wholeresults,results)
      } 
    }
  
  write.table( wholeresults,file=paste("data-raw/doni2/prevalence_alldepth.txt",sep=""),quote=F,sep="\t")
  
  #Plotting
  
  cat("\n Plotting")
  true_prevalence<-c(1,0.9,0.4,0.15,0.2)
  n_simulations=n_replicates
  n_d=n_depth
  for(imutation in 1:5){
    result_mutation=wholeresults[wholeresults$mutation_id==imutation,]
    prevalence_estimates = matrix(nrow=n_simulations, ncol=n_d)
    for(idepth in 1:length(depth)){
      for (irep in 1:n_replicates){
        prevalence_estimates[irep, idepth] = wholeresults[wholeresults$mutation_id==imutation &
                                                            wholeresults$depth==depth[idepth] & 
                                                            wholeresults$replicate==irep,"prevalence"]
      }
    }
    
    prevalence0=true_prevalence[imutation]
    
    #print the result
    error_estimates = prevalence_estimates - prevalence0
    mse_prediction=apply(error_estimates,2,function(x) mean(x*x,na.rm=T))
    names(mse_prediction) = depth
    
    names(prevalence_estimates) = depth
    
    #   df=melt(mse_prediction)
    #   df["Coverage"] = as.factor(depth)
    #   df$value = as.numeric(df$value)
    #   p1=ggplot(data=df,aes(x=Coverage,y=value,group = 1)) + geom_line(colour=jBrewColors[4],size=1.1, stat="identity")+
    #     xlab("") + ylab("Mean Square Error (MSE)")+
    #     theme(axis.text=element_text(size=12,face="bold"),  axis.title=element_text(size=12,face="bold"))
    #   
    prevalence_estimates_df=as.data.frame((prevalence_estimates))
    names(prevalence_estimates_df) =depth
    
    if(nrow(na.omit(prevalence_estimates_df))==0){
      prevalence_estimates_df[is.na(prevalence_estimates_df)]<-0
    }

    
    df=melt(prevalence_estimates_df) 
    names(df) = c("Coverage","Prevalence")
    
    means <- aggregate(Prevalence ~  Coverage, df, mean)
    p2=ggplot(data=df, aes(x=Coverage, y=Prevalence)) + geom_boxplot(fill=jBrewColors[7])+
      stat_summary(fun.y=mean, colour="darkred", geom="point",  shape=18, size=3) + 
      geom_hline(yintercept=prevalence0,color=jBrewColors[4],size=1,linetype =2)+ ylim(c(0,1)) + 
      geom_text(data = means, aes(label = format(Prevalence,digits=2), x=Coverage,y = 1.04*Prevalence, size=1),show.legend = FALSE)+
      xlab("Coverage") + ylab("Prevalence")+
      theme(axis.text=element_text(size=12,face="bold"),  axis.title=element_text(size=12,face="bold"))
    
    
    Pprev[[imutation]]<-p2
    
    
    
    resolved_df=as.data.frame(apply(prevalence_estimates_df,2,function(x) sum(!is.na(x))))
    colnames(resolved_df) = c("Resolved")
    resolved_df["Not_Resolved"] = nrow(prevalence_estimates_df) - resolved_df$Resolved
    resolved_df["Coverage"]=rownames(resolved_df)
    resolved_df$Coverage <- factor(resolved_df$Coverage, levels = resolved_df$Coverage)
    meltd<- melt(resolved_df,id.vars=3)
    p1= ggplot(meltd, aes(x=Coverage, y=value, fill=variable)) +
      geom_bar(stat="identity", position = "fill") +
      scale_y_continuous(labels = percent_format())
   
    
    Presolv[[imutation]]<-p1
    
  
    
    
    
    
    
  }
  
  
  
  list(P1=Pprev,P2=Presolv)
}



P0<-getAllMutationPlot("SNVOnly")
P0<-P0$P1
P1<-getAllMutationPlot("PhasedSNP")
P1<-P1$P1
P2<-getAllMutationPlot("FlankingSNP")
P3<-P2$P2
P2<-P2$P1


grid.arrange(P0[[1]],P1[[1]],P0[[2]],P1[[2]],P0[[3]],P1[[3]],P0[[4]],P1[[4]],P0[[5]],P1[[5]],ncol=2, 
             top=paste("SNVOnly mode                                                                                    PhasedSNP mode"))

grid.arrange(P3[[1]],P2[[1]],P1[[1]],P3[[2]],P2[[2]],P1[[2]],P3[[3]],P2[[3]],P1[[3]],P3[[4]],P2[[4]],P1[[4]],P3[[5]],P2[[5]],P1[[5]],ncol=3, 
             top=paste("Percentage resolved                                                                                     FlankingSNP                                                                                    PhasedSNP mode"))









dev.copy2pdf(file =paste("data-raw/doni2/prevalence_mutation_",imutation,".pdf",sep=""))










Exact=F

if(Exact){
  
  cat("\n We compute the exact prevalences without simulation")
  
  
  casestudy_list<-list(
    "Case Study 1"=list(lambda_S=100,mu_S=100,major_cn=1,minor_cn=1, lambda_G=100, mu_G=100,mode=mymode),
    "Case Study 2"=list(lambda_S=90,mu_S=110,major_cn=1,minor_cn=1, lambda_G=100, mu_G=100,mode=mymode),
    "Case Study 3"=list(lambda_S=70,mu_S=160,major_cn=2,minor_cn=1,lambda_G=130, mu_G=100,mode=mymode),
    "Case Study 4"=list(lambda_S=15,mu_S=185,major_cn=1,minor_cn=1,lambda_G=100, mu_G=100,mode=mymode),
    "Case Study 5"=list(lambda_S=40,mu_S=160,major_cn=2,minor_cn=0, lambda_G=120, mu_G=80,mode=mymode)
  )
  
  for(ics in 1:length(casestudy_list) ){
    
    cat("\n\n\t CASE STUDY : ", names(casestudy_list[ics]) ,"\n\n")
    
    prevalence0 =do.call(getPrevalence, casestudy_list[[ics]])
    
    print(prevalence0)
    
  }
  
  
}




