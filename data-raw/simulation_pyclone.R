
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


mymethod="PhasedSNPGeneral"
#mymethod="PhasedSNP"


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
      
      InputData = build_casestudy(lambda_G, mu_G,lambda_S,mu_S,major_cn,minor_cn)
      prevalence = getPrevalence(InputData$snp_allelecount_df, InputData$ref_allelecount_df, InputData$phasing_association_df, InputData$major_copynumber_df,InputData$minor_copynumber_df, min_cells=4, min_alleles=1,method=mymethod)
      

      
      if(mymethod=="PhasedSNPGeneral"){
        myprevalence = prevalence$Tumour1
      }else {
        myprevalence=as.numeric(strsplit(prevalence$Tumour1,":")[[1]][1])
      }
    
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
  
  df=melt(mse_prediction)
  df["Coverage"] = as.factor(depth)
  df$value = as.numeric(df$value)
  p1=ggplot(data=df,aes(x=Coverage,y=value,group = 1)) + geom_line(colour=jBrewColors[4],size=1.1, stat="identity")+
    xlab("") + ylab("Mean Square Error (MSE)")+
    theme(axis.text=element_text(size=12,face="bold"),  axis.title=element_text(size=12,face="bold"))
  
  
  prevalence_estimates_df=as.data.frame((prevalence_estimates))
  names(prevalence_estimates_df) =depth
  df=melt(prevalence_estimates_df) 
  names(df) = c("Coverage","Prevalence")
  means <- aggregate(Prevalence ~  Coverage, df, mean)
  p2=ggplot(data=df, aes(x=Coverage, y=Prevalence)) + geom_boxplot(fill=jBrewColors[7])+
    stat_summary(fun.y=mean, colour="darkred", geom="point",  shape=18, size=3) + 
    geom_hline(yintercept=prevalence0,color=jBrewColors[4],size=1,linetype =2)+ ylim(c(0,1)) + 
    geom_text(data = means, aes(label = format(Prevalence,digits=2), x=Coverage,y = 1.04*Prevalence, size=1),show.legend = FALSE)+
    xlab("Coverage") + ylab("Prevalence")+
    theme(axis.text=element_text(size=12,face="bold"),  axis.title=element_text(size=12,face="bold"))
  
  
  grid.arrange(p1,p2,ncol=1, top=paste("Mutation: ", imutation))
  
  dev.copy2pdf(file =paste("data-raw/doni2/prevalence_mutation_",imutation,".pdf",sep=""))
  

  
  
}


cat("\n We compute the exact prevalences without simulation")


casestudy_list<-list(
  "Case Study 1"=list(lambda_G=100, mu_G=100,lambda_S=100,mu_S=100,phi_G=0,major_cn=1,minor_cn=1),
  "Case Study 2"=list(lambda_G=100, mu_G=100,lambda_S=90,mu_S=110,phi_G=0,major_cn=1,minor_cn=1),
  "Case Study 3"=list(lambda_G=130, mu_G=100,lambda_S=70,mu_S=160,phi_G=0.3,major_cn=2,minor_cn=1),
  "Case Study 4"=list(lambda_G=100, mu_G=100,lambda_S=15,mu_S=185,phi_G=0,major_cn=1,minor_cn=1),
  "Case Study 5"=list(lambda_G=120, mu_G=80,lambda_S=40,mu_S=160,phi_G=0.2,major_cn=2,minor_cn=0)
)

for(ics in 1:length(casestudy_list) ){
  
  cat("\n\n\t CASE STUDY : ", names(casestudy_list[ics]) ,"\n\n")
  
  
  prevalence_estimates = matrix(nrow=n_simulations, ncol=n_d)
  cs0 = do.call(build_casestudy, casestudy_list[[ics]])
  prevalence0 = getPrevalence(cs0$snp_allelecount_df, cs0$ref_allelecount_df, cs0$phasing_association_df, cs0$major_copynumber_df,cs0$minor_copynumber_df, cs0$normalfraction_df)
  print(prevalence0)
  
}



