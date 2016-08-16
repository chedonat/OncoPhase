library(OncoPhase)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(reshape2)


jBrewColors <- brewer.pal(n = 8, name = "Set2")


casestudy_list<-list(
  "Case Study 1"=list(lambda_G=8, mu_G=5,lambda_S=3,mu_S=10,cnv_fraction=3/5,major_cn=2,minor_cn=1),
  "Case Study 2"=list(lambda_G=20, mu_G=8,lambda_S=3,mu_S=25,cnv_fraction=6/8,major_cn=3,minor_cn=1),
  "Case Study 3"=list(lambda_G=5, mu_G=3,lambda_S=3,mu_S=5,cnv_fraction=2/5,major_cn=1,minor_cn=0),
  "Case Study 4"=list(lambda_G=9, mu_G=6,lambda_S=7,mu_S=8,cnv_fraction=3/6,major_cn=2,minor_cn=1),
  "Case Study 5"=list(lambda_G=4, mu_G=4,lambda_S=2,mu_S=6,cnv_fraction=0,major_cn=1,minor_cn=1),
  "Case Study 6"=list(lambda_G=16, mu_G=8,lambda_S=14,mu_S=10,cnv_fraction=4/8,major_cn=3,minor_cn=1),
  "Case Study 7"=list(lambda_G=8, mu_G=0,lambda_S=1,mu_S=7,cnv_fraction=1,major_cn=2,minor_cn=0),
  "Case Study 8"=list(lambda_G=6, mu_G=2,lambda_S=3,mu_S=5,cnv_fraction=4/6,major_cn=1,minor_cn=0),
  "Case Study 9"=list(lambda_G=7, mu_G=1,lambda_S=6,mu_S=2,cnv_fraction=3/4,major_cn=2,minor_cn=0),
  "Case Study A"=list(lambda_G=8, mu_G=6,lambda_S=6,mu_S=8,cnv_fraction=2/6,major_cn=2,minor_cn=1),
  "Case Study B"=list(lambda_G=8, mu_G=4,lambda_S=4,mu_S=8,cnv_fraction=2/6,major_cn=2,minor_cn=0),
  "Case Study C"=list(lambda_G=8, mu_G=12,lambda_S=6,mu_S=14,cnv_fraction=4/8,major_cn=2,minor_cn=1)
  
)

casestudy_list2<-list(
  "Case Study 1"=list(lambda_G=8, mu_G=5,lambda_S=3,mu_S=10,major_cn=2,minor_cn=1),
  "Case Study 2"=list(lambda_G=20, mu_G=8,lambda_S=3,mu_S=25,major_cn=3,minor_cn=1),
  "Case Study 3"=list(lambda_G=5, mu_G=3,lambda_S=3,mu_S=5,major_cn=1,minor_cn=0),
  "Case Study 4"=list(lambda_G=9, mu_G=6,lambda_S=7,mu_S=8,major_cn=2,minor_cn=1),
  "Case Study 5"=list(lambda_G=4, mu_G=4,lambda_S=2,mu_S=6,major_cn=1,minor_cn=1),
  "Case Study 6"=list(lambda_G=16, mu_G=8,lambda_S=14,mu_S=10,major_cn=3,minor_cn=1),
  "Case Study 7"=list(lambda_G=8, mu_G=0,lambda_S=1,mu_S=7,major_cn=2,minor_cn=0),
  "Case Study 8"=list(lambda_G=6, mu_G=2,lambda_S=3,mu_S=5,major_cn=1,minor_cn=0),
  "Case Study 9"=list(lambda_G=7, mu_G=1,lambda_S=6,mu_S=2,major_cn=2,minor_cn=0),
  "Case Study A"=list(lambda_G=8, mu_G=6,lambda_S=6,mu_S=8,major_cn=2,minor_cn=1),
  "Case Study B"=list(lambda_G=8, mu_G=4,lambda_S=4,mu_S=8,major_cn=2,minor_cn=0),
  "Case Study C"=list(lambda_G=8, mu_G=12,lambda_S=6,mu_S=14,major_cn=2,minor_cn=1)
  
)



d_range = c( 10,  20, 30, 60,  120, 300, 500, 1000)
n_simulations = 200
n_d = length(d_range)



mymethod="PhasedSNPGeneral"
mymethod="PhasedSNP"

cstudy = casestudy_list2
for(ics in 1:length(cstudy) ){
#for(ics in 11:11 ){
  
  cat("\n\n\t CASE STUDY : ", names(cstudy[ics]) ,"\n\n")
  
  
  prevalence_estimates = matrix(nrow=n_simulations, ncol=n_d)
  cs0 = do.call(build_casestudy, cstudy[[ics]])
  prevalence0 = getPrevalence(cs0$snp_allelecount_df, cs0$ref_allelecount_df, cs0$phasing_association_df, cs0$major_copynumber_df,cs0$minor_copynumber_df, method=mymethod)
 
  
  
 
  if (mymethod!="PhasedSNPGeneral")
  {
    prevalence0=as.numeric(strsplit(prevalence0$Tumour1,":")[[1]][1])
  }else{
    prevalence0= prevalence0$Tumour1
  }
  
#  print( prevalence0)

 # next
  
  for ( di in 1:n_d ) {
    
    cat("\n\n\t Coverage : ", d_range[di], "\n\t\t Simulations : " )
    
    for ( si in 1:n_simulations ) {
      
      if(si%%100==0)   cat(" ", si )
      
      # generate data
      cs = do.call(build_casestudy,append(cstudy[[ics]],list(depthOfCoverage=d_range[di])))
      
      # Run the case
      prevalence = getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df, cs$major_copynumber_df,cs$minor_copynumber_df,   method=mymethod ,min_cells=4, min_alleles=0,
                                 #cnv_fraction_df =cs$cnv_fraction_df
                                 )
      
    #  print(prevalence)
      if (mymethod!="PhasedSNPGeneral")
      {
        if(!is.na(prevalence$Tumour1)){
          prevalence=as.numeric(strsplit(prevalence$Tumour1,":")[[1]][1])
        }else{
          prevalence= prevalence$Tumour1
        }

      }else{
        prevalence= prevalence$Tumour1
      }
      
     # if(prevalence>0.8 )
     #   stop(1)
      
      
      if(!is.na(prevalence))
      if(prevalence>1.0000001 || prevalence< -0.0000001)
        stop()
      #print the result
      prevalence_estimates[si, di] = prevalence
      
    }
    
  #  stop(1)
    
  }
  
  
  
  #print the result
  error_estimates = prevalence_estimates - prevalence0
  mse_prediction=apply(error_estimates,2,function(x) mean(x*x,na.rm=T))
  names(mse_prediction) = d_range
  
  
  names(prevalence_estimates) = d_range
  
  df=melt(mse_prediction)
  df["Coverage"] = as.factor(d_range)
  df$value = as.numeric(df$value)
  p1=ggplot(data=df,aes(x=Coverage,y=value,group = 1)) + geom_line(colour=jBrewColors[4],size=1.1, stat="identity")+
    xlab("") + ylab("Mean Square Error (MSE)")+
    theme(axis.text=element_text(size=12,face="bold"),  axis.title=element_text(size=12,face="bold"))
  

  
  prevalence_estimates_df=as.data.frame((prevalence_estimates))
  names(prevalence_estimates_df) =d_range
  df=melt(prevalence_estimates_df) 
  names(df) = c("Coverage","Prevalence")
  means <- aggregate(Prevalence ~  Coverage, df, mean)
  p2=ggplot(data=df, aes(x=Coverage, y=Prevalence)) + geom_boxplot(fill=jBrewColors[7])+
    stat_summary(fun.y=mean, colour="darkred", geom="point",  shape=18, size=3) + 
    geom_hline(yintercept=prevalence0,color=jBrewColors[4],size=1,linetype =2)+ ylim(c(0,1)) + 
    geom_text(data = means, aes(label = format(Prevalence,digits=2), x=Coverage,y = 1.04*Prevalence, size=1),show.legend = FALSE)+
    xlab("Coverage") + ylab("Prevalence")+
    theme(axis.text=element_text(size=12,face="bold"),  axis.title=element_text(size=12,face="bold"))
  
  
  grid.arrange(p1,p2,ncol=1, top=names(cstudy[ics]))
  
  dev.copy2pdf(file =paste("data-raw/Plots2/",names(cstudy[ics]),".pdf",sep=""))
  


  
}











