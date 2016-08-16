library(OncoPhase)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)


jBrewColors <- brewer.pal(n = 8, name = "Set2")


d_range = c( 10, 15, 20, 30, 60,  120, 300, 500, 1000,2000)
n_simulations = 3000
n_d = length(d_range)

prevalence_estimates = matrix(nrow=n_simulations, ncol=n_d)

cs0 = build_casestudy(lambda_G=16, mu_G=8, lambda_S=14, mu_S=10, phi_G=4/8, major_cn=3, minor_cn=1)
cs0= build_casestudy(lambda_G=8, mu_G=12,lambda_S=6,mu_S=14,phi_G=4/8,major_cn=2,minor_cn=1)
prevalence0 = getPrevalence(cs0$snp_allelecount_df, cs0$ref_allelecount_df, cs0$phasing_association_df, cs0$major_copynumber_df,cs0$minor_copynumber_df, cs0$normalfraction_df)

for ( di in n_d:n_d ) {
  
  cat("\n\n\t Coverage : ", d_range[di], "\n\t\t Simulations : " )
  
  for ( si in 1:n_simulations ) {
    
    if(si%%100==0)   cat(" ", si )
    
    # generate data
    cs = build_casestudy(lambda_G=16, mu_G=8, lambda_S=14, mu_S=10, phi_G=4/8, major_cn=3, minor_cn=1, depthOfCoverage=d_range[di])
    cs = build_casestudy(lambda_G=8, mu_G=12,lambda_S=6,mu_S=14,phi_G=4/8,major_cn=2,minor_cn=1#, depthOfCoverage=d_range[di]
                        ,depthOfCoverage=20
                        )
  
    # Run the case
    prevalence = getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df, cs$major_copynumber_df,cs$minor_copynumber_df, cs$normalfraction_df)
  
   # cat(" ", prevalence$Tumour1)
    
    #print the result
    prevalence_estimates[si, di] = prevalence$Tumour1
  
  }

}




#print the result
error_estimates = prevalence_estimates - prevalence0$Tumour1
mse_prediction=apply(error_estimates,2,function(x) mean(x*x,na.rm=T))
names(mse_prediction) = d_range


names(prevalence_estimates) = d_range

df=melt(mse_prediction)
df["Coverage"] = d_range
p1=ggplot(data=df,aes(x=Coverage,y=value)) + geom_line(colour=jBrewColors[4],size=1.5)


prevalence_estimates_df=as.data.frame((prevalence_estimates))
names(prevalence_estimates_df) =d_range
df=melt(prevalence_estimates_df) 
names(df) = c("Coverage","Prevalence")
means <- aggregate(Prevalence ~  Coverage, df, mean)
p2=ggplot(data=df, aes(x=Coverage, y=Prevalence)) + geom_boxplot(fill=jBrewColors[7],show.legend = FALSE)+
  stat_summary(fun.y=mean, colour="darkred", geom="point",  shape=18, size=5,show.legend = FALSE) + 
  geom_hline(yintercept=prevalence0$Tumour1,color=jBrewColors[4],size=1,show.legend = FALSE)+ ylim(c(0,1)) + 
  geom_text(data = means, aes(label = format(Prevalence,digits=2), x=Coverage,y = 1.04*Prevalence, size=2),show.legend = FALSE)


grid.arrange(p1,p2,ncol=1)







