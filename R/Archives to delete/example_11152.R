#library(OncoPhase.R)

#Take the sample list from a csv files, can also be entered manually
data_directory =  "../../OCCA/11152/"
Patient="11152"
samples_information_df<-read.csv(paste(data_directory,"samples_11152.csv",sep=""),header=T)
sample_list<- initsSamples(samples_information_df)



chrom = "chr22"



cat("\n\n Load the Phasing Association")
file_phase=paste(data_directory,
                 "Phasing/Chromosomes/phasingassociation_somatic_germline_1_", Patient,"_",chrom,".paf",sep="")
phasingassociation_chr<- read.table( file = file_phase,header=T,row.names=1, sep="\t")

cat("\n\n Load the master snp wellcount")
filesnpwellcount<-paste(data_directory,"Master/Chromosomes/AllelesCount_mastermatrix_",Patient,"_",chrom,".tsv.bz2",sep="")
mastersnpwellcountmatrix_chr<- read.table( file = filesnpwellcount,header=T,row.names=1,sep="\t")


cat("\n\n Load the master ref wellcount")
filerefwellcount<-paste(data_directory,"Master/Chromosomes/RefAllelesCount_mastermatrix_",Patient,"_",chrom,".tsv.bz2",sep="")
masterrefwellcountmatrix_chr<- read.table( file = filerefwellcount,header=T,row.names=1,sep="\t")

cat("\n\n Load the master well fraction")
file_wellfraction<-paste(data_directory,"Master/Chromosomes/AllelesFraction_mastermatrix_",Patient,"_",chrom,".tsv.bz2",sep="")
masterwellfractionmatrix_chr<- read.table( file = file_wellfraction,header=T,row.names=1,sep="\t")




#Removing the normal samples.
cat("\n\n Removong the Normal samples....")
Nsamples = intersect(colnames(mastersnpwellcountmatrix_chr[cifs:ncol(mastersnpwellcountmatrix_chr)]), normal_samples)
if(length(Nsamples )>0)
{
  for(nsample in Nsamples){
    mastersnpwellcountmatrix_chr[nsample] = NULL
    masterwellfractionmatrix_chr[nsample] = NULL
    masterrefwellcountmatrix_chr[nsample] = NULL
    
  }
  
}

tumour_samples = colnames(mastersnpwellcountmatrix_chr[cifs:ncol(mastersnpwellcountmatrix_chr)])
print(head(mastersnpwellcountmatrix_chr))
cat("\n")
print(tumour_samples)






cat("\n\n Load the copynumbers")

CopyNumber_df = data.frame()
MajorCopyNumber_df=data.frame()
MinorCopyNumber_df= data.frame()
NormalFraction_df= data.frame()
print(copynumber_groups)
for(cngroup in copynumber_groups)
{
  filecngroup<-  paste(data_directory,"CopyNumber/cn_",cngroup,"_",chrom,".csv",sep="")
  if (!file.exists(filecngroup))
    next
  
  newcopynumber<- read.csv( file = filecngroup,header=T,row.names=1)
  #We need to change the rownames so that to be only the chromosome and the position.// We remove this point of the code in the next versions.
  #make.unique since for sure we should have the same copy numbe rprofile at a particular position of the genome
  rownames(newcopynumber) = make.unique(unlist(lapply(rownames(newcopynumber), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse="_"))))
  
  
  if(ncol(CopyNumber_df)==0)
  {
    CopyNumber_df=newcopynumber[,"CopyNumber",drop=F]
    MajorCopyNumber_df=newcopynumber[,"MajorCopyNumber",drop=F]
    MinorCopyNumber_df=newcopynumber[,"MinorCopyNumber",drop=F]
    NormalFraction_df=newcopynumber[,"NormalFraction",drop=F]
  }else
  {
    CopyNumber_df=cbind(CopyNumber_df,newcopynumber[,"CopyNumber",drop=F])
    MajorCopyNumber_df=cbind(MajorCopyNumber_df, newcopynumber[,"MajorCopyNumber",drop=F])
    MinorCopyNumber_df=cbind(MinorCopyNumber_df, newcopynumber[,"MinorCopyNumber",drop=F])
    NormalFraction_df=cbind(NormalFraction_df,newcopynumber[,"NormalFraction",drop=F])
  }
}
names(CopyNumber_df) = copynumber_groups
names(MajorCopyNumber_df) = copynumber_groups
names(MinorCopyNumber_df) = copynumber_groups
names(NormalFraction_df) = copynumber_groups



nbcolumns_wellfraction=ncol(masterwellfractionmatrix_chr)

#We select only heterozygote mutations
if(OnlyHeterozygote)
{
  list_heterozygotemutations=rownames(masterrefwellcountmatrix_chr[rowSums(masterrefwellcountmatrix_chr[cifs:nbcolumns_wellfraction],na.rm=T)>0,])
  cat("\n\n\t ", length(list_heterozygotemutations), " Heterozygote mutations")
  masterwellfractionmatrix_chr=masterwellfractionmatrix_chr[list_heterozygotemutations,]
  #phasingassociation_chr=phasingassociation_chr[list_heterozygotemutations,]
  mastersnpwellcountmatrix_chr=mastersnpwellcountmatrix_chr[list_heterozygotemutations,]
  masterrefwellcountmatrix_chr=masterrefwellcountmatrix_chr[list_heterozygotemutations,]
}



#We restrict the rows
masterwellfractionmatrix_chr =masterwellfractionmatrix_chr[rowSums(masterwellfractionmatrix_chr[cifs:nbcolumns_wellfraction],na.rm=T)>0,]
#phasingassociation_chr =phasingassociation_chr[rownames(masterwellfractionmatrix_chr),]
mastersnpwellcountmatrix_chr=mastersnpwellcountmatrix_chr[rownames(masterwellfractionmatrix_chr),colnames(masterwellfractionmatrix_chr[cifs:nbcolumns_wellfraction])]
masterrefwellcountmatrix_chr=masterrefwellcountmatrix_chr[rownames(masterwellfractionmatrix_chr),colnames(masterwellfractionmatrix_chr[cifs:nbcolumns_wellfraction])]




#Important, we extract the submatrix only for somatic
somaticwellfractionmatrix_chr  = masterwellfractionmatrix_chr[masterwellfractionmatrix_chr$IsGermline==0, ]




cat("\n\t*******************\n Extract position to run \n\t *********************")

start_position_node=1
end_position_node = nrow(somaticwellfractionmatrix_chr)
if (number_node !=0)
{
  
  cat("\n\n\t We extract the ", number_node,"th part out of ",totalparts_node,"  parts")
  full_length=nrow(somaticwellfractionmatrix_chr)
  part_length=ceiling(full_length/totalparts_node)
  start_position_node= (number_node-1) * part_length + 1
  end_position_node = start_position_node + part_length -1
  if (end_position_node >nrow(somaticwellfractionmatrix_chr) )
    end_position_node = nrow(somaticwellfractionmatrix_chr)
  
  
}else
{
  cat("\n\t No subpart extracted for this node, the whole chromosome considered...")
}

cat("\n\nThis node will run from position ",start_position_node," to  position", end_position_node)


startmut=start_position_node
endmut = end_position_node


imutation=0


masterprevalence<-matrix(nrow=endmut-startmut+1,ncol=ncol(somaticwellfractionmatrix_chr))
masterprevalence<-as.data.frame(masterprevalence)
colnames(masterprevalence) <- colnames(somaticwellfractionmatrix_chr)
rownames(masterprevalence) <- rownames(somaticwellfractionmatrix_chr[startmut:endmut,])
masterprevalence[,1:nbfirstcolumns] = somaticwellfractionmatrix_chr[startmut:endmut,1:nbfirstcolumns]




