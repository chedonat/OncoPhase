#library(OncoPhase.R)

#Chromosomes
##########
#Chromsomes list
chrom_list<-paste("chr",c(1:22,"X","Y","M"),sep="")
#Chromosomes sizes
hg19_dfsize<-list(chr1=249250621,chr2=243199373,chr3=198022430, chr4=191154276, chr5=180915260, chr6=171115067, 
                   chr7=159138663, chr8=146364022, chr9=141213431, chr10=135534747, chr11=135006516, chr12=133851895, 
                   chr13=115169878, chr14=107349540,chr15=102531392, chr16=90354753, chr17=81195210, chr18=78077248,  
                   chr19=59128983, chr20=63025520,chr21=48129895,chr22=51304566, chrX=155270560, chrY=59373566,chrM=16571)


#Take the sample list from a csv files, can also be entered manually

Patient="OP1019"

data_directory =  paste("../../OCCA/",Patient,"/",sep="")

samples_information_df<-read.csv(paste(data_directory,"samples_",Patient,".csv",sep=""),header=T)
sample_list<- initsSamples(samples_information_df)

region="chr22:1000-5000000"

chrom="chr15"



cat("\n\n Load the Phasing Association")
file_phase=paste(data_directory,
                 "Phasing/phasingassociation_somatic_germline_1_", Patient,".paf",sep="")
phasingassociation_df<- read.table( file = file_phase,header=T,row.names=1, sep="\t")
phasingassociation_df <- phasingassociation_df[phasingassociation_df$Chrom==chrom,]



cat("\n\n Load the master snp wellcount")
filesnpwellcount<-paste(data_directory,"Master/Chromosomes/AllelesCount_mastermatrix_",Patient,"_",chrom,".tsv.bz2",sep="")
mastersnpwellcountmatrix_df<- read.table( file = filesnpwellcount,header=T,row.names=1,sep="\t")


cat("\n\n Load the master ref wellcount")
filerefwellcount<-paste(data_directory,"Master/Chromosomes/RefAllelesCount_mastermatrix_",Patient,"_",chrom,".tsv.bz2",sep="")
masterrefwellcountmatrix_df<- read.table( file = filerefwellcount,header=T,row.names=1,sep="\t")

# cat("\n\n Load the master well fraction")
# file_wellfraction<-paste(data_directory,"Master/Chromosomes/AllelesFraction_mastermatrix_",Patient,"_",chrom,".tsv.bz2",sep="")
# masterwellfractionmatrix_df<- read.table( file = file_wellfraction,header=T,row.names=1,sep="\t")
# 



#Keeping only the tumour samples with LFR wells data
tumoursamples= as.character(samples_information_df[samples_information_df$alleletype=="LFR_wells","name"])

columns_list=c("Chrom","Start","End","Vartype","Ref","All","IsGermline", tumoursamples)

mastersnpwellcountmatrix_df=mastersnpwellcountmatrix_df[columns_list]
masterrefwellcountmatrix_df=masterrefwellcountmatrix_df[columns_list]





cat("\n\n Load the copynumbers")

#CopyNumber_df = data.frame()
MajorCopyNumber_df=data.frame()
MinorCopyNumber_df= data.frame()
NormalFraction_df= data.frame()

copynumber_groups = unique(as.character(samples_information_df[samples_information_df$alleletype=="LFR_wells","cnvprofile"]))
print(copynumber_groups)
for(cngroup in copynumber_groups)
{
  filecngroup<-  paste(data_directory,"CopyNumber/cn_",cngroup,"_",chrom,".csv",sep="")
  if (!file.exists(filecngroup))
  {
    cat("\n\n\t The file ",filecngroup, " do not exists" )
    next
  }

  
  newcopynumber<- read.csv( file = filecngroup,header=T,row.names=1)
  #We need to change the rownames so that to be only the chromosome and the position.// We remove this point of the code in the next versions.
  #make.unique since for sure we should have the same copy numbe rprofile at a particular position of the genome
  rownames(newcopynumber) = make.unique(unlist(lapply(rownames(newcopynumber), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse="_"))))
  
  
  if(ncol(MajorCopyNumber_df)==0)
  {
   # CopyNumber_df=newcopynumber[,"CopyNumber",drop=F]
    MajorCopyNumber_df=newcopynumber[,"MajorCopyNumber",drop=F]
    MinorCopyNumber_df=newcopynumber[,"MinorCopyNumber",drop=F]
    NormalFraction_df=newcopynumber[,"NormalFraction",drop=F]
  }else
  {
   # CopyNumber_df=cbind(CopyNumber_df,newcopynumber[,"CopyNumber",drop=F])
    MajorCopyNumber_df=cbind(MajorCopyNumber_df, newcopynumber[,"MajorCopyNumber",drop=F])
    MinorCopyNumber_df=cbind(MinorCopyNumber_df, newcopynumber[,"MinorCopyNumber",drop=F])
    NormalFraction_df=cbind(NormalFraction_df,newcopynumber[,"NormalFraction",drop=F])
  }
}
#names(CopyNumber_df) = copynumber_groups
names(MajorCopyNumber_df) = copynumber_groups
names(MinorCopyNumber_df) = copynumber_groups
names(NormalFraction_df) = copynumber_groups


# Expanding the copy number matrices

Major_cn_sample_G_df=as.data.frame(matrix(ncol=length(tumoursamples),nrow=nrow(MajorCopyNumber_df)))
colnames(Major_cn_sample_G_df) = tumoursamples
rownames(Major_cn_sample_G_df) = rownames(MajorCopyNumber_df)
normalfraction_cn_sample_G_df= Major_cn_sample_G_df
Minor_cn_sample_G_df= Major_cn_sample_G_df

#We fill the matrice with the apropriate copy number
for(sample in tumoursamples)
{
  cnv_profile= as.character(samples_information_df[samples_information_df$name==sample,"cnvprofile"])
  normalfraction_cn_sample_G_df[sample]=NormalFraction_df[cnv_profile]
  Major_cn_sample_G_df[sample]=MajorCopyNumber_df[cnv_profile]
  Minor_cn_sample_G_df[sample]=MinorCopyNumber_df[cnv_profile]
}



# Uniformizing the list of mutations, The allelecounts and the copynumbe rtables shouls have the same list of mutations.

masterrefwellcountmatrix_df["chrom_pos"]=paste(masterrefwellcountmatrix_df$Chrom,masterrefwellcountmatrix_df$End,sep="_")
masterrefwellcountmatrix_df["rownames"] = rownames(masterrefwellcountmatrix_df)


chrom_pos_alleleinfo=as.character(unlist(masterrefwellcountmatrix_df["chrom_pos"]))
chrom_pos_copynumber = rownames(Major_cn_sample_G_df)
chrom_pos_intersect=intersect(chrom_pos_alleleinfo, chrom_pos_copynumber)

rownames(masterrefwellcountmatrix_df) = make.unique(chrom_pos_alleleinfo)
rownames(mastersnpwellcountmatrix_df) = make.unique(chrom_pos_alleleinfo)


Major_cn_sample_G_df=Major_cn_sample_G_df[chrom_pos_intersect,]
Minor_cn_sample_G_df=Minor_cn_sample_G_df[chrom_pos_intersect,]
normalfraction_cn_sample_G_df=  normalfraction_cn_sample_G_df[chrom_pos_intersect,]
masterrefwellcountmatrix_df=masterrefwellcountmatrix_df[chrom_pos_intersect,]
mastersnpwellcountmatrix_df=mastersnpwellcountmatrix_df[chrom_pos_intersect,]

rownames(masterrefwellcountmatrix_df) = paste(masterrefwellcountmatrix_df$Chrom,masterrefwellcountmatrix_df$End,masterrefwellcountmatrix_df$Ref,masterrefwellcountmatrix_df$All,sep="_")
rownames(mastersnpwellcountmatrix_df) = rownames(masterrefwellcountmatrix_df) 
rownames(Major_cn_sample_G_df) = rownames(masterrefwellcountmatrix_df) 
rownames(Minor_cn_sample_G_df) = rownames(masterrefwellcountmatrix_df) 
rownames(normalfraction_cn_sample_G_df) = rownames(masterrefwellcountmatrix_df) 

mastersnpwellcountmatrix_df["chrom_pos"] = NULL
masterrefwellcountmatrix_df["chrom_pos"] = NULL

phasing_association_df=phasingassociation_df["Phased_List"]
names(phasing_association_df) = "PhasedGermlines"



snp_allelecount_df=mastersnpwellcountmatrix_df
ref_allelecount_df=masterrefwellcountmatrix_df
major_copynumber_df = Major_cn_sample_G_df
minor_copynumber_df = Minor_cn_sample_G_df
normalfraction_df = normalfraction_cn_sample_G_df
nbFirstColumns=6
region = chrom


masterprevalence_df=getPrevalence(snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df,nbFirstColumns,tumoursamples, region)
  
  



