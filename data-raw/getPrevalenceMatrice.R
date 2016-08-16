#' @export
#getPrevalence_Matrice<-function( snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,mode="PhasedSNP",cnv_fraction=NULL, phasing_association_df=NULL,NormalcellContamination_df=NULL,tumoursamples=NULL,  nbFirstColumns=3, region=NULL,detail=TRUE, LocusRadius = 10000,NoPrevalence.action="Skip",Trace=FALSE,SameTumour=TRUE)
{
  
  
  #Extraction of the list of somatic mutations the cellular prevalence will be computed.
  somatic_snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$IsGermline==0, ]
  # If a region is provided, a restriction is performed on the given region
  if(!is.null(region)){
    region_parts= unlist(strsplit(region,":"))
    
    chrom = region_parts[1]
    startPosition = 1
    endPosition = hg19_dfsize[chrom]
    
    somatic_snp_allelecount_df = somatic_snp_allelecount_df[somatic_snp_allelecount_df$Chrom == chrom , ]
    
    if(length(region_parts)>1){
      coordinates = unlist(strsplit(region_parts[2],"-"))
      startPosition = as.numeric(coordinates[1])
      endPosition = as.numeric(coordinates[2])
      somatic_snp_allelecount_df = somatic_snp_allelecount_df[somatic_snp_allelecount_df$Chrom == chrom & somatic_snp_allelecount_df$Start >= startPosition & somatic_snp_allelecount_df$End <= endPosition,]
    }
  }
  
  
  #set the mode if numeric, 0=SNVOnly, 1 = PhasedSNP, 2=FlankingSNP, 3 = OptimalSNP
  numeric_mode=c("SNVOnly", "PhasedSNP","FlankingSNP","OptimalSNP")
  if(is.numeric(mode))
  {
    if(mode %in% c(0,1,2,3))
    {
      mode = numeric_mode[mode +1 ]
    }else{
      stop("\n\n Mode parameter, if numeric,  should be either 0, 1,  2 or 3")
    }
  }
  
  
  
  #We then  retrieve for each somatic mutation the list of germline mutations to consider for the prevalence computation
  # a) if PhasedSNP mode, then the considered germline are the germline mutations phased to the somatic mutation and located within the same locus
  # b) if FlankingSNP mode then the considered germline are the close germlines located within LocusRadius distance from the somatic mutation.
  #c) if OptimalSNP mode tho columns are provided, the first for the phased germline, and if only there is not phasing information for this mutation, then the second column contains the close germlines located within LocusRadius
  
  
  # LinkedGermlineMutation=getLocusGermlineMutations(somatic_snp_allelecount_df, snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,cnv_fraction,phasing_association_df,  tumoursamples,mode,  LocusRadius)
  
  
  
  #Preparing the data frame to contains the somatic mutations cellular prevalences.
  masterprevalence<-matrix(nrow=nrow(somatic_snp_allelecount_df),ncol=nbFirstColumns + length(tumoursamples))
  masterprevalence<-as.data.frame(masterprevalence)
  colnames(masterprevalence) <- c(colnames(somatic_snp_allelecount_df[1:nbFirstColumns]),tumoursamples)
  rownames(masterprevalence) <- rownames(somatic_snp_allelecount_df)
  masterprevalence[1:nbFirstColumns] = somatic_snp_allelecount_df[1:nbFirstColumns]
  
  #  for (imut in 1:nrow(masterprevalence))
  for (imut in 1:5)
  {
    
    #Mutation name and mutation position
    mut <- rownames(masterprevalence[imut,]); 
    mut_pos=as.numeric(masterprevalence[imut,"End"])
    
    #For each mutation, we need to extract one value or one vector (if multiple samples)  of :
    # - lambda_G and mu_G : Respectively Variant and reference coverage/count of the phased/nearby Germline Mutations
    # - lambda_S and mu_S : Respectively Variant and reference coverage/count of the somatic mutations
    # - major_cn and minor_cn : major and minor copy number at the locus of the somatic mutations 
    # - CNV_fraction : If provided, fraction of cells affected by the CNV
    # - NormalCell_fraction : If provided, fraction of normal cell contamination.
    
    
    
    mode_locus=mode  # Mode to consider while computing the mutation prevalence and retrievning the germline on the same locus with the mutation
    if(mode=="OptimalSNP") # If Mode = Optimal, then we  willtry first PhasedSNP, latter if no germline found, we  will try FlankingSNP
      mode_locus="PhasedSNP"
    
    
    ############################
    #####Source of information 1: The list of linked germline mutations 
    ##############
    
    if(mode!="SNVOnly")
    {
      
      linked_germline_df=getLocusGermlineMutations(somatic_snp_allelecount_df[mut,], snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,cnv_fraction,phasing_association_df,  tumoursamples,mode_locus,  LocusRadius)
      linked_germline=as.character(linked_germline_df[mut,"LinkedGermlines"])
      
      
      if(is.null(linked_germline)|| is.na(linked_germline))
        if((mode_locus=="PhasedSNP")&&(mode=="OptimalSNP")){
          mode_locus="FlankingSNP"
          linked_germline_df=getLocusGermlineMutations(somatic_snp_allelecount_df[mut,], snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,cnv_fraction,phasing_association_df,  tumoursamples,mode_locus,  LocusRadius)
          linked_germline=as.character(linked_germline_df[mut,"LinkedGermlines"])
          
        }
      
      
      #       cat("\n\n Linked for mutation ", mut,"\n")
      #       print(linked_germline)
      #       cat("\n\n mode locus ", mut,"\n")
      #       print(mode_locus)
      
      if(is.null(linked_germline)|| is.na(linked_germline))
        next
      
      
      
      
      #List of germline mutations linked to the considered somatic mutation
      linkedGermlines_list<-unlist(strsplit(linked_germline,":"))
      if (length(unlist(linkedGermlines_list))==0)
        next
      
      
    }
    
    
    ############################
    #####Source of information 2: The Copy Number 
    ##############
    
    #For the somatic mutation
    if(!is.null(cnv_fraction)) phi_cn_sample_somatic = cnv_fraction[mut,tumoursamples,drop=F]
    if(!is.null(NormalcellContamination_df))  NC_sample_somatic = NormalCellContamination_df[mut,tumoursamples,drop=F]
    major_cn_sample_somatic = major_copynumber_df[mut,tumoursamples, drop=F]
    minor_cn_sample_somatic = minor_copynumber_df[mut,tumoursamples,drop=F]
    NormalContamination_sample_somatic = minor_copynumber_df[mut,tumoursamples,drop=F]
    
    ############################
    #####Source of information 3: The Allele Count 
    
    #Somatic
    #wellfraction_somatic =snp_allelecount_df[mut,cifs:nbcolumns_wellfraction]
    snpwellcount_somatic =  data.matrix(snp_allelecount_df[mut,tumoursamples,drop=F])
    refwellcount_somatic = data.matrix(ref_allelecount_df[mut,tumoursamples,drop=F])
    
    #Germline
    # a) if mode_locus=PhasedSNP, according to the MLE (EM) estimation, at each sample, we consider the average counts of the linked germline.
    # b) if mode_locus=FlankingSNP, at each sample, we consider the closest germline mutation having all the required information (Alleles Count and Copy number Information)
    
    if(mode_locus=="PhasedSNP")
    {
      snpwellcount_germlines=  colMeans(data.matrix(snp_allelecount_df[linkedGermlines_list,tumoursamples, drop=F]), 
                                        na.rm=T)
      refwellcount_germlines =colMeans( data.matrix(ref_allelecount_df[linkedGermlines_list,tumoursamples, drop=F]), 
                                        na.rm=T)
      
    }else  if(mode_locus=="FlankingSNP")
    {
      snpwellcount_germlines=  vector("numeric",length=length(tumoursamples))
      refwellcount_germlines = vector("numeric",length=length(tumoursamples))
      
      names(snpwellcount_germlines) = tumoursamples
      names(refwellcount_germlines) = tumoursamples
      
      
      #For the Linked germline mutations
      major_cn_sample_LinkedGermline_df = major_copynumber_df[linkedGermlines_list,tumoursamples,drop=F ]
      snp_allelecount_LinkedGermline_df= data.matrix(snp_allelecount_df[linkedGermlines_list,tumoursamples, drop=F])
      ref_allelecount_LinkedGermline_df= data.matrix(ref_allelecount_df[linkedGermlines_list,tumoursamples, drop=F])
      #To fix a bug when only one germline, 
      snp_allelecount_LinkedGermline_df=as.data.frame(snp_allelecount_LinkedGermline_df)
      ref_allelecount_LinkedGermline_df=as.data.frame(ref_allelecount_LinkedGermline_df)
      
      
      for(sample in tumoursamples)
      {
        
        selectedlinkedGermlines_list = rownames(major_cn_sample_LinkedGermline_df[!is.na(major_cn_sample_LinkedGermline_df[,sample]) ,sample,drop=F])
        selectedlinkedGermlines_list = intersect(selectedlinkedGermlines_list,rownames(snp_allelecount_LinkedGermline_df[!is.na(snp_allelecount_LinkedGermline_df[,sample]),sample,drop=F] ))
        selectedlinkedGermlines_list =  intersect(selectedlinkedGermlines_list,rownames(ref_allelecount_LinkedGermline_df[!is.na(ref_allelecount_LinkedGermline_df[,sample]),sample,drop=F] ))
        
        #Compute the distance to the somatic mutation
        distance_germlines= snp_allelecount_df[selectedlinkedGermlines_list,"End", drop=F]
        distance_germlines["distance"] = abs(as.numeric(distance_germlines$End) - mut_pos)
        #select the closest germline
        closest_germline = rownames(distance_germlines[distance_germlines$distance== min(distance_germlines$distance, na.rm=T), ])
        
        if(length(closest_germline)==0)
          next
        
        #retrieve its allele counts
        snpwellcount_germlines[sample] = snp_allelecount_LinkedGermline_df[closest_germline, sample]
        refwellcount_germlines[sample] = ref_allelecount_LinkedGermline_df[closest_germline, sample]
      }
      
    }else{ #mode_locus=SNVOnly
      snpwellcount_germlines=NA
      refwellcount_germlines=NA
    }
    
    
    # wellfraction_germlines=snp_allelecount_df[linkedGermlines_list,cifs:nbcolumns_wellfraction]
    
    ##### Preparing the input for the formulw
    
    ##Allele Counts
    lambda_somatic=snpwellcount_somatic[mut, tumoursamples,drop=F] # Somatic variant counts 
    mu_somatic=refwellcount_somatic[mut, tumoursamples,drop=F] # Somatic reference counts 
    lambda_LinkedGermline<-snpwellcount_germlines #  Germline variant counts 
    mu_LinkedGermline<- refwellcount_germlines # Germline reference counts 
    
    #Normalising the somatic count to the germline count, 
    #germline and somatic are set to have the same total count of Alleles
    #Total_allele_count=lambda_LinkedGermline + mu_LinkedGermline
    #lambda_somatic = (lambda_somatic/(lambda_somatic+mu_somatic)) * Total_allele_count
    #mu_somatic= Total_allele_count - lambda_somatic
    
    ###Copy number profile (simply the one of the somatic mutation locus)
    if(!is.null(cnv_fraction)) {
      phi_cn= unlist(phi_cn_sample_somatic)
    }else{
      phi_cn=NULL
    }
    major_cn= unlist(major_cn_sample_somatic)
    minor_cn = unlist(minor_cn_sample_somatic)
    
    #Summarising the inputs
    # stop(30)
    Trace=F
    if(Trace){
      cat("\n\n The inputs are : ")
      cat("\n\t lambda_somatic :\n");print( lambda_somatic)
      cat("\n\t mu_somatic  :\n");print(  mu_somatic )
      cat("\n\t  lambda_LinkedGermline :\n");print(lambda_LinkedGermline )
      #stop(20)
      cat("\n\t  mu_LinkedGermline :\n");print(mu_LinkedGermline  )
      if(!is.null(cnv_fraction)) {cat("\n\t  phi_cn :\n");print( phi_cn )}
      cat("\n\t  major_cn :\n");print( major_cn )
      cat("\n\t minor_cn  :\n");print( minor_cn  )
    }
    
    # stop(10)
    ###Calling the prevalence quantification
    if(detail)
      detail=2
    
    prev_somatic=getPrevalence(lambda_somatic,mu_somatic,major_cn,minor_cn, lambda_LinkedGermline , mu_LinkedGermline,  detail,mode_locus,Trace,SameTumour)
    
    #  print(prev_somatic)
    # cat("\n")
    masterprevalence[mut,names(prev_somatic)]
    
    
    
    masterprevalence[mut,names(prev_somatic)] = prev_somatic
  }
  
  
  masterprevalence
}

