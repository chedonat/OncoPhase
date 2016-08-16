



#' @export
#getLocusGermlineMutations<-function(somatic_snp_allelecount_df, snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,cnv_fraction,phasing_association_df,  tumoursamples,mode="PhasedSNP", formula="matrix", LocusRadius)
{
  
  
  #Preparing the data frame to contains the Germlines linked to each somatic mutation.
  LinkedGermlines<-matrix(nrow=nrow(somatic_snp_allelecount_df),ncol=4)
  LinkedGermlines<-as.data.frame(LinkedGermlines)
  colnames(LinkedGermlines) <- c(colnames(somatic_snp_allelecount_df[1:3]),"LinkedGermlines")
  rownames(LinkedGermlines) <- rownames(somatic_snp_allelecount_df)
  LinkedGermlines[1:3] = somatic_snp_allelecount_df[1:3]
  
  imutation=0 #Simply an iterator
  
  germline_snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$IsGermline==1, ]
  
  for (imut in 1:nrow(LinkedGermlines))
  {
    #Mutation name and mutation position
    mut <- rownames(LinkedGermlines[imut,]); 
    mut_pos=as.numeric(LinkedGermlines[imut,"End"])
    
    
    
    #If mode=PhasedSNP, the candidate germline mutations are all the SNP phased to the somatic mutation 
    #If modeFlankingSNP, the candidate germline mutations are all the SNP located within LocusRadius distance of the somatic Mutation.
    
    if(mode=="PhasedSNP"){
      if(is.null(phasing_association_df))
        stop("\n\n if mode=PhasedSNP the phasing association matrice should be provided")
      CandidateGermlines = as.character(phasing_association_df[mut,"PhasedMutations"])
      if(is.null(CandidateGermlines))
        next
      CandidateGermlines<-unlist(strsplit(CandidateGermlines,":"));  
      if (length(unlist(CandidateGermlines))==0)
        next
    }else if(mode=="FlankingSNP"){
      CandidateGermlines =rownames(germline_snp_allelecount_df[germline_snp_allelecount_df$End >= mut_pos - LocusRadius & germline_snp_allelecount_df$End <= mut_pos + LocusRadius, ])
    }else{
      stop("\n\n The mode parameters shuld be either PhasedSNP either FlankingSNP")
    }
    
    
    #Now, Among the candidate germline, we want to keep only those present on the same locus with the somatic mutation
    
    submatrix_CandidateGermlines = snp_allelecount_df[CandidateGermlines, 1:3 , drop=F] #retrieve a submatrix of the CandidateGermlines with the chromosome and the position.
    
    submatrix_CandidateGermlines_leftside = submatrix_CandidateGermlines[submatrix_CandidateGermlines$End <mut_pos,1:3 , drop=F ] #germlines at the left of the mutation 
    submatrix_CandidateGermlines_rightside = submatrix_CandidateGermlines[submatrix_CandidateGermlines$End >mut_pos,1:3 , drop=F ] # germline at the right
    
    #order the leftside from highest to smallest position, leave the rightside from smallest to highest
    orders_phase=order(submatrix_CandidateGermlines_leftside["End"],decreasing=T)
    submatrix_CandidateGermlines_leftside=submatrix_CandidateGermlines_leftside[orders_phase,]
    
    
    at_least_one_good_germline=F # is there atleast one good germline kept?
    
    at_least_one_good_germline=F # is there atleast one good germline kept?
    germline_to_exclude=c()
    for(sample in tumoursamples)
    {
      samelocus_germline=c()
      if(!is.null(cnv_fraction))  
        phi_som=cnv_fraction[mut,sample]
      major_som=major_copynumber_df[mut,sample]
      minor_som=major_copynumber_df[mut,sample]
      # allelecount_som=somatic_snp_allelecount_df[mut,sample]
      #  if (is.na(allelecount_som)) #Nothing wont be done on this sample anyway withut the snp count of the somatic.
      #    next
      
      absence_copynumberprofile=c()
      count_lower_than_somatic<-c() # For more accuracy in case of abundance of germline, someone can choose to exclude germline having an allele count less than the somatic allele count.
      
      for (germ in rownames(submatrix_CandidateGermlines_leftside))
      {
        if(!is.null(cnv_fraction))   
          phi_germ=cnv_fraction[germ , sample]
        major_germ=major_copynumber_df[germ , sample]
        minor_germ=major_copynumber_df[germ, sample]
        
        if (  (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) &&
              (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          break
        }
        
        
      }
      
      for (germ in rownames(submatrix_CandidateGermlines_rightside))
      {
        if(!is.null(cnv_fraction))   
          phi_germ=cnv_fraction[germ , sample]
        major_germ=major_copynumber_df[germ , sample]
        minor_germ=major_copynumber_df[germ, sample]
        
        
        if ( (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) && 
             (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          break
        }
        
        # if( is.na(major_germ))
        #   absence_copynumberprofile=c(absence_copynumberprofile,germ)
        
      }
      
      notsamelocus_germline = setdiff(CandidateGermlines, samelocus_germline)
      #germline_to_exclude_at_this_sample=c(notsamelocus_germline,absence_copynumberprofile)
      
      germline_to_exclude=unique(c(germline_to_exclude,notsamelocus_germline))
    }
    
    # if(!  at_least_one_good_germline) next # there is no good gerline
    
    samelocus_germlines= setdiff(CandidateGermlines,germline_to_exclude)
    LinkedGermlines[mut,"LinkedGermlines"] = paste(samelocus_germlines, collapse=":")
    
  }
  
  LinkedGermlines 
  
}


