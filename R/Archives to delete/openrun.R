#' Compute somatic mutations cellular prevalence in Cancer
#' 
#' @param x 
#getPrevalence<-function(snp_allelecount_df, ref_allelecount_df,phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df,nbFirstColumns,tumoursamples, region=NULL)
{
  
  
  
  
  # Extract the somatic mutations 
  
  
  compulsory_columns=c("Chrom","End","IsGermline")
  
  if (length(setdiff(compulsory_columns,colnames(snp_allelecount_df)))>0)
    stop(" The allele count master matrices should have at least the following headers columns : ", compulsory_columns)
  
  if (length(setdiff(compulsory_columns,colnames(ref_allelecount_df)))>0)
    stop(" The allele count master matrices should have at least the following headers columns : ", compulsory_columns)
  
  
  somatic_snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$IsGermline==0, ]
  
  
  
  # Region to compute the prevalence
  
  if(!is.null(region)){
    region_parts= unlist(strsplit(region,":"))
    
    chrom = region_parts[1]
    startPosition = 1
    endPosition = hg19_dfsize[chrom]
    
    somatic_snp_allelecount_df = somatic_snp_allelecount_df[somatic_snp_allelecount_df$Chrom == chrom , ]
    
    
    if(length(region_parts)>1){
      coordinates = unlist(strsplit(region_parts[2],"-"))
      startPosition = coordinates[1]
      endPosition = coordinates[2]
      somatic_snp_allelecount_df = somatic_snp_allelecount_df[somatic_snp_allelecount_df$Chrom == chrom & somatic_snp_allelecount_df$Start >= startPosition & somatic_snp_allelecount_df$End <= endPosition,]
      
      }
    
    }
  
  
  
  
  
  
  imutation=0
  
  
  masterprevalence<-matrix(nrow=nrow(somatic_snp_allelecount_df),ncol=nbFirstColumns + length(tumoursamples))
  masterprevalence<-as.data.frame(masterprevalence)
  colnames(masterprevalence) <- c(colnames(somatic_snp_allelecount_df[1:nbFirstColumns]),tumoursamples)
  rownames(masterprevalence) <- rownames(somatic_snp_allelecount_df)
  masterprevalence[1:nbFirstColumns] = somatic_snp_allelecount_df[1:nbFirstColumns]
  
  for (imut in 1:nrow(masterprevalence))
  {
    
    
    imutation =imutation+1; 
    mut <- rownames(masterprevalence[imut,]); mut_pos=as.numeric(masterprevalence[imut,"End"])
    
    
    #  output=(((imutation-1) %%100)==0)
    
    
    # if(output)
    #   cat (" \n\n\n\n  mut ",imutation,"/ ",nrow(masterprevalence)," (", mut,") ")
    
    
    
    ############################
    #####Source of information 1: The list of phased germline
    ##############
    
    # Retrieving the phased germline, if not phased germline skip
    
    phased_germline=""
    phased_germline=  as.character(phasing_association_df[mut,"PhasedGermlines"])
    
    if(is.null(phased_germline))
      next
    phased_list<-unlist(strsplit(phased_germline,":"));  
    if (length(unlist(phased_list))==0)
      next
    
    #Keep only the germline having allele counts
    phased_list=intersect(phased_list, rownames(snp_allelecount_df[snp_allelecount_df$IsGermline==1,]))
    
    if (length(unlist(phased_list))==0)
      next
    
    #  if(output) 
    #   cat ( " :  ", length(phased_list), "Phased present in the Alelle matrice")  
    
    
    
    
    
    ############################
    #####Source of information 2: The Copy Number 
    ##############
    
    
    #For the somatic mutation
    phi_cn_sample_somatic_list = 1-normalfraction_df[mut,]
    Major_cn_sample_somatic_list = major_copynumber_df[mut, ]
    Minor_cn_sample_somatic_list = minor_copynumber_df[mut,]
    Total_cn_sample_somatic_list = Major_cn_sample_somatic_list + Minor_cn_sample_somatic_list 
    
    
    #For the phased germline mutations
    phi_cn_sample_PhasedGermline_df = 1-normalfraction_df[phased_list,]
    Major_cn_sample_PhasedGermline_df = major_copynumber_df[phased_list, ]
    Minor_cn_sample_PhasedGermline_df = minor_copynumber_df[phased_list,]
    Total_cn_sample_PhasedGermline_df = Major_cn_sample_PhasedGermline_df + Minor_cn_sample_PhasedGermline_df
    
    
    
    ############################
    #####Source of information 3: The Allele Count 
    ##############
    ############################
    
    
    #Somatic
    #wellfraction_somatic =snp_allelecount_df[mut,cifs:nbcolumns_wellfraction]
    snpwellcount_somatic =  snp_allelecount_df[mut,tumoursamples]
    refwellcount_somatic = ref_allelecount_df[mut,tumoursamples]
    
    #Germline
    # wellfraction_germlines=snp_allelecount_df[phased_list,cifs:nbcolumns_wellfraction]
    snpwellcount_germlines=  snp_allelecount_df[phased_list,tumoursamples]
    refwellcount_germlines = ref_allelecount_df[phased_list,tumoursamples]
    
    
    
    ############################
    ##### Mutations locus restriction
    ##############
    ############################
    
    # We need to pass trough each samples,  and assign NA to the entries of any phased germline appearing to not be on the same locus with the somatic mutation ( not having same phi, major_cn and minor_cn as the somatic and with no  cn breakpoint between the germline and the somatic.)
    
    submatrix_phased = snp_allelecount_df[phased_list, 1:3 , drop=F] #retrieve a submatrix of the phased with the chromosome and the position.
    
    submatrix_phased_leftside = submatrix_phased[submatrix_phased$End <mut_pos,1:3 , drop=F ] #germlines at the left of the mutation 
    submatrix_phased_rightside = submatrix_phased[submatrix_phased$End >mut_pos,1:3 , drop=F ] # germline at the right
    
    #order the leftside from highest to smallest position, leave the rightside from smallest to highest
    orders_phase=order(submatrix_phased_leftside["End"],decreasing=T)
    submatrix_phased_leftside=submatrix_phased_leftside[orders_phase,]
    
    
    #Germline
    #wellfraction_samelocus_germlines=wellfraction_germlines
    snpwellcount_samelocus_germlines= snpwellcount_germlines
    refwellcount_samelocus_germlines = refwellcount_germlines 
    
    
    at_least_one_good_germline=F # is there atleast one good germline kept?
    
    for(sample in tumoursamples)
    {
      samelocus_germline=c()
      phi_som=phi_cn_sample_somatic_list[sample]
      major_som=Major_cn_sample_somatic_list[sample]
      minor_som=Minor_cn_sample_somatic_list[sample]
      allelecount_som=snpwellcount_somatic[sample]
      if (is.na(allelecount_som)) #Nothing wont be done on this sample anyway withut the snp count of the somatic.
        next
      
      absence_copynumberprofile=c()
      count_lower_than_somatic<-c() # For more accuracy in case of abundance of germline, someone can choose to exclude germline having an allele count less than the somatic allele count.
      
      for (germ in rownames(submatrix_phased_leftside))
      {
        phi_germ=phi_cn_sample_PhasedGermline_df[germ , sample]
        major_germ=Major_cn_sample_PhasedGermline_df[germ , sample]
        minor_germ=Minor_cn_sample_PhasedGermline_df[germ, sample]
        allelecount_germ= snpwellcount_samelocus_germlines[germ,sample]
        
        if ( (phi_germ==phi_som || is.na(phi_germ)|| is.na(phi_som))&& 
             (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) &&
             (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          next
        }
        
        if(is.na(phi_germ) || is.na(major_germ))# We dont want NA to be considered as a breakpoint, but if the copy number information is absent we need to not consider this germline at that particular sample
          absence_copynumberprofile=c(absence_copynumberprofile,germ)
        
        if(allelecount_germ< 1 * allelecount_som & !is.na(allelecount_germ))
          count_lower_than_somatic<-c( count_lower_than_somatic, germ)
        
        
      }
      
      for (germ in rownames(submatrix_phased_rightside))
      {
        phi_germ=phi_cn_sample_PhasedGermline_df[germ , sample]
        major_germ=Major_cn_sample_PhasedGermline_df[germ , sample]
        minor_germ=Minor_cn_sample_PhasedGermline_df[germ, sample]
        allelecount_germ= snpwellcount_samelocus_germlines[germ,sample]
        
        if ( (phi_germ==phi_som || is.na(phi_germ)|| is.na(phi_som))&& 
             (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) && 
             (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          next
        }
        
        if(is.na(phi_germ) || is.na(major_germ))
          absence_copynumberprofile=c(absence_copynumberprofile,germ)
        
        if(allelecount_germ< 1 * allelecount_som & !is.na(allelecount_germ))
          count_lower_than_somatic<-c( count_lower_than_somatic, germ)
        
        
      }
      
      #if we assign NA to the alelle count information of germline mutation not on the same somatic mutation, we are removing them from the list at the considered sample
      notsamelocus_germline = setdiff(phased_list, samelocus_germline)
      germline_to_exclude_at_this_sample=c(notsamelocus_germline,absence_copynumberprofile
                                           #     ,count_lower_than_somatic
      )
      
      if(!at_least_one_good_germline)
        at_least_one_good_germline=(length(germline_to_exclude_at_this_sample)< length(phased_list))
      
      
      if(length(germline_to_exclude_at_this_sample) ==0) next
      
      
      
      #  wellfraction_samelocus_germlines[germline_to_exclude_at_this_sample,sample]= rep(NA, length(germline_to_exclude_at_this_sample))
      snpwellcount_samelocus_germlines[germline_to_exclude_at_this_sample,sample]= rep(NA, length(germline_to_exclude_at_this_sample))
      refwellcount_samelocus_germlines [germline_to_exclude_at_this_sample,sample]= rep(NA, length(germline_to_exclude_at_this_sample))
      
    }
    
    if(!  at_least_one_good_germline) next # there is no good gerline
    
    ############################
    ##### We prepare the input for the formulw
    ##############
    ############################
    
    #Allele Count
    ############
    # The somatic mutation
    lambda_somatic_list=snpwellcount_somatic[mut, ] # \lambda(S) across the samples (see paper)
    
    # The germline  mutation, approximation according to poisson distribution, see paper
    lambda_PhasedGermline_list<-colMeans(as.matrix(snpwellcount_samelocus_germlines), na.rm=T)  #  for \lambda(G)
    mu_PhasedGermline_list<-colMeans(as.matrix(refwellcount_samelocus_germlines), na.rm=T) # for \mu(G)
    omega_PhasedGermline_mean_list<-lambda_PhasedGermline_list/(lambda_PhasedGermline_list+mu_PhasedGermline_list) # for \omega(G)
    
    #Copy number profile (simply the one of the somatic mutation locus)
    ############
    phi_cn_list= unlist(phi_cn_sample_somatic_list)
    Major_cn_list= unlist(Major_cn_sample_somatic_list)
    Minor_cn_list = unlist(Minor_cn_sample_somatic_list)
    
    
    
    
    
    ##########
    #THE Allele count MODEL
    #########
    ########
    #########
    #####
    
    
    #First, we separate the samples according to their copy number profile. 
    CNV_groups<-list()
    for (phival in unique(phi_cn_list)) 
      for(majorval in unique(Major_cn_list)) 
        for (minorval in unique(Minor_cn_list))
        {
          
          if(!is.na(phival) && !is.na(majorval) && !is.na(minorval))
          {
            samplecnvgroups=c()
            
            for (sample in tumoursamples)
              if(!is.na(phi_cn_list[sample]) && !is.na(Major_cn_list[sample]) && !is.na(Minor_cn_list[sample]))
                if ((phi_cn_list[sample] ==phival) &&(Major_cn_list[sample]==majorval) && (Minor_cn_list[sample]== minorval ))
                  samplecnvgroups=c(samplecnvgroups,sample)
                if(length(samplecnvgroups)!=0)
                  CNV_groups[length(CNV_groups) +1] = list(samplecnvgroups)
          }
        }
    
    
    
    
    
    sigma_PhasedGermline_list=Major_cn_list-Major_cn_list
    tau_PhasedGermline_list= phi_cn_list * sigma_PhasedGermline_list + (1-phi_cn_list) # for \tau(G)
    hatlambda_somatic_list=lambda_somatic_list-lambda_somatic_list
    hatlambda_PhasedGermline_list=lambda_PhasedGermline_list-lambda_PhasedGermline_list
    condition_list=rep(NA,length(lambda_somatic_list))
    names(condition_list)= names(lambda_somatic_list)
    
    
    
    
    #We run the model on each group
    for (icnvgroup in 1 : length(CNV_groups))
    {
      
      sample_cnvgroup=CNV_groups[[icnvgroup]]
      
      
      # Estimate of sigma. See the paper
      #################################
      A=0
      B=0
      Asnp=0
      Bref=0
      for (sample_i in sample_cnvgroup)
      {
        
        
        for ( germlinemut in rownames(snpwellcount_samelocus_germlines) )
        {
          if (!is.na(snpwellcount_samelocus_germlines[germlinemut,sample_i]) #&& (snpwellcount_samelocus_germlines[germlinemut,sample_i]>0)
          )
            Asnp= Asnp + snpwellcount_samelocus_germlines[germlinemut,sample_i]
          if (!is.na(refwellcount_samelocus_germlines[germlinemut,sample_i]) #&& (refwellcount_samelocus_germlines[germlinemut,sample_i]>0)
          )
            Bref=Bref+refwellcount_samelocus_germlines[germlinemut,sample_i]
          
          if(!is.na(snpwellcount_samelocus_germlines[germlinemut,sample_i]) #&& (snpwellcount_samelocus_germlines[germlinemut,sample_i]>0) 
             &&
             !is.na(refwellcount_samelocus_germlines[germlinemut,sample_i]) #&& (refwellcount_samelocus_germlines[germlinemut,sample_i]>0)
          )
          {
            A= A + snpwellcount_samelocus_germlines[germlinemut,sample_i]
            B=B+refwellcount_samelocus_germlines[germlinemut,sample_i]
          }
        }
      }
      
      if((is.na(A)) || (is.na(B)))
      {
        A=Asnp
        B=Bref
      }
      if(A>=B)
        sigma_PhasedGermline_list[sample_cnvgroup] = Major_cn_list[sample_cnvgroup] 
      if(A<B)
        sigma_PhasedGermline_list[sample_cnvgroup] = Minor_cn_list[sample_cnvgroup] 
      
      
      tau_PhasedGermline_list[sample_cnvgroup]= phi_cn_list[sample_cnvgroup] * sigma_PhasedGermline_list[sample_cnvgroup] + (1-phi_cn_list[sample_cnvgroup]) # for \tau(G)
      
      
      #Allele_Count
      alpha_list=1/tau_PhasedGermline_list[sample_cnvgroup]
      beta_list=(phi_cn_list[sample_cnvgroup]  * sigma_PhasedGermline_list[sample_cnvgroup] )/ tau_PhasedGermline_list[sample_cnvgroup]
      
      #condition_list2=lambda_somatic_list - alpha_condition * lambda_PhasedGermline_list# The condition, if >=0 then case 1 else case 2
      lambda_S=as.numeric(lambda_somatic_list[sample_cnvgroup])
      lambda_G=as.matrix(snpwellcount_samelocus_germlines[sample_cnvgroup])
      alpha=as.numeric(alpha_list[sample_cnvgroup])
      beta=as.numeric(beta_list[sample_cnvgroup])
      
      EM_parameters=bestAllele(lambda_S, lambda_G, alpha, beta)
      
      
      hatlambda_somatic_list[sample_cnvgroup] =EM_parameters$hatlambda_S
      hatlambda_PhasedGermline_list[sample_cnvgroup] =EM_parameters$hatlambda_G
      condition=EM_parameters$bestC
      condition_list[sample_cnvgroup]=rep(condition, length(lambda_S))
      
      
    }
    
    
    ############################
    #####   ############################
    ##### Intermediary computation
    ##############
    ############################
    
    u_list =phi_cn_list * hatlambda_PhasedGermline_list/tau_PhasedGermline_list # for u
    v_list=(1-phi_cn_list) * hatlambda_PhasedGermline_list/tau_PhasedGermline_list # for v
    
    sum_cells_list<-u_list+v_list # estimation of the number of cells if well counts and the K * the number of cells of reads count with k = amplification factor
    sum_allele_list= hatlambda_PhasedGermline_list +  mu_PhasedGermline_list # estimation of the number of allele
    
    
    
    ############################
    ##### Prevalence computation
    ##############
    ############################
    
    #`The user should specify a minimum count of cells and allele for the prevalence to be computed. By default we take 6 cells and 8 alleles.
    min_cells=4
    min_alleles=6
    
    
    
    prev_somatic_list =vector("numeric", length=length(hatlambda_somatic_list))
    names(prev_somatic_list) = names(hatlambda_somatic_list)
    prev_somatic_list[prev_somatic_list==0]<-NA
    
    
    for(sample in tumoursamples)
    {
      
      if (is.na(condition_list[sample]))
        next
      
      if(is.na(sum_cells_list[sample]) || is.na(sum_allele_list[sample]))
        next
      if (sum_cells_list[sample]< min_cells)
        next
      if(sum_allele_list[sample]< min_alleles)
        next
      
      if (condition_list[sample]=="C2")#CONTEXTE 2
        prev_somatic_list[sample] = as.numeric(unlist(   phi_cn_list[sample] + (1-phi_cn_list[sample]) * ((hatlambda_somatic_list[sample] - u_list[sample] * sigma_PhasedGermline_list[sample])/v_list[sample])   ))
      else if (condition_list[sample]=="C1")
        prev_somatic_list[sample] = as.numeric(unlist( (hatlambda_somatic_list[sample] / hatlambda_PhasedGermline_list[sample])* tau_PhasedGermline_list[sample]    ))
      
      #if(sample=="O13_A_ABpre")
      #   exit()
      if(!is.na(prev_somatic_list[sample]) && prev_somatic_list[sample] > 1.00000001 ) exit()
      if(!is.na(prev_somatic_list[sample]) && prev_somatic_list[sample] < 0 ) exit()
      
    }
    
    
    masterprevalence[mut,names(prev_somatic_list)] = prev_somatic_list
    
    #     
    #     if (output)
    #     {
    #       #cat("\nTotal_cn_Glist :", Total_cn_Glist, " Major_cn_Glist : ", Major_cn_Glist, "  Minor_cn_Glist : ", Minor_cn_Glist, " phi_Glist : ", phi_Glist)
    #       cat(" \n\nphi_G : ", phi_cn_list)
    #       cat("\n\nsigma_G", sigma_PhasedGermline_list)
    #       cat("\n lambda_PhasedGermline_mean: ", lambda_PhasedGermline_list, " \n lambda_somatic_list: " )
    #       cat(unlist(lambda_somatic_list))
    #       cat("\n hatlambda_PhasedGermline_mean: ", hatlambda_PhasedGermline_list, " \n hatlambda_somatic_list: " )
    #       cat(unlist(hatlambda_somatic_list))
    #       cat("\n Aggregate intermediary values :")
    #       cat("\n tau_G: ", tau_PhasedGermline_list, "sigma_G:", sigma_PhasedGermline_list)
    #       cat("\n u_list: ", u_list, " \n v_list:", v_list)
    #       cat("\n Contexte :", unlist(condition_list) )
    #       # cat("\n Alpha :",alpha_condition)
    #       
    #       # cat("\n\n Prev_somatic_list2 :", prev_somatic_list2 )
    #       
    #       cat("\n prev_somatic_list :", prev_somatic_list )
    #       
    #     }
    
    
    
    masterprevalence[mut,names(prev_somatic_list)] = prev_somatic_list
    
    list_prev=unlist(prev_somatic_list)
    list_prev=list_prev[!is.na(list_prev)]
    
    
    
    
  }
  
  
  masterprevalence
  
  
}


