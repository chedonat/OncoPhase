
#' @export
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y=y[x > (qnt[1] - H)] 
  y=y[x < (qnt[2] + H)] 
  y
}


#' Somatic mutations cellular prevalence using haplotype phasing on a multi sample study.
#' 
#' This is a generic function to compute the cellular prevalence of somatic mutations in
#'  cancer using haplotype phasing.  The function applies the model to a range of mutations located at a given genomic region or at the whole genome scale. The model computes the prevalence of a somatic
#'   mutation relatively to close and eventually phased germline mutations. It uses three sources
#'    of information as input : The allelic counts, the phasing information and the 
#'    copy number alteration.  Multiple tumor samples can be provided for the prevalence computation.
#' 
#' @param snp_allelecount_df A data frame containing for each mutation the  allelic 
#' counts of the variant at each tumor samples. The data frame should contains at least the following three columns among its firsts columns: Chrom (The mutation
#'  chromosome) , End (The mutation position) and IsGermline (is the mutation a germline
#'   or somatic mutation).
#' @param ref_allelecount_df A data frame containing for each mutation the allelic count
#'  of the reference at each tumor sample. The data frame should contains at least the following three columns among its firsts columns:  Chrom (The mutation
#’ chromosome) , End (The mutation position) and IsGermline (is the mutation a Germline
#'    or Somatic mutation)  
#' @param phasing_association_df A data frame containing for each somatic mutation, 
#' a colon separated list of germline SNP phased to it.
#' @param major_copynumber_df A data frame containing for each mutation, its  major
#’ chromosomal copy number at each tumor samples.
#' @param minor_copynumber_df A data frame containing for each mutation the minor
#'  chromosomal copy number at each tumor samples.
#' @param CNVfraction_df, If provided, represents a data frame containing for each mutation,  the fraction of
#'  cells affected by the copy number alteration. Used only if the mode is "PhasedSNP" and formula is "General".
#' @param nbFirstColumns Number of first columns in snp_allelecount_df to reproduce in
#'  the output dataframe e.g: Chrom, Pos, Vartype. Columns from  nbFirstColumns +1 to the last column should contains the information needed for the prevalence computation at each tumour sample
#' @param region The region of the genome to consider for the prevalence computation  in the format chrom:start-end 
#' e.g "chr22:179800-98767 
#' @param tumoursamples : The list of tumor samples to consider for the prevalence
#’ computation.  This samples should be present as column header in the data frame
#'   snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df
#'   and  CNVfraction_df. If not provided, the headers from nbFirstColumns + 1 to 
#'   the last column of snp_allelecount_df is retrieved and its intersection with the
#' other inputted data frames headers is considered.
#' @param mode The mode under which the prevalence is computed  (default : PhasedSNP , alternatives methods  are FlankingSNP, OptimalSNP,and SNVOnly).  Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=FlankingSNP and 3 = OptimalSNP
#' #@param formula The formula used to compute the prevalence. can be either "matrix" for the linear equations or "General" for the exact allele count cases. Default : Matrix
#' @param min_cells Minimum number of cells (default 2). In case the estimated number of cells sequenced at the locus of the mutation is less than min_cells, NA is returned.
#' @param min_alleles Minimum number of alleles. (default 4). In case the estimated number of alleles sequenced at the locus of the mutation is less than min_alleles, NA is returned.
#' @param detail when set to TRUE, a detailed output is generated containing, the context and the detailed prevalence for each group of cells (germline cells, cells affected by one of the two genomic alterations (SNV or CNV) but not both, cells affected by  by both copynumber alteration and SNV ). Default : TRUE.
#' @param LocusCoverage when set to true, lambda_S and mu_S might be adjusted if necessary so that they meet the rules lambda_S <= lambda_G. mu_S >= mu_G and lambda_S + mu_S = lambda_G + mu_G. Not used if mode=SNVOnly,  Default = FALSE 
#' @param SomaticCountAdjust when set to 1, lambda_S and mu_S might be adjusted if necessary so that they meet the rules lambda_S <= lambda_G and  mu_S >= mu_G. If set to 2, in addition to the previous adjustment,  the criteria lambda_S + mu_S ~ lambda_G + mu_G should also be met. mu_S is then adjusted to lambda_G + mu_G  - mu_S when its tested to not follow Poiss(lambda_G + mu_G  - mu_S). Not used if mode=SNVOnly,  Default = 0.
#' 
#' 
#' @return A data frame containing :
#'  \describe{
#'        \item{}{Column 1 to NbFirstcolumn of the input data frame snp_allelecount_df. 
#'        This will generally include the chromosome and the position of the mutation plus
#'        any other columns to report in the prevalence dataframe (e.g REF, ALL, ...) }
#'         \item{}{One column per tumour sample reporting the prevalence of the mutation 
#'         at each samples}
#'      }
#'      
#' @examples
#' 
#' #Example 1: Loading a simple example data set with two somatic mutations, 5 germlines SNP, and 3 tumor samples
#' data(simpleExample2)
#' se=simpleExample2
#' prevalence_df=getPrevalenceMultiSamples(se$snp_allelecount_df, se$ref_allelecount_df, se$major_copynumber_df,se$minor_copynumber_df,phasing_association_df=se$phasing_association_df, )
#' print(prevalence_df)
#' 
#' #Chrom     End IsGermline  Tumour1        Tumour2        Tumour3
#' #mutation2  chr2 3003000          0 C2:0|0|1 C2:0.15|0|0.85 C2:0.12|0|0.88
#' #mutation6  chr2 4008000          0 C1:1|0|0       C1:1|0|0 C2:0|0.24|0.76
#' 
#' #Example 2 : Computing somatic mutation cellular prevalence on chromosome 15 of  patient 11152 (data retrieved from a parallel study)
#' 
#' data("chr22_1152")
#' ds=chr22_11152
#' masterprevalence_df=getPrevalenceMultiSamples(ds$snp_allelecount_df, ds$ref_allelecount_df,  ds$major_copynumber_df,ds$minor_copynumber_df,phasing_association_df = ds$phasing_association_df, cnv_fraction=ds$CNVFraction_df,nbFirstColumns=6,detail=FALSE)
#' print(head(masterprevalence_df))
#' 
#' data("chr10_11152")
#' df=chr10_11152
#' masterprevalence_df=getPrevalenceMultiSamples(df$snp_allelecount_df, df$ref_allelecount_df, df$major_copynumber_df,df$minor_copynumber_df,phasing_association_df=df$phasing_association_df, cnv_fraction=df$CNVFraction_df,nbFirstColumns=6, region="chr10:50000000-180000000")
#' print(head(masterprevalence_df))
#' 
#' 
#' 
#'@seealso \code{\link{getPrevalence}}
#' @export
getPrevalenceMultiSamples<-function(snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,mode="PhasedSNP",cnv_fraction=NULL, phasing_association_df=NULL,NormalcellContamination_df=NULL,tumoursamples=NULL,  nbFirstColumns=3, region=NULL,detail=TRUE,  LocusRadius = 10000,SameTumour=TRUE,LocusCoverage=1,ProgressOutputs=T,SomaticCountAdjust=0)
{
  
  
  
  # Extract the somatic mutations 
  
  cat("\n\n Cellular prevalence computation using OncoPhase")
  
  
  compulsory_columns=c("Chrom","End","IsGermline")
  
  if (length(setdiff(compulsory_columns,colnames(snp_allelecount_df)))>0){
    stop(" The allele count master matrices should have at least the following headers
         columns : ")
    print(compulsory_columns)
  }
  
  
  if (length(setdiff(compulsory_columns,colnames(ref_allelecount_df)))>0){
    stop(" The allele count master matrices should have at least the following 
         headers columns : ")
    print(compulsory_columns)
  }
  
  #Restriction to the region
  
  
  # If a region is provided, a restriction is performed on the given region
  if(!is.null(region)){
    
    cat("\n\n Restriction of the matrices within the region ", region,"...\n\n")
    region_parts= unlist(strsplit(region,":"))
    
    chrom = region_parts[1]
    startPosition = 1
    endPosition = hg19_dfsize[chrom]
    
    #snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,mode="PhasedSNP",cnv_fraction=NULL, phasing_association_df=NULL,NormalcellContamination_df=NULL
    snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$Chrom == chrom , ]
    ref_allelecount_df = ref_allelecount_df[ref_allelecount_df$Chrom == chrom , ]
    major_copynumber_df = major_copynumber_df[major_copynumber_df$Chrom == chrom , ]
    minor_copynumber_df = minor_copynumber_df[minor_copynumber_df$Chrom == chrom , ]
    if(!is.null(cnv_fraction))
      cnv_fraction = cnv_fraction[cnv_fraction$Chrom == chrom , ]
    if(!is.null(phasing_association_df))
      phasing_association_df = phasing_association_df[phasing_association_df$Chrom == chrom , ]
    if(!is.null(NormalcellContamination_df))
      NormalcellContamination_df = NormalcellContamination_df[NormalcellContamination_df$Chrom == chrom , ]
    
    if(length(region_parts)>1){
      coordinates = unlist(strsplit(region_parts[2],"-"))
      startPosition = as.numeric(coordinates[1])
      endPosition = as.numeric(coordinates[2])
      
      snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$Start >= startPosition & snp_allelecount_df$End <= endPosition,]
      ref_allelecount_df = ref_allelecount_df[ref_allelecount_df$Start >= startPosition & ref_allelecount_df$End <= endPosition,]
      major_copynumber_df = major_copynumber_df[major_copynumber_df$Start >= startPosition & major_copynumber_df$End <= endPosition,]
      minor_copynumber_df = minor_copynumber_df[minor_copynumber_df$Start >= startPosition & minor_copynumber_df$End <= endPosition,]
      
      if(!is.null(cnv_fraction))
        cnv_fraction = cnv_fraction[cnv_fraction$Start >= startPosition & cnv_fraction$End <= endPosition,]
      if(!is.null(phasing_association_df))
        phasing_association_df= phasing_association_df[phasing_association_df$Start >= startPosition & phasing_association_df$End <= endPosition,]
      if(!is.null(NormalcellContamination_df))
        NormalcellContamination_df = NormalcellContamination_df[NormalcellContamination_df$Start >= startPosition & NormalcellContamination_df$End <= endPosition,]
      
    }
  }
  
  
  
  
  
  
  
  if (is.null(tumoursamples)){
    tumoursamples = colnames(snp_allelecount_df[(nbFirstColumns+1):ncol(snp_allelecount_df)])
  }
  
  #print(colnames(snp_allelecount_df))
  
  
  tumoursamples = Reduce(intersect,list(tumoursamples,colnames(snp_allelecount_df),
                                        colnames(ref_allelecount_df),
                                        colnames(major_copynumber_df), 
                                        colnames(minor_copynumber_df)
  ))
  
  if(!is.null(cnv_fraction))
    tumoursamples =intersect(tumoursamples,colnames(cnv_fraction) )
  
  
  if(length(tumoursamples) ==0)
  {
    stop(" None of the tumour samples provided is present in the five  matrices :
         snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df, cnv_fraction")
  }
  
  snp_allelecount_df=numeric_column(snp_allelecount_df,tumoursamples)
  ref_allelecount_df=numeric_column(ref_allelecount_df,tumoursamples)  
  major_copynumber_df=numeric_column(major_copynumber_df,tumoursamples)
  minor_copynumber_df=numeric_column(minor_copynumber_df,tumoursamples)
  if(!is.null(cnv_fraction)) cnv_fraction=numeric_column(cnv_fraction,tumoursamples)
  
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
  
 
  
    
    #Extraction of the list of somatic mutations the cellular prevalence will be computed.
    somatic_snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$IsGermline==0, ]

    
    
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
    
    
    #Initialising the masterprevalence matrice, setting the headers
    prevalence_columns=tumoursamples
    if(detail ==1)
    {
      prevalence_columns=c()
      for(sample in tumoursamples){
        prevalence_columns=c(prevalence_columns,paste(sample,c("Prevalence","Context","Germ","Alt","Both","ResidualNorm","InputValues"),sep="_"))
      }
    }
    masterprevalence<-matrix(nrow=nrow(somatic_snp_allelecount_df),ncol=nbFirstColumns + length(prevalence_columns))
    masterprevalence<-as.data.frame(masterprevalence)
    colnames(masterprevalence) <- c(colnames(somatic_snp_allelecount_df[1:nbFirstColumns]),prevalence_columns)
    rownames(masterprevalence) <- rownames(somatic_snp_allelecount_df)
    masterprevalence[1:nbFirstColumns] = somatic_snp_allelecount_df[1:nbFirstColumns]
    
    
    
    
    # masterprevalence=masterprevalence[listover_estimated,]
    #  mut=""
    Nbmutations = nrow(masterprevalence)
    
    
    cat("\n\n Number of mutations : ", Nbmutations)
    
    if(ProgressOutputs )
    {
      TraceProgress=T
      cat(" \n A progression message giving information about the mutation under processing will be displayed each  (N/100)th mutation. Set ProgressOutputs to FALSE to turn off this. ")
      trace_step= ceiling(Nbmutations/100)
    }
    
    for (imut in 1:nrow(masterprevalence))
      #for (imut in 1:5)
    {
      
      
      #Mutation name and mutation position
      mut <- rownames(masterprevalence[imut,]); 
      mut_pos=as.numeric(masterprevalence[imut,"End"])
      
      if(imut%%trace_step==0){
        cat("\n Mutation ", mut, " ", imut, "/",Nbmutations,sep="" )
      }
      
      
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
        
        #Id mode=Optimal then in case no germline is found with PhasedSNP the search is relaunched with FlankingSNP
        
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
        
      }else   if(mode_locus=="FlankingSNP")
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
        
      }else  if(mode=="SNVOnly"){ 
        snpwellcount_germlines=NA
        refwellcount_germlines=NA
      }else{
        stop("\n\n mode shpu;d be one of SNVOnly, PhasedSNP, FlankingSNP or OptimalSNP")
      }
      
      
      # wellfraction_germlines=snp_allelecount_df[linkedGermlines_list,cifs:nbcolumns_wellfraction]
      
      ##### Preparing the input for the formulw
      
      ##Allele Counts
      lambda_somatic=snpwellcount_somatic[mut, tumoursamples,drop=F] # Somatic variant counts 
      mu_somatic=refwellcount_somatic[mut, tumoursamples,drop=F] # Somatic reference counts 
      lambda_LinkedGermline<-snpwellcount_germlines #  Germline variant counts 
      mu_LinkedGermline<- refwellcount_germlines # Germline reference counts 
      
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
      if(Trace ){
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
     # if(detail)
       # detail=2
      
      
      Trace=F
      prev_somatic=getPrevalence(lambda_somatic,mu_somatic,major_cn,minor_cn, lambda_LinkedGermline , mu_LinkedGermline,  detail ,mode_locus,Trace,SameTumour,LocusCoverage,SomaticCountAdjust=SomaticCountAdjust)
      
      
      
      if (detail!=1)
        {
        masterprevalence[mut,names(prev_somatic)] = prev_somatic
        }else{
        
        
        for(sample in tumoursamples)
        {
          prevalence=prev_somatic[[sample]]  
          if(is.null(prevalence) || is.na(prevalence)){
            masterprevalence[mut,paste(sample,"Prevalence",sep="_")] = NA
            next()
          }

          
          masterprevalence[mut,paste(sample,"Prevalence",sep="_")] = prevalence$Prevalence
          masterprevalence[mut,paste(sample,"Context",sep="_")] = prevalence$Context
          masterprevalence[mut,paste(sample,"Germ",sep="_")] = prevalence$DetailedPrevalence["Germ"]
          masterprevalence[mut,paste(sample,"Alt",sep="_")] = prevalence$DetailedPrevalence["Alt"]
          masterprevalence[mut,paste(sample,"Both",sep="_")] = prevalence$DetailedPrevalence["Both"]
          masterprevalence[mut,paste(sample,"ResidualNorm",sep="_")] = prevalence$ResidualNorm
          masterprevalence[mut,paste(sample,"InputValues",sep="_")] = prevalence$InputValues
          
          
        }
        
        
      }
      
      
      
    }
    
 
  masterprevalence
  
  }






#' Somatic mutations cellular prevalence on a single sample.
#' 
#' This is a generic function to compute the cellular prevalence of somatic mutations in
#'  cancer.  The function applies the model to a range of mutations located at a given genomic region or at the whole genome scale. The model computes the prevalence of a somatic
#'   mutation relatively to close and eventually phased germline mutations. It uses three sources
#'    of information as input : The allelic counts, the phasing information and the 
#'    copy number alteration.  Multiple tumor samples can be provided for the prevalence computation.
#' 
#' @param input_df A data frame containing for each mutations :
#' \describe{
#'        \item{lambda_S}{Alelle counts supporting the SNV}
#'        \item{mu_S}{Alelle counts supporting the reference at the SNV locus}
#'        \item{major_cn}{Major copy number at the SNV locus}
#'        \item{minor_cn}{Minor copy number at the SNV locus}
#'        \item{lambda_G}{Alelle counts supporting the SNP}
#'        \item{mu_G}{Alelle counts supporting the reference at the SNP}
#'      }
#' @param nbFirstColumns Number of first columns in snp_allelecount_df to reproduce in
#'  the output dataframe e.g: Chrom, Pos, Vartype. Columns from  nbFirstColumns +1 to the last column should contains the information needed for the prevalence computation at each tumour sample
#' @param region The region of the genome to consider for the prevalence computation  in the format chrom:start-end 
#' e.g "chr22:179800-98767 
#' @param mode The mode under which the prevalence is computed  (default : PhasedSNP , alternatives methods  are FlankingSNP, OptimalSNP,and SNVOnly).  Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=FlankingSNP and 3 = OptimalSNP
#' #@param formula The formula used to compute the prevalence. can be either "matrix" for the linear equations or "General" for the exact allele count cases. Default : Matrix
#' @param detail when set to TRUE, a detailed output is generated containing, the context and the detailed prevalence for each group of cells (germline cells, cells affected by one of the two genomic alterations (SNV or CNV) but not both, cells affected by  by both copynumber alteration and SNV ). Default : TRUE.
#' @param LocusCoverage when set to true, lambda_S and mu_S might be adjusted if necessary so that they meet the rules lambda_S <= lambda_G. mu_S >= mu_G and lambda_S + mu_S = lambda_G + mu_G. Not used if mode=SNVOnly,  Default = FALSE.
#' 
#' @param SomaticCountAdjust when set to 1, lambda_S and mu_S might be adjusted if necessary so that they meet the rules lambda_S <= lambda_G and  mu_S >= mu_G. If set to 2, in addition to the previous adjustment,  the criteria lambda_S + mu_S ~ lambda_G + mu_G should also be met. mu_S is then adjusted to lambda_G + mu_G  - mu_S when its tested to not follow Poiss(lambda_G + mu_G  - mu_S). Not used if mode=SNVOnly,  Default = 0.
#' 
#' @return A data frame containing :
#'  \describe{
#'        \item{}{Column 1 to NbFirstcolumn of the input data frame snp_allelecount_df. 
#'        This will generally include the chromosome and the position of the mutation plus
#'        any other columns to report in the prevalence dataframe (e.g REF, ALL, ...) }
#'         \item{}{and the following information}
#'   \describe{
#'        \item{Prev}{The Cellular Prevalence of the mutation}
#'        \item{Germ}{The proportion of cells with a normal genotype}
#'        \item{Alt}{The proportion of cells with only the CNA if the context C=C1 or with only the SNV if the context C=C2}
#'        \item{Both}{The proportion of cells with both the SNV and the SCNA}
#'        \item{Context}{Context at the mutation. If C1 then the SNV occured after the SCNA, if C=c2 then the SNV occured before the SCNA}
#'        \item{residual}{Residual after limSolveapproximation.}
#'      }
#'      }
#'      
#' @examples
#' 
#' #Example 1: 
#' 
#' input_file=system.file("extdata","phylogeny1_d300_n80.tsv", package = "OncoPhase")
#' input_df<-read.table(input_file,header=TRUE)
#' rownames(input_df) = input_df$mutation_id
#' print(input_df)
#' #  mut_id lambda_S mu_S major_cn minor_cn lambda_G mu_G
#' #a      a      151  152        1        1      151  135
#' #b      b      123  176        1        1      161  150
#' #c      c       94  209        2        1      176  134
#' #d      d       23  283        1        1      155  144
#' #e      e       60  228        2        0      174  125
#' 
#' prevalence_df=getPrevalenceSingleSample(input_df,nbFirstColumns = 1)
#' 
#' print(prevalence_df)
#' #  mut_id   Prev   Germ    Alt   Both Residual Context
#' #  a      a 0.9967 0.0017 0.0017 0.9967  3.1e-03      C1
#' #  b      b 0.8230 0.0890 0.0890 0.8230  1.3e-03      C1
#' #  c      c 0.4010 0.6000 0.0910 0.3100  3.9e-33      C2
#' #  d      d 0.1500 0.4200 0.4200 0.1500  1.4e-03      C1
#' #  e      e 0.2490 0.7500 0.0890 0.1600  5.1e-31      C2
#' 
#' 
#'@seealso \code{\link{getPrevalence}}
#' @export
getPrevalenceSingleSample<-function(input_df,mode="PhasedSNP",  nbFirstColumns=0, region=NULL,detail=TRUE,  LocusCoverage=1,SomaticCountAdjust=0)
{
  
#Check the compulsory columns
  compulsory_columns=c("varcounts_snv","refcounts_snv","major_cn","minor_cn")
  
  if (length(setdiff(compulsory_columns,colnames(input_df)))>0){
    stop(" The allele count master matrices should have at least the following headers
         columns : ", compulsory_columns)
   # print(compulsory_columns)
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
  
  
  #Prepare the master matrices for the prevalence.
  masterprevalence=as.data.frame(matrix(nrow=nrow(input_df), ncol=nbFirstColumns+6))
  rownames(masterprevalence) = rownames(input_df)
  if(nbFirstColumns>0){
    masterprevalence[1:nbFirstColumns] = input_df[1:nbFirstColumns]
    colnames(masterprevalence) = c(colnames(input_df[1:nbFirstColumns]),"Prevalence", "Germ","Alt","Both","Residual","Context")
  }else{
    colnames(masterprevalence) = c("Prevalence", "Germ","Alt","Both","Residual","Context")
  }
  
  
  for(imut in 1: nrow(input_df))
  {
    lambda_S = input_df[imut,"varcounts_snv"]
    mu_S = input_df[imut,"refcounts_snv"]
    minor_cn=input_df[imut,"minor_cn"]
    major_cn=input_df[imut,"major_cn"]
    lambda_G=input_df[imut,"varcounts_snp"]
    mu_G=input_df[imut,"refcounts_snp"]
    
    InputValues=paste(lambda_S,mu_S,minor_cn,major_cn,lambda_G,mu_G,sep=":")
    prevalence= getPrevalence(lambda_S, mu_S, major_cn, minor_cn, lambda_G, mu_G, detail=T, mode=mode ,LocusCoverage=LocusCoverage, SomaticCountAdjust=SomaticCountAdjust)
    
    if(is.null(prevalence) ||  length(prevalence)==0 || is.na(prevalence[[1]]) )
      next
    
    masterprevalence[imut,"Prevalence"] = prevalence[[1]]$Prevalence
    masterprevalence[imut,"Germ"] = prevalence[[1]]$DetailedPrevalence["Germ"]
    masterprevalence[imut,"Alt"] = prevalence[[1]]$DetailedPrevalence["Alt"]
    masterprevalence[imut,"Both"] = prevalence[[1]]$DetailedPrevalence["Both"]
    masterprevalence[imut,"Context"] = prevalence[[1]]$Context
    masterprevalence[imut,"Residual"] = prevalence[[1]]$ResidualNorm
    masterprevalence[imut,"InputValues"] = InputValues
  }
  
  
  
  masterprevalence
  
}







#' Compute cellular prevalence at a single mutation point
#' 
#' This is a generic function to compute the cellular prevalence of a somatic mutation point (SNV) using one of the following mode :
#'  \describe{
#'        \item{PhasedSNP}{ The prevalence is computed relatively to a Phased SNP}
#'        \item{FlankingSNP}{ The prevalence is computed relatively to a neighbour SNP located on the same locus with the somatic SNV. NA is returned if the prevalence can not be resolved without knowing the phasing information between the SNP and the SNV.}
#'        \item{SNVOnly} {The prevalence is computed using only the SNV information without the usage of any nearby SNP}
#'     }
#'     The methods is  particularly highly accurate when the locus of the SNP is also affected by a somatic copy number alteration (SCNA). The method detect the temporal relationship between the two alterations (C1: SNV occured afterthe SCNA ; C2: SNV occured before the SCNA)  and  compute the detailed prevalence of each of the following group of cells (if detail is set to TRUE) : 
#'     #'  \describe{
#'        \item{Germ}{ Cells having a germline genotype  at the locus of the SNV. That is No SNV, no SCNA}
#'        \item{Alt}{ Cells having one alternative of  the two somatic alteration. That is either the SCNA, either the SNV not both.}
#'        \item{Germ}{ Cells having both somatic alterations. That is the SNV and the SCNA}
#'     }
#' 
#' @param lambda_S : A count (or a vector of counts  if multiple samples ) of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param mu_S : A count (or a vector of counts  if multiple samples ) of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn  major copy number (or a vector  if multiple samples ) 
#' at the locus of the mutation
#' @param minor_cn : minor copy number (or a vector  if multiple samples)  
#'  at the locus of the mutation 
#' @param lambda_G  A count (or a vector of counts  if multiple samples) of 
#' alleles supporting the variant sequence of  the Germline SNP
#' @param mu_G  A count (or a vector of counts  if multiple samples) of 
#' alleles  supporting the reference sequence  of  the Germline SNP
#' @param mode The mode under which the prevalence is computed  (default : PhasedSNP , alternatives methods  are FlankingSNP and SNVOnly). Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=FlankingSNP
#' @param detail when set to 0, the function simply output the cellular prevalence of the somatic mutation. if set to 1,  a detailed output is generated containing:
#'  \describe{
#'        \item{Context}{ The associated context (C=1 or C=2) }
#'        \item{Prevalence}{The computed cellular prevalence}
#'        \item{DetailedPrevvalence}{the detailed prevalence for each group of cells (germline cells (Germ), cells affected by one of the two genomic alterations (Alt), cells affected  by both genomic alterations (Both)}
#'        \item{ResidualNorm}{Norm of the residuals when an approximation is performed}
#'        \item{CondensedPrevalence}{A colon separated list of the above fields (Context, Prevalence, Detailedprevalence and ResidualNorm). The detailed prevalence are separated by "|"}
#' }
#'   if detail is set to 2, the function outputs the CondensedPrevalence field above.     
#' @param Trace if set to TRUE, print the trace of the computation.    
#'  
#' @return   The cellular prevalence if detail =0, a detailed output if detail = 1, and a condensed output if detail =2. See the usage of the parameter detail above.
#' @param LocusCoverage when set to true, lambda_S and mu_S might be adjusted if necessary so that they meet the rules lambda_S <= lambda_G. mu_S >= mu_G and lambda_S + mu_S = lambda_G + mu_G. Not used if mode=SNVOnly,  Default = FALSE 
#' @param SomaticCountAdjust when set to 1, lambda_S and mu_S might be adjusted if necessary so that they meet the rules lambda_S <= lambda_G and  mu_S >= mu_G. If set to 2, in addition to the previous adjustment,  the criteria lambda_S + mu_S ~ lambda_G + mu_G should also be met. mu_S is then adjusted to lambda_G + mu_G  - mu_S when its tested to not follow Poiss(lambda_G + mu_G  - mu_S). Not used if mode=SNVOnly,  Default = 0.
#'      
#'     
#'      
#' @examples
#' #Example 1
#' prevalence=getPrevalence(lambda_S=14,mu_S=10,major_cn=3,minor_cn=1,lambda_G=16, mu_G=8  )
#' print(prevalence)
#' #Sample_1 
#' #0.75   
#' #The above example  gives the same  prevalence with modeFlankingSNP but not with mode SNVOnly
#' prevalence=getPrevalence(lambda_S=14,mu_S=10,major_cn=3,minor_cn=1,lambda_G=16, mu_G=8 ,mode="FlankingSNP")
#' print(prevalence)
#' #Sample_1 
#' #0.75 
#' prevalence=getPrevalence(lambda_S=14,mu_S=10,major_cn=3,minor_cn=1,lambda_G=16, mu_G=8 ,mode="SNVOnly")
#' print(prevalence)
#' #Sample_1 
#' #0.79
#'  
#' #Example 2 Case Study A (see paper)
#' prevalence = getPrevalence(lambda_S=6,mu_S=8,major_cn=2,minor_cn=1,lambda_G=8, mu_G=6, detail=1)
#' print(prevalence)
#' #$Sample_1
#' #$Sample_1$Context
#' #[1] "C2"
#' #
#' #$Sample_1$Prevalence
#' #[1] 0.66
#' #
#' #$Sample_1$DetailedPrevalence
#' #Germ  Alt Both 
#' #0.33 0.33 0.33 
#' #
#' #$Sample_1$ResidualNorm
#' #residual 
#' #2e-31 
#' #
#' #$Sample_1$CondensedPrevalence
#' #[1] "C2:0.66:0.33|0.33|0.33:2e-31"
#'   
#' ## The above example is not resolvable without phasing information
#'  prevalence =  getPrevalence(lambda_S=6,mu_S=8,major_cn=2,minor_cn=1,lambda_G=8, mu_G=6, mode="FlankingSNP",detail=TRUE)
#'  print(prevalence)
#' #Warning message:
#' #  In getFlankingSNPPrevalence(lambda_S, mu_S, major_cn, minor_cn,  :
#' #                                The prevalence is not resolved without the knowledge of the Phased Germline
#' #NA
#'## it gives an inaccurate  prevalence under mode "SNVOnly"    
#' prevalence= getPrevalence(lambda_S=6,mu_S=8,major_cn=2,minor_cn=1,lambda_G=8, mu_G=6, mode="SNVOnly",detail=0)
#' print(prevalence)
#' #Sample_1 
#' #1 
#' #Example 3 Case Study B (see paper) Not resolvable without phasing information
#' prevalence=getPrevalence(lambda_S=4,mu_S=8,major_cn=2,minor_cn=0,lambda_G=8, mu_G=4, detail=2)
#' print(prevalence)
#' #Sample_1 
#' # "C2:0.33:0.67|0|0.Prevalence(lambda_S=6,mu_S=8,major_cn=2,minor_cn=1,lambda_G=8, mu_G=6, detail=TRUE,Trace=TRUE )33:5.2e-32" 
#' 
#' #Example 4 Case Study A (see paper) Not resolvable without phasing information
#' prevalence = getPrevalence(lambda_S=6,mu_S=8,major_cn=2,minor_cn=1,lambda_G=8, mu_G=6 )
#' print(prevalence)      
#' #Sample_1 
#' #0.66  
#' 
#' #Example 5 We group case study A, B and C above to form a multisample case
#' prevalence=getPrevalence(lambda_S=c(6,4,6),mu_S=c(8,8,14),major_cn=c(2,2,2),minor_cn=c(1,0,1),lambda_G=c(8,8,8), mu_G=c(6,4,12) )                   
#' print(prevalence)
#' #Sample_1 Sample_2 Sample_3 
#' #0.66     0.33     0.75   
#' 
#' @seealso \code{\link{getPrevalence}},  \code{\link{getPhasedSNPPrevalenceGeneral}},   \code{\link{getPrevalenceLinear}}, \code{\link{getPrevalenceSNVOnly}}                                                 
#' @export
getPrevalence<-function(lambda_S,mu_S,major_cn,minor_cn, lambda_G=NULL, mu_G=NULL,  detail=0, mode="PhasedSNP",Trace=FALSE,SameTumour=TRUE,LocusCoverage=1,SomaticCountAdjust=0)
  {
  

  
  
  N=length(lambda_S) # Number of samples
  
  if((length(mu_S)!=N) || 
     (!is.null(lambda_G) && (mode !="SNVOnly") &&((length(lambda_G) !=N) || (length(mu_G) !=N))) ||
     (length(major_cn) !=N) || (length(minor_cn)!=N) )
    stop("\n\nThe vectors passed as input should have the same size\n\n")
  
  #set the mode if numeric, 0=SNVOnly, 1 = PhasedSNP, 2=FlankingSNP
  numeric_mode=c("SNVOnly", "PhasedSNP","FlankingSNP")
  if(is.numeric(mode))
  {
    if(mode %in% c(0,1,2))
    {
      mode = numeric_mode[mode +1 ]
    }else{
      stop("\n\n Mode parameter, if numeric,  should be either 0, 1 or 2")
    }
  }
  
  

  
 
  if(mode=="PhasedSNP"){
    prev_somatic=getPhasedSNPPrevalence( lambda_S,mu_S , major_cn,minor_cn, lambda_G , mu_G,detail,Trace=Trace,LocusCoverage,SomaticCountAdjust)
  }else if(mode=="FlankingSNP"){
    prev_somatic=getFlankingSNPPrevalence( lambda_S,mu_S ,  major_cn,minor_cn,lambda_G , mu_G, detail,Trace=Trace,SameTumour,LocusCoverage,SomaticCountAdjust)
  }else if(mode=="SNVOnly"){
    prev_somatic=getSNVOnlyPrevalence(lambda_S,mu_S ,major_cn,minor_cn, detail, Trace=Trace )
  }else {
    stop("parameter mode should be either FlankingSNP, PhasedSNP or SNVOnly")
  }
  
  
  
  
  
  prev_somatic
}

#' @export
getPhasedSNPPrevalence<-function( lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G, detail=0,Trace=FALSE,LocusCoverage=1,SomaticCountAdjust=0)
  {
  
  
  
    #if(length(lambda_G)>1)
    {
      tumoursamples= names(lambda_G)
      if (is.null(tumoursamples)){
        tumoursamples = paste("Sample",c(1:length(lambda_G)),sep="_")
      }
      
      #To avoid some side error, we transform them into vector
      lambda_G=as.vector(lambda_G)
      mu_G=as.vector(mu_G)
      lambda_S=as.vector(lambda_S)
      mu_S = as.vector(mu_S)
      major_cn = as.vector(major_cn)
      minor_cn=as.vector(minor_cn)
      
      names(lambda_G) =  tumoursamples
      names(mu_G) =  tumoursamples
      names(lambda_S) =  tumoursamples
      names(mu_S) =  tumoursamples
      names(major_cn) =  tumoursamples
      names(minor_cn) =  tumoursamples
      
    }
  

  #Initialisation of prev_S
  if(detail==1)
  { prev_S = list()
  }else{
    prev_S =vector("numeric", length=length(lambda_G))
    names(prev_S) = names(lambda_G)
    prev_S[prev_S==0]<-NA  
    }

    
    for(sample in tumoursamples)
    {
      
    #  Trace=(sample=="ABpre")
     # if(Trace) cat("\n\n\n Computing the prevalence on sample :", sample)
      args_list=list(
        lambda_S=lambda_S[sample],
        mu_S=mu_S[sample],
        major_cn=major_cn[sample],
        minor_cn=minor_cn[sample],
        lambda_G=lambda_G[sample],
        mu_G=mu_G[sample],#/ omega_G[sample] - lambda_G[sample],
        detail=1,
        Trace=Trace,
        LocusCoverage=LocusCoverage,
        SomaticCountAdjust=SomaticCountAdjust)
      
      if(anyNA(args_list))
        next
      if(lambda_G[sample]+mu_G[sample] ==0)
        next
      if(lambda_S[sample]+ mu_S[sample] ==0)
        next
      
      prevalence=do.call(getPhasedSNPPrevalence_on_singlemutation, args_list)
      
     # print(prevalence)
      
      #  if detail, the context and  tree type of prevalence are collapsed else only the prevalence is outputed
      
      prev_S[sample] =NA
      
     # if(length(prevalence) > 0 && !is.na(prevalence))
     #   {
        if(detail==2){
          if(length(prevalence) > 0)
          prev_S[sample] = prevalence$CondensedPrevalence
        } else if(detail==1){
          if(length(prevalence) > 0)
          prev_S[sample] =list(prevalence)
        }else{
          if(length(prevalence) > 0)
          prev_S[sample] =prevalence$Prevalence
        }
     # }else{
     #   prev_S[sample] =NA
     # }

    }
    
    if(Trace) {cat("\n\n\n Computed prevalences are  :\n")
    print(prev_S)}
    
    
    prev_S
  }
  
  

#' @export
getFlankingSNPPrevalence<-function( lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G, detail=FALSE,Trace=False,SameTumour=TRUE,LocusCoverage=1,SomaticCountAdjust=0)
  {
  
  #For each case, we compute the prevalence twice. By considering the somatic to be phased to the germline SNP or phased with the alternative chromosome harboring the reference of the Germline. The latter is achieved just by 
 
    prevalence_phasedSNP = getPhasedSNPPrevalence(lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G, detail=1,Trace=0,LocusCoverage=LocusCoverage,SomaticCountAdjust=SomaticCountAdjust)
     prevalence_phasedREF= getPhasedSNPPrevalence(lambda_S,mu_S,major_cn,minor_cn, mu_G,lambda_G, detail=1, Trace=0,LocusCoverage=LocusCoverage,SomaticCountAdjust=SomaticCountAdjust)
     
     if(Trace){
       cat("\n\n Prevalence case phased with SNP \n ")
       print( prevalence_phasedSNP)
       cat("\n Prevalence case phased with REF \n ")
       print( prevalence_phasedREF)
     }
     
     
     
     #Initialisation
     prevalence_flankingSNP=prevalence_phasedSNP
     
     
     # we now have choose the prevalence having the highest residual, if Same tumour we sum the residual across all the sample and we compare, else We go sample by sample.
     #if both have the same residual, we check into sample prevalence if they have the same prevalence, in that case the value of that prevalence will be assigned. Otherwise NA will be assigned.
     
     
     if(SameTumour)
     {
       #We sum the residuals and we choose the prevalence with the lower residuals.
       SNP_residual=sum(as.numeric(unlist(lapply(prevalence_phasedSNP, function(x)  if(!is.na(x[["ResidualNorm"]])) as.numeric(x[["ResidualNorm"]]) else NA ))) ,na.rm=T)
       REF_residual=sum(as.numeric(unlist(lapply( prevalence_phasedREF, function(x)  if(!is.na(prevalence_phasedREF[1])) as.numeric(x[["ResidualNorm"]]) else NA ))),na.rm=T)
       
       if(SNP_residual<1e-8)
         SNP_residual=0
       if(REF_residual<1e-8)
         REF_residual=0
       
       if(SNP_residual<REF_residual){
         prevalence_flankingSNP=prevalence_phasedSNP
       }else if(SNP_residual>REF_residual) {
         prevalence_flankingSNP=prevalence_phasedREF
       } else{
         for(sample in names(prevalence_flankingSNP))
           if (prevalence_phasedSNP[[sample]]$Prevalence != prevalence_phasedREF[[sample]]$Prevalence)
             prevalence_flankingSNP[[sample]] =NA 
       }
       
       
     }else{
       for(sample in names(prevalence_flankingSNP))
       {
         #First any residual less than 1e-8 is considered to be 0.
         
         if(prevalence_phasedSNP[[sample]]$ResidualNorm <1e-8)
           prevalence_phasedSNP[[sample]]$ResidualNorm=0
         if(prevalence_phasedREF[[sample]]$ResidualNorm <1e-8)
           prevalence_phasedREF[[sample]]$ResidualNorm=0  
         
         if(prevalence_phasedSNP[[sample]]$ResidualNorm < prevalence_phasedREF[[sample]]$ResidualNorm){
           prevalence_flankingSNP[[sample]] = prevalence_phasedSNP[[sample]]  
         } else if (prevalence_phasedSNP[[sample]]$ResidualNorm > prevalence_phasedREF[[sample]]$ResidualNorm) {
           prevalence_flankingSNP[[sample]] = prevalence_phasedREF[[sample]]  
         }else if (prevalence_phasedSNP[[sample]]$Prevalence != prevalence_phasedREF[[sample]]$Prevalence){
           prevalence_flankingSNP[[sample]] =NA
           
           warning("The prevalence is not resolved without the knowledge of the Phased Germline at this sample")
         }
       }
     }
     
     
     
     
     
     
     
       if(detail==1){
         prev_S= prevalence_flankingSNP
       }else{
         prev_S =vector("numeric", length=length(prevalence_flankingSNP))
         names(prev_S) = names(prevalence_flankingSNP)
         prev_S[prev_S==0]<-NA 
         for(sample in names(prevalence_flankingSNP))
         {
           if(length(prevalence_flankingSNP[[sample]])==1)
             if(is.na(prevalence_flankingSNP[[sample]]))
             next
           
           if(detail==2){
             prev_S[sample] = prevalence_flankingSNP[[sample]]$CondensedPrevalence             
           }else{
             prev_S[sample] = prevalence_flankingSNP[[sample]]$Prevalence               
           }
         }
       }
       
  if(Trace){
        cat("\n Final prevalence \n ")
    print( prev_S)
  }
  
 prev_S
}



#' @export
getPhasedSNPPrevalence_on_singlemutation<-function(lambda_S,mu_S,major_cn,minor_cn, lambda_G, mu_G, detail=FALSE,Trace=FALSE,LocusCoverage=1,SomaticCountAdjust=0)
  {
  
  
 
  

  Prevalence=NA
  DetailedPrevvalence=NA
  
  if(SomaticCountAdjust>0){
      if(lambda_S > lambda_G)
        lambda_S=lambda_G
      if(mu_S < mu_G)
        mu_S = mu_G 
    }
  
  if(SomaticCountAdjust==2){
    if(mu_S < qpois(0.001,max(0,mu_G + lambda_G - lambda_S)))
      mu_S=mu_G + lambda_G - lambda_S
  }

    
    
    
    #We compute the prevalence for the two contexts and we choose the one with the less residual
    if(Trace) cat("\n\n\n Context : C1 (C=0)  SNV after CNA \n **********")
    PrevalenceCond_C1 = getPrevalenceLinear(lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G,"C1",Trace,LocusCoverage)
    if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA \n **********")
    PrevalenceCond_C2 = getPrevalenceLinear(lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G,"C2",Trace,LocusCoverage)
    
    
    
    PrevalenceCond = PrevalenceCond_C1
    context="C1"
    if(as.numeric(PrevalenceCond_C2["residual"]) < as.numeric(PrevalenceCond_C1["residual"])){
      PrevalenceCond = PrevalenceCond_C2
      context="C2"
    }
    
    
    PrevalenceCond=as.numeric(format(PrevalenceCond,digits=2))
    names(PrevalenceCond) = names(PrevalenceCond_C2)
    
    
    
    AllPrevalences=PrevalenceCond[1:3]
    if(context=="C2"){
      Prevalence=sum(PrevalenceCond["Alt"],PrevalenceCond["Both"],na.rm=T)
    }
    if(context=="C1"){
      Prevalence=PrevalenceCond["Both"]
    }
    
    residualNorm = PrevalenceCond["residual"]
    
    condensedPrevalence=paste( context,Prevalence,paste(AllPrevalences,collapse="|"),residualNorm, sep=":")
    lambda_G=as.numeric(format(round(lambda_G, 2), nsmall = 2))
    mu_G=as.numeric(format(round(mu_G, 2), nsmall = 2))
    input_values = paste(lambda_S,mu_S,major_cn,minor_cn, lambda_G, mu_G,sep=":")
    
    if(detail){
      Prevalence_output = list(Context=context,Prevalence=Prevalence,DetailedPrevalence=AllPrevalences,ResidualNorm=residualNorm, CondensedPrevalence = condensedPrevalence,InputValues=input_values)
    }else{
      Prevalence_output = Prevalence
    }
    
    if(Trace) 
    {
      cat("\n\n\n\t\t ***** Final Prevalence is \n")
      print(Prevalence_output)
    }

  
  
  
  
  
  Prevalence_output
  
}



#' @export
getSNVOnlyPrevalence<-function(lambda_S,mu_S,major_cn,minor_cn, detail=FALSE,Trace=FALSE)
{
  
    tumoursamples= colnames(lambda_S)
  if(is.null(tumoursamples))
    tumoursamples= names(lambda_S)   


  if (is.null(tumoursamples))
    tumoursamples = paste("Sample",c(1:length(lambda_S)),sep="_")
 
  
  
  
  #To avoid some side error, we transform them into vector
  lambda_S=as.vector(lambda_S)
  mu_S = as.vector(mu_S)
  major_cn = as.vector(major_cn)
  minor_cn=as.vector(minor_cn)
  
  names(lambda_S) =  tumoursamples
  names(mu_S) =  tumoursamples
  names(major_cn) =  tumoursamples
  names(minor_cn) =  tumoursamples

  
  
  
  if(detail==1)
  { prev_S = list()
  }else{
    prev_S =vector("numeric", length=length(lambda_S))
    names(prev_S) = names(lambda_S)
    prev_S[prev_S==0]<-NA  
  }
  
  
  
  for(sample in tumoursamples)
  {
    if(Trace) cat("\n\n\n Computing the prevalence on sample :", sample)
    args_list=list(
      lambda_S=lambda_S[sample],
      mu_S=mu_S[sample],
      major_cn=major_cn[sample],
      minor_cn=minor_cn[sample],
      detail=1,
      Trace=Trace)
    
    if(anyNA(args_list))
      next
    if(lambda_S[sample]+ mu_S[sample] ==0)
      next
    
    prevalence=do.call(getSNVOnlyPrevalence_on_singlemutation, args_list)
 
    # print(prevalence)
    
    #  if detail, the context and  tree type of prevalence are collapsed else only the prevalence is outputed
    if(detail==2){
      prev_S[sample] = prevalence$CondensedPrevalence
    } else if(detail==1){
      prev_S[sample] =list(prevalence)
    }else{
      prev_S[sample] =prevalence$Prevalence
    }
  }
  
  if(Trace) {cat("\n\n\n Computed prevalences are  :\n")
    print(prev_S)}
  
  
  prev_S
}


#' @export
getSNVOnlyPrevalence_on_singlemutation<-function(lambda_S,mu_S,major_cn,minor_cn, detail=FALSE,Trace=FALSE)
  {
  

  Prevalence=NA
  DetailedPrevvalence=NA
  sigma=1
  if(Trace) cat("\n The parameters are : ", c(lambda_S,mu_S,major_cn,minor_cn))
  
  
  #We compute the prevalence for the two contexts and we choose the one with the less residual
  if(Trace) cat("\n\n\n Context : C1 (C=0)  SNV after CNA \n **********")
  PrevalenceCond_C1 = getPrevalenceSNVOnly(lambda_S,mu_S,major_cn,minor_cn,"C1")
  if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA and sigma=major copy number \n **********")
  PrevalenceCond_C2_major = getPrevalenceSNVOnly(lambda_S,mu_S,major_cn,minor_cn,"C2",major_cn)
  if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA and sigma=minor copy number \n **********")
  PrevalenceCond_C2_minor = getPrevalenceSNVOnly(lambda_S,mu_S,major_cn,minor_cn,"C2",minor_cn)
  
  PrevalenceCond = PrevalenceCond_C1
  context="C1"
  if(min(as.numeric(PrevalenceCond_C2_major["residual"]),as.numeric(PrevalenceCond_C2_minor["residual"]) ) < as.numeric(PrevalenceCond_C1["residual"])){

    context="C2"
    PrevalenceCond = PrevalenceCond_C2_major
    sigma=major_cn
    if(as.numeric(PrevalenceCond_C2_minor["residual"]) < as.numeric(PrevalenceCond_C2_major["residual"])){
      PrevalenceCond = PrevalenceCond_C2_minor
      sigma=minor_cn
    }
    
  }
  
  
  PrevalenceCond=as.numeric(format(PrevalenceCond,digits=2))
  names(PrevalenceCond) = names(PrevalenceCond_C1)

  if(Trace){
    cat("\n Prevalence: \n")
    print(PrevalenceCond)
  }
  
  
  AllPrevalences=PrevalenceCond[1:3]
  
  
  if(context=="C2"){
    Prevalence=sum(PrevalenceCond["Alt"],PrevalenceCond["Both"],na.rm=T)
  }
  if(context=="C1"){
    Prevalence=PrevalenceCond["Both"]
  }
  
  residualNorm = PrevalenceCond["residual"]
  
  
  condensedPrevalence=paste( context,Prevalence,paste(AllPrevalences,collapse="|"),residualNorm, sep=":")
  
  
  input_values = paste(lambda_S,mu_S,major_cn,minor_cn,sep=":")
  
  if(detail){
    Prevalence_output = list(Context=context,Prevalence=Prevalence,DetailedPrevalence=AllPrevalences,ResidualNorm=residualNorm, CondensedPrevalence = condensedPrevalence,InputValues=input_values)
  }else{
    Prevalence_output = Prevalence
  }
  
  
  if(Trace) {cat("\n\n\n\t\t ***** Final Prevalence is \n")
  print(Prevalence_output)
  }
  Prevalence_output
  
  
  
  
}





#' Generate the matrices C, W and M from a set of parameters.
#' 
#' This is a generic function to generate the matrices of the linear system (see the paper) from the allele counts and the copy number information.
#' 
#' 
#' @param lambda_S : A count of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param mu_S : A count of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn: major copy number  at the locus of the mutation
#' @param minor_cn : minor copy number   at the locus of the mutation 
#' @param lambda_G  A count of 
#' alleles supporting the variant sequence of  the Germline SNP
#' @param mu_G : A count of 
#' alleles  supporting the reference sequence  of  the Germline SNP
#' @param context represents either the situation of a mutation which occurred after the CNV ("C1") or the context of a mutation which occurred before the CNV ("C2"). If not provided, the right context will be estimated from the input
#'  
#'  
#' @return   the matrices W, C and M for the linear system of prevalence computation.
#'      
#'     
#'      
#' @examples
#' 
#' Matrices = getMatrices(8, 5,2,1,3,10,"C1")
#'  
#'  print(Matrices)
#' # $context
#' # [1] "C1"
#' # 
#' # $W
#' # SNP       SNV
#' # SNP 0.6153846 0.0000000
#' # SNV 0.0000000 0.2307692
#' # 
#' # $M
#' # Germ Alt Both
#' # SNP    1   2    2
#' # SNV    0   0    1
#' # 
#' # $C
#' # Germ Alt Both
#' # SNP    2   3    3
#' # SNV    2   3    3
#'  
#' @export
getMatrices<-function(lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G,context,LocusCoverage=1){
  total_cn = major_cn + minor_cn
  
  if((LocusCoverage>0)){
    if(LocusCoverage==1){
      Locus_Coverage= mu_G+lambda_G
    }else if (LocusCoverage==2){
      Locus_Coverage=(mu_G+lambda_G + mu_S +lambda_S) /2
    }

    omega_G = min(1, lambda_G/Locus_Coverage)
    omega_S= min(1, lambda_S/Locus_Coverage)
  }else {
    if (LocusCoverage!=0)
      warnings("\n\n LocusCoverage should be one of 0, 1 or 2. 0 will be considered")
    omega_G = lambda_G/(mu_G+lambda_G)
    omega_S= lambda_S/(mu_S +lambda_S) 
  }
  
  W=matrix(c(omega_G,0,0,omega_S),ncol=2,nrow=2)
  colnames(W)= c("SNP","SNV")
  rownames(W) = c("SNP","SNV")
  
  C=matrix(nrow=2,ncol=3)
  colnames(C) = c("Germ","Alt","Both")
  rownames(C) = c("SNP","SNV")
  M=C
  
  if(major_cn<minor_cn)
  {
    stop("\n The major copy number ", major_cn, " can not be less than the minor copy number", minor_cn)
  }
  #if the germline VAF is < 0.5 then sigma = major_Cn else minor_CN
  if(omega_G > 0.5) sigma =major_cn
  if(omega_G <= 0.5) sigma =minor_cn
  
  if(context=="C1"){
    M["SNP",] = c(1,sigma,sigma)
    M["SNV",] = c(0,0,1)
    C["SNP",] = c(2,total_cn,total_cn)
    C["SNV",] = c(2,total_cn,total_cn)
  }else{
    M["SNP",] = c(1,1,sigma)
    M["SNV",] = c(0,1,sigma)
    C["SNP",] = c(2,2,total_cn)
    C["SNV",] = c(2,2,total_cn)
  }
  
  list(context=context,W=W,M=M,C=C)
}




#' Generate the matrices C, W and M from a set of parameters under the mode "SNVOnly.
#' 
#' This is a generic function to generate the matrices of the linear system (see the paper) from the allele counts and the copy number information under the SNVOnly mode.
#' 
#' 
#' @param lambda_S : A count of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param mu_S : A count of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn: major copy number  at the locus of the mutation
#' @param minor_cn : minor copy number   at the locus of the mutation 
#' @param context represents either the situation of a mutation which occurred after the CNV ("C1") or the context of a mutation which occurred before the CNV ("C2"). If not provided, the right context will be estimated from the input
#'  
#'  
#' @return   the matrices W, C and M for the linear system of prevalence computation.
#'      
#'     
#'      
#' @examples
#' 
#' Matrices = getMatricesSNVOnly(3,10,2,1,"C1")
#'  
#'  print(Matrices)
#' #$context
#' #[1] "C1"
#' #
#' #$W
#' #SNV
#' #SNV 0.2307692
#' #
#' #$M
#' #Germ Alt Both
#' #SNV    0   0    1
#' #
#' #$C
#' #Germ Alt Both
#' #SNV    2   3    3
#'  
#' @export
getMatricesSNVOnly<-function(lambda_S,mu_S,major_cn,minor_cn,context, sigma=NULL){
  total_cn = major_cn + minor_cn
  #omega_G = lambda_G/(mu_G+lambda_G)
  omega_S= lambda_S/(mu_S +lambda_S) 
  W=matrix(c(omega_S),ncol=1,nrow=1)
  colnames(W)= c("SNV")
  rownames(W) = c("SNV")
  
  C=matrix(nrow=1,ncol=3)
  colnames(C) = c("Germ","Alt","Both")
  rownames(C) = c("SNV")
  M=C
  
  if(major_cn<minor_cn)
  {
    stop("\n The major copy number ", major_cn, " can not be less than the minor copy number", minor_cn)
  }
  
  
  #if the germline VAF is < 0.5 then sigma = major_Cn else minor_CN
 # if(omega_G > 0.5) sigma =major_cn
 # if(omega_G <= 0.5) sigma =minor_cn
  
  if(context=="C1"){
   # M["SNP",] = c(1,sigma,sigma)
    M["SNV",] = c(0,0,1)
    #C["SNP",] = c(2,total_cn,total_cn)
    C["SNV",] = c(2,total_cn,total_cn)
  }else{
    if (is.null(sigma)){
      stop("\n\n If the context is C_2, the value of sigma is required. Please provide it or set the context to C1\n\n")
    }else{
      if((sigma!= major_cn) && (sigma != minor_cn))
      {
        stop("\n\n sigma should be either the major copy number either the minor copy number. Please check your parameters and try again\n\n")
      }
    }

    
    #M["SNP",] = c(1,1,sigma)
    M["SNV",] = c(0,1,sigma)
   # C["SNP",] = c(2,2,total_cn)
    C["SNV",] = c(2,2,total_cn)
  }
  
  list(context=context,W=W,M=M,C=C)
}




#' Compute the cellular prevalence of each group of cells
#' 
#' This is a generic function to compute the detailed prevalence of a single mutation using the linear system of  the model.
#' 
#' 
#' @param lambda_G : A count of 
#' alleles supporting the variant sequence of  the Germline SNP
#' @param mu_G : A count of 
#' alleles  supporting the reference sequence  of  the Germline SNP
#' @param lambda_S : A count of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param mu_S : A count of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn: major copy number at the locus of the mutation
#' @param minor_cn : minor copy number (or a vector of copy number if multiple tumor samples)
#’ at the locus of the mutation 
#' @param context represents either the situation of a mutation which occurred after the CNV ("C1") or the context of a mutation which occurred before the CNV ("C2"). If not provided, the right context will be estimated from the input
#'  
#'  
#' @return   A list of the three cellular prevalence of each of the three groups of cells
#'      
#' @examples
#' 
#' Prevalences = getPrevalenceLinear(8, 5,2,1,3,10,"C1")
#'  
#'   print(Prevalences)
#' # Germ  Alt Both 
#' # 0.4  0.0  0.6 
#' 
#' @seealso \code{\link{getPrevalence}},   \code{\link{getMatrices}}
#' @export
getPrevalenceLinear<-function(lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G,context,Trace=FALSE,LocusCoverage=1)
  {
  

  
  if(is.na(mu_S) || is.na(mu_G)){
    Prevalence_output=NA
  }else{
    
    
    
  }
  
  
  
  if(Trace){
    cat("\n\n The input :\n ")
    cat(" lambda_S :", lambda_S," mu_S :", mu_S," major_cn :", major_cn," minor_cn :", minor_cn," lambda_G :", lambda_G, " mu_G :", mu_G," context :", context)
  }

   
  matrix=getMatrices(lambda_S,mu_S,major_cn,minor_cn,lambda_G, mu_G,context,LocusCoverage=LocusCoverage)
  if(Trace){
    cat("\n\n The matrices :\n ")
    print(matrix)
  }
  
  #We solve a system A*X =B
 
  
  
  #matrix A
  A=matrix(ncol=3,nrow=3)
  A[1,] = c(1,1,1)
  A[c(2,3),] = matrix$W %*% matrix$C-matrix$M
  B=c(1,0,0)
  

  
  # equality constraints
  e = matrix( c(1, 1, 1), nrow=1, byrow=TRUE) # theta_G + theta_A + theta_B = 1
  f = 1
  
  # bounds (first 3 rows theta_G, theta_A, theta_B > 0, second 3 rows -theta_G, -theta_A, -theta_B > -1 (equiv. theta_G, theta_A, theta_b < 1))
  g = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, -1), nrow=6, byrow=TRUE)
  h = matrix(c(0, 0, 0, -1, -1, -1), nrow=6, byrow=TRUE)
  
  #   cat("\n A \n") 
  #   print( A)
  #   cat("\n B \n") 
  #   print( B)
  #   cat("\n E \n") 
  #   print( e)
  #   cat("\n F \n") 
  #   print( f)
  #   cat("\n G \n") 
  #   print( g)
  #   cat("\n H \n") 
  #   print( h)
  # use constrained linear solver
  
  if(Trace){
    
    
    cat("\n \n matrices in  system A*X = B\n")
    cat("\n A \n")
    print(A)
    cat("\n B\n")
    print(B)
   # cat("\n ** Result obtained with solve() of  A*X=B\n")
    #P=solve(A,B)
    #names(P) = c("Germ","Alt","Both")
    #print(P)
    
    
    #cat("\n\n ** Result obtained from lsei without bounds\n") 
   # print(  lsei( A = A, B = B, E = e, F = f))
    
    cat("\n\n ** Result obtained from lsei with bounds\n") 
    print(  lsei( A = A, B = B, E = e, F = f, G = g, H = h))
    
  }

  
  #coutput_with_bounds 
  Prevalence= lsei( A = A, B = B, E = e, F = f, G = g, H = h)
  #print(as.numeric(format(Prevalence$X,digits=2)))
  
  Prevalence=c(as.numeric(format(Prevalence$X,digits=2)),Prevalence$solutionNorm)
  #print(Prevalence)
  Prevalence=as.numeric(format(Prevalence,digits=2))
  names(Prevalence) = c("Germ","Alt","Both","residual")
  
  if(Trace){
    cat("\n\n The prevalences (from lsei() of  limSolve with bounds) :\n ")
    print(Prevalence)
  }
  

  
  Prevalence
  
  
  
}





#' Compute the cellular prevalence of each group of cells in case of SNVOnly mode
#' 
#' This is a generic function to compute the detailed prevalence of a single mutation using the linear system making the model.
#' 
#' 
#' @param lambda_S : A count of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param mu_S : A count of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn: major copy number at the locus of the mutation
#' @param minor_cn : minor copy number (or a vector of copy number if multiple tumor samples)
#’ at the locus of the mutation 
#' @param context represents either the situation of a mutation which occurred after the CNV ("C1") or the context of a mutation which occurred before the CNV ("C2"). If not provided, the right context will be estimated from the input
#' @param sigma The parental copy number of  the chromosome harboring the mutation locus. Only needed if the context = C2. Should be either the major copy number either minor copy number
#'  
#'  
#' @return   A list of the three cellular prevalence of each of the three groups of cells
#'      
#' @examples
#' 
#' Prevalences = getPrevalenceSNVOnly(3,10,2,1,"C2",2)
#'  
#' print(Prevalences)
#' #Germ      Alt     Both residual 
#' #0.60     0.31     0.09     0.00 
#' 
#' @seealso \code{\link{getPrevalence}},   \code{\link{getMatricesSNVonly}}
#' @export
getPrevalenceSNVOnly<-function(lambda_S,mu_S,major_cn,minor_cn,context,sigma=NULL,Trace=FALSE)
  {
  
   if(Trace){
     cat("\n\n The input :\n ")
     cat(" lambda_S :", lambda_S," mu_S :", mu_S," major_cn :", major_cn," minor_cn :", minor_cn,"sigma", sigma, " context :", context)
   }
   
   
  matrix=getMatricesSNVOnly(lambda_S,mu_S,major_cn,minor_cn,context,sigma)
  
   if(Trace){
     cat("\n\n The matrices :\n ")
     print(matrix)
   }
   
  #We solve a system A*X =B
  Trace=FALSE
  
  
  if(Trace){
    print(  matrix)
  }
  #matrix A
  A=matrix(ncol=3,nrow=2)
  A[1,] = c(1,1,1)
  A[2,] = matrix$W %*% matrix$C-matrix$M
  B=c(1,0)
  
  # equality constraints
  e = matrix( c(1, 1, 1), nrow=1, byrow=TRUE) # theta_G + theta_A + theta_B = 1
  f = 1
  
  # bounds (first 3 rows theta_G, theta_A, theta_B > 0, second 3 rows -theta_G, -theta_A, -theta_B > -1 (equiv. theta_G, theta_A, theta_b < 1))
  g = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, -1), nrow=6, byrow=TRUE)
  h = matrix(c(0, 0, 0, -1, -1, -1), nrow=6, byrow=TRUE)
  
#   cat("\n A \n") 
#   print( A)
#   cat("\n B \n") 
#   print( B)
#   cat("\n E \n") 
#   print( e)
#   cat("\n F \n") 
#   print( f)
#   cat("\n G \n") 
#   print( g)
#   cat("\n H \n") 
#   print( h)
  # use constrained linear solver
  #coutput_with_bounds 
  Prevalence= lsei( A = A, B = B, E = e, F = f, G = g, H = h)
  Prevalence=c(as.numeric(format(Prevalence$X,digits=2)),Prevalence$solutionNorm)
  names(Prevalence) = c("Germ","Alt","Both","residual")
  
  Prevalence
  
  #cat("\n\n\n")
 # print(Prevalence)
}











#' @export
getLocusGermlineMutations<-function(somatic_snp_allelecount_df, snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,cnv_fraction,phasing_association_df,  tumoursamples,mode="PhasedSNP",  LocusRadius)
{
  
  #We then  retrieve for each somatic mutation the list of germline mutations to consider for the prevalence computation
  # a) if PhasedSNP mode, then the considered germline are the germline mutations phased to the somatic mutation and located within the same locus
  # b) if FlankingSNP mode then the considered germline are the close germlines located within LocusRadius distance from the somatic mutation.
  #c) if OptimalSNP mode tho columns are provided, the first for the phased germline, and if only there is not phasing information for this mutation, then the second column contains the close germlines located within LocusRadius
  
  
  # LinkedGermlineMutation=getLocusGermlineMutations(somatic_snp_allelecount_df, snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,cnv_fraction,phasing_association_df,  samples_to_pool,mode,  LocusRadius)
  
  
  
  Association=("Phased_List" %in% colnames(phasing_association_df)) 
  if(!Association)
  PhasingCode=TRUE
  
  
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
    
    
    
    if(mode=="PhasedSNP"){
      
      if(Association){
        
        if(is.null(phasing_association_df))
          stop("\n\n if mode=PhasedSNP the phasing association matrice should be provided")
        CandidateGermlines = as.character(phasing_association_df[mut,"Phased_List"])
        if(is.null(CandidateGermlines))
          next
        CandidateGermlines<-unlist(strsplit(CandidateGermlines,":"));  
        if (length(unlist(CandidateGermlines))==0)
          next
        
      }else if(PhasingCode){
        
        PhasedSamples=intersect(colnames(snp_allelecount_df[4:ncol(snp_allelecount_df)]), colnames(phasing_association_df))
        if(length(PhasedSamples)==0){
          stop("\n\n\t No samples in common between the columns of the Allele Count and the phasing Association.")
        }
        CandidateGermlines=c()
        for(sample in PhasedSamples)
        {
          mut_phasingcode=as.character(phasing_association_df[mut,sample])
          if(length(mut_phasingcode)==0 || is.na(mut_phasingcode))
            next
          
          CandidateGermlines = unique(c(CandidateGermlines,rownames(phasing_association_df[!is.na(as.character(unlist(phasing_association_df[sample]))) & as.character(unlist(phasing_association_df[sample]))==mut_phasingcode,])))
          
        }
        
      }else{
        stop("\n'n Unknown Phasing Information format. Format should be an Association or a matrice of phasing codes. Check your phasing information matrice and try again")
      }
      
    }else if(mode=="FlankingSNP"){
      CandidateGermlines =rownames(germline_snp_allelecount_df[germline_snp_allelecount_df$End >= mut_pos - LocusRadius & germline_snp_allelecount_df$End <= mut_pos + LocusRadius, ])
    }else if(mode=="SNVOnly"){
      CandidateGermlines =NA
    }else{
      stop("\n\n The mode parameters shuld be either PhasedSNP either FlankingSNP")
    }
    
    CandidateGermlines=setdiff(CandidateGermlines, mut)
    
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
      minor_som=minor_copynumber_df[mut,sample]
      
      absence_copynumberprofile=c()
      count_lower_than_somatic<-c() # For more accuracy in case of abundance of germline, someone can choose to exclude germline having an allele count less than the somatic allele count.
      
      for (germ in rownames(submatrix_CandidateGermlines_leftside))
      {
        if(!is.null(cnv_fraction))   
          phi_germ=cnv_fraction[germ , sample]
        major_germ=major_copynumber_df[germ , sample]
        minor_germ=minor_copynumber_df[germ, sample]
        
        
        #         #To check zygocity, if minor_cn>0 and refcount_cnp=0, then homozygote
        #         varcount_germ= snp_allelecount_df[germ,sample]
        #         refcount_germ= ref_allelecount_df[germ,sample]
        #         if(!is.na(minor_germ) && !is.na(refcount_germ))
        #         if((minor_germ>0) && !is.na(refcount_germ) && (refcount_germ==0)){
        #           #stop()
        #           # cat("\n We stop here : ",refcount_germ, " for mut ", germ )
        #           next        
        #         }
        #         
        #         
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
        minor_germ=minor_copynumber_df[germ, sample]
        
        
        #         #To check zygocity, if minor_cn>0 and refcount_cnp=0, then homozygote
        #         varcount_germ= snp_allelecount_df[germ,sample]
        #         refcount_germ= ref_allelecount_df[germ,sample]
        #         if(!is.na(minor_germ) && !is.na(refcount_germ))
        #         if((minor_germ>0) && !is.na(refcount_germ) && (refcount_germ==0)){
        #           #stop()
        #           # cat("\n We stop here : ",refcount_germ, " for mut ", germ )
        #           next        
        #         }
        #         
        
        if ( (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) && 
             (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          break
        }
      }
      
      notsamelocus_germline = setdiff(CandidateGermlines, samelocus_germline)
      germline_to_exclude=unique(c(germline_to_exclude,notsamelocus_germline))
    }
    samelocus_germlines= setdiff(CandidateGermlines,germline_to_exclude)
    LinkedGermlines[mut,"LinkedGermlines"] = paste(samelocus_germlines, collapse=":")
    
  }
  
  LinkedGermlines 
  
}




#### UNCOMMENT TO INCLUDE THE GENERAL FORMULA
#############################################
#############################################


# 
# 
# 
# 
# #' Compute detailed prevalence at a single mutation point under the Phased SNP mode using the General (exact) formula
# #' 
# #' This is a generic function to compute the prevalence at a single somatic mutation point using a phased Germline SNP.
# #' 
# #' @param lambda_S : A count (or a vector of counts) of 
# #' alleles supporting the variant sequence of  the somatic mutation 
# #' @param mu_S : A count (or a vector of counts) of 
# #' alleles supporting the reference sequence  of  the somatic mutation
# #' @param major_cn  major copy number (or a vector ) 
# #' at the locus of the mutation
# #' @param minor_cn : minor copy number (or a vector )  
# #'  at the locus of the mutation 
# #' @param lambda_G  A count (or a vector of counts if multiple samples) of 
# #' alleles supporting the variant sequence of  the Germline SNP
# #' @param mu_G  A count (or a vector of counts) of 
# #' alleles  supporting the reference sequence  of  the Germline SNP
# #' @param cnv_fraction If provided, represents the fraction of cells affected by the copy number alteration. This value, if not provided, is computed from the allelic count information and copy number information. Default NULL
# #' @param formula Can be either "Matrix" either "General", specify if the prevalence should be computed using the linear form formula or the General form formula. Default "Matrix"
# #' @param detail In case form="Matrix", when set to TRUE, a detailed output is generated containing, the context and the detailed prevalence for each group of cells (germline cells (Germ), cells affected by one of the two genomic alterations (Alt), cells affected by  by both genomic alterations (Both) ).
# #' @param SameTumour In the case of a multiple sample computation, specify if the samples are originating from the same tumour. 
# #'  
# #' @return   \describe{
# #'        \item{}{ if form="general", the function return a numerical value representing the prevalence at the somatic mutation.}
# #'        \item{}{ if form="matrix", the function return a list containing the following data frames:
# #'  \describe{
# #'        \item{Context}{ The associated context}
# #'        \item{Prevalence}{The computed prevalence}
# #'        \item{DetailedPrevvalence}{Detailed prevalence for each of the three genotype groups separated by "|". The three groups are Germline mutations, mutations harboring one of the two alterations (CNV or SNP) mutations harboring both alterations }
# #'        
# #'      }
# #'       }
# #'       }
# #'      
# #'     
# #'      
# #' @examples
# #' 
# #' # We reproduce here the case study No 6 of the paper
# #' #prevalence=getPhasedSNPPrevalenceGeneral(lambda_S=14,mu_S=10,major_cn=3,minor_cn=1, lambda_G=16, mu_G=8, cnv_fraction=4/8 )
# #'  
# #' @seealso \code{\link{getPrevalence}} 
# #' @export
# getPhasedSNPPrevalenceGeneral<-function(lambda_S,mu_S,major_cn,minor_cn, lambda_G, mu_G, cnv_fraction=NULL,min_cells=1, min_alleles=1,SameTumour=TRUE)
#   
# {
#   tumoursamples= names(lambda_G)
#   
#   if (is.null(tumoursamples)){
#     tumoursamples = paste("Sample",c(1:length(lambda_G)),sep="_")
#     names(lambda_G) =  tumoursamples
#     names(mu_G) =  tumoursamples
#     names(lambda_S) =  tumoursamples
#     names(mu_S) =  tumoursamples
#     names(major_cn) =  tumoursamples
#     names(minor_cn) =  tumoursamples
#     
#   }
#   
#   #   mut=rownames(lambda_S)[1]
#   #   
#   #   print(mut)
#   #   if(is.null(mut))
#   #   {
#   #     mut="somatic"
#   #     rownames(lambda_S) = mut
#   #   }
#   #   
#   #First, we splitthe samples according to their copy number profile. 
#   #If they come for the same tumour, then samples having the same copy numbe rprofile are grouped together.
#   #If they dont come from the same tumour, then each sample is put in its own group.
#   CNV_groups<-list()
#   if(SameTumour)
#   {
#     # for (phival in unique(phi_cn)) 
#     for(majorval in unique(major_cn)) 
#       for (minorval in unique(minor_cn))
#       {
#         
#         if( !is.na(majorval) && !is.na(minorval))
#         {
#           samplecnvgroups=c()
#           
#           for (sample in tumoursamples)
#             if(!is.na(major_cn[sample]) && !is.na(minor_cn[sample]))
#               if ((major_cn[sample]==majorval) && (minor_cn[sample]== minorval ))
#                 samplecnvgroups=c(samplecnvgroups,sample)
#               
#               if(length(samplecnvgroups)!=0)
#                 CNV_groups[length(CNV_groups) +1] = list(samplecnvgroups)
#         }
#       }
#     
#   }else{
#     
#     for (sample in tumoursamples)
#       CNV_groups[length(CNV_groups) +1] = list(sample)
#   }
#   
#   
#   sigma_G=major_cn-major_cn
#   rho_G=major_cn-major_cn
#   
#   phi_cn=cnv_fraction
#   if(is.null(cnv_fraction)) phi_cn=major_cn-major_cn
#   Ncells_list=major_cn-major_cn
#   #tau_G= phi_cn * sigma_G + (1-phi_cn) # for \tau(G)
#   tau_G=major_cn-major_cn
#   hatlambda_S=lambda_S-lambda_S
#   hatlambda_G=lambda_G-lambda_G
#   context_list=rep(NA,length(lambda_S))
#   names(context_list)= names(lambda_S)
#   
#   
#   #We run the model on each group
#   for (icnvgroup in 1 : length(CNV_groups))
#   {
#     sample_cnvgroup=CNV_groups[[icnvgroup]]
#     
#     # Estimate of sigma. See the paper
#     #################################
#     A=0
#     B=0
#     Asnp=0
#     Bref=0
#     for (sample_i in sample_cnvgroup)
#     {
#       if (!is.na(lambda_G[sample_i])  )
#         Asnp= Asnp + lambda_G[sample_i]
#       if (!is.na(mu_G[sample_i]))
#         Bref=Bref+mu_G[sample_i]
#       
#       if(!is.na(lambda_G[sample_i])&& !is.na(mu_G[sample_i]) )
#       {
#         A= A + lambda_G[sample_i]
#         B=B+mu_G[sample_i]
#       }
#     }
#     
#     if((is.na(A)) || (is.na(B)))
#     {
#       A=Asnp
#       B=Bref
#     }
#     if(A>=B)
#       sigma_G[sample_cnvgroup] = major_cn[sample_cnvgroup] 
#     if(A<B)
#       sigma_G[sample_cnvgroup] = minor_cn[sample_cnvgroup] 
#     
#     rho_G[sample_cnvgroup] = major_cn[sample_cnvgroup]  + minor_cn[sample_cnvgroup]  - sigma_G[sample_cnvgroup]
#     
#     if(is.null(cnv_fraction)) {
#       A_equation = lambda_G[sample_cnvgroup] - mu_G[sample_cnvgroup]
#       B_equation= mu_G[sample_cnvgroup] * (sigma_G[sample_cnvgroup] - 1) -  lambda_G[sample_cnvgroup] * (rho_G[sample_cnvgroup] - 1)
#       
#       phi_cn[sample_cnvgroup] = A_equation / B_equation
#       phi_cn[phi_cn>1] = 1
#       
#       Ncells_list[sample_cnvgroup] =  lambda_G[sample_cnvgroup]  * B_equation / ( B_equation + A_equation * (sigma_G[sample_cnvgroup] - 1)  )
#       
#       noCNV_cases = which(A_equation==0 && B_equation==0)
#       noCNV_cases = which( B_equation==0)
#       if(length(noCNV_cases )>0){
#         Ncells_list[noCNV_cases] = lambda_G[noCNV_cases]
#         phi_cn[noCNV_cases] = rep(0,length(noCNV_cases))
#       }
#     }
#     
#     tau_G[sample_cnvgroup]= phi_cn[sample_cnvgroup] * sigma_G[sample_cnvgroup] + (1-phi_cn[sample_cnvgroup]) # for \tau(G)
#     
#     alpha_list=1/tau_G[sample_cnvgroup]
#     beta_list=(phi_cn[sample_cnvgroup]  * sigma_G[sample_cnvgroup] )/ tau_G[sample_cnvgroup]
#     
#     hasTau0 = which(tau_G==0)
#     if(length(hasTau0)>0){
#       alpha_list[hasTau0] = rep(0,length(hasTau0))
#       beta_list[hasTau0] = rep(0,length(hasTau0))
#     }
#     
#     lambda_S=as.numeric(lambda_S[sample_cnvgroup])
#     lambda_G=lambda_G[sample_cnvgroup,drop=F]
#     alpha=as.numeric(alpha_list[sample_cnvgroup])
#     beta=as.numeric(beta_list[sample_cnvgroup])
#     
#     EM_parameters=bestAllele(lambda_S, lambda_G, alpha, beta)
#     
#     print(EM_parameters)
#     
#     hatlambda_S[ sample_cnvgroup] =EM_parameters$hatlambda_S
#     hatlambda_G[sample_cnvgroup] =EM_parameters$hatlambda_G
#     context=EM_parameters$bestC
#     context_list[sample_cnvgroup]=rep(context, length(lambda_S))
#   }
#   
#   ##### Intermediary computation
#   u_list =phi_cn * hatlambda_G/tau_G # for u
#   v_list=(1-phi_cn) * hatlambda_G/tau_G # for v
#   
#   if(length(hasTau0)>0){
#     u_list[hasTau0] = rep(0,length(hasTau0))
#     v_list[hasTau0] = rep(0,length(hasTau0))
#   }
#   
#   sum_cells_list<-u_list+v_list # estimation of the number of cells if well counts and the K * the number of cells of reads count with k = amplification factor
#   sum_allele_list= hatlambda_G +  mu_G # estimation of the number of allele
#   
#   ##### Prevalence computation
#   
#   prev_S =vector("numeric", length=length(lambda_G))
#   names(prev_S) = names(lambda_G)
#   prev_S[prev_S==0]<-NA
#   
#   for(sample in tumoursamples)
#   {
#     
#     if (is.na(context_list[sample]))
#       next
#     
#     if(is.na(sum_cells_list[sample]) || is.na(sum_allele_list[sample]))
#       next
#     if (sum_cells_list[sample]< min_cells)
#       next
#     if(sum_allele_list[sample]< min_alleles)
#       next
#     
#     if (context_list[sample]=="C2")#CONTEXTE 2{
#     {
#       if(v_list[sample]!=0){
#         prev_S[sample] = as.numeric(unlist(   phi_cn[sample] + (1-phi_cn[sample]) * ((hatlambda_S[sample] - u_list[sample] * sigma_G[sample])/v_list[sample])   ))
#       }else{
#         prev_S[sample] = as.numeric(unlist(   phi_cn[sample] ))
#       }
#     } else if (context_list[sample]=="C1"){
#       prev_S[sample] = as.numeric(unlist( (hatlambda_S[sample] / hatlambda_G[sample])* tau_G[sample]    ))
#     }
#     if(tau_G[sample]==0)
#       prev_S[sample] =0
#   }
#   prev_S
#   
# }
# 






