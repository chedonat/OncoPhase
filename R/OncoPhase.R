############################################################
#
# OncoPhase : files Organization 
# data.R - file to generate the Roxygen help documents for the data
# OncoPhase.R -  this file is the main file of the package
# OncoPhase_methods.R - This files contains the main implementations of the prevalence computation
# OncoPhase_EMmodel.R - This files contains the expecattion-maximisation model




#Roxigen help for the package (top help)

#' OncoPhase package for somatic mutations cellular prevalence quantification using haplotype phasing
#' 
#' The main function for somatic mutation cellular prevalence computation are    \code{\link{getPrevalence}} and  \code{\link{getPrevalenceMultiSamples}} . \code{\link{getPrevalence}}  compute the cellular prevalence at a particular mutations under 4 distincts modes : PhasedSNP, FlankingSNP, OptimalSNP and SNVOnly.  \code{\link{getPrevalenceMultiSamples}}  computes the cellular prevalence of a list of mutations located at a given region of the genome. It can also work on a whole genome scale.
#' See the manual and examples at \code{\link{getPrevalence}} and  \code{\link{getPrevalenceMultiSamples}}  for more details.
#' 
#'  The function \code{\link{getPrevalenceLinear}} compute the prevalence of a given mutation by directly solving the linear system associated to the model.
#' 
#' Input data for simple case studies can be generated with the function \code{\link{build_casestudy}}.
#' 
#' The package include experimental data for chromosome 22 from a  patients retrieved from a parallel clinical study.(see \code{\link{chr22_11152}} )
#' 
#' For more detailed information on usage, see the package vignette, by typing
#' \code{vignette("OncoPhase")}. All support questions should be emailed to the authors.
#'
#' @references
#'
#' OncoPhase reference:
#' 
#' OncoPhase: A package for computing Somatic Mutation cellular Prevalence in cancer using haplotype phasing. Bioinformatics 2016. Submitted
#'
#' OncoPhase reference:
#' 
#' ” Ovarian cancer haplotype sequencing reveals ubiquitous SOX2 overexpression in the premalignant fallopian tube epithelium”
#'
#' @author Donatien Chedom-Fotso, Ahmed Ahmed, Christopher Yau.
#' 
#' @docType package
#' @name aOncoPhase
#' @aliases OncoPhase-package
#' @keywords package
#' @import limSolve
NULL



#'  @export
hg19_dfsize<-list(chr1=249250621,chr2=243199373,chr3=198022430, chr4=191154276, 
                  chr5=180915260, chr6=171115067, chr7=159138663, chr8=146364022, 
                  chr9=141213431, chr10=135534747, chr11=135006516, chr12=133851895, 
                  chr13=115169878, chr14=107349540,chr15=102531392, chr16=90354753,
                  chr17=81195210, chr18=78077248, chr19=59128983, chr20=63025520,
                  chr21=48129895,chr22=51304566, chrX=155270560, chrY=59373566,
                  chrM=16571)





# Welcome message
#################

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\n************\nWelcome to OncoPhase package\n************")
}

# Set up some custom options
############################

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.OncoPhase <- list(
    OncoPhase.path = "~/R",
    OncoPhase.install.args = "",
    OncoPhase.name = "Donatien Chedom-Fotso",
    OncoPhase.desc.author = '"Donatien Chedom-Fotso <donatien.chedom@gmail.com> 
    [aut, cre]"',
    OncoPhase.desc.license = "GPL-2",
    OncoPhase.desc.suggests = NULL,
    OncoPhase.desc = list(),
    OncoPhase.PatientName = "NoName",
    OncoPhase.WD = getwd()
  )
  toset <- !(names(op.OncoPhase) %in% names(op))
  if(any(toset)) options(op.OncoPhase[toset])
  
  invisible()
}



# Main prevalence function
##########################
##########################
##########################



setClass("SampleList" )

setClass("SampleCnvAssociation" )



Sample <- setClass(
  #Set the name of the class
  "Sample",
  
  #Define the slots
  slots = c(
    name = "character",
    ID = "character",
    type = "character",
    cnvprofile = "character",
    group = "character",
    alleletype = "character",
    variantfiletype = "character",
    variantfilename = "character",
    phasinganalysis = "logical",
    prevalenceanalysis = "logical"
  ),
  
  #Set the defaults values for the slots
  prototype = list(
    name = "",
    type = "tumour",
    variantfiletype = "varfile",
    alleletype = "wells",
    phasinganalysis = TRUE,
    prevalenceanalysis = TRUE
  ),
  
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  
  validity = function(object)
  {
    if(! object@type %in% c("tumour", "normal")) {
      return(" Incorrect sample type assigned. The sample type should be one of ")
    }
    if(! object@variantfiletype %in% c("varfile", "mastervarfile", "vcf")) {
      return(" Incorrect sample type assigned. The sample type should be one of varfile, mastervarfile, vcf")
    }
    
  }
)


initsSamples<-function(samples_df)
{
  sampleList<-list()
  
  # Control of the entries, check the match between the list of properties and the 
  #column 
  # of the samples data frame
  slotList = slotNames("Sample")
  colList =  colnames(samples_df)
  
  presentSlots = intersect(slotList, colList)
  absentSlots = setdiff(slotList, colList)
  toignoreFields = setdiff(colList, slotList)
  
  if (length(presentSlots)>0){
    cat("\n\n\t\t The following properties of samples wil be retrieved : \n\t\t")
    print(presentSlots)
  }
  
  if (length(absentSlots)>0){
    cat("\n\n\t\t The following properties of samples are not present in the fields 
        provided : \n\t\t")
    print(absentSlots)
  }
  
  if (length(toignoreFields)>0){
    cat("\n\n\t\t The following fields provided ar not part of a sample property : \n\t\t")
    print(toignoreFields)
  }
  
  if (length(presentSlots)==0){
    stop("\n\n\t\t Error, No sample property  provided in the data frame, check your input")
  }
  
  
  
  
  
  
  
  for (isample in 1:nrow(samples_df)){
    newSample=Sample(
      name = as.character(samples_df[isample,"name"]),
      ID =  as.character(samples_df[isample,"ID"]),
      type =  as.character(samples_df[isample,"type"]),
      cnvprofile =  as.character(samples_df[isample,"cnvprofile"]),
      group =  as.character(samples_df[isample,"group"]),
      alleletype =  as.character(samples_df[isample,"alleletype"]),
      variantfiletype =  as.character(samples_df[isample,"variantfiletype"]),
      variantfilename =  as.character(samples_df[isample,"variantfilename"]),
      phasinganalysis =  as.logical(samples_df[isample,"phasinganalysis"]),
      prevalenceanalysis =  as.logical(samples_df[isample,"prevalenceanalysis"])
    )
    
    sampleList<-append(sampleList, newSample)
  }
  
  return(sampleList)
}


setGeneric("setProperty",
           def=function(theObject, xname, xvalue){
             standardGeneric("setProperty")
           }
)

setMethod(f="setProperty",
          signature = "Sample",
          definition = function(theObject, xname, xvalue){
            slot(theObject,xname) = xvalue
            return(theObject)
          }
)

#' @export
numeric_column<-function(df,tumoursamples)
{
  for (col in tumoursamples)
  {
    df[col] =  as.numeric(unlist(df[col]))
  }
  df
}







