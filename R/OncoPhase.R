############################################################
#
# OncoPhase : Files Organization 
# data.R - file to generate the Roxygen help documents for the data
# OncoPhase.R -  Main file of the package (This file)
# OncoPhase_methods.R  File containing the methods implementation

#Roxigen help for the package (top help)

#' 
#' OncoPhase: An R package for somatic mutations cellular prevalence quantification using haplotype phasing.
#' 
#' OncoPhase uses haplotype phase information when required to accurately compute mutational cellular prevalence. OncoPhase utilizes three sources of information: the phasing information, the copy number variation, and the allele counts.  It takes as input a combination of phased SNV and SNP allele-specific sequence read counts and local allele-specific copy numbers to determine the proportion of cells harboring the SNV and compute specific and detailed mutation cellular prevalence for each of the following groups of cells: 
#'  \describe{
#'       \item{Germ}{ Germline cells having a normal genotype with no mutations and no copy number alteration at the considered locus. }
#'       \item{Alt}{ Cells harboring one  alternative between the two somatic alterations. That is either only the SNV if C=1 (SNV occurred before SCNA) or only the SCNA if C=0 (SNV occurred after the SCNA).}
#'       \item{Both}{ Cells harboring both somatic alterations. That is the SNV and the SCNA}
#'    }
#' 
#' OncoPhase can also compute the mutation cellular prevalence without requiring any nearby phased SNP when no phasing information is  available or when explicitely specified by the choice of the mode.  OncoPhase can be run  under three different modes : 
#'  \describe{
#'        \item{PhasedSNP}{ Phasing information is required. The prevalence is computed relatively to a nearby Phased SNP whose allelic counts should be provided}
#'        \item{SNVOnly}{The prevalence is computed using only the SNV information without the usage of any nearby SNP}
#'        \item{Ultimate}{ This is the default mode. For a given mutation, the method checks if the phasing information is required to compute an accurate cellular prevalence. If it is not, the SNVOnly mode  is used. If instead the phasing information is required the mode is then set to PhasedSNP if allelic counts of a phased nearby SNP are provided.  This is done by first computing the prevalence under the SNVOnly mode. If the data do not fit into this mode (hiogh residual of the linear model), then the prevalence is computed using PhasedSNP mode.
#'     }
#'     }
#'     OncoPhase also infer the context establishing the temporal relationship between the SNV and the copy number lateration affecting the mutation locus. Two context exists :
#'  \describe{
#'       \item{C1}{ The SNV occured after the copy number alteration }
#'       \item{C2}{ The SNV occured before the copy number alteration}
#'    }  
#' 
#' 
#' 
#' The main functions of OncoPhase package  are    \code{\link{getPrevalence}}  and   \code{\link{getSamplePrevalence}}.   For more detailed information on usage, see the package vignette, by typing
#' \code{vignette("OncoPhase")}. All support questions should be emailed to the authors.
#' 
#' @references
#' 
#' OncoPhase reference:
#' 
#' Chedom-Fotso Donatien, Ahmed Ashour Ahmed, and Christopher Yau. "OncoPhase: Quantification of somatic mutation cellular prevalence using phase information." bioRxiv (2016): 046631.
#' 
#' 
#' @author Donatien Chedom-Fotso, Ahmed Ahmed, Christopher Yau.
#' 
#' @docType package
#' @name A-OncoPhase
#' @aliases OncoPhase-package
#' @keywords package
#' @import limSolve
NULL




#' @export
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
  
  
  #require(limSolve)
  
  invisible()
}


