#' Convert VCF data into long-format loci sequences
#'
#' This function reconstructs full-length sequences for each individual
#' and contig from a VCF file and a reference genome. It outputs sequences
#' in a long format suitable for downstream phylogenetic or population
#' genetic analyses. Kept sites is in vcftools format (https://vcftools.github.io/man_latest.html) and reference in dDocent fasta format (https://ddocent.com//)
#'
#' @param wd Character. Working directory. Default is \code{getwd}.
#' @param vcf.name Character. Path to the VCF file.
#' @param kept.sites Character. File containing retained SNP positions.
#' @param reference Character. Path to a reference FASTA file.
#' @param rm.initial.seq Character. Sequence to remove from the beginning
#'   of each reconstructed sequence. For example, to delete EcoRI sequence "NTAATTC" (default = \code{""}).
#' @param rm.final.seq Character. Sequence to remove from the end of each
#'   reconstructed sequence (default = \code{""}).
#' @param phased.vcf Logical. Whether the VCF is phased (default = \code{TRUE}).
#' @param out.name Character. Output file name.
#'
#' @return No return value. Writes reconstructed sequences to file.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Reads SNP positions and VCF genotype data.
#'   \item Reconstructs sequences using the reference genome.
#'   \item Handles phased and unphased data differently:
#'   \itemize{
#'     \item Phased: reconstructs two haplotypes per individual.
#'     \item Unphased: encodes heterozygous sites using IUPAC ambiguity codes.
#'   }
#'   \item Replaces missing data ('.') with '?'.
#'   \item Outputs sequences in a simple text format.
#' }
#'
#' @examples
#' \dontrun{
#' vcf2longLoci(
#'   vcf.name = "data.vcf",
#'   kept.sites = "sites.txt",
#'   reference = "ref.fasta",
#'   phased.vcf = TRUE,
#'   out.name = "output.txt"
#' )
#' }
#'
#' @importFrom ape read.FASTA
#' @importFrom vcfR read.vcfR
#'
#' @export
vcf2longLoci <- function(wd = getwd, vcf.name, kept.sites, reference, 
                         rm.initial.seq="", rm.final.seq="",
                         phased.vcf=T,
                         out.name){
  sites <- read.table(kept.sites, sep="\t", header=T)
  vcf <- read.vcf(vcf.name, to=nrow(sites))
  chr <- unique(sites$CHROM)
  ref <- read.FASTA(reference)
  
  # filter reference to match kept.sites file
  keep <- labels(ref)[labels(ref) %in% chr]
  ref <- ref[keep]
  
  ind <- rownames(vcf)
  cat("A total of:\n\t", 
      nrow(sites), "were read in vcf\n\t",
      length(chr), "contigs will be reconstructed\n\t",
      length(ind), "of individuals found in vcf\n")
  matrix <-NULL
  CHR.count <- 0
  for(i in 1:length(chr)){
    cat("Working on contig ", chr[i], "...", sep="")
    chr.snps <- sites[sites$CHROM == chr[i],]
    pos <- chr.snps$POS
    ref.chr <- unlist(as.character(ref[i]))
    matrix <- c(matrix, paste(length(ind), 
                              length(ref.chr)-length(unlist(strsplit(rm.initial.seq, split="")))-length(unlist(strsplit(rm.final.seq, split="")))))
    for(k in 1:length(ind)){
      if(phased.vcf){ # if phased vcf, reconstruct two alleles
        cat("Phased vcf provided. Two alleles will be reconstructed\n")
        ind.seq_0 <- ref.chr
        ind.seq_1 <- ref.chr
        SNP.count <- 1+CHR.count
        for(j in 1:length(pos)){
          SNP <- as.data.frame(vcf[k,SNP.count])
          SNP <- unlist(strsplit(as.character(SNP$.), split=""))[c(1,3)]
          ind.seq_0[pos[j]] <- SNP[1]
          ind.seq_1[pos[j]] <- SNP[2]
          SNP.count <- SNP.count +1
        }
        #ind 0
        ind.seq_0 <-paste(ind.seq_0, collapse="")
        ind.seq_0 <- toupper(ind.seq_0)
        ind.seq_0 <- gsub(pattern=paste("^", rm.initial.seq, sep=""), replacement="", ind.seq_0)
        ind.seq_0<- gsub(pattern=paste(rm.final.seq,"$", sep=""), replacement="", ind.seq_0)
        ind.seq_0<- gsub(pattern="\\.", replacement="?", ind.seq_0)
        new.ind_0<-paste(ind[k], "_0      ",ind.seq_0, sep="")
        
        #ind 1
        ind.seq_1 <-paste(ind.seq_1, collapse="")
        ind.seq_1 <- toupper(ind.seq_1)
        ind.seq_1 <- gsub(pattern=paste("^", rm.initial.seq, sep=""), replacement="", ind.seq_1)
        ind.seq_1<- gsub(pattern=paste(rm.final.seq,"$", sep=""), replacement="", ind.seq_1)
        ind.seq_1<- gsub(pattern="\\.", replacement="?", ind.seq_1)
        new.ind_1<-paste(ind[k], "_1      ",ind.seq_1, sep="")
        matrix <- c(matrix, new.ind_0, new.ind_1)
        
      }else{ #code heterozygote sites as ambiguities
        ind.seq <- ref.chr
        SNP.count <- 1+CHR.count
        for(j in 1:length(pos)){
          SNP <- as.data.frame(vcf[k,SNP.count])
          SNP <- unlist(strsplit(as.character(SNP$.), split=""))[c(1,3)]
          if(SNP[1] != SNP[2]){
            if((SNP[1] == "C" & SNP[2] == "T") | (SNP[2] == "C" & SNP[1] == "T")){
              new.SNP <- "Y"
            }else if((SNP[1] == "A" & SNP[2] == "G") | (SNP[2] == "A" & SNP[1] == "G")){
              new.SNP <- "R"
            }else if((SNP[1] == "A" & SNP[2] == "T") | (SNP[2] == "A" & SNP[1] == "T")){
              new.SNP <- "W"
            }else if((SNP[1] == "G" & SNP[2] == "C") | (SNP[2] == "G" & SNP[1] == "C")){
              new.SNP <- "S"
            }else if((SNP[1] == "T" & SNP[2] == "G") | (SNP[2] == "T" & SNP[1] == "G")){
              new.SNP <- "K"
            }else if((SNP[1] == "C" & SNP[2] == "A") | (SNP[2] == "C" & SNP[1] == "A")){
              new.SNP <- "M"
            }
          }else{
            new.SNP <- SNP[1]
          }
          ind.seq[pos[j]] <- new.SNP
          SNP.count <- SNP.count +1
        }
        ind.seq<-paste(ind.seq, collapse="")
        ind.seq <- toupper(ind.seq)
        ind.seq<- gsub(pattern=paste("^", rm.initial.seq, sep=""), replacement="", ind.seq)
        ind.seq<- gsub(pattern=paste(rm.final.seq,"$", sep=""), replacement="", ind.seq)
        ind.seq<- gsub(pattern="\\.", replacement="?", ind.seq)
        new.ind<-paste(ind[k], "     ",ind.seq)
        matrix <- c(matrix, new.ind)
      }
    }
    CHR.count <- SNP.count-1
    cat("done!\n")
  }
  
  write(matrix, out.name)
}
