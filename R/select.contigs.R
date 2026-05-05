#' Filter and sample contigs based on SNP density and missing data
#'
#' This function filters contigs (loci) based on SNP density and/or missing data,
#' and optionally subsamples a fixed number of loci. The filtered dataset is
#' written to a tab-delimited file. Kept sites, SNP density and missing data input are in vcftools format (https://vcftools.github.io/man_latest.html)
#'
#' @param wd Character. Working directory. Default is \code{getwd()}.
#' @param kept.sites Character. Path to a tab-delimited file containing retained sites (in vcftools format).
#' @param n.loci Integer or \code{"all"}. Number of loci to retain after filtering.
#' @param snpden Character or NULL. Path to a file with SNP counts per contig (in vcftools format).
#'   If \code{NULL}, SNP density filtering is skipped.
#' @param SNP.density.range Numeric vector of length 2. Quantile range used to
#'   filter contigs based on SNP density (default = \code{c(0.05, 0.95)}).
#' @param lmiss Character or NULL. Path to file with missing data per SNP (in vcftools format).
#'   If \code{NULL}, missing data filtering is skipped.
#' @param min.lmiss Numeric. Maximum allowed missing data per SNP.
#' @param out.name Character. Output file name.
#'
#' @return No return value. Writes filtered loci to file.
#'
#' @details
#' The function applies the following filters:
#' \itemize{
#'   \item Removes contigs outside the specified SNP density quantile range.
#'   \item Removes contigs with SNPs exceeding the allowed missing data threshold.
#'   \item Optionally subsamples a fixed number of loci.
#' }
#'
#' @examples
#' \dontrun{
#' select.contigs(
#'   kept.sites = "sites.txt",
#'   n.loci = 100,
#'   snpden = "snp_density.txt",
#'   lmiss = "missing.txt",
#'   min.lmiss = 0.2
#' )
#' }
#'
#' @export
select.contigs <- function(wd=getwd(), kept.sites, n.loci, 
                           snpden, SNP.density.range=c(0.05, 0.95),
                           lmiss, min.lmiss, 
                           out.name=paste(kept.sites, ".out", sep="")){
  setwd(wd)
  sites <- read.table(kept.sites, sep="\t", header=T)
  
  if(is.null(snpden) == F){
    snp <- read.table(snpden, sep="\t", header=T)
    ci <-quantile(snp$SNP_COUNT, SNP.density.range) 
    cat("Total SNPs per contig within confidence interval = from", ci[1], "to", ci[2],"\nContigs out of this range will be removed\n")
    snp<-snp[snp$SNP_COUNT >= ci[1] & snp$SNP_COUNT <= ci[2],]
    sites <- sites[sites$CHROM %in% snp$CHROM,]
  }
  if(is.null(lmiss) == F){
    missing <- read.table(lmiss, sep="\t", header=T)
    chr.list <- unique(sites$CHR)
    keep.chr <- NULL
    for(i in 1:length(chr.list)){
      subtable <- missing[missing$CHR == chr.list[i],]
      cat("Working on contig", i,  "/", length(chr.list), "\n")
      if(any(subtable$F_MISS > min.lmiss) == F){ # keep chr only if all SNPs are above the min missing allowed
        keep.chr <- c(keep.chr, unique(subtable$CHR))
      }
    }
    
    sites <- sites[sites$CHROM %in% keep.chr ,]
  }
  if(n.loci != "all" & unique(sites$CHROM) >= n.loci){
    sampled <- sample(unique(sites$CHROM), n.loci, replace=F)
    sites <- sites[sites$CHROM %in% sampled,]
  }
  
  write.table(sites, out.name, sep="\t", col.names=T, quote=F, row.names = F)
  
  cat(unique(sites$CHROM), "loci recovered and saved\n")
}
