#' Compute genealogical divergence index (gdi) from MCMC output
#'
#' This function calculates the genealogical divergence index (gdi)
#' for each lineage using posterior samples of divergence time (\code{tau})
#' and population size (\code{theta}) obtained from an MCMC analysis.
#' The gdi is computed as \eqn{gdi = 1 - e^{-2\tau/\theta}}.
#'
#' The genealogical divergence index was originally proposed by
#' Jackson et al. (2017) and later formalized for species delimitation
#' by Leaché et al. (2019). It provides a quantitative measure of the
#' strength of divergence between pairs of taxa and is commonly used
#' in species delimitation analyses performed in \code{BPP}
#' (Yang and Rannala 2010) or \code{iBPP}
#' (Solís-Lemus et al. 2015).
#'
#' Summary statistics for \code{tau}, \code{theta}, and \code{gdi}
#' are written to a tab-delimited output file.
#'
#' @param wd Character. Working directory. Default is \code{getwd()}.
#' @param mcmc Data frame or matrix containing MCMC samples. The first column
#'   must correspond to iteration/state, and parameter columns must include
#'   names starting with \code{"tau_"} and \code{"theta_"}.
#' @param burnin Numeric. Proportion of samples to discard as burn-in
#'   (default = 0.1).
#' @param out.name Character. Name of the output file
#'   (default = \code{"gdi.output.txt"}).
#'
#' @return No return value. A summary table is written to file.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Removes burn-in samples based on the specified proportion.
#'   \item Extracts matching \code{tau} and \code{theta} parameters.
#'   \item Computes gdi for each posterior sample.
#'   \item Outputs summary statistics (Min, 1st Qu., Median, Mean, 3rd Qu., Max).
#' }
#'
#' @references
#' Jackson, N. D., Morales, A. E., Carstens, B. C., \& O'Meara, B. C. (2017).
#' PHRAPL: phylogeographic inference using approximate likelihoods.
#' Systematic Biology, 66(6), 1045–1053.
#'
#' Leaché, A. D., Zhu, T., Rannala, B., \& Yang, Z. (2019).
#' The spectre of too many species.
#' Systematic Biology, 68(1), 168–181.
#'
#' Yang, Z., \& Rannala, B. (2010).
#' Bayesian species delimitation using multilocus sequence data.
#' Proceedings of the National Academy of Sciences, 107(20), 9264–9269.
#'
#' Solís-Lemus, C., Knowles, L. L., \& Ané, C. (2015).
#' Bayesian species delimitation combining multiple genes and traits
#' in a unified framework.
#' Evolution, 69(2), 492–507. https://doi.org/10.1111/evo.12582 :contentReference[oaicite:0]{index=0}
#'
#' @examples
#' \dontrun{
#' get.gdi(mcmc = my_mcmc, burnin = 0.2)
#' }
#'
#' @export
get.gdi <- function(wd=getwd(), mcmc, burnin=0.1, out.name="gdi.output.txt"){
  tau.col <- grep(pattern="tau_", value=T, colnames(mcmc))
  node.names <- gsub(pattern="tau_", replacement="", tau.col)
  
  burnin.state <- mcmc[nrow(mcmc),1]*burnin
  postbunin.mcmc <- mcmc[mcmc[,1] > burnin.state,]
  tau.summary.table <- NULL
  theta.summary.table <- NULL
  gdi.summary.table <- NULL
  for(i in 1:length(tau.col)){
    #get theta
    theta.col <- paste("theta_", node.names[i],sep="")
    theta <-postbunin.mcmc[,theta.col]
    #get tau
    tau <- postbunin.mcmc[,tau.col[i]]
    #calculate gdi #gdi=1−e^(−2tau/theta)
    gdi <- 1-(exp((-2*tau)/theta))
    mean(gdi)
    tau.summary.table <- rbind(tau.summary.table, summary(tau));
    theta.summary.table <- rbind(theta.summary.table, summary(theta));
    gdi.summary.table <- rbind(gdi.summary.table, summary(gdi));
  }
  colnames(tau.summary.table) <- paste("tau_", colnames(tau.summary.table), sep="")
  colnames(theta.summary.table) <- paste("theta_", colnames(theta.summary.table), sep="")
  colnames(gdi.summary.table) <- paste("gdi_", colnames(gdi.summary.table), sep="")
  summary.table <- cbind(node.names, gdi.summary.table, tau.summary.table, theta.summary.table)
  colnames(summary.table)[1] <- "lineage"
  write.table(summary.table, out.name, sep="\t", col.names = T, row.names = F, quote = F)
  
}
