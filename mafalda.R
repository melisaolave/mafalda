v0.0


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
