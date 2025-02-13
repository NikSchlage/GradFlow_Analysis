##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##

#--------------------------------------------------------------------------------------------------------#
#   LOAD PACKAGES AND CODE FOR ANALYSIS                                                                  #
#--------------------------------------------------------------------------------------------------------#
library(tidyverse) ## for function "alpha(col, a)"
library(tikzDevice)
library(stringr)   ## required for str_pad()
source("/home/nikolas/Documents/NPV/njjn_analysis/npv_analysis_jan_2023/npv.R")
source("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/R_Code/mixing_probe_gf.R")

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   LOAD WORKSPACE IMAGE                                                                                 #
#                                                                                                        #
#      * containing the cf ratio ratio_1_g5 = Im(<X1g5>) / Re(<Xmx_g5>)                                  #
#      * containing the cf ratio ratio_2_g5 = Im(<X2g5>) / Re(<Xmx_g5>)                                  #
#      * containing the cf ratio ratio_3_g5 = Im(<X3g5mulmin1>) / Re(<Xmx_g5>)                           #
#                                                                                                        #
#      * containing the cf ratio ratio_4_id = Re(<X1idmulmin1>) / Re(<Xmx_id_mulmin1>)                   #
#      * containing the cf ratio ratio_5_id = Re(<X2idmulmin1>) / Re(<Xmx_id_mulmin1>)                   #
#      * containing the cf ratio ratio_6_id = Re(<X3id>) / Re(<Xmx_id_mulmin1>)                          #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
my_data_path <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/01_boot_cf_data/NLoopsamples_01/01_boot_cf_data_Nconfs1126_4qFlavorChange/bootstrap_cf_mx_prb_oet_gf_notimedilution_Nconfs1126_Nsamples1_4qFlavorChange_part3.RData"
load(my_data_path)

#--------------------------------------------------------------------------------------------------------#
#   SET PLOTTING PARAMETERS                                                                              #
#--------------------------------------------------------------------------------------------------------#
Nconfs <- length(my_confs)
Nsrcs <- 1
Nsamples <- 1


##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##

##===================##
##                   ##
       x_u <- 6
##                   ##
##===================##

fit_res_tmp2 <- list()
fit_res <- list()

for ( i in 1:21 ) {
  for ( j in 0:(x_u - 1) ) {

    fit_res_tmp <- fit.const( y = ratio_6_id_ud_du[[i]]$cf.tsboot$t0,
                              bsamples = ratio_6_id_ud_du[[i]]$cf.tsboot$t,
                              x = seq( j, length(ratio_6_id_ud_du[[i]]$cf.tsboot$t0 - 1), 1 ),
                              x.lower = j,
                              x.upper = x_u,
                              start.par = c(1),
                              useCov = TRUE,
                              boot.fit = TRUE,
                              fit.method = "lm",
                              autoproceed = FALSE,
                              every,
                              cov_fn = cov,
                              error = sd )

    fit_res_tmp2[[j+1]] <- fit_res_tmp
  }
  fit_res[[i]] <- fit_res_tmp2
}

my_res                        <- list()
ratio_val                     <- c()
ratio_err                     <- c()
redchisqr                     <- c()
pval                          <- c()
res_ratio_6_id_ud_du <- list()

for ( i in 1:21 ) {
  my_res <- fit_res[[i]]
  for ( j in 0:(x_u - 1) ) {
    ratio_val[j+1] <- my_res[[j+1]]$t0
    ratio_err[j+1] <- my_res[[j+1]]$se
    redchisqr[j+1] <- my_res[[j+1]]$chisqr / my_res[[j+1]]$dof
    pval[j+1]      <- my_res[[j+1]]$Qval
  }

  res_ratio_6_id_ud_du[[i]] <- cbind( ratio_val, ratio_err, redchisqr, pval )

  write.table(x=res_ratio_6_id_ud_du[[i]], file=paste0("fitres_tu", x_u , "_ratio_6_id_ud_du_tgf", my_gft[i] , ".txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
}




mycols <- c( "black", "green", "blue", "red", "antiquewhite3",
             "aquamarine3", "azure3", "blue3", "bisque3", "blueviolet",
             "red3", "burlywood3", "brown3", "cadetblue3", "chartreuse3",
             "chocolate3", "cornflowerblue", "coral3", "cyan3", "darkblue",
             "darkgoldenrod3" )

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,3],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,2),
                       xlab="$t_l$",
                       ylab="$\\chi^2\\, /\\, \\mathrm{ndof}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,8.), y=c(1,1), se=0.025, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,3],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,4],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,1),
                       xlab="$t_l$",
                       ylab="$p\\,\\,\\mathrm{value}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,12.), y=c(0.5,0.5), se=0.025/6, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,4],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), standAlone = TRUE, width = 5, height = 2,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[1]][,1],
                       dy=res_ratio_6_id_ud_du[[1]][,2],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0.002,0.01),
                       xlab="$t_l$",
                       ylab="$\\langle O_{aa}^{qb}\\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf})$",
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[1] )

for ( i in 2:21 ) {
hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[i]][,1],
                       dy=res_ratio_6_id_ud_du[[i]][,2],
                       rep=TRUE,
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

#############################################################################
#############################################################################
#############################################################################

##===================##
##                   ##
       x_u <- 5
##                   ##
##===================##

fit_res_tmp2 <- list()
fit_res <- list()

for ( i in 1:21 ) {
  for ( j in 0:(x_u - 1) ) {

    fit_res_tmp <- fit.const( y = ratio_6_id_ud_du[[i]]$cf.tsboot$t0,
                              bsamples = ratio_6_id_ud_du[[i]]$cf.tsboot$t,
                              x = seq( j, length(ratio_6_id_ud_du[[i]]$cf.tsboot$t0 - 1), 1 ),
                              x.lower = j,
                              x.upper = x_u,
                              start.par = c(1),
                              useCov = TRUE,
                              boot.fit = TRUE,
                              fit.method = "lm",
                              autoproceed = FALSE,
                              every,
                              cov_fn = cov,
                              error = sd )

    fit_res_tmp2[[j+1]] <- fit_res_tmp
  }
  fit_res[[i]] <- fit_res_tmp2
}

my_res                        <- list()
ratio_val                     <- c()
ratio_err                     <- c()
redchisqr                     <- c()
pval                          <- c()
res_ratio_6_id_ud_du <- list()

for ( i in 1:21 ) {
  my_res <- fit_res[[i]]
  for ( j in 0:(x_u - 1) ) {
    ratio_val[j+1] <- my_res[[j+1]]$t0
    ratio_err[j+1] <- my_res[[j+1]]$se
    redchisqr[j+1] <- my_res[[j+1]]$chisqr / my_res[[j+1]]$dof
    pval[j+1]      <- my_res[[j+1]]$Qval
  }

  res_ratio_6_id_ud_du[[i]] <- cbind( ratio_val, ratio_err, redchisqr, pval )

  write.table(x=res_ratio_6_id_ud_du[[i]], file=paste0("fitres_tu", x_u , "_ratio_6_id_ud_du_tgf", my_gft[i] , ".txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
}




mycols <- c( "black", "green", "blue", "red", "antiquewhite3",
             "aquamarine3", "azure3", "blue3", "bisque3", "blueviolet",
             "red3", "burlywood3", "brown3", "cadetblue3", "chartreuse3",
             "chocolate3", "cornflowerblue", "coral3", "cyan3", "darkblue",
             "darkgoldenrod3" )

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,3],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,2),
                       xlab="$t_l$",
                       ylab="$\\chi^2\\, /\\, \\mathrm{ndof}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,8.), y=c(1,1), se=0.025, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,3],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,4],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,1),
                       xlab="$t_l$",
                       ylab="$p\\,\\,\\mathrm{value}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,12.), y=c(0.5,0.5), se=0.025/6, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,4],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), standAlone = TRUE, width = 5, height = 2,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[1]][,1],
                       dy=res_ratio_6_id_ud_du[[1]][,2],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0.002,0.01),
                       xlab="$t_l$",
                       ylab="$\\langle O_{aa}^{qb}\\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf})$",
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[1] )

for ( i in 2:21 ) {
hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[i]][,1],
                       dy=res_ratio_6_id_ud_du[[i]][,2],
                       rep=TRUE,
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

#############################################################################
#############################################################################
#############################################################################

##===================##
##                   ##
       x_u <- 4
##                   ##
##===================##

fit_res_tmp2 <- list()
fit_res <- list()

for ( i in 1:21 ) {
  for ( j in 0:(x_u - 1) ) {

    fit_res_tmp <- fit.const( y = ratio_6_id_ud_du[[i]]$cf.tsboot$t0,
                              bsamples = ratio_6_id_ud_du[[i]]$cf.tsboot$t,
                              x = seq( j, length(ratio_6_id_ud_du[[i]]$cf.tsboot$t0 - 1), 1 ),
                              x.lower = j,
                              x.upper = x_u,
                              start.par = c(1),
                              useCov = TRUE,
                              boot.fit = TRUE,
                              fit.method = "lm",
                              autoproceed = FALSE,
                              every,
                              cov_fn = cov,
                              error = sd )

    fit_res_tmp2[[j+1]] <- fit_res_tmp
  }
  fit_res[[i]] <- fit_res_tmp2
}

my_res                        <- list()
ratio_val                     <- c()
ratio_err                     <- c()
redchisqr                     <- c()
pval                          <- c()
res_ratio_6_id_ud_du <- list()

for ( i in 1:21 ) {
  my_res <- fit_res[[i]]
  for ( j in 0:(x_u - 1) ) {
    ratio_val[j+1] <- my_res[[j+1]]$t0
    ratio_err[j+1] <- my_res[[j+1]]$se
    redchisqr[j+1] <- my_res[[j+1]]$chisqr / my_res[[j+1]]$dof
    pval[j+1]      <- my_res[[j+1]]$Qval
  }

  res_ratio_6_id_ud_du[[i]] <- cbind( ratio_val, ratio_err, redchisqr, pval )

  write.table(x=res_ratio_6_id_ud_du[[i]], file=paste0("fitres_tu", x_u , "_ratio_6_id_ud_du_tgf", my_gft[i] , ".txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
}




mycols <- c( "black", "green", "blue", "red", "antiquewhite3",
             "aquamarine3", "azure3", "blue3", "bisque3", "blueviolet",
             "red3", "burlywood3", "brown3", "cadetblue3", "chartreuse3",
             "chocolate3", "cornflowerblue", "coral3", "cyan3", "darkblue",
             "darkgoldenrod3" )

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,3],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,2),
                       xlab="$t_l$",
                       ylab="$\\chi^2\\, /\\, \\mathrm{ndof}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,8.), y=c(1,1), se=0.025, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,3],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,4],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,1),
                       xlab="$t_l$",
                       ylab="$p\\,\\,\\mathrm{value}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,12.), y=c(0.5,0.5), se=0.025/6, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,4],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), standAlone = TRUE, width = 5, height = 2,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[1]][,1],
                       dy=res_ratio_6_id_ud_du[[1]][,2],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0.002,0.01),
                       xlab="$t_l$",
                       ylab="$\\langle O_{aa}^{qb}\\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf})$",
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[1] )

for ( i in 2:21 ) {
hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[i]][,1],
                       dy=res_ratio_6_id_ud_du[[i]][,2],
                       rep=TRUE,
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

#############################################################################
#############################################################################
#############################################################################

##===================##
##                   ##
       x_u <- 3
##                   ##
##===================##

fit_res_tmp2 <- list()
fit_res <- list()

for ( i in 1:21 ) {
  for ( j in 0:(x_u - 1) ) {

    fit_res_tmp <- fit.const( y = ratio_6_id_ud_du[[i]]$cf.tsboot$t0,
                              bsamples = ratio_6_id_ud_du[[i]]$cf.tsboot$t,
                              x = seq( j, length(ratio_6_id_ud_du[[i]]$cf.tsboot$t0 - 1), 1 ),
                              x.lower = j,
                              x.upper = x_u,
                              start.par = c(1),
                              useCov = TRUE,
                              boot.fit = TRUE,
                              fit.method = "lm",
                              autoproceed = FALSE,
                              every,
                              cov_fn = cov,
                              error = sd )

    fit_res_tmp2[[j+1]] <- fit_res_tmp
  }
  fit_res[[i]] <- fit_res_tmp2
}

my_res                        <- list()
ratio_val                     <- c()
ratio_err                     <- c()
redchisqr                     <- c()
pval                          <- c()
res_ratio_6_id_ud_du <- list()

for ( i in 1:21 ) {
  my_res <- fit_res[[i]]
  for ( j in 0:(x_u - 1) ) {
    ratio_val[j+1] <- my_res[[j+1]]$t0
    ratio_err[j+1] <- my_res[[j+1]]$se
    redchisqr[j+1] <- my_res[[j+1]]$chisqr / my_res[[j+1]]$dof
    pval[j+1]      <- my_res[[j+1]]$Qval
  }

  res_ratio_6_id_ud_du[[i]] <- cbind( ratio_val, ratio_err, redchisqr, pval )

  write.table(x=res_ratio_6_id_ud_du[[i]], file=paste0("fitres_tu", x_u , "_ratio_6_id_ud_du_tgf", my_gft[i] , ".txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
}




mycols <- c( "black", "green", "blue", "red", "antiquewhite3",
             "aquamarine3", "azure3", "blue3", "bisque3", "blueviolet",
             "red3", "burlywood3", "brown3", "cadetblue3", "chartreuse3",
             "chocolate3", "cornflowerblue", "coral3", "cyan3", "darkblue",
             "darkgoldenrod3" )

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,3],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,2),
                       xlab="$t_l$",
                       ylab="$\\chi^2\\, /\\, \\mathrm{ndof}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,8.), y=c(1,1), se=0.025, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,3],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,4],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,1),
                       xlab="$t_l$",
                       ylab="$p\\,\\,\\mathrm{value}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,12.), y=c(0.5,0.5), se=0.025/6, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,4],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), standAlone = TRUE, width = 5, height = 2,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[1]][,1],
                       dy=res_ratio_6_id_ud_du[[1]][,2],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0.002,0.01),
                       xlab="$t_l$",
                       ylab="$\\langle O_{aa}^{qb}\\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf})$",
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[1] )

for ( i in 2:21 ) {
hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[i]][,1],
                       dy=res_ratio_6_id_ud_du[[i]][,2],
                       rep=TRUE,
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

#############################################################################
#############################################################################
#############################################################################

##===================##
##                   ##
       x_u <- 2
##                   ##
##===================##

fit_res_tmp2 <- list()
fit_res <- list()

for ( i in 1:21 ) {
  for ( j in 0:(x_u - 1) ) {

    fit_res_tmp <- fit.const( y = ratio_6_id_ud_du[[i]]$cf.tsboot$t0,
                              bsamples = ratio_6_id_ud_du[[i]]$cf.tsboot$t,
                              x = seq( j, length(ratio_6_id_ud_du[[i]]$cf.tsboot$t0 - 1), 1 ),
                              x.lower = j,
                              x.upper = x_u,
                              start.par = c(1),
                              useCov = TRUE,
                              boot.fit = TRUE,
                              fit.method = "lm",
                              autoproceed = FALSE,
                              every,
                              cov_fn = cov,
                              error = sd )

    fit_res_tmp2[[j+1]] <- fit_res_tmp
  }
  fit_res[[i]] <- fit_res_tmp2
}

my_res                        <- list()
ratio_val                     <- c()
ratio_err                     <- c()
redchisqr                     <- c()
pval                          <- c()
res_ratio_6_id_ud_du <- list()

for ( i in 1:21 ) {
  my_res <- fit_res[[i]]
  for ( j in 0:(x_u - 1) ) {
    ratio_val[j+1] <- my_res[[j+1]]$t0
    ratio_err[j+1] <- my_res[[j+1]]$se
    redchisqr[j+1] <- my_res[[j+1]]$chisqr / my_res[[j+1]]$dof
    pval[j+1]      <- my_res[[j+1]]$Qval
  }

  res_ratio_6_id_ud_du[[i]] <- cbind( ratio_val, ratio_err, redchisqr, pval )

  write.table(x=res_ratio_6_id_ud_du[[i]], file=paste0("fitres_tu", x_u , "_ratio_6_id_ud_du_tgf", my_gft[i] , ".txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
}




mycols <- c( "black", "green", "blue", "red", "antiquewhite3",
             "aquamarine3", "azure3", "blue3", "bisque3", "blueviolet",
             "red3", "burlywood3", "brown3", "cadetblue3", "chartreuse3",
             "chocolate3", "cornflowerblue", "coral3", "cyan3", "darkblue",
             "darkgoldenrod3" )

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,3],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,2),
                       xlab="$t_l$",
                       ylab="$\\chi^2\\, /\\, \\mathrm{ndof}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,8.), y=c(1,1), se=0.025, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,3],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_redchisqr.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), standAlone = TRUE, width = 5, height = 2)
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0-0.4, (x_u - 1)-0.4,1),
                       y=res_ratio_6_id_ud_du[[1]][,4],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0,1),
                       xlab="$t_l$",
                       ylab="$p\\,\\,\\mathrm{value}$",
                       pch=19, cex=0.5, col=mycols[1] )

plt.errband( x=c(-1.,12.), y=c(0.5,0.5), se=0.025/6, lwd = 0., col="red", a=0.5 )

for ( i in 1:21 ) {
hadron::plotwitherror( x=seq(0+(i-2)*0.04-0.4, (x_u - 1)+(i-2)*0.04-0.4,1),
                       y=res_ratio_6_id_ud_du[[i]][,4],
                       rep=TRUE,
                       pch=19, cex=0.5, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du_pval.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#

tikzDevice::tikz(file = paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), standAlone = TRUE, width = 5, height = 2,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.4, 5, 1.5, 0.5))

hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[1]][,1],
                       dy=res_ratio_6_id_ud_du[[1]][,2],
                       xlim=c(-0.5,x_u - 0.5),
                       ylim=c(0.002,0.01),
                       xlab="$t_l$",
                       ylab="$\\langle O_{aa}^{qb}\\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf})$",
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[1] )

for ( i in 2:21 ) {
hadron::plotwitherror( x=seq(0,x_u - 1,1),
                       y=res_ratio_6_id_ud_du[[i]][,1],
                       dy=res_ratio_6_id_ud_du[[i]][,2],
                       rep=TRUE,
                       pch=19, cex=0.5, cex.lab=0.9, col=mycols[i] )
}

dev.off()
# Compile the tex file
tools::texi2dvi(paste0("fitres_tu", x_u, "_ratio_6_id_ud_du.tex"), pdf=TRUE)

#---------------------------------------------------------------------------#