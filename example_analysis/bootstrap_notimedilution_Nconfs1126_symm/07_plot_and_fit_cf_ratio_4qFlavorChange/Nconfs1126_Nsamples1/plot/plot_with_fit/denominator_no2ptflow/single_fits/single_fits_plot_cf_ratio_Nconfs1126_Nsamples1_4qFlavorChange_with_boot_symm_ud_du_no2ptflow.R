##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##

#--------------------------------------------------------------------------------------------------------#
#   LOAD PACKAGES AND CODE FOR ANALYSIS                                                                  #
#--------------------------------------------------------------------------------------------------------#
library(stringr)   ## required for str_pad()
source("/home/nikolas/Documents/NPV/njjn_analysis/npv_analysis_jan_2023/npv.R")

library(hadron)
library(devtools)
library(latex2exp)
library(tidyverse)
library(tikzDevice)

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   LOAD WORKSPACE IMAGE                                                                                 #
#                                                                                                        #
#      * containing the cf ratio ratio_1_g5_ud = Im(<X1g5_ud>) / Re(<Xmx_g5_ud>)                         #
#      * containing the cf ratio ratio_2_g5_ud = Im(<X2g5_ud>) / Re(<Xmx_g5_ud>)                         #
#      * containing the cf ratio ratio_3_g5_ud = Im(<X3g5_ud_mulmin1>) / Re(<Xmx_g5_ud>)                 #
#                                                                                                        #
#      * containing the cf ratio ratio_4_id_ud = Re(<X1id_ud_mulmin1>) / Re(<Xmx_id_ud_mulmin1>)         #
#      * containing the cf ratio ratio_5_id_ud = Re(<X2id_ud_mulmin1>) / Re(<Xmx_id_ud_mulmin1>)         #
#      * containing the cf ratio ratio_6_id_ud = Re(<X3id_ud>) / Re(<Xmx_id_ud_mulmin1>)                 #
#                                                                                                        #
#      * containing the cf ratio ratio_1_g5_du = Im(<X1g5_du>) / Re(<Xmx_g5_du>)                         #
#      * containing the cf ratio ratio_2_g5_du = Im(<X2g5_du>) / Re(<Xmx_g5_du>)                         #
#      * containing the cf ratio ratio_3_g5_du = Im(<X3g5_du_mulmin1>) / Re(<Xmx_g5_du>)                 #
#                                                                                                        #
#      * containing the cf ratio ratio_4_id_du = Re(<X1id_du_mulmin1>) / Re(<Xmx_id_du_mulmin1>)         #
#      * containing the cf ratio ratio_5_id_du = Re(<X2id_du_mulmin1>) / Re(<Xmx_id_du_mulmin1>)         #
#      * containing the cf ratio ratio_6_id_du = Re(<X3id_du>) / Re(<Xmx_id_du_mulmin1>)                 #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
my_data_path <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/01_boot_cf_data/NLoopsamples_01/01_boot_cf_data_Nconfs1126_4qFlavorChange/bootstrap_cf_mx_prb_oet_gf_notimedilution_Nconfs1126_Nsamples1_4qFlavorChange_part3.RData"
load(my_data_path)

print(ratio_4_id_ud_du_no2ptflow[[1]]$bsamples)

stop()

#--------------------------------------------------------------------------------------------------------#
#   SET PLOTTING PARAMETERS                                                                              #
#--------------------------------------------------------------------------------------------------------#
my_cols     <- c( "black", "green", "blue", "red", "antiquewhite3",
                  "aquamarine3", "azure3", "blue3", "bisque3", "blueviolet",
                  "red3", "burlywood3", "brown3", "cadetblue3", "chartreuse3",
                  "chocolate3", "cornflowerblue", "coral3", "cyan3", "darkblue",
                  "darkgoldenrod3" )

my_legend   <- c( "$t_\\mathrm{gf} = 0.00$", "$t_\\mathrm{gf} = 0.01$", "$t_\\mathrm{gf} = 0.02$",
                  "$t_\\mathrm{gf} = 0.03$", "$t_\\mathrm{gf} = 0.04$", "$t_\\mathrm{gf} = 0.05$",
                  "$t_\\mathrm{gf} = 0.06$", "$t_\\mathrm{gf} = 0.07$", "$t_\\mathrm{gf} = 0.08$",
                  "$t_\\mathrm{gf} = 0.09$", "$t_\\mathrm{gf} = 0.10$", "$t_\\mathrm{gf} = 0.11$",
                  "$t_\\mathrm{gf} = 0.12$", "$t_\\mathrm{gf} = 0.13$", "$t_\\mathrm{gf} = 0.14$",
                  "$t_\\mathrm{gf} = 0.15$", "$t_\\mathrm{gf} = 0.16$", "$t_\\mathrm{gf} = 0.17$",
                  "$t_\\mathrm{gf} = 0.18$", "$t_\\mathrm{gf} = 0.19$", "$t_\\mathrm{gf} = 0.20$" )

Nconfs <- length(my_confs)
Nsrcs <- 1
Nsamples <- 1


##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
#                                                                                                                                                 #
# PLOTS OF RATIOS WITHOUT!!! FLOWED 2PT-FUNCTION IN THE DENOMINATOR:                                                                              #
#                                                                                                                                                 #
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

print_g5 <- FALSE

if ( print_g5 == TRUE ) {

tikzDevice::tikz( file = "cf_ratio1_vvqb_g5_ud_du_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_1_g5_ud_du_no2ptflow[[1]]$cf0, dy=ratio_1_g5_ud_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{vv}^{qb} \\, O_{\\gamma_5} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud+du} (0)$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=0.95 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_1_g5_ud_du_no2ptflow[[i]]$cf0, dy=ratio_1_g5_ud_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=0.95, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio1_vvqb_g5_ud_du_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio2_vvcc_g5_ud_du_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_2_g5_ud_du_no2ptflow[[1]]$cf0, dy=ratio_2_g5_ud_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-g5}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{vv}^{cc} \\, O_{\\gamma_5} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud+du} (0)$",
               xlim=c(0,40),
               ylim=c(0,7),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=0.95 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_2_g5_ud_du_no2ptflow[[i]]$cf0, dy=ratio_2_g5_ud_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=0.95, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio2_vvcc_g5_ud_du_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio3_aaqb_g5_ud_du_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_3_g5_ud_du_no2ptflow[[1]]$cf0, dy=ratio_3_g5_ud_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-aa-f1-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{aa}^{qb} \\, O_{\\gamma_5} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud+du} (0)$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=0.95 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_3_g5_ud_du_no2ptflow[[i]]$cf0, dy=ratio_3_g5_ud_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=0.95, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio3_aaqb_g5_ud_du_no2ptflow.tex", pdf=TRUE )

}

###################################################################################################################################################
###################################################################################################################################################

plt.errband1 <- function (x, y, se1, se2, lwd = 1.5, col="gray", a=0.25) {
  polyx <- c( x, rev(x) )
  polyval <- c( (y + se2), rev(y - se1) )
  ## for a = 0.25 the opacity is 25%
  polycol <- alpha( col, a )
  
  polygon( x=polyx, y=polyval, col=polycol, border=FALSE )
  lines( x, y, col=col, lwd=lwd )
}

plt.errband2 <- function (x, y, se1, se2, lwd = 1.5, col="gray", a=0.25) {
  polyx <- c( (x + se2), rev(x - se1) )
  polyval <- c( y, rev(y) )
  ## for a = 0.25 the opacity is 25%
  polycol <- alpha( col, a )

  polygon( x=polyx, y=polyval, col=polycol, border=FALSE )
  lines( x, y, col=col, lwd=lwd )
}

plt.errband3 <- function (x, y, se1, se2, lwd = 1.5, col="gray", a=0.25) {
  polyx <- c( x, rev(x) )
  polyval1 <- c( (y + se1), rev(y - se1) )
  polyval2 <- c( (y + se2), rev(y - se2) )
  ## for a = 0.25 the opacity is 25%
  polycol1 <- alpha( col, a )
  polycol2 <- alpha( col, a*1.2 )
  
  polygon( x=polyx, y=polyval1, col=polycol1, border=FALSE )
  polygon( x=polyx, y=polyval2, col=polycol2, border=FALSE )
  lines( x, y, col=col, lwd=lwd )
}


tgf_1 <- seq(4,19,3)
tgf_2 <- seq(5,20,3)
tgf_3 <- seq(6,21,3)

shift <- 0.1

my_cols_new <- c( rep("blue",3),rep("red",3),rep("forestgreen",3),
                  rep("gold",3),rep("magenta4",3),rep("darkorange1",3),rep("black",3) )

res1 <- read.table("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/07_plot_and_fit_cf_ratio_4qFlavorChange/Nconfs1126_Nsamples1/fit/cf_ratio4_vvqb_id_ud_du_no2ptflow/AIC_fit_results_cf_ratio_4_vvqb_id_ud_du_no2ptflow.txt")






fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_part1 <- list()
fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_part2 <- list()
fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_part3 <- list()
fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_part4 <- list()
fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_part5 <- list()

for ( i in 1:21 ) {

  fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_xlower0[[i]] <- fit.const( y = ratio_4_id_ud_du_no2ptflow[[i]]$cf0,
                                                                   bsamples = ratio_4_id_ud_du_no2ptflow[[i]],
                                                                   x = seq(0,10,1),
                                                                   x.lower = 0, x.upper = 6,
                                                                   start.par = 1.,
                                                                   useCov = TRUE,
                                                                   boot.fit = TRUE,
                                                                   fit.method = "optim",
                                                                   cov_fn = cov,
                                                                   error = sd )

  fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_xlower1[[i]] <- fit.const( y = ratio_4_id_ud_du_no2ptflow[[i]]$cf0,
                                                                   bsamples = ratio_4_id_ud_du_no2ptflow[[i]],
                                                                   x = seq(0,10,1),
                                                                   x.lower = 1, x.upper = 6,
                                                                   start.par = 1.,
                                                                   useCov = TRUE,
                                                                   boot.fit = TRUE,
                                                                   fit.method = "optim",
                                                                   cov_fn = cov,
                                                                   error = sd )

  fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_xlower2[[i]] <- fit.const( y = ratio_4_id_ud_du_no2ptflow[[i]]$cf0,
                                                                   bsamples = ratio_4_id_ud_du_no2ptflow[[i]],
                                                                   x = seq(0,10,1),
                                                                   x.lower = 2, x.upper = 6,
                                                                   start.par = 1.,
                                                                   useCov = TRUE,
                                                                   boot.fit = TRUE,
                                                                   fit.method = "optim",
                                                                   cov_fn = cov,
                                                                   error = sd )

  fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_xlower3[[i]] <- fit.const( y = ratio_4_id_ud_du_no2ptflow[[i]]$cf0,
                                                                   bsamples = ratio_4_id_ud_du_no2ptflow[[i]],
                                                                   x = seq(0,10,1),
                                                                   x.lower = 3, x.upper = 6,
                                                                   start.par = 1.,
                                                                   useCov = TRUE,
                                                                   boot.fit = TRUE,
                                                                   fit.method = "optim",
                                                                   cov_fn = cov,
                                                                   error = sd )

  fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_xlower4[[i]] <- fit.const( y = ratio_4_id_ud_du_no2ptflow[[i]]$cf0,
                                                                   bsamples = ratio_4_id_ud_du_no2ptflow[[i]],
                                                                   x = seq(0,10,1),
                                                                   x.lower = 4, x.upper = 6,
                                                                   start.par = 1.,
                                                                   useCov = TRUE,
                                                                   boot.fit = TRUE,
                                                                   fit.method = "optim",
                                                                   cov_fn = cov,
                                                                   error = sd )

}






print(fit_cf_ratio4_vvqb_id_ud_du_no2ptflow_xlower4)








stop()









tikzDevice::tikz( file = "cf_ratio4_vvqb_id_ud_du_no2ptflow_part1.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_4_id_ud_du_no2ptflow[[1]]$cf0[1:11], dy=ratio_4_id_ud_du_no2ptflow[[1]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{vv}^{qb} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.001,0.009),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,6+shift*3,1), y=rep(res1[1,2],5), se1=rep(res1[1,3],5), se2=rep(res1[1,4],5), lwd = 1.5, col=my_cols_new[1])
j <- 2
for ( i in tgf_1 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_4_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_4_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,6+shift*j,1), y=rep(res1[i,2],5), se1=rep(res1[i,3],5), se2=rep(res1[i,4],5), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[1],my_legend[tgf_1]), col=c(my_cols_new[1],my_cols_new[tgf_1]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio4_vvqb_id_ud_du_no2ptflow_part1.tex", pdf=TRUE )


tikzDevice::tikz( file = "cf_ratio4_vvqb_id_ud_du_no2ptflow_part2.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_4_id_ud_du_no2ptflow[[2]]$cf0[1:11], dy=ratio_4_id_ud_du_no2ptflow[[2]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{vv}^{qb} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.001,0.0085),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,6+shift*3,1), y=rep(res1[2,2],5), se1=rep(res1[2,3],5), se2=rep(res1[2,4],5), lwd = 1.5, col=my_cols_new[2])
j <- 2
for ( i in tgf_2 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_4_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_4_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,6+shift*j,1), y=rep(res1[i,2],5), se1=rep(res1[i,3],5), se2=rep(res1[i,4],5), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[2],my_legend[tgf_2]), col=c(my_cols_new[2],my_cols_new[tgf_2]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio4_vvqb_id_ud_du_no2ptflow_part2.tex", pdf=TRUE )


tikzDevice::tikz( file = "cf_ratio4_vvqb_id_ud_du_no2ptflow_part3.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_4_id_ud_du_no2ptflow[[3]]$cf0[1:11], dy=ratio_4_id_ud_du_no2ptflow[[3]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{vv}^{qb} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.001,0.0075),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,6+shift*3,1), y=rep(res1[3,2],5), se1=rep(res1[3,3],5), se2=rep(res1[3,4],5), lwd = 1.5, col=my_cols_new[3])
j <- 2
for ( i in tgf_3 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_4_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_4_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,6+shift*j,1), y=rep(res1[i,2],5), se1=rep(res1[i,3],5), se2=rep(res1[i,4],5), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[3],my_legend[tgf_3]), col=c(my_cols_new[3],my_cols_new[tgf_3]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio4_vvqb_id_ud_du_no2ptflow_part3.tex", pdf=TRUE )

###################################################################################################################################################

res2 <- read.table("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/07_plot_and_fit_cf_ratio_4qFlavorChange/Nconfs1126_Nsamples1/fit/cf_ratio5_vvcc_id_ud_du_no2ptflow/AIC_fit_results_cf_ratio_5_vvcc_id_ud_du_no2ptflow.txt")

tikzDevice::tikz( file = "cf_ratio5_vvcc_id_ud_du_no2ptflow_part1.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_5_id_ud_du_no2ptflow[[1]]$cf0[1:11], dy=ratio_5_id_ud_du_no2ptflow[[1]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{vv}^{cc} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.003,0.028),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,6+shift*3,1), y=rep(res2[1,2],5), se1=rep(res2[1,3],5), se2=rep(res2[1,4],5), lwd = 1.5, col=my_cols_new[1])
j <- 2
for ( i in tgf_1 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_5_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_5_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,6+shift*j,1), y=rep(res2[i,2],5), se1=rep(res2[i,3],5), se2=rep(res2[i,4],5), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[1],my_legend[tgf_1]), col=c(my_cols_new[1],my_cols_new[tgf_1]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio5_vvcc_id_ud_du_no2ptflow_part1.tex", pdf=TRUE )


tikzDevice::tikz( file = "cf_ratio5_vvcc_id_ud_du_no2ptflow_part2.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_5_id_ud_du_no2ptflow[[2]]$cf0[1:11], dy=ratio_5_id_ud_du_no2ptflow[[2]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{vv}^{cc} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.004,0.024),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,6+shift*3,1), y=rep(res2[2,2],5), se1=rep(res2[2,3],5), se2=rep(res2[2,4],5), lwd = 1.5, col=my_cols_new[2])
j <- 2
for ( i in tgf_2 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_5_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_5_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,6+shift*j,1), y=rep(res2[i,2],5), se1=rep(res2[i,3],5), se2=rep(res2[i,4],5), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[2],my_legend[tgf_2]), col=c(my_cols_new[2],my_cols_new[tgf_2]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio5_vvcc_id_ud_du_no2ptflow_part2.tex", pdf=TRUE )


tikzDevice::tikz( file = "cf_ratio5_vvcc_id_ud_du_no2ptflow_part3.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_5_id_ud_du_no2ptflow[[3]]$cf0[1:11], dy=ratio_5_id_ud_du_no2ptflow[[3]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{vv}^{cc} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.003,0.023),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,6+shift*3,1), y=rep(res2[3,2],5), se1=rep(res2[3,3],5), se2=rep(res2[3,4],5), lwd = 1.5, col=my_cols_new[3])
j <- 2
for ( i in tgf_3 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_5_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_5_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,6+shift*j,1), y=rep(res2[i,2],5), se1=rep(res2[i,3],5), se2=rep(res2[i,4],5), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[3],my_legend[tgf_3]), col=c(my_cols_new[3],my_cols_new[tgf_3]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio5_vvcc_id_ud_du_no2ptflow_part3.tex", pdf=TRUE )

###################################################################################################################################################

res3 <- read.table("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/07_plot_and_fit_cf_ratio_4qFlavorChange/Nconfs1126_Nsamples1/fit/cf_ratio6_aaqb_id_ud_du_no2ptflow/AIC_fit_results_cf_ratio_6_aaqb_id_ud_du_no2ptflow.txt")

tikzDevice::tikz( file = "cf_ratio6_aaqb_id_ud_du_no2ptflow_part1.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_6_id_ud_du_no2ptflow[[1]]$cf0[1:11], dy=ratio_6_id_ud_du_no2ptflow[[1]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-aa-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{aa}^{qb} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.001,0.009),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,8+shift*3,1), y=rep(res3[1,2],7), se1=rep(res3[1,3],7), se2=rep(res3[1,4],7), lwd = 1.5, col=my_cols_new[1])
j <- 2
for ( i in tgf_1 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_6_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_6_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,8+shift*j,1), y=rep(res3[i,2],7), se1=rep(res3[i,3],7), se2=rep(res3[i,4],7), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[1],my_legend[tgf_1]), col=c(my_cols_new[1],my_cols_new[tgf_1]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio6_aaqb_id_ud_du_no2ptflow_part1.tex", pdf=TRUE )


tikzDevice::tikz( file = "cf_ratio6_aaqb_id_ud_du_no2ptflow_part2.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_6_id_ud_du_no2ptflow[[2]]$cf0[1:11], dy=ratio_6_id_ud_du_no2ptflow[[2]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-aa-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{aa}^{qb} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.001,0.0085),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,8+shift*3,1), y=rep(res3[2,2],7), se1=rep(res3[2,3],7), se2=rep(res3[2,4],7), lwd = 1.5, col=my_cols_new[2])
j <- 2
for ( i in tgf_2 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_6_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_6_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,8+shift*j,1), y=rep(res3[i,2],7), se1=rep(res3[i,3],7), se2=rep(res3[i,4],7), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[2],my_legend[tgf_2]), col=c(my_cols_new[2],my_cols_new[tgf_2]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio6_aaqb_id_ud_du_no2ptflow_part2.tex", pdf=TRUE )


tikzDevice::tikz( file = "cf_ratio6_aaqb_id_ud_du_no2ptflow_part3.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0+shift*3,10+shift*3,1), y=ratio_6_id_ud_du_no2ptflow[[3]]$cf0[1:11], dy=ratio_6_id_ud_du_no2ptflow[[3]]$tsboot.se[1:11],
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-aa-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\frac{1}{2}\\langle O_{aa}^{qb} \\, O_\\mathbb{1} \\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud+du} (0)$",
               xlim=c(0,12.5),
               ylim=c(0.001,0.0075),
               pch=19, cex=0.5, col=my_cols_new[1], cex.axis=1., cex.lab=1. )
plt.errband1(x=seq(2+shift*3,8+shift*3,1), y=rep(res3[3,2],7), se1=rep(res3[3,3],7), se2=rep(res3[3,4],7), lwd = 1.5, col=my_cols_new[3])
j <- 2
for ( i in tgf_3 ){
  plotwitherror( x=seq(0+shift*j,10+shift*j,1), y=ratio_6_id_ud_du_no2ptflow[[i]]$cf0[1:11], dy=ratio_6_id_ud_du_no2ptflow[[i]]$tsboot.se[1:11],
                 pch=19, cex=0.5, col=my_cols_new[i], cex.axis=1., cex.lab=1., rep="TRUE" )
  plt.errband1(x=seq(2+shift*j,8+shift*j,1), y=rep(res3[i,2],7), se1=rep(res3[i,3],7), se2=rep(res3[i,4],7), lwd = 1.5, col=my_cols_new[i])
  j <- j - 1
}
legend( "bottomright", legend=c(my_legend[3],my_legend[tgf_3]), col=c(my_cols_new[3],my_cols_new[tgf_3]), pch=19, cex=0.7 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.7)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio6_aaqb_id_ud_du_no2ptflow_part3.tex", pdf=TRUE )

##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
