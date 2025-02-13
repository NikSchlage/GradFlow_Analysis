##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##

#--------------------------------------------------------------------------------------------------------#
#   LOAD PACKAGES AND CODE FOR ANALYSIS                                                                  #
#--------------------------------------------------------------------------------------------------------#
library(stringr)   ## required for str_pad()
source("/home/nikolas/Documents/NPV/njjn_analysis/npv_analysis_jan_2023/npv.R")

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
my_data_path <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/01_boot_cf_data/NLoopsamples_01/01_boot_cf_data_Nconfs1126_4qFlavorChange/bootstrap_cf_mx_prb_oet_gf_notimedilution_Nconfs1126_Nsamples1_4qFlavorChange_part2.RData"
load(my_data_path)

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
# PLOTS OF RATIOS WITH!!! FLOWED 2PT-FUNCTION IN THE DENOMINATOR:                                                                                 #
#                                                                                                                                                 #
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

#========================================================================================================#
#--------------------------------------------------------------------------------------------------------#
#   Plot the cf ratios for keys "f0-[...]-f1-[...]"                                                      #
#--------------------------------------------------------------------------------------------------------#
#========================================================================================================#
tikzDevice::tikz( file = "cf_ratio1_f0-vv-f1-g5_qb_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_1_g5_ud[[1]]$cf0, dy=ratio_1_g5_ud[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{vv}^{qb} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_1_g5_ud[[i]]$cf0, dy=ratio_1_g5_ud[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio1_f0-vv-f1-g5_qb_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio2_f0-vv-f1-g5_cc_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_2_g5_ud[[1]]$cf0, dy=ratio_2_g5_ud[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-g5}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{vv}^{cc} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,7),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_2_g5_ud[[i]]$cf0, dy=ratio_2_g5_ud[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio2_f0-vv-f1-g5_cc_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio3_f0-aa-f1-g5_qb_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_3_g5_ud[[1]]$cf0, dy=ratio_3_g5_ud[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-aa-f1-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{aa}^{qb} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_3_g5_ud[[i]]$cf0, dy=ratio_3_g5_ud[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio3_f0-aa-f1-g5_qb_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio4_f0-vv-f1-id_qb_Nsamples8_nolog.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_4_id_ud[[1]]$cf0, dy=ratio_4_id_ud[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{vv}^{qb} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,0.01),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_4_id_ud[[i]]$cf0, dy=ratio_4_id_ud[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio4_f0-vv-f1-id_qb_Nsamples8_nolog.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio5_f0-vv-f1-id_cc_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_5_id_ud[[1]]$cf0, dy=ratio_5_id_ud[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{vv}^{cc} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,0.04),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_5_id_ud[[i]]$cf0, dy=ratio_5_id_ud[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio5_f0-vv-f1-id_cc_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio6_f0-aa-f1-id_qb_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_6_id_ud[[1]]$cf0, dy=ratio_6_id_ud[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-aa-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{aa}^{qb} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,0.01),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_6_id_ud[[i]]$cf0, dy=ratio_6_id_ud[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio6_f0-aa-f1-id_qb_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

#========================================================================================================#
#--------------------------------------------------------------------------------------------------------#
#   Plot the cf ratios for keys "f1-[...]-f0-[...]"                                                      #
#--------------------------------------------------------------------------------------------------------#
#========================================================================================================#
tikzDevice::tikz( file = "cf_ratio7_f1-vv-f0-g5_qb_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_1_g5_du[[1]]$cf0, dy=ratio_1_g5_du[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-vv-f0-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f1-g5-f0-g5}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{vv}^{qb} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_1_g5_du[[i]]$cf0, dy=ratio_1_g5_du[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio7_f1-vv-f0-g5_qb_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio8_f1-vv-f0-g5_cc_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_2_g5_du[[1]]$cf0, dy=ratio_2_g5_du[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-vv-f0-g5}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f1-g5-f0-g5}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{vv}^{cc} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,7),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_2_g5_du[[i]]$cf0, dy=ratio_2_g5_du[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio8_f1-vv-f0-g5_cc_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio9_f1-aa-f0-g5_qb_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_3_g5_du[[1]]$cf0, dy=ratio_3_g5_du[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-aa-f0-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f1-g5-f0-g5}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{aa}^{qb} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_3_g5_du[[i]]$cf0, dy=ratio_3_g5_du[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio9_f1-aa-f0-g5_qb_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio10_f1-vv-f0-id_qb_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_4_id_du[[1]]$cf0, dy=ratio_4_id_du[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-vv-f0-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f1-id-f0-id}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{vv}^{qb} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,0.01),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_4_id_du[[i]]$cf0, dy=ratio_4_id_du[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio10_f1-vv-f0-id_qb_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio11_f1-vv-f0-id_cc_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_5_id_du[[1]]$cf0, dy=ratio_5_id_du[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-vv-f0-id}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f1-id-f0-id}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{vv}^{cc} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,0.04),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_5_id_du[[i]]$cf0, dy=ratio_5_id_du[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio11_f1-vv-f0-id_cc_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio12_f1-aa-f0-id_qb_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_6_id_du[[1]]$cf0, dy=ratio_6_id_du[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-aa-f0-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f1-id-f0-id}^\\mathrm{2pt}(t_\\mathrm{gf}) $",
               ylab="$\\langle O_{aa}^{qb} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf})$",
               xlim=c(0,40),
               ylim=c(0,0.01),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_6_id_du[[i]]$cf0, dy=ratio_6_id_du[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio12_f1-aa-f0-id_qb_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################


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

#========================================================================================================#
#--------------------------------------------------------------------------------------------------------#
#   Plot the cf ratios for keys "f0-[...]-f1-[...]"                                                      #
#--------------------------------------------------------------------------------------------------------#
#========================================================================================================#
tikzDevice::tikz( file = "cf_ratio1_f0-vv-f1-g5_qb_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_1_g5_ud_no2ptflow[[1]]$cf0, dy=ratio_1_g5_ud_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{vv}^{qb} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud} (0)$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_1_g5_ud_no2ptflow[[i]]$cf0, dy=ratio_1_g5_ud_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio1_f0-vv-f1-g5_qb_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio2_f0-vv-f1-g5_cc_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_2_g5_ud_no2ptflow[[1]]$cf0, dy=ratio_2_g5_ud_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-g5}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{vv}^{cc} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud} (0)$",
               xlim=c(0,40),
               ylim=c(0,7),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_2_g5_ud_no2ptflow[[i]]$cf0, dy=ratio_2_g5_ud_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio2_f0-vv-f1-g5_cc_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio3_f0-aa-f1-g5_qb_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_3_g5_ud_no2ptflow[[1]]$cf0, dy=ratio_3_g5_ud_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-aa-f1-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-g5-f1-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{aa}^{qb} \\, O_{\\gamma_5} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{ud} (0)$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_3_g5_ud_no2ptflow[[i]]$cf0, dy=ratio_3_g5_ud_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio3_f0-aa-f1-g5_qb_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio4_f0-vv-f1-id_qb_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_4_id_ud_no2ptflow[[1]]$cf0, dy=ratio_4_id_ud_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{vv}^{qb} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud} (0)$",
               xlim=c(0,40),
               ylim=c(0,0.01),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_4_id_ud_no2ptflow[[i]]$cf0, dy=ratio_4_id_ud_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio4_f0-vv-f1-id_qb_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio5_f0-vv-f1-id_cc_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_5_id_ud_no2ptflow[[1]]$cf0, dy=ratio_5_id_ud_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-vv-f1-id}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{vv}^{cc} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud} (0)$",
               xlim=c(0,40),
               ylim=c(0,0.04),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_5_id_ud_no2ptflow[[i]]$cf0, dy=ratio_5_id_ud_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio5_f0-vv-f1-id_cc_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio6_f0-aa-f1-id_qb_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_6_id_ud_no2ptflow[[1]]$cf0, dy=ratio_6_id_ud_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f0-aa-f1-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f0-id-f1-id}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{aa}^{qb} \\, O_\\mathbb{1} \\rangle_{ud} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{ud} (0)$",
               xlim=c(0,40),
               ylim=c(0,0.01),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_6_id_ud_no2ptflow[[i]]$cf0, dy=ratio_6_id_ud_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio6_f0-aa-f1-id_qb_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

#========================================================================================================#
#--------------------------------------------------------------------------------------------------------#
#   Plot the cf ratios for keys "f1-[...]-f0-[...]"                                                      #
#--------------------------------------------------------------------------------------------------------#
#========================================================================================================#
tikzDevice::tikz( file = "cf_ratio7_f1-vv-f0-g5_qb_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_1_g5_du_no2ptflow[[1]]$cf0, dy=ratio_1_g5_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-vv-f0-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f1-g5-f0-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{vv}^{qb} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{du} (0)$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_1_g5_du_no2ptflow[[i]]$cf0, dy=ratio_1_g5_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio7_f1-vv-f0-g5_qb_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio8_f1-vv-f0-g5_cc_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_2_g5_du_no2ptflow[[1]]$cf0, dy=ratio_2_g5_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-vv-f0-g5}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f1-g5-f0-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{vv}^{cc} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{du} (0)$",
               xlim=c(0,40),
               ylim=c(0,7),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_2_g5_du_no2ptflow[[i]]$cf0, dy=ratio_2_g5_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio8_f1-vv-f0-g5_cc_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio9_f1-aa-f0-g5_qb_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_3_g5_du_no2ptflow[[1]]$cf0, dy=ratio_3_g5_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-aa-f0-g5}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f1-g5-f0-g5}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{aa}^{qb} \\, O_{\\gamma_5} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_{\\gamma_5} \\, O_{\\gamma_5} \\rangle_{du} (0)$",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_3_g5_du_no2ptflow[[i]]$cf0, dy=ratio_3_g5_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio9_f1-aa-f0-g5_qb_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio10_f1-vv-f0-id_qb_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_4_id_du_no2ptflow[[1]]$cf0, dy=ratio_4_id_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-vv-f0-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f1-id-f0-id}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{vv}^{qb} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{du} (0)$",
               xlim=c(0,40),
               ylim=c(0,0.01),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_4_id_du_no2ptflow[[i]]$cf0, dy=ratio_4_id_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio10_f1-vv-f0-id_qb_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio11_f1-vv-f0-id_cc_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_5_id_du_no2ptflow[[1]]$cf0, dy=ratio_5_id_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-vv-f0-id}^\\mathrm{2pt,cc}(t_\\mathrm{gf}) / C_\\mathrm{f1-id-f0-id}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{vv}^{cc} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{du} (0)$",
               xlim=c(0,40),
               ylim=c(0,0.04),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_5_id_du_no2ptflow[[i]]$cf0, dy=ratio_5_id_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio11_f1-vv-f0-id_cc_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio12_f1-aa-f0-id_qb_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_6_id_du_no2ptflow[[1]]$cf0, dy=ratio_6_id_du_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$C_\\mathrm{f1-aa-f0-id}^\\mathrm{2pt,qb}(t_\\mathrm{gf}) / C_\\mathrm{f1-id-f0-id}^\\mathrm{2pt}(0) $",
               ylab="$\\langle O_{aa}^{qb} \\, O_\\mathbb{1} \\rangle_{du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1} \\rangle_{du} (0)$",
               xlim=c(0,40),
               ylim=c(0,0.01),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_6_id_du_no2ptflow[[i]]$cf0, dy=ratio_6_id_du_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio12_f1-aa-f0-id_qb_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
