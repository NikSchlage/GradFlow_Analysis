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
#      * containing the cf ratio ratio_1_g5 = Im(<X1g5>) / Re(<Xmx_g5>)                                  #
#      * containing the cf ratio ratio_2_g5 = Im(<X2g5>) / Re(<Xmx_g5>)                                  #
#      * containing the cf ratio ratio_3_g5 = Im(<X3g5mulmin1>) / Re(<Xmx_g5>)                           #
#                                                                                                        #
#      * containing the cf ratio ratio_4_id = Re(<X1idmulmin1>) / Re(<Xmx_id_mulmin1>)                   #
#      * containing the cf ratio ratio_5_id = Re(<X2idmulmin1>) / Re(<Xmx_id_mulmin1>)                   #
#      * containing the cf ratio ratio_6_id = Re(<X3id>) / Re(<Xmx_id_mulmin1>)                          #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
my_data_path <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/01_boot_cf_data/NLoopsamples_01/01_boot_cf_data_Nconfs1126/bootstrap_cf_mx_prb_oet_gf_notimedilution_Nconfs1126_Nsamples1.RData"
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

tikzDevice::tikz( file = "cf_ratio_4_vvqb_id_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(12,25,1), y=ratio_4_id[[1]]$cf0[13:26], dy=ratio_4_id[[1]]$tsboot.se[13:26],
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{qb}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) $",
               xlim=c(12,25),
               #ylim=c(-2,10),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(12,25,1), y=ratio_4_id[[i]]$cf0[13:26], dy=ratio_4_id[[i]]$tsboot.se[13:26],
                 xlab="$t / a$",
                 ylab="$\\langle O_{vv}^{qb}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) $",
                 xlim=c(12,25),
                 #ylim=c(-2,10),
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2) #, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_4_vvqb_id_Nsamples8.tex", pdf=TRUE )


stop()

boot <- FALSE

if ( boot == TRUE )
{
#--------------------------------------------------------------------------------------------------------#
#   Bootstrap the flavor-conserving cf ratios C^{3pt}(t_{gf}) / C^{2pt}(t_{gf}):                         #
#--------------------------------------------------------------------------------------------------------#
  ratio_1_g5 <- ratio_cf( cf2pt=Xmx_g5,
                          cf3pt=X1g5,
                          gamma_mix="g5",
                          T=my_T,
                          conflist=my_confs,
                          gradflow_t=my_gft,
                          flow_c2pt=TRUE )

  ratio_2_g5 <- ratio_cf( cf2pt=Xmx_g5,
                          cf3pt=X2g5,
                          gamma_mix="g5",
                          T=my_T,
                          conflist=my_confs,
                          gradflow_t=my_gft,
                          flow_c2pt=TRUE )

  ratio_3_g5 <- ratio_cf( cf2pt=Xmx_g5,
                          cf3pt=X3g5mulmin1,
                          gamma_mix="g5",
                          T=my_T,
                          conflist=my_confs,
                          gradflow_t=my_gft,
                          flow_c2pt=TRUE )

  ratio_4_id <- ratio_cf( cf2pt=Xmx_id_mulmin1,
                          cf3pt=X1idmulmin1,
                          gamma_mix="id",
                          T=my_T,
                          conflist=my_confs,
                          gradflow_t=my_gft,
                          flow_c2pt=TRUE )

  ratio_5_id <- ratio_cf( cf2pt=Xmx_id_mulmin1,
                          cf3pt=X2idmulmin1,
                          gamma_mix="id",
                          T=my_T,
                          conflist=my_confs,
                          gradflow_t=my_gft,
                          flow_c2pt=TRUE )

  ratio_6_id <- ratio_cf( cf2pt=Xmx_id_mulmin1,
                          cf3pt=X3id,
                          gamma_mix="id",
                          T=my_T,
                          conflist=my_confs,
                          gradflow_t=my_gft,
                          flow_c2pt=TRUE )

#--------------------------------------------------------------------------------------------------------#
#   Bootstrap the flavor-conserving cf ratios C^{3pt}(t_{gf}) / C^{2pt}(t_{gf}=0):                       #
#--------------------------------------------------------------------------------------------------------#
  ratio_1_g5_no2ptflow <- ratio_cf( cf2pt=Xmx_g5,
                                    cf3pt=X1g5,
                                    gamma_mix="g5",
                                    T=my_T,
                                    conflist=my_confs,
                                    gradflow_t=my_gft,
                                    flow_c2pt=FALSE )

  ratio_2_g5_no2ptflow <- ratio_cf( cf2pt=Xmx_g5,
                                    cf3pt=X2g5,
                                    gamma_mix="g5",
                                    T=my_T,
                                    conflist=my_confs,
                                    gradflow_t=my_gft,
                                    flow_c2pt=FALSE )

  ratio_3_g5_no2ptflow <- ratio_cf( cf2pt=Xmx_g5,
                                    cf3pt=X3g5mulmin1,
                                    gamma_mix="g5",
                                    T=my_T,
                                    conflist=my_confs,
                                    gradflow_t=my_gft,
                                    flow_c2pt=FALSE )

  ratio_4_id_no2ptflow <- ratio_cf( cf2pt=Xmx_id_mulmin1,
                                    cf3pt=X1idmulmin1,
                                    gamma_mix="id",
                                    T=my_T,
                                    conflist=my_confs,
                                    gradflow_t=my_gft,
                                    flow_c2pt=FALSE )

  ratio_5_id_no2ptflow <- ratio_cf( cf2pt=Xmx_id_mulmin1,
                                    cf3pt=X2idmulmin1,
                                    gamma_mix="id",
                                    T=my_T,
                                    conflist=my_confs,
                                    gradflow_t=my_gft,
                                    flow_c2pt=FALSE )

  ratio_6_id_no2ptflow <- ratio_cf( cf2pt=Xmx_id_mulmin1,
                                    cf3pt=X3id,
                                    gamma_mix="id",
                                    T=my_T,
                                    conflist=my_confs,
                                    gradflow_t=my_gft,
                                    flow_c2pt=FALSE )
}

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

tikzDevice::tikz( file = "cf_ratio_1_vvqb_g5_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_1_g5[[1]]$cf0, dy=ratio_1_g5[[1]]$tsboot.se,
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{qb}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) / \\langle O_{\\gamma_5}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) $",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_1_g5[[i]]$cf0, dy=ratio_1_g5[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_1_vvqb_g5_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_2_vvcc_g5_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_2_g5[[1]]$cf0, dy=ratio_2_g5[[1]]$tsboot.se,
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{cc}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) / \\langle O_{\\gamma_5}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) $",
               xlim=c(0,40),
               ylim=c(0,7),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_2_g5[[i]]$cf0, dy=ratio_2_g5[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_2_vvcc_g5_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_3_aaqb_g5_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_3_g5[[1]]$cf0, dy=ratio_3_g5[[1]]$tsboot.se,
               xlab="$t / a$",
               ylab="$\\langle O_{aa}^{qb}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) / \\langle O_{\\gamma_5}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) $",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_3_g5[[i]]$cf0, dy=ratio_3_g5[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_3_aaqb_g5_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_4_vvqb_id_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(12,26,1), y=ratio_4_id[[1]]$cf0[13:27], dy=ratio_4_id[[1]]$tsboot.se[13:27],
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{qb}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) $",
               xlim=c(12,26),
               ylim=c(0,6),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(12,26,1), y=ratio_4_id[[i]]$cf0[13:27], dy=ratio_4_id[[i]]$tsboot.se[13:27],
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2) #, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_4_vvqb_id_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_5_vvcc_id_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(18,22,1), y=ratio_5_id[[1]]$cf0[19:23], dy=ratio_5_id[[1]]$tsboot.se[19:23],
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{cc}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) $",
               xlim=c(18,24),
               ylim=c(0,15),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(18,22,1), y=ratio_5_id[[i]]$cf0[19:23], dy=ratio_5_id[[i]]$tsboot.se[19:23],
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2) #, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_5_vvcc_id_Nsamples8.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_6_aaqb_id_Nsamples8.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(18,22,1), y=ratio_6_id[[1]]$cf0[19:23], dy=ratio_6_id[[1]]$tsboot.se[19:23],
               xlab="$t / a$",
               ylab="$\\langle O_{aa}^{qb}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) $",
               xlim=c(18,24),
               ylim=c(0,6),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(18,22,1), y=ratio_6_id[[i]]$cf0[19:23], dy=ratio_6_id[[i]]$tsboot.se[19:23],
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2) #, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_6_aaqb_id_Nsamples8.tex", pdf=TRUE )




###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
#                                                                                                                                                 #
# PLOTS OF RATIOS WITHOUT!!! FLOWED 2PT-FUNCTION IN THE DENOMINATOR:                                                                              #
#                                                                                                                                                 #
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_1_vvqb_g5_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_1_g5_no2ptflow[[1]]$cf0, dy=ratio_1_g5_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{qb}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) / \\langle O_{\\gamma_5}\\, O_{\\gamma_5} \\rangle(0) $",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_1_g5_no2ptflow[[i]]$cf0, dy=ratio_1_g5_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_1_vvqb_g5_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_2_vvcc_g5_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_2_g5_no2ptflow[[1]]$cf0, dy=ratio_2_g5_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{cc}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) / \\langle O_{\\gamma_5}\\, O_{\\gamma_5} \\rangle(0) $",
               xlim=c(0,40),
               ylim=c(0,7),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_2_g5_no2ptflow[[i]]$cf0, dy=ratio_2_g5_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_2_vvcc_g5_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_3_aaqb_g5_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=ratio_3_g5_no2ptflow[[1]]$cf0, dy=ratio_3_g5_no2ptflow[[1]]$tsboot.se,
               xlab="$t / a$",
               ylab="$\\langle O_{aa}^{qb}\\, O_{\\gamma_5} \\rangle(t_\\mathrm{gf}) / \\langle O_{\\gamma_5}\\, O_{\\gamma_5} \\rangle(0) $",
               xlim=c(0,40),
               ylim=c(0,2.5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=ratio_3_g5_no2ptflow[[i]]$cf0, dy=ratio_3_g5_no2ptflow[[i]]$tsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_3_aaqb_g5_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_4_vvqb_id_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(18,22,1), y=ratio_4_id_no2ptflow[[1]]$cf0[19:23], dy=ratio_4_id_no2ptflow[[1]]$tsboot.se[19:23],
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{qb}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle(0) $",
               xlim=c(18,24),
               ylim=c(0,6),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(18,22,1), y=ratio_4_id_no2ptflow[[i]]$cf0[19:23], dy=ratio_4_id_no2ptflow[[i]]$tsboot.se[19:23],
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2) #, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_4_vvqb_id_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_5_vvcc_id_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(18,22,1), y=ratio_5_id_no2ptflow[[1]]$cf0[19:23], dy=ratio_5_id_no2ptflow[[1]]$tsboot.se[19:23],
               xlab="$t / a$",
               ylab="$\\langle O_{vv}^{cc}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle(0) $",
               xlim=c(18,24),
               ylim=c(0,15),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(18,22,1), y=ratio_5_id_no2ptflow[[i]]$cf0[19:23], dy=ratio_5_id_no2ptflow[[i]]$tsboot.se[19:23],
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2) #, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_5_vvcc_id_Nsamples8_no2ptflow.tex", pdf=TRUE )

###################################################################################################################################################

tikzDevice::tikz( file = "cf_ratio_6_aaqb_id_Nsamples8_no2ptflow.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(18,22,1), y=ratio_6_id_no2ptflow[[1]]$cf0[19:23], dy=ratio_6_id_no2ptflow[[1]]$tsboot.se[19:23],
               xlab="$t / a$",
               ylab="$\\langle O_{aa}^{qb}\\, O_\\mathbb{1} \\rangle(t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle(0) $",
               xlim=c(18,24),
               ylim=c(0,6),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )
for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(18,22,1), y=ratio_6_id_no2ptflow[[i]]$cf0[19:23], dy=ratio_6_id_no2ptflow[[i]]$tsboot.se[19:23],
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2) #, rep="TRUE" )
}
legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{confs}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=", Nsamples, "$"), line=0.6, cex=0.8)

dev.off()
# Compile the tex file
tools::texi2dvi( "cf_ratio_6_aaqb_id_Nsamples8_no2ptflow.tex", pdf=TRUE )
