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
#   LOAD WORKSPACE IMAGE                                                                                 #
#      * containing the unsymmetrized c2pt bootstrap samples   c2pt_mx_gamma5_cf_boot[[i]]               #
#      * containing the unsymmetrized c2pt bootstrap samples   c2pt_mx_id_cf_boot[[i]]                   #
#      * containing the unsymmetrized c2pt bootstrap samples   c2pt_mx_id_mulmin1_cf_boot[[i]]           #
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

##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##


#========================================================================================================#
#--------------------------------------------------------------------------------------------------------#
#   Plot the c2pt mixing_probe_src_oet_gf functions                                                      #
#--------------------------------------------------------------------------------------------------------#
#========================================================================================================#
tikzDevice::tikz( file = "c2pt_mx_gamma5_re.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=c2pt_mx_gamma5_cf_boot[[1]]$cf0, dy=c2pt_mx_gamma5_cf_boot[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$\\mathrm{Re}\\big[\\langle O_{\\gamma_5}(t_f)\\, O_{\\gamma_5}(t_i) \\rangle_u\\big]$",
               ylab="$\\langle O_{\\gamma_5}\\, O_{\\gamma_5} \\rangle (t_\\mathrm{gf},t)$",
               log="y",
               xlim=c(0,40),
               ylim=c(2,6000),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )

for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=c2pt_mx_gamma5_cf_boot[[i]]$cf0, dy=c2pt_mx_gamma5_cf_boot[[i]]$tsboot.se,
                 log="y",
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}

legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{conf}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=1$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "c2pt_mx_gamma5_re.tex", pdf=TRUE )

#--------------------------------------------------------------------------------------------------------#

tikzDevice::tikz( file = "c2pt_mx_gamma5_im.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=c2pt_mx_gamma5_cf_boot[[1]]$icf0, dy=c2pt_mx_gamma5_cf_boot[[1]]$itsboot.se,
               xlab="$t / a$",
               ylab="$\\mathrm{Im}\\big[\\langle O_{\\gamma_5}(t_f)\\, O_{\\gamma_5}(t_i) \\rangle_u\\big]$",
               xlim=c(0,40),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )

for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=c2pt_mx_gamma5_cf_boot[[i]]$icf0, dy=c2pt_mx_gamma5_cf_boot[[i]]$itsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}

legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{conf}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=1$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "c2pt_mx_gamma5_im.tex", pdf=TRUE )

#--------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------#

tikzDevice::tikz( file = "c2pt_mx_id_re.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=c2pt_mx_id_mulmin1_cf_boot[[1]]$cf0, dy=c2pt_mx_id_mulmin1_cf_boot[[1]]$tsboot.se,
               xlab="$t / a$",
               #ylab="$\\mathrm{Re}\\big[-\\langle O_\\mathbb{1}(t_f)\\, O_\\mathbb{1}(t_i) \\rangle_u\\big]$",
               ylab="$\\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle (t_\\mathrm{gf},t)$",
               log="y",
               xlim=c(0,40),
               ylim=c(1e-1,1e5),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )

for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=c2pt_mx_id_mulmin1_cf_boot[[i]]$cf0, dy=c2pt_mx_id_mulmin1_cf_boot[[i]]$tsboot.se,
                 log="y",
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}

legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{conf}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=1$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "c2pt_mx_id_re.tex", pdf=TRUE )

#--------------------------------------------------------------------------------------------------------#

tikzDevice::tikz( file = "c2pt_mx_id_im.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=seq(0,32,1), y=c2pt_mx_id_mulmin1_cf_boot[[1]]$icf0, dy=c2pt_mx_id_mulmin1_cf_boot[[1]]$itsboot.se,
               xlab="$t / a$",
               ylab="$\\mathrm{Im}\\big[-\\langle O_\\mathbb{1}(t_f)\\, O_\\mathbb{1}(t_i) \\rangle_u\\big]$",
               xlim=c(0,40),
               pch=19, cex=0.5, col=my_cols[1], cex.axis=1., cex.lab=1.2 )

for ( i in 2:length(my_gft) ){
  plotwitherror( x=seq(0,32,1), y=c2pt_mx_id_mulmin1_cf_boot[[i]]$icf0, dy=c2pt_mx_id_mulmin1_cf_boot[[i]]$itsboot.se,
                 pch=19, cex=0.5, col=my_cols[i], cex.axis=1., cex.lab=1.2, rep="TRUE" )
}

legend( "topright", legend=my_legend, col=my_cols, pch=19, cex=0.555 )
title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{conf}=", Nconfs, ",$ $N_\\mathrm{srcs}=", Nsrcs, ",$ $N_\\mathrm{samples}=1$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "c2pt_mx_id_im.tex", pdf=TRUE )


##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##