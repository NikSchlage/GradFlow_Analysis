#--------------------------------------------------------------------------------------------------------#
#   LOAD PACKAGES AND CODE FOR ANALYSIS                                                                  #
#--------------------------------------------------------------------------------------------------------#
library(dplyr)                                                     ## required for filter()
devtools::load_all("/home/nikolas/Software/hadron")                ## required for new_matrixfit()
source("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/R_Code/kchi_and_c2pt_gf.R") ## load R functions required for gf-Zchi analysis
source("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/R_Code/mixing_probe_gf.R")


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
my_kinop_treelevel_path <- read.table( "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/08_plot_treelevel_kinetic_operator/kchi_treelevel_data.txt",
                                       header=TRUE )

my_Zchi_path <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/03_bootstrap_notimedilution_NLoopsamples8/12_boot_Zchi/bootstrap_Zchi_Nconfs350_Nsamples10.RData"
load(my_Zchi_path)

# res_ratio <- read.table("../AIC_fit_results_cf_ratio2_c2ptc2pt0_id_ud_du.txt")

# print(res_ratio[,2])
# print(res_ratio[,1])

Nconfs <- 1126 #length(my_confs)
Nsrcs <- 1
Nsamples <- 1


# Zchi-central-value: Zchi_MeanVal[[i]][-1]
# Zchi-error: Zchi_ErrorVal[[i]][-1]

ratio_msbar      <- c()
ratio_err1_msbar <- c()
ratio_err2_msbar <- c()

ratio_msbar_tgfsqr      <- c()
ratio_err1_msbar_tgfsqr <- c()
ratio_err2_msbar_tgfsqr <- c()

for( i in 2:length(my_gft) ){

  ratio_msbar      <- c( ratio_msbar, my_kinop_treelevel_path[i,2] )
  ratio_err1_msbar <- c( ratio_err1_msbar, my_kinop_treelevel_path[i,3] )
  ratio_err2_msbar <- c( ratio_err2_msbar, my_kinop_treelevel_path[i,4] )

  ratio_msbar_tgfsqr      <- c( ratio_msbar_tgfsqr, my_kinop_treelevel_path[i,2] * my_gft[i]^2 )
  ratio_err1_msbar_tgfsqr <- c( ratio_err1_msbar_tgfsqr, my_kinop_treelevel_path[i,3] * my_gft[i]^2 )
  ratio_err2_msbar_tgfsqr <- c( ratio_err2_msbar_tgfsqr, my_kinop_treelevel_path[i,4] * my_gft[i]^2 )

}

print(ratio_msbar)

tikzDevice::tikz( file = "kinop_treelevel.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=my_gft[-1], y=ratio_msbar, # dy=ratio_err1_msbar, mdy=ratio_err2_msbar,
               col="dodgerblue4",
               cex=1.2, pch=20,
               xlab="$t_\\mathrm{gf} / a^2$",
               ylab="$\\langle K(t_\\mathrm{gf}) \\rangle^{(0)}$")

#title(paste0("$\\mathrm{cA211a.30.32}$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "kinop_treelevel.tex", pdf=TRUE )




tikzDevice::tikz( file = "kinop_treelevel_tgfsqr.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=my_gft[-1]/32^2, y=ratio_msbar_tgfsqr, # dy=ratio_err1_msbar_tgfsqr, mdy=ratio_err2_msbar_tgfsqr,
               col="dodgerblue4",
               cex=1.2, pch=20,
               xlab="$t_\\mathrm{gf} / L^2$",
               ylab="$t_\\mathrm{gf}^2 \\cdot \\langle K(t_\\mathrm{gf}) \\rangle^{(0)}$")

#title(paste0("$\\mathrm{cA211a.30.32}$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "kinop_treelevel_tgfsqr.tex", pdf=TRUE )