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
my_data_path <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/04_bootstrap_notimedilution_NLoopsamples8_symm/01_boot_cf_data/NLoopsamples_01/01_boot_cf_data_Nconfs563/bootstrap_cf_mx_prb_oet_gf_notimedilution_Nconfs563_Nsamples1.RData"
load(my_data_path)

my_Zchi_path <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/03_bootstrap_notimedilution_NLoopsamples8/12_boot_Zchi/bootstrap_Zchi_Nconfs350_Nsamples10.RData"
load(my_Zchi_path)

res_ratio <- read.table("../AIC_fit_results_cf_ratio_5_vvcc_id.txt")

print(res_ratio[,2])
print(res_ratio[,1])

Nconfs <- length(my_confs)
Nsrcs <- 1
Nsamples <- 1


# Zchi-central-value: Zchi_MeanVal[[i]][-1]
# Zchi-error: Zchi_ErrorVal[[i]][-1]

ratio_msbar      <- c()
ratio_err1_msbar <- c()
ratio_err2_msbar <- c()

for( i in 2:length(my_gft) ){

  xichi <- convert_to_msbar( tgf=res_ratio[i,1], mubar=3*5.067731*0.097, Nc=3, Nf=3, D=5, alphaS=0.249 )

  ratio_msbar      <- c( ratio_msbar, res_ratio[i,2] * Zchi_MeanVal[[1]][i] / xichi[[1]] )
  ratio_err1_msbar <- c( ratio_err1_msbar, res_ratio[i,4] * Zchi_MeanVal[[1]][i] / xichi[[1]] )
  ratio_err2_msbar <- c( ratio_err2_msbar, res_ratio[i,3] * Zchi_MeanVal[[1]][i] / xichi[[1]] )

  xichi <- c()

}

print(ratio_msbar)

tikzDevice::tikz( file = "AIC_fit_results_cf_ratio_5_vvcc_id_Nconfs563_Nsamples1_msbar.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=my_gft[-1], y=ratio_msbar, dy=ratio_err1_msbar, mdy=ratio_err2_msbar,
               col="dodgerblue4",
               cex=1.2,
               xlab="$t_\\mathrm{gf}$",
               ylab="$(Z_\\chi / \\xi_\\chi) \\cdot \\langle O_{vv}^{cc}\\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf})$", log="xy" )

title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{samples}=1,$ $N_\\mathrm{confs}=", Nconfs, "$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "AIC_fit_results_cf_ratio_5_vvcc_id_Nconfs563_Nsamples1_msbar.tex", pdf=TRUE )