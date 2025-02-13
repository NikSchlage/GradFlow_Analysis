library(dplyr)

source("/home/nikolas/Documents/NPV/njjn_analysis/npv_analysis_jan_2023/npv.R")
source("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/R_Code/mixing_probe_gf.R")


Nconfs <- 1126

my_tgf <- c( 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
             0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20 )

res1 <- read.table("./AIC_fit_results_cf_ratio_4_vvqb_id_ud_du.txt")
res2 <- get_best_redchisqr_pval_result_new( iratio=4, tu1=2, tu2=5, my_tgf=my_tgf, print_txt=TRUE )
res3 <- get_best_redchisqr_pval_result_new( iratio=4, tu1=2, tu2=6, my_tgf=my_tgf, rm_last_row_for_tu=6, print_txt=TRUE )

print(res2[,2])
print(res2[,1])


tikzDevice::tikz( file = "AIC_fit_results_cf_ratio_4_vvqb_id_ud_du_Nconfs1126_Nsamples1.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=res1[,1], y=res1[,2], dy=res1[,4], mdy=res1[,3],
               col="dodgerblue4",
               cex=0.5, pch=19,
               xlab="$t_\\mathrm{gf} / a^2$",
               ylab="$\\langle O_{vv}^{qb}\\, O_\\mathbb{1}\\rangle_{ud+du} (t_\\mathrm{gf}) / \\langle O_\\mathbb{1} \\, O_\\mathbb{1}\\rangle_{ud+du} (t_\\mathrm{gf})$" )

plotwitherror( x=res2[,1]+0.003, y=res2[,4], dy=res2[,5], mdy=res2[,5],
               col="red3",
               cex=0.5, pch=19, rep=TRUE )

plotwitherror( x=res3[,1]+0.006, y=res3[,4], dy=res3[,5], mdy=res3[,5],
               col="forestgreen",
               cex=0.5, pch=19, rep=TRUE )

title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{samples}=1,$ $N_\\mathrm{confs}=", Nconfs, "$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "AIC_fit_results_cf_ratio_4_vvqb_id_ud_du_Nconfs1126_Nsamples1.tex", pdf=TRUE )