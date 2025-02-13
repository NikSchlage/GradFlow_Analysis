
Nconfs <- 563

res <- read.table("./AIC_fit_results_cf_ratio_5_vvcc_id.txt")

print(res[,2])
print(res[,1])


tikzDevice::tikz( file = "AIC_fit_results_cf_ratio_5_vvcc_id_Nconfs563_Nsamples1.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
plotwitherror( x=res[,1], y=res[,2], dy=res[,4], mdy=res[,3],
               col="dodgerblue4",
               cex=1.2,
               log="xy",
               xlab="$t_\\mathrm{gf}$",
               ylab="$\\langle O_{vv}^{cc}\\, O_\\mathbb{1}\\rangle (t_\\mathrm{gf}) / \\langle O_\\mathbb{1}\\, O_\\mathbb{1} \\rangle (t_\\mathrm{gf})$" )

title(paste0("$\\mathrm{cA211a.30.32,}$ $N_\\mathrm{samples}=1,$ $N_\\mathrm{confs}=", Nconfs, "$"), line=0.6)

dev.off()
# Compile the tex file
tools::texi2dvi( "AIC_fit_results_cf_ratio_5_vvcc_id_Nconfs563_Nsamples1.tex", pdf=TRUE )