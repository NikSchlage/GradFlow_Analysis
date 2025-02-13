##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##

#--------------------------------------------------------------------------------------------------------#
#   LOAD PACKAGES AND CODE FOR ANALYSIS                                                                  #
#--------------------------------------------------------------------------------------------------------#
my_data_path <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/05_bootstrap_notimedilution_Nconfs1126_symm/01_boot_cf_data/Kchi/02_boot_kchi_data_treelevel/bootstrap_Kchi_Nconfs1_Nsamples10_treelevel.RData"
load(my_data_path)

mp_data_path <- read.table( "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_DATA/06_qbig_kchi_data_Nsamples10_treelevel/idgauge_kinetic_op_gf_treelevel_mp/kinop.u.L032.mu0.003000.tl",
                            header=FALSE )


print(Re(Kchi_treelevel_data))

my_legend_1 <- c( "$\\text{Niko:}\\,\\,\\,\\,\\,\\,\\,\\,\\,\\,\\, L = 32,\\,\\,\\,\\, a \\mu = 0.003,\\,\\,\\,\\ \\kappa = 0.125,\\,\\,\\,\\ N_\\mathrm{samples}=10$",
                  "$\\text{Marcus:}\\,\\,\\,\\, L = 32,\\,\\,\\,\\, a \\mu = 0.003,\\,\\,\\,\\ \\kappa = 0.125,\\,\\,\\,\\ N_\\mathrm{samples}=1$" )

my_legend_2 <- c( "$\\text{Niko:}\\,\\,\\,\\, L = 32,\\,\\,\\,\\, a \\mu = 0.003,\\,\\,\\,\\ \\kappa = 0.125,\\,\\,\\,\\ N_\\mathrm{samples}=10$" )

kchi_treelevel_data <- cbind( my_gft, Re(Kchi_treelevel_data) )
colnames(kchi_treelevel_data) <- c( "t_{gf}", "kchi-treelevel" )

write.table( x=kchi_treelevel_data, file="kchi_treelevel_data.txt", append = FALSE, sep = " ", dec = ".",
             row.names = FALSE, col.names = TRUE )

tikzDevice::tikz( file = "kchi_treelevel.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{slashed, amsmath}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

hadron::plotwitherror( x=my_gft,
                       y=Re(Kchi_treelevel_data),
                       xlab="$t_\\mathrm{gf} / a^2$",
                       ylab="$-\\langle \\bar{\\chi} (t_\\mathrm{gf},x)\\, \\frac{1}{2} \\overleftrightarrow{\\slashed{D}}\\, \\chi (t_\\mathrm{gf},x) \\rangle_\\text{vol-avg}^{(0)}$",
                       pch=19, cex=0.6, cex.axis=1., cex.lab=1.1, col="red" )

hadron::plotwitherror( x=my_gft+0.003,
                       y=mp_data_path[,2],
                       rep=TRUE,
                       pch=19, cex=0.6, cex.axis=1., cex.lab=1.1, col="blue" )

legend( "topright", legend=my_legend_1, col=c("red","blue"), lwd=2., cex=0.8 )

dev.off()
# Compile the tex file
tools::texi2dvi( "kchi_treelevel.tex", pdf=TRUE )




Kchi_treelevel_data_tgfsqr <- c()

for( i in 1:length(my_gft) ){
  Kchi_treelevel_data_tgfsqr <- c( Kchi_treelevel_data_tgfsqr, Kchi_treelevel_data[i] * my_gft[i]^2 )
}

tikzDevice::tikz( file = "kchi_treelevel_tgfsqr.tex", standAlone = TRUE, width = 4.5, height = 3,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{slashed, amsmath}") )
par(mar = c(4.5, 4.5, 1.4, 0.25))

plot( x=my_gft/32^2,
      y=Re(Kchi_treelevel_data_tgfsqr),
      xlab="$t_\\mathrm{gf} / L^2$",
      ylab="$-t_\\mathrm{gf}^2 \\cdot \\langle \\bar{\\chi} (t_\\mathrm{gf},x)\\, \\frac{1}{2} \\overleftrightarrow{\\slashed{D}}\\, \\chi (t_\\mathrm{gf},x) \\rangle_\\text{vol-avg}^{(0)}$",
      xlim=c(0., 0.000205),
      ylim=c(0., 0.012),
      pch=19, cex=0.6, cex.axis=1., cex.lab=1.1, col="red", xaxt="n" )

axis(1,at=c(0.00000, 0.00005, 0.00010, 0.00015, 0.00020),labels=format(c(0.00000, 0.00005, 0.00010, 0.00015, 0.00020),scientific=TRUE))

legend( "bottomright", legend=my_legend_2, col="red", lwd=2., cex=0.8 )

dev.off()
# Compile the tex file
tools::texi2dvi( "kchi_treelevel_tgfsqr.tex", pdf=TRUE )