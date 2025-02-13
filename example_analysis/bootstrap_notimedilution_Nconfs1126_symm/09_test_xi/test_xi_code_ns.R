#--------------------------------------------------------------------------------------------------------#
#   LOAD PACKAGES AND CODE FOR ANALYSIS                                                                  #
#--------------------------------------------------------------------------------------------------------#
source("/home/nikolas/Documents/NPV/njjn_analysis/npv_analysis_jan_2023/npv.R")

library(tidyverse)
library(dplyr)                                                     ## required for filter()
devtools::load_all("/home/nikolas/Software/hadron")                ## required for new_matrixfit()
#source("./xi_code_ns.R")
source("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/R_Code/kchi_and_c2pt_gf.R")

#----------------------------------------------------#


gftmax <- 0.01 # GeV^2
gftnum <- 1000
my_gft <- (1:gftnum)  * gftmax / gftnum

alpha  <- 1.77^2 / ( 4 * pi )
dalpha <- 0.05 * alpha

set.seed(12354) 
my_alphaS <- rnorm (n = 1000, mean = alpha, sd = dalpha )

my_mubar <- 3.
my_Nc <- 3
my_Nf <- 3


# Note: 1 GeV = 5.067731 fm^{-1} and fm^{-1} = 0.097 a^{-1}


res <- array( dim=c(gftnum, 3) )

fout <- "xi_t.txt"

cat( "# Nf       = ", my_Nf, "\n", file=fout, append=F )
cat( "# mubar    = ", my_mubar, "\n", file=fout, append=T )
cat( "# alpha_S  = ", alpha, " +/- ", dalpha, "\n", file=fout, append=T )

for ( i in 1:gftnum )
{
  my_gft_tmp <- my_gft[i]

  res_tmp <- convert_to_msbar( tgf=my_gft_tmp, mubar=my_mubar, Nc=my_Nc, Nf=my_Nf, D=4, alphaS=my_alphaS )

  res[i,] <- c( my_gft_tmp, mean(res_tmp[[1]]) , sqrt( var(res_tmp[[1]]) ) )

  cat ( formatC( res[i,], width=25, digits=16, format="e" ), "\n", file=fout, append=T )

}


#----------------------------------------------------#

my_legend <- c( "$N_c=3,\\,\\,\\,\\, N_f=3,\\,\\,\\,\\, \\bar{\\mu}=3\\,\\mathrm{GeV},\\,\\,\\,\\, \\alpha^{N_f=3}(\\bar{\\mu}=3\\,\\mathrm{GeV}) = 0.249(12)$" )
my_cols   <- c("red")

tikzDevice::tikz( file = "xi_t.tex", standAlone = TRUE, width = 4.5, height = 3)
par(mar = c(4.5, 4.5, 1.4, 0.25))

devtools::load_all("/home/nikolas/Software/hadron")   # required for plot.cf()
hadron::plotwitherror( x=res[,1],
                       y=res[,2],
                       col="red",
                       type="l",
                       lwd=2.0,
                       xlab="$t_\\mathrm{gf} / \\mathrm{GeV}^{-2}$",
                       ylab="$\\xi_\\chi (t_\\mathrm{gf}, \\bar{\\mu})$" )
plt.errband(x=res[,1], y=res[,2], se=res[,3], lwd = 1.5, col="red")

legend( "topright", legend=my_legend, col=my_cols, lwd=2., cex=0.8 )

dev.off()
# Compile the tex file
tools::texi2dvi( "xi_t.tex", pdf=TRUE )

#----------------------------------------------------#
