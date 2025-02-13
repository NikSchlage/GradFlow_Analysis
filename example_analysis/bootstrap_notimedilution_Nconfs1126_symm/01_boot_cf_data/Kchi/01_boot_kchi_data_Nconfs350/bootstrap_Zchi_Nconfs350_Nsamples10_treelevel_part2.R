##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##

#--------------------------------------------------------------------------------------------------------#
#   LOAD PACKAGES AND CODE FOR ANALYSIS                                                                  #
#--------------------------------------------------------------------------------------------------------#
library(dplyr)                                                     ## required for filter()
devtools::load_all("/home/nikolas/Software/hadron")                ## required for new_matrixfit()
source("/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_ANALYSIS/R_Code/kchi_and_c2pt_gf.R") ## load R functions required for this
                                                                   ## gf-Zchi and gf-Zpcp analysis

## TEST:

test <- FALSE

if ( test == TRUE ) {
  mypath <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_DATA/04_qbig_kchi_data_Nsamples10"
  gradflow_t <- c(0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20)
  Nsamples <- 10

  mydata <- read.h5.gf.new( path=mypath, basename_start="test_gf.c", basename_end=".h5",
                            conflist=c(0),
                            tvec=gradflow_t,
                            Nsamples=Nsamples )

  myoutput <- analysis.gf.new( tvec=gradflow_t, x=mydata, Nsamples=Nsamples, icplx=1, err="", format="h5" )

  print(mydata)
  print(myoutput)
  stop()
}


#--------------------------------------------------------------------------------------------------------#
#   SET BOOTSTRAP PARAMETERS                                                                             #
#--------------------------------------------------------------------------------------------------------#
## bootstrapping parameters:
boot.R <- 400   ## number of bootstrap samples
boot.l <- 2     ## block length
seed <- 1234

##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##

#========================================================================================================#
#                                                                                                        #
#   BOOTSTRAP FOR ANALYSIS:                                                                              #
#                                                                                                        #
#   * Extract gradflow time from K_\chi data files in h5 format                                          #
#   * Read K_\chi data files in h5 format                                                                #
#   * Calculate Z_\chi = (-2*Nc*Nf) / K_\chi                                                             #
#   * Bootstrap Z_\chi data                                                                              #
#                                                                                                        #
#========================================================================================================#

Nconfigs <- 350
conflist <- seq(from=0, to=(Nconfigs-1)*4, by=4)
samplelist <- seq(from=6, to=8, by=1)

VOLUME <- 64 * 32^3
Nc <- 3   ## number of colors (r, g, b)
Nf <- 1   ## number of flavors 1 because u or d

mypath <- "/home/nikolas/Documents/NPV/njjn_analysis/GRADIENT_FLOW/mixing_probe_src_oet_GRADFLOW/GRADFLOW_DATA/04_qbig_kchi_data_Nsamples10"
my_gft <- c(0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20)

##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##


#--------------------------------------------------------------------------------------------------------#
#    Bootstrap Analysis for K_\chi, Z_\chi and (Z_p * Z_\chi)                                            #
#--------------------------------------------------------------------------------------------------------#
res_Kchi_bootfit_gft        <- list()
Kchi_MeanVal                <- list()
Kchi_ErrorVal               <- list()
Kchi_DErrorVal              <- list()
Kchi_TauintVal              <- list()
Kchi_DTauintVal             <- list()
res_gftsqr_Kchi_bootfit_gft <- list()
gftsqr_Kchi_MeanVal         <- list()
gftsqr_Kchi_ErrorVal        <- list()
gftsqr_Kchi_DErrorVal       <- list()
gftsqr_Kchi_TauintVal       <- list()
gftsqr_Kchi_DTauintVal      <- list()
res_Zchi_bootfit_gft        <- list()
Zchi_MeanVal                <- list()
Zchi_ErrorVal               <- list()
Zchi_DErrorVal              <- list()
Zchi_TauintVal              <- list()
Zchi_DTauintVal             <- list()
res_gftsqr_Zchi_bootfit_gft <- list()
gftsqr_Zchi_MeanVal         <- list()
gftsqr_Zchi_ErrorVal        <- list()
gftsqr_Zchi_DErrorVal       <- list()
gftsqr_Zchi_TauintVal       <- list()
gftsqr_Zchi_DTauintVal      <- list()
boot_Zpcp_res               <- list()


for ( isample in 1:length(samplelist) ) {
  print(paste0("Start bootstrap procedure for Nsamples = ", samplelist[isample]))

  #--------------------------------------------------------------------------------------------------------#
  #   Read (VOLUMEN * K_\chi) data files in h5 format                                                      #
  #--------------------------------------------------------------------------------------------------------#
  print("Read (VOLUMEN * Kchi) data files in h5 format")

  VOLUME_Kchi_data    <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  VOLUME_Kchi_data_up <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  VOLUME_Kchi_data_dn <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  Kchi                <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  gftsqr_Kchi         <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  Kchi_up             <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  Kchi_dn             <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  Zchi                <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  gftsqr_Zchi         <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  Zchi_up             <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  Zchi_dn             <- array( rep(NA, length(conflist)*length(my_gft)), dim=c(length(my_gft),length(conflist)) )
  
  for ( iconfig in 1:Nconfigs ) {
    gfdata <- read.h5.gf.new( path=mypath,
                              basename_start="test_gf.c",
                              basename_end=".h5",
                              conflist=conflist[iconfig],
                              tvec=my_gft,
                              Nsamples = samplelist[isample] )

    VOLUME_Kchi_data_re <- analysis.gf.new( tvec=my_gft,
                                            x=gfdata,
                                            Nsamples=samplelist[isample],
                                            icplx=1,
                                            err="",
                                            format="h5" )
    VOLUME_Kchi_data_im <- analysis.gf.new( tvec=my_gft,
                                            x=gfdata,
                                            Nsamples=samplelist[isample],
                                            icplx=2,
                                            err="",
                                            format="h5" )
  
    ## only consider K_\chi data "VOLUME_Kchi_data_re[[4]]" and "VOLUME_Kchi_data_im[[4]]"
    ## where we have averaged over
    ##    * samples,
    ##    * mu-values,
    ##    * flavors
    ## Thus, VOLUME_Kchi_data_re[[4]] contains one single vector of length "length(my_gft)"
    ## These results are stored as complex numbers in the array "VOLUME_Kchi_data"
    VOLUME_Kchi_data[,iconfig] <- complex( real = VOLUME_Kchi_data_re[[4]],
                                           imaginary = VOLUME_Kchi_data_im[[4]] )

    ## now also consider K_\chi data "VOLUME_Kchi_data_re[[3]]" and "VOLUME_Kchi_data_im[[3]]"
    ## where we have averaged over
    ##    * samples,
    ##    * mu-values.
    ## Thus, VOLUME_Kchi_data_re[[3]] contains two vectors of length "length(my_gft)", namely
    ##    * "VOLUME_Kchi_data_re[[3]][,1,1]" for flavor "up" and
    ##    * "VOLUME_Kchi_data_re[[3]][,2,1]" for flavor "down".
    ## These results are stored as complex numbers in the arrays
    ##    * "VOLUME_Kchi_data_up" and
    ##    * "VOLUME_Kchi_data_dn"
    VOLUME_Kchi_data_up[,iconfig] <- complex( real = VOLUME_Kchi_data_re[[3]][,1,1],
                                              imaginary = VOLUME_Kchi_data_im[[3]][,1,1] )

    VOLUME_Kchi_data_dn[,iconfig] <- complex( real = VOLUME_Kchi_data_re[[3]][,2,1],
                                              imaginary = VOLUME_Kchi_data_im[[3]][,2,1] )
  }

  #--------------------------------------------------------------------------------------------------------#
  #   Store                VOLUME x K_\chi = VOLUME * Re(K_\chi)                                           #
  #   Calculate            K_\chi = VOLUME x Re(K_\chi) / VOLUME                                           #
  #   Calculate t_{gf}^2 * K_\chi                                                                          #
  #   Calculate            Z_\chi =         ((-2) * Nc * Nf * VOLUME) / (VOLUME x K_\chi)                  #
  #   Calculate t_{gf}^2 * Z_\chi = gft^2 * ((-2) * Nc * Nf * VOLUME) / (VOLUME x K_\chi)                  #
  #                                                                                                        #
  #   Store                K_\chi_up = Re(K_\chi_up)                                                       #
  #   Calculate            Z_\chi_up = ((-2) * Nc * Nf * VOLUME) / (VOLUME x K_\chi_up)                    #
  #                                                                                                        #
  #   Store                K_\chi_dn = Re(K_\chi_dn)                                                       #
  #   Calculate            Z_\chi_dn = ((-2) * Nc * Nf * VOLUME) / (VOLUME x K_\chi_dn)                    #
  #--------------------------------------------------------------------------------------------------------#
  VOLUME_Kchi_data    <- Re(VOLUME_Kchi_data)
  VOLUME_Kchi_data_up <- Re(VOLUME_Kchi_data_up)
  VOLUME_Kchi_data_dn <- Re(VOLUME_Kchi_data_dn)

  for ( igft in 1:length(my_gft) ) {
    for ( iconfig in 1:Nconfigs ) {
      Kchi[igft,iconfig]        <- VOLUME_Kchi_data[igft,iconfig] / VOLUME
      gftsqr_Kchi[igft,iconfig] <- (my_gft[igft]^2 * VOLUME_Kchi_data[igft,iconfig]) / VOLUME
      Kchi_up[igft,iconfig]     <- VOLUME_Kchi_data_up[igft,iconfig] / VOLUME
      Kchi_dn[igft,iconfig]     <- VOLUME_Kchi_data_dn[igft,iconfig] / VOLUME

      Zchi[igft,iconfig]        <- ((-2) * Nc * Nf * VOLUME) / VOLUME_Kchi_data[igft,iconfig]
      gftsqr_Zchi[igft,iconfig] <- my_gft[igft]^2 * ((-2) * Nc * Nf * VOLUME) / VOLUME_Kchi_data[igft,iconfig]
      Zchi_up[igft,iconfig]     <- ((-2) * Nc * Nf * VOLUME) / VOLUME_Kchi_data_up[igft,iconfig]
      Zchi_dn[igft,iconfig]     <- ((-2) * Nc * Nf * VOLUME) / VOLUME_Kchi_data_dn[igft,iconfig]
    }
  }

  #--------------------------------------------------------------------------------------------------------#
  #   Bootstrap K_\chi data                                                                                #
  #--------------------------------------------------------------------------------------------------------#
  print("Bootstrap Kchi data")

  res_Kchi_bootfit_gft[[isample]] <- apply(X=Kchi, MARGIN=c(1), FUN=new.bootstrap.analysis, boot.R=400)

  Kchi_MeanVal_tmp    <- c()
  Kchi_ErrorVal_tmp   <- c()
  Kchi_DErrorVal_tmp  <- c()
  Kchi_TauintVal_tmp  <- c()
  Kchi_DTauintVal_tmp <- c()

  for ( igft in 1:length(my_gft) ){
    Kchi_MeanVal_tmp    <- c( Kchi_MeanVal_tmp,    as.numeric(res_Kchi_bootfit_gft[[isample]][[igft]][4]) )
    Kchi_ErrorVal_tmp   <- c( Kchi_ErrorVal_tmp,   as.numeric(res_Kchi_bootfit_gft[[isample]][[igft]][5]) )
    Kchi_DErrorVal_tmp  <- c( Kchi_DErrorVal_tmp,  as.numeric(res_Kchi_bootfit_gft[[isample]][[igft]][6]) )
    Kchi_TauintVal_tmp  <- c( Kchi_TauintVal_tmp,  as.numeric(res_Kchi_bootfit_gft[[isample]][[igft]][7]) )
    Kchi_DTauintVal_tmp <- c( Kchi_DTauintVal_tmp, as.numeric(res_Kchi_bootfit_gft[[isample]][[igft]][8]) )
  }

  Kchi_MeanVal[[isample]]    <- Kchi_MeanVal_tmp
  Kchi_ErrorVal[[isample]]   <- Kchi_ErrorVal_tmp
  Kchi_DErrorVal[[isample]]  <- Kchi_DErrorVal_tmp
  Kchi_TauintVal[[isample]]  <- Kchi_TauintVal_tmp
  Kchi_DTauintVal[[isample]] <- Kchi_DTauintVal_tmp

  #--------------------------------------------------------------------------------------------------------#
  #   Bootstrap t_{gf}^2 * K_\chi data                                                                     #
  #--------------------------------------------------------------------------------------------------------#
  print("Bootstrap tgf**2 * Kchi data")

  res_gftsqr_Kchi_bootfit_gft[[isample]] <- apply(X=gftsqr_Kchi, MARGIN=c(1), FUN=new.bootstrap.analysis, boot.R=400)

  gftsqr_Kchi_MeanVal_tmp    <- c()
  gftsqr_Kchi_ErrorVal_tmp   <- c()
  gftsqr_Kchi_DErrorVal_tmp  <- c()
  gftsqr_Kchi_TauintVal_tmp  <- c()
  gftsqr_Kchi_DTauintVal_tmp <- c()

  for ( igft in 1:length(my_gft) ){
    gftsqr_Kchi_MeanVal_tmp    <- c( gftsqr_Kchi_MeanVal_tmp,    as.numeric(res_gftsqr_Kchi_bootfit_gft[[isample]][[igft]][4]) )
    gftsqr_Kchi_ErrorVal_tmp   <- c( gftsqr_Kchi_ErrorVal_tmp,   as.numeric(res_gftsqr_Kchi_bootfit_gft[[isample]][[igft]][5]) )
    gftsqr_Kchi_DErrorVal_tmp  <- c( gftsqr_Kchi_DErrorVal_tmp,  as.numeric(res_gftsqr_Kchi_bootfit_gft[[isample]][[igft]][6]) )
    gftsqr_Kchi_TauintVal_tmp  <- c( gftsqr_Kchi_TauintVal_tmp,  as.numeric(res_gftsqr_Kchi_bootfit_gft[[isample]][[igft]][7]) )
    gftsqr_Kchi_DTauintVal_tmp <- c( gftsqr_Kchi_DTauintVal_tmp, as.numeric(res_gftsqr_Kchi_bootfit_gft[[isample]][[igft]][8]) )
  }

  gftsqr_Kchi_MeanVal[[isample]]    <- gftsqr_Kchi_MeanVal_tmp
  gftsqr_Kchi_ErrorVal[[isample]]   <- gftsqr_Kchi_ErrorVal_tmp
  gftsqr_Kchi_DErrorVal[[isample]]  <- gftsqr_Kchi_DErrorVal_tmp
  gftsqr_Kchi_TauintVal[[isample]]  <- gftsqr_Kchi_TauintVal_tmp
  gftsqr_Kchi_DTauintVal[[isample]] <- gftsqr_Kchi_DTauintVal_tmp

  #--------------------------------------------------------------------------------------------------------#
  #   Bootstrap Z_\chi data                                                                                #
  #--------------------------------------------------------------------------------------------------------#
  print("Bootstrap Zchi data")

  res_Zchi_bootfit_gft[[isample]] <- apply(X=Zchi, MARGIN=c(1), FUN=new.bootstrap.analysis, boot.R=400)

  Zchi_MeanVal_tmp    <- c()
  Zchi_ErrorVal_tmp   <- c()
  Zchi_DErrorVal_tmp  <- c()
  Zchi_TauintVal_tmp  <- c()
  Zchi_DTauintVal_tmp <- c()

  for ( igft in 1:length(my_gft) ){
    Zchi_MeanVal_tmp    <- c( Zchi_MeanVal_tmp,    as.numeric(res_Zchi_bootfit_gft[[isample]][[igft]][4]) )
    Zchi_ErrorVal_tmp   <- c( Zchi_ErrorVal_tmp,   as.numeric(res_Zchi_bootfit_gft[[isample]][[igft]][5]) )
    Zchi_DErrorVal_tmp  <- c( Zchi_DErrorVal_tmp,  as.numeric(res_Zchi_bootfit_gft[[isample]][[igft]][6]) )
    Zchi_TauintVal_tmp  <- c( Zchi_TauintVal_tmp,  as.numeric(res_Zchi_bootfit_gft[[isample]][[igft]][7]) )
    Zchi_DTauintVal_tmp <- c( Zchi_DTauintVal_tmp, as.numeric(res_Zchi_bootfit_gft[[isample]][[igft]][8]) )
  }

  Zchi_MeanVal[[isample]]    <- Zchi_MeanVal_tmp
  Zchi_ErrorVal[[isample]]   <- Zchi_ErrorVal_tmp
  Zchi_DErrorVal[[isample]]  <- Zchi_DErrorVal_tmp
  Zchi_TauintVal[[isample]]  <- Zchi_TauintVal_tmp
  Zchi_DTauintVal[[isample]] <- Zchi_DTauintVal_tmp

  #--------------------------------------------------------------------------------------------------------#
  #   Bootstrap t_{gf}^2 * Z_\chi data                                                                     #
  #--------------------------------------------------------------------------------------------------------#
  print("Bootstrap tgf**2 * Zchi data")

  res_gftsqr_Zchi_bootfit_gft[[isample]] <- apply(X=gftsqr_Zchi, MARGIN=c(1), FUN=new.bootstrap.analysis, boot.R=400)

  gftsqr_Zchi_MeanVal_tmp    <- c()
  gftsqr_Zchi_ErrorVal_tmp   <- c()
  gftsqr_Zchi_DErrorVal_tmp  <- c()
  gftsqr_Zchi_TauintVal_tmp  <- c()
  gftsqr_Zchi_DTauintVal_tmp <- c()

  for ( igft in 1:length(my_gft) ){
    gftsqr_Zchi_MeanVal_tmp    <- c( gftsqr_Zchi_MeanVal_tmp,    as.numeric(res_gftsqr_Zchi_bootfit_gft[[isample]][[igft]][4]) )
    gftsqr_Zchi_ErrorVal_tmp   <- c( gftsqr_Zchi_ErrorVal_tmp,   as.numeric(res_gftsqr_Zchi_bootfit_gft[[isample]][[igft]][5]) )
    gftsqr_Zchi_DErrorVal_tmp  <- c( gftsqr_Zchi_DErrorVal_tmp,  as.numeric(res_gftsqr_Zchi_bootfit_gft[[isample]][[igft]][6]) )
    gftsqr_Zchi_TauintVal_tmp  <- c( gftsqr_Zchi_TauintVal_tmp,  as.numeric(res_gftsqr_Zchi_bootfit_gft[[isample]][[igft]][7]) )
    gftsqr_Zchi_DTauintVal_tmp <- c( gftsqr_Zchi_DTauintVal_tmp, as.numeric(res_gftsqr_Zchi_bootfit_gft[[isample]][[igft]][8]) )
  }

  gftsqr_Zchi_MeanVal[[isample]]    <- gftsqr_Zchi_MeanVal_tmp
  gftsqr_Zchi_ErrorVal[[isample]]   <- gftsqr_Zchi_ErrorVal_tmp
  gftsqr_Zchi_DErrorVal[[isample]]  <- gftsqr_Zchi_DErrorVal_tmp
  gftsqr_Zchi_TauintVal[[isample]]  <- gftsqr_Zchi_TauintVal_tmp
  gftsqr_Zchi_DTauintVal[[isample]] <- gftsqr_Zchi_DTauintVal_tmp

  #--------------------------------------------------------------------------------------------------------#
  #   Compute Z_p * c_p bootstrap data                                                                     #
  #--------------------------------------------------------------------------------------------------------#
  #print("Compute Zp * cp bootstrap data")
  #
  #boot_Zpcp_res[[isample]] <- boot_Zpcp_gradflow( res_c2pt_cf_bootfit_gft=res_c2pt_cf_bootfit_gft,
  #                                                res_zchi_bootfit_gft=res_Zchi_bootfit_gft[[isample]],
  #                                                gft=my_gft, Nc=Nc, Nf=Nf )
}


##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
print("Start saving the workspace image")

#--------------------------------------------------------------------------------------------------------#
#   SAVE WORKSPACE IMAGE                                                                                 #
#--------------------------------------------------------------------------------------------------------#
save.image( file = paste0("bootstrap_Zchi_Nconfs", Nconfigs, "_Nsamples6to8.RData") )

##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##
##==================================================================================================================================================================================================================================##