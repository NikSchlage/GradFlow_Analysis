#========================================================================================================#
#                                                                                                        #
#   R FUNCTIONS FOR EFFMASS AND C2PT ANALYSIS:                                                           #
#                                                                                                        #
#   * Create meson c2pt keys                                                                             #
#   * Read c2pt h5 files                                                                                 #
#   * Bootstrap c2pt data                                                                                #
#   * Bootstrap effective mass                                                                           #
#                                                                                                        #
#========================================================================================================#

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Create meson c2pt keys                                                                               #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
correlators_key_meson_2pt_gradflow <- function( bwd_flav, fwd_flav,
                                                snk_gamma, src_gamma,
                                                src_ts,
                                                id_sample,
                                                snk_p, src_p,
                                                gradflow_t )
{
  stopifnot( length(snk_p) == 3 )
  stopifnot( length(src_p) == 3 )
  stopifnot( is.character(fwd_flav) )
  stopifnot( is.character(bwd_flav) )
  stopifnot( is.integer(src_ts) )
  stopifnot( is.integer(id_sample) )
  stopifnot( is.integer(snk_gamma) )
  stopifnot( is.integer(src_gamma) )
  stopifnot( is.double(gradflow_t) )

  sprintf("/%s-gf-%s-gi/gf%d/gi%d/t%d/s%d/pix%dpiy%dpiz%d/gft%f/px%dpy%dpz%d",
          bwd_flav,
          fwd_flav,
          snk_gamma,
          src_gamma,
          src_ts,
          id_sample,
          src_p[1],src_p[2],src_p[3],
          gradflow_t,
          snk_p[1],snk_p[2],snk_p[3])
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Read c2pt h5 files                                                                                   #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
correlators_read_h5_meson_2pt_gradflow <- function ( path,
                                                     basename_start, basename_end,
                                                     T,
                                                     conflist,
                                                     bwd_flav, fwd_flav,
                                                     snk_gamma, src_gamma,
                                                     src_ts,
                                                     id_sample,
                                                     snk_p, src_p,
                                                     gradflow_t )
{
  x <- array( rep(NA, T*length(conflist)*length(gradflow_t)), dim=c(length(conflist), length(gradflow_t), T) )
  
  for( i in 1:length(conflist) ){
    h5file <- rhdf5::H5Fopen( paste0( path, "/", basename_start, str_pad(conflist[i], 4, pad = "0"), basename_end ),
                              flags = "H5F_ACC_RDONLY" )

    for( j in 1:length(gradflow_t) ) {

      h5key <- correlators_key_meson_2pt_gradflow( bwd_flav=bwd_flav, fwd_flav=fwd_flav,
                                                   snk_gamma=snk_gamma, src_gamma=src_gamma,
                                                   src_ts=src_ts,
                                                   id_sample=id_sample,
                                                   snk_p=snk_p, src_p=src_p,
                                                   gradflow_t=gradflow_t[j] )

      h5dataset <- hadron::h5_get_dataset( h5f=h5file, key=h5key )

    x[i,j,] <- complex( real = h5dataset[c(TRUE, FALSE)], imaginary = h5dataset[c(FALSE, TRUE)])
    }
    rhdf5::H5Fclose(h5file)
  }
  
  return(x)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Bootstrap c2pt data                                                                                  #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
boot_correlators_meson_2pt_gradflow <- function ( x, gradflow_t, boot.R=400, boot.l=2, seed=1234 )
{
  c2pt_cf_boot <- list()

  for ( i in 1:length(gradflow_t) ) {
    c2pt_cf <- corr_to_cf(x[,i,])
    c2pt_cf_boot[[i]] <- hadron::bootstrap.cf(c2pt_cf, boot.R = boot.R, boot.l = boot.l, seed = seed)
  }

  return(c2pt_cf_boot)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Bootstrap effective mass                                                                             #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
boot_effmass_meson_2pt_gradflow <- function ( x, gradflow_t )
{
  c2pt_effmass_boot <- list()

  for ( i in 1:length(gradflow_t) ) {
    cf_effmass_symm <- hadron::symmetrise.cf( cf=x[[i]] )
    c2pt_effmass_boot[[i]] <- bootstrap.effectivemass(cf=cf_effmass_symm, type="solve")
  }

  return(c2pt_effmass_boot)
}


#========================================================================================================#
#                                                                                                        #
#   R FUNCTIONS FOR ZCHI AND ZPCP ANALYSIS:                                                              #
#                                                                                                        #
#   * Read K_\chi data in ascii format                                                                   #
#   * Determine gradient-flow time values from K_\chi data in ascii format                               #
#   * Process gradient-flow K_\chi data: average over samples/ average over gauge configs/               #
#                                        sum over mu/ average over flavors                               #
#   * Compute ratio of "ringed" to MSbar quark-field renormalization constant                            #
#       ( \xi_\chi := \mathring{Z}_\chi / Z_\chi ) which is required for conversion to MSbar scheme      #
#   * Bootstrap data                                                                                     #
#   * Calculate Z_p * c_p using bootstrap samples                                                        #
#                                                                                                        #
#========================================================================================================#

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Read K_\chi data in ascii format                                                                     #
#                                                                                                        #
#   Takes                                                                                                #
#      * Path constituents                                        "path",                                #
#                                                                 "basename_start",                      #
#                                                                 "basename_end"                         #
#      * List of gauge configs for which the files are to be read "conflist"                             #
#        (for single gauge config this list has only one element)                                        #
#      * Total number of samples                                  "Nsamples"                             #
#                                                                                                        #
#                                                                                                        #
#   Returns                                                                                              #
#      * K_\chi data array                                        "x"                                    #
#                                                                                                        #
#   NOTE: This function is used to read in the gradient-flowed K_\chi data in ascii format               #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
read.ascii.gf.new <- function ( path=path, basename_start=basename_start, basename_end=basename_end,
                                conflist=conflist,
                                Nsamples=Nsamples ) {  
  x <- c()
  
  for( i in 1:length(conflist) ){
    x_tmp2 <- read.table( paste0( path, "/", basename_start, conflist[i], basename_end ),
                          header = FALSE, sep = " ", dec = ".",
                          colClasses=c("integer","character","numeric", "integer","numeric","numeric") )
    x_tmp1 <- x_tmp2[x_tmp2$V1 < Nsamples, ]
    x <- rbind( x, x_tmp1 )
  }
  
  return(x)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Determine gradient-flow time values from K_\chi data in ascii format                                 #
#                                                                                                        #
#   Takes                                                                                                #
#      * Path constituents                                        "path",                                #
#                                                                 "basename_start",                      #
#                                                                 "basename_end"                         #
#      * List of gauge configs for which the files are to be read "conflist"                             #
#        (for single gauge config this list has only one element)                                        #
#      * Total number of samples                                  "Nsamples"                             #
#                                                                                                        #
#                                                                                                        #
#   Returns                                                                                              #
#      * gardient-flow time vector                                "tvec"                                 #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
tvec.gf.new <- function ( path=path, basename_start=basename_start, basename_end=basename_end,
                          conflist=conflist,
                          Nsamples=Nsamples ) {
  muvec <- c(0,1,2,3)
  Nflavors <- 2
  
  x_aux_raw <- read.table( paste0( path, "/", basename_start, conflist[1], basename_end ),
                           header = FALSE, sep = " ", dec = "." )
  x_aux <- x_aux_raw[x_aux_raw$V1 < Nsamples,]
  
  length_tvec <- length(x_aux$V3) / ( Nsamples * Nflavors * length(muvec) )
  tvec <- array( x_aux$V3, dim=c(length(muvec), length_tvec) )[1,]

  return(tvec)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Create Zchi keys for reading h5 files                                                                #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#

zchi_key_gradflow <- function( isample,
                               flavor,
                               gft,
                               imu )
{

  my_key <- sprintf( "/s%d/%s/t%.6f/mu%d",
                     isample,
                     flavor,
                     gft,
                     imu )
  return(my_key)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Read K_\chi data in h5 format                                                                        #
#                                                                                                        #
#   Takes                                                                                                #
#      * Path constituents                                        "path",                                #
#                                                                 "basename_start",                      #
#                                                                 "basename_end"                         #
#      * List of gauge configs for which the files are to be read "conflist"                             #
#        (for single gauge config this list has only one element)                                        #
#      * Total number of samples                                  "Nsamples"                             #
#                                                                                                        #
#                                                                                                        #
#   Returns                                                                                              #
#      * K_\chi data array                                        "x"                                    #
#                                                                                                        #
#   NOTE: This function is used to read in the gradient-flowed K_\chi data in h5 format                  #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
read.h5.gf.new <- function ( path=path,
                             basename_start=basename_start,
                             basename_end=basename_end,
                             conflist=conflist,
                             tvec=tvec,
                             Nsamples=Nsamples ) {
  muvec <- c(0,1,2,3)
  Nflavors <- 2
  x <- c()

  for ( iconf in 1:length(conflist) ) {
    for ( isample in 1:Nsamples ) {
      for ( iflav in 1:Nflavors ) {
        if( (iflav-1) == 0 ) {  # 0 means flavor tag "up"
          flav <- "up"
        } else {                # 1 means flavor tag "down"
          flav <- "dn"
        }
        for ( igft in 1:length(tvec) ) {
          for ( imu in 1:length(muvec) ) {
            h5file <- rhdf5::H5Fopen( paste0( path, "/", basename_start, conflist[iconf], basename_end ), flags = "H5F_ACC_RDONLY" )
            h5key <- zchi_key_gradflow(isample-1, flav, tvec[igft], muvec[imu])
            h5dataset <- hadron::h5_get_dataset( h5f=h5file, key=h5key ) 
            x <- rbind(x, c(isample-1, flav, tvec[igft], imu-1, h5dataset))
          }
        }
      }
    }
  }

  return(x)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Process gradient-flow K_\chi data: average over samples/ average over gauge configs/                 #
#                                      sum over mu/ average over flavors                                 #
#                                                                                                        #
#   Takes                                                                                                #
#      * gradient-flow time vector                               "tvec"                                  #
#      * gradient flow data                                      "x"                                     #
#        (use function "read.ascii.gf.new()"                                                             #
#         or  function "read.h5.gf.new()"                                                                #
#         to get data in appropriate format)                                                             #
#      * Total number of samples                                 "Nsamples"                              #
#      * Complex-number index                                    "icplx"                                 #
#           ** Re: icplx = 1                                                                             #
#           ** Im: icplx = 2                                                                             #
#      * error function to be applied                            "err"                                   #
#        (must be sd, uwerr or left empty)                                                               #
#      * format of K_\chi data files                             "format"                                #
#        (must be either "ascii" or "h5")                                                                #
#                                                                                                        #
#   Returns                                                                                              #
#      * if Nconfs == 1:                                                                                 #
#           ** average over samples (sample mean)                "res_sample_mean"                       #
#           ** error after averaging over samples (sample error) "res_sample_sd"                         #
#              (only if err=="sd")                                                                       #
#           ** gradient-flow time vector                         "tvec"                                  #
#           ** average over samples/ sum over mu                 "res_sample_mean_mu_sum"                #
#           ** average over samples/ sum over mu/                                                        #
#              average over flavors                              "res_sample_mean_mu_sum_flavor_mean"    #
#      * if Nconfs > 1:                                                                                  #
#           ** average over samples/ average over gauge configs  "res_sample_mean_config_mean"           #
#           ** error after averaging over samples/ gauge configs "res_config_sd" or "res_config_uwerr"   #
#              (only if err is specified)                                                                #
#           ** gradient-flow time vector                         "tvec"                                  #
#           ** average over samples/ average over gauge configs/ "res_sample_mean_config_mean_mu_sum"    #
#              sum over mu                                                                               #
#           ** average over samples/ average over gauge configs/ "res_sample_mean_config_mean_           #
#              sum over mu/ average over flavors                  mu_sum_flavor_mean"                    #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#

analysis.gf.new <- function ( tvec=tvec, x=x, Nsamples=Nsamples, icplx=icplx, err="", format="ascii" )
{
  muvec <- c(0,1,2,3)
  Nflavors <- 2

  if ( format == "h5" ) {
    Nconfs <- as.numeric(sum(x[,1]==0)) / ( Nflavors*length(tvec)*length(muvec) )

    if ( icplx==1 ){
      X <- array( as.numeric(x[,5]), dim=c(length(muvec), length(tvec), Nflavors, Nsamples, Nconfs) )
    } else {
      X <- array( as.numeric(x[,6]), dim=c(length(muvec), length(tvec), Nflavors, Nsamples, Nconfs) )
    }

  } else if ( format == "ascii" ) {
    Nconfs <- as.numeric(sum(x$V1==0)) / ( Nflavors*length(tvec)*length(muvec) )

    xin <- x %>% filter( x$V1 < Nsamples )

    if ( icplx==1 ){
      X <- array( xin$V5, dim=c(length(muvec), length(tvec), Nflavors, Nsamples, Nconfs) )
    } else {
      X <- array( xin$V6, dim=c(length(muvec), length(tvec), Nflavors, Nsamples, Nconfs) )
    }

  } else {
    stop("Error: 'format' either has to be 'h5' or 'ascii'.")
  }
  
  #----------------------------------------------------#
  #   Averaging over samples                           #
  #   ( if Nconfs == 1 )                               #
  #----------------------------------------------------#
  res_sample_mean <- apply( X=X, MARGIN=c(1,2,3,5), FUN=mean )
  
  if ( Nconfs == 1 ) {
    res_sample_mean_mu_sum             <- apply( X=res_sample_mean,        MARGIN=c(2,3,4), FUN=sum  )
    res_sample_mean_mu_sum_flavor_mean <- apply( X=res_sample_mean_mu_sum, MARGIN=c(1,3),   FUN=mean )

    if ( err == "sd" ) {
      res_sample_sd <- apply( X=X, MARGIN=c(1,2,3,5), FUN=sd )
      res <- list( res_sample_mean,
                   res_sample_sd,
                   tvec,
                   res_sample_mean_mu_sum,
                   res_sample_mean_mu_sum_flavor_mean )    
    } else if ( err == "uwerr" ) {
      stop( "Error, uwerr cannot be calculated for single gauge config!" )
    } else {
      res <- list( res_sample_mean,
                   tvec,
                   res_sample_mean_mu_sum,
                   res_sample_mean_mu_sum_flavor_mean ) 
    }

  } else if ( Nconfs > 1 ) {
  #----------------------------------------------------#
  #   Additional averaging over gauge configurations   #
  #   ( if Nconfs > 1 )                                #
  #----------------------------------------------------#
    res_sample_mean_config_mean                    <- apply( X=res_sample_mean,                    MARGIN=c(1,2,3), FUN=mean )
    res_sample_mean_config_mean_mu_sum             <- apply( X=res_sample_mean_config_mean,        MARGIN=c(2,3),   FUN=sum  )
    res_sample_mean_config_mean_mu_sum_flavor_mean <- apply( X=res_sample_mean_config_mean_mu_sum, MARGIN=c(2),     FUN=mean )
  
    if ( err == "sd" ) {
      res_config_sd <- apply( X=res_sample_mean, MARGIN=c(1,2,3), FUN=sd )
      res <- list( res_sample_mean_config_mean,
                   res_config_sd,
                   tvec,
                   res_sample_mean_config_mean_mu_sum,
                   res_sample_mean_config_mean_mu_sum_flavor_mean )
    } else if ( err == "uwerr" ) { 
      res_config_uwerr <- apply( X=res_sample_mean, MARGIN=c(1,2,3), FUN=hadron::uwerrprimary )
      res <- list( res_sample_mean_config_mean,
                   res_config_uwerr,
                   tvec,
                   res_sample_mean_config_mean_mu_sum,
                   res_sample_mean_config_mean_mu_sum_flavor_mean )
    } else {
      res <- list( res_sample_mean_config_mean,
                   tvec,
                   res_sample_mean_config_mean_mu_sum,
                   res_sample_mean_config_mean_mu_sum_flavor_mean )
    }

  } else {
    stop( "Error, the number of gauge configurations must be Nconfs >= 1" )
  }
  
  return(res)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Compute ratio of "ringed" to MSbar quark-field renormalization constant                              #
#     ( \xi_\chi := \mathring{Z}_\chi / Z_\chi ) which is required for conversion to MSbar scheme        #
#                                                                                                        #
#   Takes                                                                                                #
#      * tgf    = gradient flow time provided as a vector                                                #
#      * mubar  = mass scale                                                                             #
#      * Nc     = number of quark colors                                                                 #
#      * Nf     = number of quark flavors which are active at the scale mubar                            #
#      * D      = dimension of Euclidean space-time                                                      #
#      * alphaS = QCD running coupling constant ( alphaS(mubar) = g^2(mubar^2)/(4*pi) )                  #
#                                                                                                        #
#   Returns                                                                                              #
#      * xichi  = ratio of "ringed" quark-field renormalization constant                                 #
#                          to MSbar quark-field renormalization constant                                 #
#      * Zg     = renormalization constant                                                               #
#      * g      = renormalized (physical) coupling constant at mass sacle mubar                          #
#      * g0     = bare (unphysical) coupling constant (independent of mubar)                             #
#      * tgf    = gradient flow time                                                                     #
#                                                                                                        #
#   as list elements                                                                                     #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#

 convert_to_msbar <- function ( tgf, mubar, Nc, Nf, D, alphaS )
{
  ## Euler-Mascheroni constant:
  gammaE    <- 0.57721

  ## trace normalization in the fundamental representation times number of quark flavors:
  TF        <- Nf / 2.

  ## MSbar scheme parameter introduced in [Harlander, Kluth, Lange (2019)]:
  L         <- log( 2 * mubar^2 * tgf ) + gammaE

  ## quadratic Casimir eigenvalue of the adjoint representation of the gauge group:
  CA        <- Nc

  ## quadratic Casimir eigenvalue of the fundamental representation of the gauge group:
  CF        <- ( Nc^2 - 1 ) / ( 2*Nc )

  ## three-loop constant introduced in [Harlander, Kluth, Lange (2019)]:
  C2        <- -23.8*CA*CF + 30.4*CF^2 - 3.92*CF*TF

  ## first coefficient of the beta-function expressed as a perturbative series:
  beta0     <- (11/3.)*CA - (4./3.)*TF

  ## second coefficient of the beta-function expressed as a perturbative series:
  beta1     <- (34./3.)*CA^2 - ( 4*CF + (20./3.)*CA )*TF

  ## inverse MSbar quark-field renormalization Z_\chi^{-1} coefficient [Harlander, Kluth, Lange (2019)]:
  gammachi0 <- 6. * CF

  ## inverse MSbar quark-field renormalization Z_\chi^{-1} coefficient [Harlander, Kluth, Lange (2019)]:
  gammachi1 <- CA*CF * ( (223./3.) - 16. * log(2) ) - CF^2 * ( 3. + 16. * log(2) ) - (44./3.)*CF*TF

  ## considering D-dimensional Euclidean space-time this parameter results from dimensional regularization with D=4-2*epsilon:
  epsilon <- ( 4 - D )/2

  ## renormalized (physical) coupling constant at mass sacle mubar:
  g <- sqrt( 4*pi*alphaS )

  ## renormalization constant:
  Zg <- 1 - ( g^2 / (4*pi)^2 ) * ( beta0 / (2*epsilon) ) + ( g^4 / (4*pi)^4 ) * ( (3/8)*(beta0^2 / epsilon^2) - beta1/(4*epsilon) )

  ## bare (unphysical) coupling constant (independent of mubar):
  g0 <- ( (mubar * exp( gammaE/2 ))/sqrt(4*pi) )^epsilon * Zg * g

  ## ratio of "ringed" quark-field renormalization constant to MSbar quark-field renormalization constant ( \xi_\chi := \mathring{Z}_\chi / Z_\chi ):
  xichi <- 1 + ( g^2 / (4*pi)^2 ) * ( (gammachi0 / 2.)*L - 3.*CF*log(3) - 4.*CF*log(2) ) + ( g^4 / (4*pi)^4 ) * ( (gammachi0 / 4.)*(beta0 + (gammachi0 / 2.))*L^2 + ( (gammachi1 / 2.) - (gammachi0 / 2.)*(beta0 + (gammachi0 / 2.) )*log(3) - (2./3.)*gammachi0*(beta0 + (gammachi0 / 2.))*log(2) )*L + C2 )

  res <- list( xichi, Zg, g0, tgf )
  return(res)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Bootstrap data                                                                                       #
#                                                                                                        #
#   ( can be used e.g. to bootstrap K_\chi data )                                                        #
#                                                                                                        #
#   Takes                                                                                                #
#      * data array of dimension dim=c(T, Nsamples), e.g. dim=c(length(gft) , Nconfs)                    #
#                                                                                                        #
#   Returns                                                                                              #
#      * Total number of bootstrap samples        "NBootSamples"                                         #
#      * Blocksize                                "Blocksize"                                            #
#      * Bootstrap seed                           "BootSeed"                                             #
#      * Mean of bootstrap samples                "BootMean"                                             #
#      * Error of bootstrap samples               "BootError"                                            #
#      * Error of the error of bootstrap samples  "BootDError"                                           #
#      * Integrated autocorrelation time          "BootTauint"                                           #
#      * Error on integrated autocorrelation time "BootDTauint"                                          #
#      * Bootstrap bias                           "BootBias"                                             #
#      * Bootstrap samples of data                "BootSamples"                                          #
#                                                                                                        #
#   NOTE: Complex number cannot be processed by this function so either pass real-part data              #
#         or imaginary-part data to this function                                                        #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
new.bootstrap.analysis <- function( data, skip=0, boot.R=400, boot.l=2,
                                    seed=1234, sim="geom", pl=FALSE )
{
  stopifnot(boot.l >= 1)
  stopifnot(boot.l <= nrow(data))
  stopifnot(boot.R >= 1)

  data <- data[skip:length(data)]
  data.mean = mean(data)
  error.naive = sd(data)/sqrt(length(data))

  NBootSamples <- numeric()
  Blocksize <-  numeric()
  BootSeed <- numeric()
  BootMean <- numeric()
  BootError <- numeric()
  BootDError <- numeric()
  BootTauint <- numeric()
  BootDTauint <- numeric()
  BootBias <- numeric()
  Blocksize[1] <- 1
  BootMean[1] <- 0.
  BootError[1] <- 0.
  BootDError[1] <- 0.
  BootTauint[1] <- 0.
  BootDTauint[1] <- 0.
  BootBias[1] <- 0.

  ndata <- block.ts(data, l=boot.l)

  ## now we bootstrap the data
  data.tsboot <- boot::boot(data=ndata, statistic=meanindexed, R=boot.R)

  ## use the same seed for data.sdboot as for data.tsboot
  set.seed(data.tsboot$seed)
  data.sdboot <- boot::boot(data=ndata, statistic=sd.index, R=boot.R)

  NBootSamples <- boot.R
  Blocksize <- boot.l
  BootSeed <- seed
  BootMean <- data.tsboot$t0[1]
  BootError <- sd(data.tsboot$t[,1])
  BootDError <- sd(data.sdboot$t[,1])/sqrt(length(ndata))
  BootTauint <- sd(data.tsboot$t[,1])^2/error.naive^2/2
  BootDTauint <- ( sd(data.sdboot$t[,1])^2/error.naive^2/2 )/sqrt(length(ndata))
  BootBias <- data.tsboot$t0[1] - mean(data.tsboot$t[,1])
  BootSamples <- data.tsboot$t[,1]

  res <- list( NBootSamples, Blocksize, BootSeed, BootMean, BootError, BootDError, BootTauint, BootDTauint, BootBias, BootSamples )

  return(res)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Calculate Z_p * c_p using bootstrap samples                                                          #
#                                                                                                        #
#   Takes                                                                                                #
#      * (symmetrized) c2pt bootstrap fit results "res_c2pt_cf_bootfit_gft" including the ampitude,      #
#         e.g. can be loaded via "load("../03_fit_and_plot_c2pt/c2pt_bootfit_Nconfs_120_gf.RData")"      #
#      * bootstrapped Z_chi data "res_zchi_bootfit_gft"                                                  #
#      * gardient flow time vector "gft"                                                                 #
#      * Number of colors "Nc"                                                                           #
#      * Number of flavors "Nf"                                                                          #
#                                                                                                        #
#   Returns                                                                                              #
#      * Total number of bootstrap samples "boot.R"                                                      #
#      * Blocksize                         "boot.l"                                                      #
#      * Bootstrap seed                    "seed"                                                        #
#      * Mean of bootstrap samples         "boot_Zpcp_gft_mean"                                          #
#      * Error of bootstrap samples        "boot_Zpcp_gft_err"                                           #
#      * Bootstrap samples of Z_p * c_p    "boot_Zpcp_gft"                                               #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#

boot_Zpcp_gradflow <- function ( res_c2pt_cf_bootfit_gft=res_c2pt_cf_bootfit_gft,
                                 res_zchi_bootfit_gft=res_zchi_bootfit_gft,
                                 gft=gft,
                                 Nc=Nc, Nf=Nf )
{
  boot_amplitude_gft <- c()

  for ( i in 1:length(my_gft) ) {
   boot_amplitude_gft <- rbind( boot_amplitude_gft, res_c2pt_cf_bootfit_gft[[i]]$opt.tsboot[2,] / res_c2pt_cf_bootfit_gft[[1]]$opt.tsboot[2,] )
  }

  boot_Zpcp_gft <- array( rep(NA, boot.R*length(gft)), dim=c(length(gft), boot.R) )

  for ( i in 1:length(gft) ) {
    for ( j in 1:boot.R ) {
      boot_Zpcp_gft[i,j] <- (res_zchi_bootfit_gft[[i]][[10]][j])^2 * boot_amplitude_gft[i,j]
    }
  }

  boot_Zpcp_gft_mean <- apply( X=boot_Zpcp_gft, MARGIN=1, FUN=mean)
  boot_Zpcp_gft_err <- apply( X=boot_Zpcp_gft, MARGIN=1, FUN=sd)

  res <- list( boot.R, boot.l, seed, boot_Zpcp_gft_mean, boot_Zpcp_gft_err, boot_Zpcp_gft ) 
  return(res)
}