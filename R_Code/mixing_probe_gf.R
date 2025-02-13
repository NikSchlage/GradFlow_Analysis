
#========================================================================================================#
#                                                                                                        #
#   R FUNCTIONS FOR MIXING PROBE SRC OET ANALYSIS:                                                       #
#                                                                                                        #
#   * Shift a vector forwards ("fwd") or backwards ("bwd")                                               #
#   * Create keys for 2-point mixing functions <O_X O_X>, X = \gamma_5, id                               #
#   * Create keys for 2-point mixing functions <O_X1 O_X2>, X1, X2 = \gamma_5, id                        #
#   * Create keys for 2-point mixing functions <O_X O_X>, X = \gamma_5, id                               #
#       for h5 files containing Nsamples data                                                            #
#   * Create keys for 2-point mixing functions <O_X1 O_X2>, X1, X2 = \gamma_5, id                        #
#     for h5 files containing Nsamples data                                                              #
#   * Create keys for 3-point mixing functions with mixing operator at src <O_{4q} O_X>, X=\gamma_5, id  #
#   * Create keys for 3-point mixing functions with mixing operator at src <O_{4q} O_X>, X=\gamma_5, id  #
#       for h5 files containing Nsamples data                                                            #
#   * Read flowed <O_X O_X> data from mixing_probe_src_oet_gf files in h5 format                         #
#   * Read flowed <O_X1 O_X2> data from mixing_probe_src_oet_gf files in h5 format                       #
#   * Read flowed <O_{4q} O_X> data from mixing_probe_src_oet_gf_c3pt files in h5 format                 #
#   * Read flowed <O_{4q} O_X> data from mixing_probe_src_oet_gf_c3pt files in h5 format                 #
#     for h5 data files containing on the one hand keys without sample IDs and on the other hand         #
#     keys with sample IDs                                                                               #
#   * Compute mean value, uwerr, integrated autocorr time for given gradflow-corr data                   #
#       preprocessed with "c2pt_read_h5_mx_oet_gradflow()" or "c3pt_read_h5_mx_prb_src_oet_gradflow()"   #
#   * Bootstrap mixing_probe_src_oet_gf cf data                                                          #
#   * Calculate (cf_1 + cf_2)/2                                                                          #
#   * Calculate and bootstrap (cf_1 + cf_2)/2                                                            #
#   * Calculate cf ratio R(t) = C3pt(t) / C2pt(t)                                                        #
#                                                                                                        #
#   * Perform a Linear Fit (needed for AIC weighted plateau Fit)                                         #
#   * Calculate the weighted sum of normal distributions for uniroot (needed for AIC weighted            #
#       plateau Fit)                                                                                     #
#   * Calculate the Akaike Information Criterion (AIC) weight                                            #
#   * Perform a AIC weighted plateau Fit                                                                 #
#                                                                                                        #
#   * Extract the best-redchisqr-pval result from given plateau-fit-result data                          #
#   * Extract the best-redchisqr-pval result from given plateau-fit-result data (version 2)              #
#                                                                                                        #
#========================================================================================================#

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Shift a vector forwards ("fwd") or backwards ("bwd")                                                 #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
shift_vec <- function(x, n, dir="fwd"){
  stopifnot( length(x) >= n )
  if( n == 0 ){
    return(x)
  }

  if( dir == "fwd" ){
    return(c(x[(length(x)-n+1):length(x)], x[1:(length(x)-n)]))
  } else if( dir == "bwd" ){
    return(c(x[(n+1):length(x)], x[1:n]))
  } else {
    print("ERROR: The direction parameter 'dir' either must be set to 'fwd' or to 'bwd'")
  }
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Create keys for 2-point mixing functions <O_X O_X>, X = \gamma_5, id                                 #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c2pt_key_mx_prb_src_oet_gradflow <- function( diagram_type,
                                              h5i_group,
                                              gradflow_t,
                                              flavor_ind_1,
                                              gamma_mix,
                                              flavor_ind_2,
                                              px, py, pz )
{

  sprintf( "/%s/tau%6.4f/f%d-%s-f%d-%s/%s/px%dpy%dpz%d",
           h5i_group,
           gradflow_t,
           flavor_ind_1,
           gamma_mix,
           flavor_ind_2,
           gamma_mix,
           diagram_type,
           px, py, pz )
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Create keys for 2-point mixing functions <O_X1 O_X2>, X1, X2 = \gamma_5, id                          #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c2pt_key_mx_prb_src_oet_gradflow_new <- function( diagram_type,
                                                  h5i_group,
                                                  gradflow_t,
                                                  flavor_ind_1,
                                                  gamma_mix_1,
                                                  flavor_ind_2,
                                                  gamma_mix_2,
                                                  px, py, pz )
{

  sprintf( "/%s/tau%6.4f/f%d-%s-f%d-%s/%s/px%dpy%dpz%d",
           h5i_group,
           gradflow_t,
           flavor_ind_1,
           gamma_mix_1,
           flavor_ind_2,
           gamma_mix_2,
           diagram_type,
           px, py, pz )
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Create keys for 2-point mixing functions <O_X O_X>, X = \gamma_5, id                                 #
#   for h5 files containing Nsamples data                                                                #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c2pt_key_mx_prb_src_oet_gradflow_nsamples <- function( diagram_type,
                                                       h5i_group,
                                                       gradflow_t,
                                                       isample,
                                                       flavor_ind_1,
                                                       gamma_mix,
                                                       flavor_ind_2,
                                                       px, py, pz )
{

  sprintf( "/%s/tau%6.4f/ns%d/f%d-%s-f%d-%s/%s/px%dpy%dpz%d",
           h5i_group,
           gradflow_t,
           isample,
           flavor_ind_1,
           gamma_mix,
           flavor_ind_2,
           gamma_mix,
           diagram_type,
           px, py, pz )
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Create keys for 2-point mixing functions <O_X1 O_X2>, X1, X2 = \gamma_5, id                          #
#   for h5 files containing Nsamples data                                                                #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c2pt_key_mx_prb_src_oet_gradflow_nsamples_new <- function( diagram_type,
                                                           h5i_group,
                                                           gradflow_t,
                                                           isample,
                                                           flavor_ind_1,
                                                           gamma_mix_1,
                                                           flavor_ind_2,
                                                           gamma_mix_2,
                                                           px, py, pz )
{

  sprintf( "/%s/tau%6.4f/ns%d/f%d-%s-f%d-%s/%s/px%dpy%dpz%d",
           h5i_group,
           gradflow_t,
           isample,
           flavor_ind_1,
           gamma_mix_1,
           flavor_ind_2,
           gamma_mix_2,
           diagram_type,
           px, py, pz )
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Create keys for 3-point mixing functions with mixing operator at src <O_{4q} O_X>, X=\gamma_5, id    #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c3pt_key_mx_prb_src_oet_gradflow <- function( diagram_type,
                                              h5i_group,
                                              gradflow_t,
                                              flavor_ind_1,
                                              gamma_tag_4q,
                                              flavor_ind_2,
                                              gamma_mix,
                                              px, py, pz )
{

  stopifnot( is.character(diagram_type) ) # "b" or "d"
  stopifnot( is.character(h5i_group) )    # "qb" or "cc"
  stopifnot( is.double(gradflow_t) )
  stopifnot( is.integer(flavor_ind_1) )   # 0="u" or 1="d"
  stopifnot( is.character(gamma_tag_4q) ) # "vv" or "aa"
  stopifnot( is.integer(flavor_ind_2) )   # 0="u" or 1="d"
  stopifnot( is.character(gamma_mix) )    # "id" or "g5"

  sprintf( "/%s/tau%6.4f/f%d-%s-f%d-%s/%s/px%dpy%dpz%d",
           h5i_group,
           gradflow_t,
           flavor_ind_1,
           gamma_tag_4q,
           flavor_ind_2,
           gamma_mix,
           diagram_type,
           px, py, pz )
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Create keys for 3-point mixing functions with mixing operator at src <O_{4q} O_X>, X=\gamma_5, id    #
#   for h5 files containing Nsamples data                                                                #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c3pt_key_mx_prb_src_oet_gradflow_nsamples <- function( diagram_type,
                                                       h5i_group,
                                                       gradflow_t,
                                                       isample,
                                                       flavor_ind_1,
                                                       gamma_tag_4q,
                                                       flavor_ind_2,
                                                       gamma_mix,
                                                       px, py, pz )
{

  stopifnot( is.character(diagram_type) ) # "b" or "d"
  stopifnot( is.character(h5i_group) )    # "qb" or "cc"
  stopifnot( is.double(gradflow_t) )
  stopifnot( is.integer(isample) )
  stopifnot( is.integer(flavor_ind_1) )   # 0="u" or 1="d"
  stopifnot( is.character(gamma_tag_4q) ) # "vv" or "aa"
  stopifnot( is.integer(flavor_ind_2) )   # 0="u" or 1="d"
  stopifnot( is.character(gamma_mix) )    # "id" or "g5"

  sprintf( "/%s/tau%6.4f/ns%d/f%d-%s-f%d-%s/%s/px%dpy%dpz%d",
           h5i_group,
           gradflow_t,
           isample,
           flavor_ind_1,
           gamma_tag_4q,
           flavor_ind_2,
           gamma_mix,
           diagram_type,
           px, py, pz )
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Read flowed <O_X O_X> data from mixing_probe_src_oet_gf files in h5 format                           #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c2pt_read_h5_mx_oet_gradflow <- function ( path,
                                           outfile_prefix,
                                           T,
                                           h5i_group,
                                           conflist,
                                           src_ts_list,
                                           isample_list,
                                           flavor_1,
                                           gamma_mix,
                                           flavor_2,
                                           px, py, pz,
                                           gradflow_t,
                                           re_im_x=FALSE,
                                           average=FALSE )
{
  if ( flavor_1 == "up" ) {
    flavor_ind_1 <- as.integer(0)
  } else if ( flavor_1 == "down" ) {
    flavor_ind_1 <- as.integer(1)
  } else {
    print("Error: 'flavor_1' either must be 'up' or 'down'")
  }

  if ( flavor_2 == "up" ) {
    flavor_ind_2 <- as.integer(0)
  } else if ( flavor_2 == "down" ) {
    flavor_ind_2 <- as.integer(1)
  } else {
    print("Error: 'flavor_2' either must be 'up' or 'down'")
  }

  if( is.null(isample_list) == TRUE ){
    if ( average == TRUE ) {
      stop("Error: h5 file contains only one sample. Averaging over 'length(isample_list)' samples is therefore not possible!")
    }

    x <- array( rep(NA, T*length(conflist)*length(gradflow_t)),
                dim=c(length(conflist), length(gradflow_t), T) )
  
    for( i in 1:length(conflist) ){
      h5file <- rhdf5::H5Fopen( paste0( path, "/", outfile_prefix, ".", conflist[i], ".t", src_ts_list[i], ".s0.h5" ),
                                flags = "H5F_ACC_RDONLY" )

      for( j in 1:length(gradflow_t) ) {

        h5key <- c2pt_key_mx_prb_src_oet_gradflow( diagram_type="m",
                                                   h5i_group=h5i_group,
                                                   gradflow_t=gradflow_t[j],
                                                   flavor_ind_1=flavor_ind_1,
                                                   gamma_mix=gamma_mix,
                                                   flavor_ind_2=flavor_ind_2,
                                                   px=px,
                                                   py=py,
                                                   pz=pz )

        h5dataset <- hadron::h5_get_dataset( h5f=h5file, key=h5key )

        if ( re_im_x == FALSE ) {
          x[i,j,] <- complex( real = h5dataset[c(TRUE, FALSE)], imaginary = h5dataset[c(FALSE, TRUE)] )
        } else {
          x[i,j,] <- complex( real = h5dataset[c(FALSE, TRUE)], imaginary = h5dataset[c(TRUE, FALSE)] )
        }
        x[i,j,] <- shift_vec( x=x[i,j,], n=src_ts_list[i], dir="bwd" )
      }
      rhdf5::H5Fclose(h5file)
    }
  } else {
    isample_list <- as.integer(isample_list)

    x <- array( rep(NA, T*length(conflist)*length(isample_list)*length(gradflow_t)),
                dim=c(length(conflist), length(isample_list), length(gradflow_t), T) )

    for( i in 1:length(conflist) ){
      for( j in 1:length(isample_list) ){
        h5file <- rhdf5::H5Fopen( paste0( path, "/", outfile_prefix, ".", conflist[i], ".t", src_ts_list[i], ".s0.h5" ),
                                  flags = "H5F_ACC_RDONLY" )

        for( k in 1:length(gradflow_t) ) {

          h5key <- c2pt_key_mx_prb_src_oet_gradflow_nsamples( diagram_type="m",
                                                              h5i_group=h5i_group,
                                                              gradflow_t=gradflow_t[k],
                                                              isample=isample_list[j],
                                                              flavor_ind_1=flavor_ind_1,
                                                              gamma_mix=gamma_mix,
                                                              flavor_ind_2=flavor_ind_2,
                                                              px=px,
                                                              py=py,
                                                              pz=pz )

          h5dataset <- hadron::h5_get_dataset( h5f=h5file, key=h5key )

          if ( re_im_x == FALSE ) {
            x[i,j,k,] <- complex( real = h5dataset[c(TRUE, FALSE)], imaginary = h5dataset[c(FALSE, TRUE)] )
          } else {
            x[i,j,k,] <- complex( real = h5dataset[c(FALSE, TRUE)], imaginary = h5dataset[c(TRUE, FALSE)] )
          }
          x[i,j,k,] <- shift_vec( x=x[i,j,k,], n=src_ts_list[i], dir="bwd" )
        }
        rhdf5::H5Fclose(h5file)
      }
    }

    if ( average==TRUE && length(isample_list) > 1 ) {
      print("Averaging over samples. For raw data use 'average=FALSE'.")
      x <- apply( X=x, MARGIN=c(1,3,4), FUN=mean )
    } else if ( average==TRUE && length(isample_list) == 1 ) {
      stop("Error: Only one sample from h5 file is selected. Averaging over 'length(isample_list)' samples is therefore not possible!")
    } else if ( average==FALSE && length(isample_list) > 1 ) {
      print("No averaging over samples. Use 'average=TRUE' for averaging.")
    } else {
      print("Note: Only one sample has been selected.")
      x <- apply( X=x, MARGIN=c(1,3,4), FUN=mean )
    }
  }
  
  return(x)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Read flowed <O_X1 O_X2> data from mixing_probe_src_oet_gf files in h5 format                         #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c2pt_read_h5_mx_oet_gradflow_new <- function ( path,
                                               outfile_prefix,
                                               T,
                                               h5i_group,
                                               conflist,
                                               src_ts_list,
                                               isample_list,
                                               flavor_1,
                                               gamma_mix_1,
                                               flavor_2,
                                               gamma_mix_2,
                                               px, py, pz,
                                               gradflow_t,
                                               re_im_x=FALSE,
                                               average=FALSE )
{
  if ( flavor_1 == "up" ) {
    flavor_ind_1 <- as.integer(0)
  } else if ( flavor_1 == "down" ) {
    flavor_ind_1 <- as.integer(1)
  } else {
    print("Error: 'flavor_1' either must be 'up' or 'down'")
  }

  if ( flavor_2 == "up" ) {
    flavor_ind_2 <- as.integer(0)
  } else if ( flavor_2 == "down" ) {
    flavor_ind_2 <- as.integer(1)
  } else {
    print("Error: 'flavor_2' either must be 'up' or 'down'")
  }

  if( is.null(isample_list) == TRUE ){
    if ( average == TRUE ) {
      stop("Error: h5 file contains only one sample. Averaging over 'length(isample_list)' samples is therefore not possible!")
    }

    x <- array( rep(NA, T*length(conflist)*length(gradflow_t)),
                dim=c(length(conflist), length(gradflow_t), T) )
  
    for( i in 1:length(conflist) ){
      h5file <- rhdf5::H5Fopen( paste0( path, "/", outfile_prefix, ".", conflist[i], ".t", src_ts_list[i], ".s0.h5" ),
                                flags = "H5F_ACC_RDONLY" )

      for( j in 1:length(gradflow_t) ) {

        h5key <- c2pt_key_mx_prb_src_oet_gradflow_new( diagram_type="m",
                                                       h5i_group=h5i_group,
                                                       gradflow_t=gradflow_t[j],
                                                       flavor_ind_1=flavor_ind_1,
                                                       gamma_mix_1=gamma_mix_1,
                                                       flavor_ind_2=flavor_ind_2,
                                                       gamma_mix_2=gamma_mix_2,
                                                       px=px,
                                                       py=py,
                                                       pz=pz )

        h5dataset <- hadron::h5_get_dataset( h5f=h5file, key=h5key )

        if ( re_im_x == FALSE ) {
          x[i,j,] <- complex( real = h5dataset[c(TRUE, FALSE)], imaginary = h5dataset[c(FALSE, TRUE)] )
        } else {
          x[i,j,] <- complex( real = h5dataset[c(FALSE, TRUE)], imaginary = h5dataset[c(TRUE, FALSE)] )
        }
        x[i,j,] <- shift_vec( x=x[i,j,], n=src_ts_list[i], dir="bwd" )
      }
      rhdf5::H5Fclose(h5file)
    }
  } else {
    isample_list <- as.integer(isample_list)

    x <- array( rep(NA, T*length(conflist)*length(isample_list)*length(gradflow_t)),
                dim=c(length(conflist), length(isample_list), length(gradflow_t), T) )

    for( i in 1:length(conflist) ){
      for( j in 1:length(isample_list) ){
        h5file <- rhdf5::H5Fopen( paste0( path, "/", outfile_prefix, ".", conflist[i], ".t", src_ts_list[i], ".s0.h5" ),
                                  flags = "H5F_ACC_RDONLY" )

        for( k in 1:length(gradflow_t) ) {

          h5key <- c2pt_key_mx_prb_src_oet_gradflow_nsamples_new( diagram_type="m",
                                                                  h5i_group=h5i_group,
                                                                  gradflow_t=gradflow_t[k],
                                                                  isample=isample_list[j],
                                                                  flavor_ind_1=flavor_ind_1,
                                                                  gamma_mix_1=gamma_mix_1,
                                                                  flavor_ind_2=flavor_ind_2,
                                                                  gamma_mix_2=gamma_mix_2,
                                                                  px=px,
                                                                  py=py,
                                                                  pz=pz )

          h5dataset <- hadron::h5_get_dataset( h5f=h5file, key=h5key )

          if ( re_im_x == FALSE ) {
            x[i,j,k,] <- complex( real = h5dataset[c(TRUE, FALSE)], imaginary = h5dataset[c(FALSE, TRUE)] )
          } else {
            x[i,j,k,] <- complex( real = h5dataset[c(FALSE, TRUE)], imaginary = h5dataset[c(TRUE, FALSE)] )
          }
          x[i,j,k,] <- shift_vec( x=x[i,j,k,], n=src_ts_list[i], dir="bwd" )
        }
        rhdf5::H5Fclose(h5file)
      }
    }

    if ( average==TRUE && length(isample_list) > 1 ) {
      print("Averaging over samples. For raw data use 'average=FALSE'.")
      x <- apply( X=x, MARGIN=c(1,3,4), FUN=mean )
    } else if ( average==TRUE && length(isample_list) == 1 ) {
      stop("Error: Only one sample from h5 file is selected. Averaging over 'length(isample_list)' samples is therefore not possible!")
    } else if ( average==FALSE && length(isample_list) > 1 ) {
      print("No averaging over samples. Use 'average=TRUE' for averaging.")
    } else {
      print("Note: Only one sample has been selected.")
      x <- apply( X=x, MARGIN=c(1,3,4), FUN=mean )
    }
  }
  
  return(x)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Read flowed <O_{4q} O_X> data from mixing_probe_src_oet_gf_c3pt files in h5 format                   #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c3pt_read_h5_mx_prb_src_oet_gradflow <- function ( path,
                                                   outfile_prefix,
                                                   T,
                                                   diagram_type,
                                                   h5i_group,
                                                   conflist,
                                                   src_ts_list,
                                                   isample_list=NULL,
                                                   flavor_1,
                                                   gamma_tag_4q,
                                                   flavor_2,
                                                   gamma_mix,
                                                   px, py, pz,
                                                   gradflow_t,
                                                   re_im_x=FALSE,
                                                   average=FALSE,
                                                   print_details=FALSE )
{
  if ( flavor_1 == "up" ) {
    flavor_ind_1 <- as.integer(0)
  } else if ( flavor_1 == "down" ) {
    flavor_ind_1 <- as.integer(1)
  } else {
    stop("Error: 'flavor_1' either must be 'up' or 'down'")
  }

  if ( flavor_2 == "up" ) {
    flavor_ind_2 <- as.integer(0)
  } else if ( flavor_2 == "down" ) {
    flavor_ind_2 <- as.integer(1)
  } else {
    stop("Error: 'flavor_2' either must be 'up' or 'down'")
  }
  
  if( is.null(isample_list) == TRUE ){
    if ( average == TRUE ) {
      stop("Error: h5 file contains only one sample. Averaging over 'length(isample_list)' samples is therefore not possible!")
    }

    x <- array( rep(NA, T*length(conflist)*length(gradflow_t)),
                dim=c(length(conflist), length(gradflow_t), T) )

    for( i in 1:length(conflist) ){
      h5file <- rhdf5::H5Fopen( paste0( path, "/", outfile_prefix, ".", conflist[i], ".t", src_ts_list[i], ".s0.h5" ),
                                flags = "H5F_ACC_RDONLY" )

      if ( print_details == TRUE ) {
        print( paste0("conf = ", conflist[i], " iconf =", i ) )
      }

      for( j in 1:length(gradflow_t) ) {

        h5key <- c3pt_key_mx_prb_src_oet_gradflow( diagram_type=diagram_type,
                                                   h5i_group=h5i_group,
                                                   gradflow_t=gradflow_t[j],
                                                   flavor_ind_1=flavor_ind_1,
                                                   gamma_tag_4q=gamma_tag_4q,
                                                   flavor_ind_2=flavor_ind_2,
                                                   gamma_mix=gamma_mix,
                                                   px=px,
                                                   py=py,
                                                   pz=pz )

        if ( print_details == TRUE ) {
          print( h5key )
        }

        h5dataset <- hadron::h5_get_dataset( h5f=h5file, key=h5key )

        if ( re_im_x == FALSE ) {
          x[i,j,] <- complex( real = h5dataset[c(TRUE, FALSE)], imaginary = h5dataset[c(FALSE, TRUE)] )
        } else {
          x[i,j,] <- complex( real = h5dataset[c(FALSE, TRUE)], imaginary = h5dataset[c(TRUE, FALSE)] )
        }
        x[i,j,] <- shift_vec( x=x[i,j,], n=src_ts_list[i], dir="bwd" )
      }
      rhdf5::H5Fclose(h5file)
    }
  } else {
    isample_list <- as.integer(isample_list)

    x <- array( rep(NA, T*length(conflist)*length(isample_list)*length(gradflow_t)),
                dim=c(length(conflist), length(isample_list), length(gradflow_t), T) )

    for( i in 1:length(conflist) ){
      for( j in 1:length(isample_list) ){
        h5file <- rhdf5::H5Fopen( paste0( path, "/", outfile_prefix, ".", conflist[i], ".t", src_ts_list[i], ".s0.h5" ),
                                  flags = "H5F_ACC_RDONLY" )

        if ( print_details == TRUE ) {
          print( paste0("conf = ", conflist[i], " iconf =", i ) )
        }

        for( k in 1:length(gradflow_t) ) {

          h5key <- c3pt_key_mx_prb_src_oet_gradflow_nsamples( diagram_type=diagram_type,
                                                              h5i_group=h5i_group,
                                                              gradflow_t=gradflow_t[k],
                                                              isample=isample_list[j],
                                                              flavor_ind_1=flavor_ind_1,
                                                              gamma_tag_4q=gamma_tag_4q,
                                                              flavor_ind_2=flavor_ind_2,
                                                              gamma_mix=gamma_mix,
                                                              px=px,
                                                              py=py,
                                                              pz=pz )

          if ( print_details == TRUE ) {
            print( h5key )
          }

          h5dataset <- hadron::h5_get_dataset( h5f=h5file, key=h5key )

          if ( re_im_x == FALSE ) {
            x[i,j,k,] <- complex( real = h5dataset[c(TRUE, FALSE)], imaginary = h5dataset[c(FALSE, TRUE)] )
          } else {
            x[i,j,k,] <- complex( real = h5dataset[c(FALSE, TRUE)], imaginary = h5dataset[c(TRUE, FALSE)] )
          }
          x[i,j,k,] <- shift_vec( x=x[i,j,k,], n=src_ts_list[i], dir="bwd" )
        }
        rhdf5::H5Fclose(h5file)
      }
    }

    if ( average==TRUE && length(isample_list) > 1 ) {
      print("Averaging over samples. For raw data use 'average=FALSE'.")
      x <- apply( X=x, MARGIN=c(1,3,4), FUN=mean )
    } else if ( average==TRUE && length(isample_list) == 1 ) {
      stop("Error: Only one sample from h5 file is selected. Averaging over 'length(isample_list)' samples is therefore not possible!")
    } else if ( average==FALSE && length(isample_list) > 1 ) {
      print("No averaging over samples. Use 'average=TRUE' for averaging.")
    } else {
      print("Note: Only one sample has been selected.")
      # Note: apply( ..., FUN=mean) is only used to reduce the dimensions of the single-sample array x and
      # y <- x
      # x <- apply( X=x, MARGIN=c(1,3,4), FUN=mean )
      # print(all(x==y[,1,,]))
      #  would print "TRUE"
      x <- apply( X=x, MARGIN=c(1,3,4), FUN=mean )
    }
  }

  return(x)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Read flowed <O_{4q} O_X> data from mixing_probe_src_oet_gf_c3pt files in h5 format                   #
#   for h5 data files containing on the one hand keys without sample IDs and on the other hand           #
#   keys with sample IDs                                                                                 #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
c3pt_read_h5_mx_prb_src_oet_gradflow_onesample_mixedkeys <- function ( path,
                                                                       outfile_prefix,
                                                                       T,
                                                                       diagram_type,
                                                                       h5i_group,
                                                                       conflist,
                                                                       src_ts_list,
                                                                       iconflist_key_with_isample,
                                                                       flavor_1,
                                                                       gamma_tag_4q,
                                                                       flavor_2,
                                                                       gamma_mix,
                                                                       px, py, pz,
                                                                       gradflow_t,
                                                                       re_im_x=FALSE,
                                                                       average=FALSE )
{
  x1 <- c3pt_read_h5_mx_prb_src_oet_gradflow( path=path,
                                              outfile_prefix=outfile_prefix,
                                              T=T,
                                              diagram_type=diagram_type,
                                              h5i_group=h5i_group,
                                              conflist=conflist[-iconflist_key_with_isample],
                                              src_ts_list=src_ts_list[-iconflist_key_with_isample],
                                              isample_list=c(NULL),
                                              flavor_1=flavor_1,
                                              gamma_tag_4q=gamma_tag_4q,
                                              flavor_2=flavor_2,
                                              gamma_mix=gamma_mix,
                                              px=px, py=py, pz=pz,
                                              gradflow_t=gradflow_t,
                                              re_im_x=re_im_x,
                                              average=FALSE )

  x2 <- c3pt_read_h5_mx_prb_src_oet_gradflow( path=path,
                                              outfile_prefix=outfile_prefix,
                                              T=T,
                                              diagram_type=diagram_type,
                                              h5i_group=h5i_group,
                                              conflist=conflist[iconflist_key_with_isample],
                                              src_ts_list=src_ts_list[iconflist_key_with_isample],
                                              isample_list=c(0),
                                              flavor_1=flavor_1,
                                              gamma_tag_4q=gamma_tag_4q,
                                              flavor_2=flavor_2,
                                              gamma_mix=gamma_mix,
                                              px=px, py=py, pz=pz,
                                              gradflow_t=gradflow_t,
                                              re_im_x=re_im_x,
                                              average=FALSE )

  x <- abind(x1, x2, along = 1)

  return(x)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Compute mean value, uwerr, integrated autocorr time for given gradflow-corr data                     #
#   preprocessed with "c2pt_read_h5_mx_oet_gradflow()" or "c3pt_read_h5_mx_prb_src_oet_gradflow()"       #
#   (gradflow-corr data for "Nsamples" loop samples and "Nconfs" gauge configurations)                   #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
analysis.uwerr.corr.gf <- function( x, Nsamples1=TRUE, return_uwerr_obj=FALSE ) {
  if( Nsamples1 == TRUE ) {
    # In this case: dim(x) = c(Nconfs, Ngradflowtpoints, T)
    res_sample_mean <- x
  } else {
    # In this case: dim(x) = c(Nconfs, Nsamples, Ngradflowtpoints, T)
    res_sample_mean <- apply( X=x, MARGIN=c(1,3,4), FUN=mean )
  }

  res_config_uwerr <- apply( X=res_sample_mean, MARGIN=c(2,3), FUN=hadron::uwerrprimary )

  value   <- array( rep(NA, dim(res_config_uwerr)[1]*dim(res_config_uwerr)[2]), dim=c(dim(res_config_uwerr)[1], dim(res_config_uwerr)[2]) )
  dvalue  <- array( rep(NA, dim(res_config_uwerr)[1]*dim(res_config_uwerr)[2]), dim=c(dim(res_config_uwerr)[1], dim(res_config_uwerr)[2]) )
  ddvalue <- array( rep(NA, dim(res_config_uwerr)[1]*dim(res_config_uwerr)[2]), dim=c(dim(res_config_uwerr)[1], dim(res_config_uwerr)[2]) )
  tauint  <- array( rep(NA, dim(res_config_uwerr)[1]*dim(res_config_uwerr)[2]), dim=c(dim(res_config_uwerr)[1], dim(res_config_uwerr)[2]) )
  dtauint <- array( rep(NA, dim(res_config_uwerr)[1]*dim(res_config_uwerr)[2]), dim=c(dim(res_config_uwerr)[1], dim(res_config_uwerr)[2]) )

  for( i in 1:dim(res_config_uwerr)[1] ) {
    for ( j in 1:dim(res_config_uwerr)[2] ) {
      value[i,j]   <- res_config_uwerr[i,j][[1]]$value
      dvalue[i,j]  <- res_config_uwerr[i,j][[1]]$dvalue
      ddvalue[i,j] <- res_config_uwerr[i,j][[1]]$ddvalue
      tauint[i,j]  <- res_config_uwerr[i,j][[1]]$tauint
      dtauint[i,j] <- res_config_uwerr[i,j][[1]]$dtauint
    }
  }

  if ( return_uwerr_obj == FALSE ) {
    res <- list(value, dvalue, ddvalue, tauint, dtauint)
    names(res) <- c("value", "dvalue", "ddvalue", "tauint", "dtauint")
  } else {
    res <- res_config_uwerr
  }

  return(res)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Bootstrap mixing_probe_src_oet_gf cf data                                                            #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
boot_cf_mixing_probe_src_oet_gf <- function ( x,
                                              gradflow_t,
                                              boot.R=400,
                                              boot.l=2,
                                              seed=1234 )
{
  my_cf_boot <- list()

  for ( i in 1:length(gradflow_t) ) {
    my_cf <- corr_to_cf(x[,i,])
    my_cf_boot[[i]] <- hadron::bootstrap.cf( my_cf,
                                             boot.R = boot.R,
                                             boot.l = boot.l,
                                             seed = seed )
  }

  return(my_cf_boot)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Calculate (cf_1 + cf_2)/2                                                                            #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
corr2cf_mixing_probe_src_oet_gf_average_cf2pts <- function ( x1,
                                                             x2,
                                                             Nconfs,
                                                             T,
                                                             gradflow_t )
{
  my_cf_avg <- array( rep(NA,Nconfs*length(gradflow_t)*T),dim=c(Nconfs,length(gradflow_t),T) )

  for ( i in 1:Nconfs ) {
    for ( j in 1:T ) {
      for ( k in 1:length(gradflow_t) ) {
        my_cf_avg[i,k,j] <- (x1[i,k,j] + x2[i,k,j]) / 2
      }
    }
  }

  return(my_cf_avg)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Calculate and bootstrap (cf_1 + cf_2)/2                                                              #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
boot_cf_mixing_probe_src_oet_gf_average_cf2pts <- function ( x1,
                                                             x2,
                                                             gradflow_t,
                                                             boot.R=400,
                                                             boot.l=2,
                                                             seed=1234 )
{
  my_cf_boot <- list()

  for ( i in 1:length(gradflow_t) ) {
    my_cf1 <- corr_to_cf(x1[,i,])
    my_cf2 <- corr_to_cf(x2[,i,])
    my_cf_avg <- mul.cf((my_cf1 + my_cf2),0.5)
    my_cf_boot[[i]] <- hadron::bootstrap.cf( my_cf_avg,
                                             boot.R = boot.R,
                                             boot.l = boot.l,
                                             seed = seed )
  }

  return(my_cf_boot)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Calculate cf ratio R(t) = C3pt(t) / C2pt(t)                                                          #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
# returns the correlation function ratio
#
# R(t) = C3pt(t) / C2pt(t)
#
# in form of a cf object.

ratio_cf <- function( cf2pt, cf3pt, gamma_mix, T, conflist, gradflow_t, flow_c2pt=TRUE, boot.R=400, boot.l=2, seed=1234, symm=FALSE ) {

  res_boot_c2pt <- list()
  res_boot_c3pt <- list()
  res_boot <- list()

  x2pt <- array( rep(NA, T*length(conflist)*length(gradflow_t)), dim=c(length(conflist), length(gradflow_t), T) )
  x3pt <- array( rep(NA, T*length(conflist)*length(gradflow_t)), dim=c(length(conflist), length(gradflow_t), T) )
  ratio <- array( rep(NA, T*length(conflist)*length(gradflow_t)), dim=c(length(conflist), length(gradflow_t), T) )

  if ( gamma_mix == "g5" ) {
    for( i in 1:length(conflist) ){
      for( j in 1:length(gradflow_t) ) {
          x2pt[i,j,] <- complex( real = Re(cf2pt[i,j,]), imaginary = Im(cf2pt[i,j,]))
          x3pt[i,j,] <- complex( real = Im(cf3pt[i,j,]), imaginary = Re(cf3pt[i,j,]))
      }
    }
  } else if ( gamma_mix == "id" ) {
    for( i in 1:length(conflist) ){
      for( j in 1:length(gradflow_t) ) {
          x2pt[i,j,] <- complex( real = Re(cf2pt[i,j,]), imaginary = Im(cf2pt[i,j,]))
          x3pt[i,j,] <- complex( real = Re(cf3pt[i,j,]), imaginary = Im(cf3pt[i,j,]))
      }
    }
  } else {
    print("Error: gamma_mix either nust be 'g5' or 'id'")
  }

  for ( i in 1:length(gradflow_t) ) {

    if ( flow_c2pt == TRUE ) {
      cf2pt_cf <- corr_to_cf(x2pt[,i,])
    } else { 
      cf2pt_cf <- corr_to_cf(x2pt[,1,])
    }
    cf3pt_cf <- corr_to_cf(x3pt[,i,])

    res_boot_c2pt[[i]] <- hadron::bootstrap.cf( cf2pt_cf,
                                                boot.R = boot.R,
                                                boot.l = boot.l,
                                                seed = seed )

    res_boot_c3pt[[i]] <- hadron::bootstrap.cf( cf3pt_cf,
                                                boot.R = boot.R,
                                                boot.l = boot.l,
                                                seed = seed )

    if ( symm == TRUE ) {
      res_boot_c2pt[[i]] <- hadron::symmetrise.cf( cf=res_boot_c2pt[[i]] )
      res_boot_c3pt[[i]] <- hadron::symmetrise.cf( cf=res_boot_c3pt[[i]] )
    }

    res_boot[[i]] <- res_boot_c3pt[[i]] / res_boot_c2pt[[i]]
  }
  
  return(res_boot)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Calculate cf ratio R(t) = C2pt_1(t) / C2pt_2(t)                                                      #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
# returns the correlation function ratio
#
# R(t) = C2pt_1(t) / C2pt_2(t)
#
# in form of a cf object.

ratio_cf_2pt <- function( cf2pt_1, cf2pt_2, T, conflist, gradflow_t, flow_c2pt_2=TRUE, boot.R=400, boot.l=2, seed=1234, symm=FALSE ) {

  res_boot_c2pt_1 <- list()
  res_boot_c2pt_2 <- list()
  res_boot <- list()

  x2pt_1 <- array( rep(NA, T*length(conflist)*length(gradflow_t)), dim=c(length(conflist), length(gradflow_t), T) )
  x2pt_2 <- array( rep(NA, T*length(conflist)*length(gradflow_t)), dim=c(length(conflist), length(gradflow_t), T) )
  ratio <- array( rep(NA, T*length(conflist)*length(gradflow_t)), dim=c(length(conflist), length(gradflow_t), T) )

  for( i in 1:length(conflist) ){
    for( j in 1:length(gradflow_t) ) {
        x2pt_1[i,j,] <- complex( real = Re(cf2pt_1[i,j,]), imaginary = Im(cf2pt_1[i,j,]))
        x2pt_2[i,j,] <- complex( real = Re(cf2pt_2[i,j,]), imaginary = Im(cf2pt_2[i,j,]))
    }
  }

  for ( i in 1:length(gradflow_t) ) {

    cf2pt_1_cf <- corr_to_cf(x2pt_1[,i,])

    if ( flow_c2pt_2 == TRUE ) {
      cf2pt_2_cf <- corr_to_cf(x2pt_2[,i,])
    } else { 
      cf2pt_2_cf <- corr_to_cf(x2pt_2[,1,])
    }

    res_boot_c2pt_1[[i]] <- hadron::bootstrap.cf( cf2pt_1_cf,
                                                  boot.R = boot.R,
                                                  boot.l = boot.l,
                                                  seed = seed )

    res_boot_c2pt_2[[i]] <- hadron::bootstrap.cf( cf2pt_2_cf,
                                                  boot.R = boot.R,
                                                  boot.l = boot.l,
                                                  seed = seed )

    if ( symm == TRUE ) {
      res_boot_c2pt_1[[i]] <- hadron::symmetrise.cf( cf=res_boot_c2pt_1[[i]] )
      res_boot_c2pt_2[[i]] <- hadron::symmetrise.cf( cf=res_boot_c2pt_2[[i]] )
    }

    res_boot[[i]] <- res_boot_c2pt_1[[i]] / res_boot_c2pt_2[[i]]
  }
  
  return(res_boot)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Perform a Linear Fit (needed for AIC weighted plateau Fit)                                           #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
fit_ratio <- function(fht_ratio, ti=0, tf=63, tl1=5, tl2=13, tu1=8, tu2=16) {
  fit_R <- c()
  fit_R_tmp <- c()
  chisq_reduced_min <- 0.
  chisq_reduced_max <- 10.
  
  for (l in tl1:tl2) {
    for ( u in tu1:tu2) {
      tryCatch({
        if (u < l+1) {
          message( "skip fit with l = ", l , " u = ", u , ", because ", u , " <= ", l )
          next
        } else if (abs(u-l) < 2) {
          message( "skip fit with l = ", l , " u = ", u , ", because abs(u-l) = ", abs(u-l), " < 2" )
          next
        } else {
          fit_R_tmp <- fit.const( y = fht_ratio$cf.tsboot$t0,
                                  bsamples = fht_ratio$cf.tsboot$t,
                                  x = c(ti:tf),
                                  x.lower = l,
                                  x.upper = u,
                                  useCov = TRUE,
                                  fit.method = "lm",
                                  start.par = 0.2e-2 )
        
          chisq <- fit_R_tmp$chisqr
          dof <- fit_R_tmp$dof
        
          if ( chisq/dof < chisq_reduced_max && chisq/dof > chisq_reduced_min )  {
            fit_R <- rbind(fit_R, cbind(l, u, fit_R_tmp$t0, fit_R_tmp$se, chisq/dof, fit_R_tmp$Qval))
          } else {
            message( "skip fit with l = ", l , " u = ", u , ", because chisq/dof = ", chisq/dof )
            next
          }
        }
      }, error = function(e){cat("Error:", conditionMessage(e), "\n")})
    }
  }
  colnames(fit_R) <- c("t_l", "t_u", "cf_ratio", "err cf_ratio", "chisqr/dof", "pval")
  return(fit_R)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Calculate the weighted sum of normal distributions for uniroot (needed for AIC weighted plateau Fit) #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
wdf <- function(x, w=w, m=m, s=s, a=0 ) {
  res <- sum(w * pnorm(x, mean=m, sd=s) ) - a
  return(res)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Calculate the Akaike Information Criterion (AIC) weight                                              #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
w_aic <- function (chisqr, Npar, Ndata) {
  res <- exp(-0.5*(chisqr + 2*Npar - Ndata))
  return(res)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Perform a AIC weighted plateau Fit                                                                   #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
AIC_fit <- function( Rdata, gradflow_t, filename, x_lab, t_i=0, t_f=63, tl_1=1, tl_2=61, tu_1=3, tu_2=63, x1=0.0015, x2=0.0030, Npoints=1000 ) {
  
  if (is.list(Rdata[[1]]) == FALSE) {
    Rdata <- list(Rdata)
  } else if (is.list(Rdata[[1]]) == TRUE && length(t_i) == 1) {
    tl_1 <- rep(tl_1, length(Rdata))
    tl_2 <- rep(tl_2, length(Rdata))
    tu_1 <- rep(tu_1, length(Rdata))
    tu_2 <- rep(tu_2, length(Rdata))
  }
  
  fit_res <- NA
  
  for (i in 1:length(Rdata)) {
    fit_tmp <- fit_ratio(fht_ratio=Rdata[[i]], ti=t_i, tf=t_f, tl1=tl_1[i], tl2=tl_2[i], tu1=tu_1[i], tu2=tu_2[i])
    fit_res <- rbind(fit_res, fit_tmp)
  }
  
  fit_res <- fit_res[-1,]
  
  if (is.array(fit_res) == FALSE) {
    m_list <- fit_res[3]
    e_list <- fit_res[4]
    chisqr <- fit_res[5]
    Ndata <- fit_res[2] - fit_res[1] + 1
  } else {
    m_list <- fit_res[,3]
    e_list <- fit_res[,4]
    chisqr <- fit_res[,5]
    Ndata <- fit_res[,2] - fit_res[,1] + 1
  }

  Npar <- rep(1,length(Ndata))
  x <- seq(from=x1,to=x2,by=(x2-x1)/(Npoints-1))

  w_list <- w_aic(chisqr, Npar, Ndata)
  w <- w_list / sum ( w_list )

  lower <- min ( m_list ) - 5*max ( e_list )
  upper <- max ( m_list ) + 5*max ( e_list )

  u_m <- uniroot( wdf, lower=lower, upper=upper , check.conv = T, tol = 1.e-12, maxiter = 100000, a=0.50, s=e_list, w=w, m=m_list )
  u_l <- uniroot( wdf, lower=lower, upper=upper , check.conv = T, tol = 1.e-12, maxiter = 100000, a=0.16, s=e_list, w=w, m=m_list )
  u_u <- uniroot( wdf, lower=lower, upper=upper , check.conv = T, tol = 1.e-12, maxiter = 100000, a=0.84, s=e_list, w=w, m=m_list )
  
  res <- list()
  res$x <- lower +  (( ( upper-lower)/1000 ) * (0:1000))
 
  for ( i in 1:length(res$x) ) 
  {
    res$y[i] <- wdf (x=res$x[i], w=w, m=m_list, s=e_list, a=0 )
    c( res$x[i], res$y[i] )
  }

  xlim_l <- u_m$root - 6*(u_m$root - u_l$root)
  xlim_u <- u_m$root + 6*(u_u$root - u_m$root)
  
  # tikzDevice::tikz(file = as.character(paste0("CDF_",filename,".tex")), standAlone = TRUE, width = 5, height = 4)
  # par(mar=c(4,5,0.5,0.5)+.75)
  tikzDevice::tikz(file = as.character(paste0("CDF_",filename,".tex")), standAlone = TRUE, width = 4.6, height = 4,
                  packages=c(getOption( "tikzLatexPackages" ),"\\usepackage{bbold}") )
  par(mar=c(4,4.5,0.6,0.85)+.75)
  plot(res,
       xlim=c(xlim_l,xlim_u),
       type="l",
       lwd=2.5,
       cex.axis=1.8,
       cex.lab=1.8,
       col="dodgerblue3",
       xlab=x_lab,
       ylab="CDF")
  abline(h=0., lwd=1.5, col="black")
  abline(v=u_m$root, lwd=2.5, col="firebrick1")
  plt.errband2(x=rep(u_m$root,3000),
               y=seq(0,2999),
               se1=rep(u_m$root - u_l$root,3000),
               se2=rep(u_u$root - u_m$root,3000),
               lwd=0,
               col="firebrick1")
  ### legend("topleft", inset=0.02, legend=c("AIC CDF","median","68\\% CI"), fill=NULL, col=c("dodgerblue3","firebrick1",alpha( col="firebrick1", a=0.25 )), lty=c(1,1,1), lwd=c(3,3,15), pch=c(NA,NA,NA), border=c(NA,NA,NA), box.col="gray38", cex=1.2)
  legend("topleft", inset=0.02, legend=c("AIC CDF","median","68\\% CI"), fill=NULL, col=c("dodgerblue3","firebrick1",alpha( col="firebrick1", a=0.25 )), lty=c(1,1,1), lwd=c(3,3,15), pch=c(NA,NA,NA), border=c(NA,NA,NA), bty="n", cex=1.35)
  title(paste0(gradflow_t), line=0.6)
  dev.off()
  tools::texi2dvi(file = as.character(paste0("CDF_",filename,".tex")), pdf=TRUE)

  res <- paste0("res_", filename, "_gft", gradflow_t )

  res <- cbind(gradflow_t, u_m$root, u_m$root - u_l$root, u_u$root - u_m$root)
  colnames(res) <- c("gradflow_t", "cf_ratio", "err_l", "err_u")
  write.table(res, file = as.character(paste0("CDF_",filename,".txt")), sep = " ", col.names = TRUE, row.names = FALSE)
 
  return ( NULL )
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Extract the best-redchisqr-pval result from given plateau-fit-result data                            #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
get_best_redchisqr_pval_result <- function ( iratio, tu1, tu2, my_tgf, no2ptflow=FALSE, print_txt=FALSE ) {
  print("function get_best_redchisqr_pval_result() requires library dplyr")

  res <- c()

  for ( itgf in my_tgf ) {
    res_out <- c()

    for ( itu in seq(tu1, tu2, 1) ) {
      if ( itgf == 0. ) {
        if ( no2ptflow == TRUE ) {
          res_in <- read.table(paste0("../fitres_tu", itu, "_ratio_", iratio, "_id_ud_du_no2ptflow_tgf0.txt"), header=TRUE)
        } else {
          res_in <- read.table(paste0("../fitres_tu", itu, "_ratio_", iratio, "_id_ud_du_tgf0.txt"), header=TRUE)
        }
      } else {
        if ( no2ptflow == TRUE ) {
          res_in <- read.table(paste0("../fitres_tu", itu, "_ratio_", iratio, "_id_ud_du_no2ptflow_tgf", itgf, ".txt"), header=TRUE)
        } else {
          res_in <- read.table(paste0("../fitres_tu", itu, "_ratio_", iratio, "_id_ud_du_tgf", itgf, ".txt"), header=TRUE)
        }
      }

      weight_vec <- c()

      for ( i in 1:length(res_in[,1]) ) {
        if ( res_in[i,3] <= 1 ) {
          redchisqr_weight <- 1 - res_in[i,3]
        } else if ( res_in[i,3] > 1 ) {
          redchisqr_weight <- res_in[i,3] - 1
        }

        if ( res_in[i,4] <= 0.5 ) {
          pval_weight <- 1 - 2 * res_in[i,4]
        } else if ( res_in[i,4] > 0.5 ) {
          pval_weight <- 2 * res_in[i,4] - 1
        }

        weight_vec[i] <- redchisqr_weight + pval_weight
      }

      itl <- which.min(weight_vec) - 1

      if ( itu - itl == 1 ) {
        weight_vec <- weight_vec[-which.min(weight_vec)]
        itl <- which.min(weight_vec) - 1
      } else {
        itl <- itl
      }

      res_out_tmp <- cbind( itgf, itl, itu, res_in[which.min(weight_vec),] )
      res_out <- rbind( res_out, res_out_tmp )
    }

    weight_vec <- c()

    for ( i in 1:length(res_out[,1]) ) {
      if ( res_out[i,6] <= 1 ) {
        redchisqr_weight <- 1 - res_out[i,6]
      } else if ( res_out[i,6] > 1 ) {
        redchisqr_weight <- res_out[i,6] - 1
      }

      if ( res_out[i,7] <= 0.5 ) {
        pval_weight <- 1 - 2 * res_out[i,7]
      } else if ( res_out[i,7] > 0.5 ) {
        pval_weight <- 2 * res_out[i,7] - 1
      }

      weight_vec[i] <- redchisqr_weight + pval_weight
    }

    res <- rbind( res, res_out[which.min(weight_vec),] )
  } # end of loop over itgf

  if ( print_txt == TRUE ) {
    if ( no2ptflow == TRUE ) {
      write.table(x=res, file=paste0("final_data_ratio", iratio, "_no2ptflow.txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
    } else {
      write.table(x=res, file=paste0("final_data_ratio", iratio, ".txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
    }
  }

  return(res)
}

#--------------------------------------------------------------------------------------------------------#
#                                                                                                        #
#   Extract the best-redchisqr-pval result from given plateau-fit-result data (version 2)                #
#                                                                                                        #
#--------------------------------------------------------------------------------------------------------#
get_best_redchisqr_pval_result_new <- function ( iratio, tu1, tu2, my_tgf, rm_last_row_for_tu=NULL, no2ptflow=FALSE, print_txt=FALSE ) {
  print("function get_best_redchisqr_pval_result() requires library dplyr")

  res <- c()

  for ( itgf in my_tgf ) {
    res_out <- c()

    for ( itu in seq(tu1, tu2, 1) ) {
      if ( itgf == 0. ) {
        if ( no2ptflow == TRUE ) {
          res_in <- read.table(paste0("../fitres_tu", itu, "_ratio_", iratio, "_id_ud_du_no2ptflow_tgf0.txt"), header=TRUE)
        } else {
          res_in <- read.table(paste0("../fitres_tu", itu, "_ratio_", iratio, "_id_ud_du_tgf0.txt"), header=TRUE)
        }
      } else {
        if ( no2ptflow == TRUE ) {
          res_in <- read.table(paste0("../fitres_tu", itu, "_ratio_", iratio, "_id_ud_du_no2ptflow_tgf", itgf, ".txt"), header=TRUE)
        } else {
          res_in <- read.table(paste0("../fitres_tu", itu, "_ratio_", iratio, "_id_ud_du_tgf", itgf, ".txt"), header=TRUE)
        }
      }

      res_in <- res_in[-length(res_in[,1]),]

      if ( !is.null(rm_last_row_for_tu) ) {
        if ( itu == rm_last_row_for_tu ) {
          res_in <- res_in[-length(res_in[,1]),]
        } else {
          res_in <- res_in
        }
      }

      weight_vec <- c()

      for ( i in 1:length(res_in[,1]) ) {
        if ( res_in[i,3] <= 1 ) {
          redchisqr_weight <- 1 - res_in[i,3]
        } else if ( res_in[i,3] > 1 ) {
          redchisqr_weight <- res_in[i,3] - 1
        }

        if ( res_in[i,4] <= 0.5 ) {
          pval_weight <- 1 - 2 * res_in[i,4]
        } else if ( res_in[i,4] > 0.5 ) {
          pval_weight <- 2 * res_in[i,4] - 1
        }

        weight_vec[i] <- redchisqr_weight + pval_weight
      }

      itl <- which.min(weight_vec) - 1

      res_out_tmp <- cbind( itgf, itl, itu, res_in[which.min(weight_vec),] )
      res_out <- rbind( res_out, res_out_tmp )
    }

    weight_vec <- c()

    for ( i in 1:length(res_out[,1]) ) {
      if ( res_out[i,6] <= 1 ) {
        redchisqr_weight <- 1 - res_out[i,6]
      } else if ( res_out[i,6] > 1 ) {
        redchisqr_weight <- res_out[i,6] - 1
      }

      if ( res_out[i,7] <= 0.5 ) {
        pval_weight <- 1 - 2 * res_out[i,7]
      } else if ( res_out[i,7] > 0.5 ) {
        pval_weight <- 2 * res_out[i,7] - 1
      }

      weight_vec[i] <- redchisqr_weight + pval_weight
    }

    res <- rbind( res, res_out[which.min(weight_vec),] )
  } # end of loop over itgf

  if ( print_txt == TRUE ) {
    if ( is.null(rm_last_row_for_tu) ) {
      if ( no2ptflow == TRUE ) {
        write.table(x=res, file=paste0("final_data_ratio", iratio, "_no2ptflow.txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
      } else {
        write.table(x=res, file=paste0("final_data_ratio", iratio, ".txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
      }
    } else {
      write.table(x=res, file=paste0("final_data_ratio", iratio, "_without_tl", (rm_last_row_for_tu - 2), "_to_tu", rm_last_row_for_tu, ".txt"), append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
    }
  }

  return(res)
}
