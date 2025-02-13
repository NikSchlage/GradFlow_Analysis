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