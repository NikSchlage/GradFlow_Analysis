Lfunc <- function(mu,t)
{
  gE <- 0.57721
  return ( log ( 2 * mu^2 * t ) + gE )
}

Nf <- 3  # number flavors, up, down, strange

Nc <- 3
CA <- Nc
CF <- (Nc^2 -1) / (2 * Nc )
TF <- 1 / 2. * Nf


b0 <- 11 * CA / 3. - 4./3. * TF

b1 <- 34./3. * CA^2 - (4 * CF + 20./3. *CA ) * TF

g0 <- 6. * CF

g1 <- CA * CF * ( 223./3. - 16. * log(2) ) - CF^2 * ( 3. + 16. * log(2) ) - 44./3. * CF * TF


C2 <-  -23.8 * CA * CF + 30.4 * CF^2 - 3.92 * CF * TF

# g^(Nf =3 ) (3 GeV) = 1.77 means
muref <- 3.
alpharef <- 1.77^2 / ( 4 * pi )

###################################################
ralpha <- function(alpha,mu)
###################################################
{
 # ... running of alphaQCD
}


###################################################
chimatch <- function(alpha,mu,t)
#  NLO, up to order alpha^3
###################################################
{
  return (
    1 + alpha/(4*pi) * ( g0/2. * Lfunc(mu,t) - 3. * CF * log(3) - 4. * CF * log(2) )
      + (alpha/(4*pi))^2 * (
    g0/4.*(b0 + g0/2.) * Lfunc(mu,t)^2 + ( g1/2. - g0/2.*(b0 + g0/2.)*log(3) -2./3.*g0*(b0+g0/2.)*log(2) )*Lfunc(mu,t) + C2 
  )
    )
}

###################################################
#
probe_matching <- function()
#
###################################################
{
  mu <- muref

  tmax <- 0.01 # GeV^2
  tnum <- 1000

  tt <- (1:tnum)  * tmax / tnum

  alpha  <- alpharef
  dalpha <- 0.05 * alpha

  alpha_s <- rnorm (n = 1000, mean = alpha, sd = dalpha )

  res <- array ( dim=c(tnum, 3) )

  fout <- "probe_matching-smallt.out"

  cat( "# nf    = ", Nf, "\n", file=fout, append=F )
  cat( "# mu    = ", mu, "\n", file=fout, append=T )
  cat( "# alpha = ", alpha, " +/- ", dalpha, "\n", file=fout, append=T )

  for ( i in 1:tnum )
  {
    t <- tt[i]

    u <- chimatch( alpha_s, mu = mu, t = t)

    res[i,] <- c( t, mean(u) , sqrt( var(u) ) )

    cat ( formatC( res[i,], width=25, digits=16, format="e" ), "\n", file=fout, append=T )

  }
 
  return(invisible(res))
}
