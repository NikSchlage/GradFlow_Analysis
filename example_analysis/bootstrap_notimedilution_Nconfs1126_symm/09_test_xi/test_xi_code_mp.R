source("./xi_code_mp.R")

mu <- muref

tmax <- 0.01 # GeV^2
tnum <- 1000

tt <- (1:tnum)  * tmax / tnum

alpha  <- alpharef
dalpha <- 0.05 * alpha

set.seed(12354)
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