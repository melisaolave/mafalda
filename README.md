
MAny Functions to Approach Large DAtasets (mafalda) in Ecology and Evolution is (going to be) an R library (sorry, still in progress)

questions to molave@mendoza-conicet.gob.ar

follow updates in Twitter @Melisa_Olave and Mastodon melisaolave@ecoevo.social


Description of functions included (v0.0):


get.gdi = a function to compute genealogical divergence index (gdi) (Jackson et al. 2017, Syst. Biol.) as proposed by Leache et al. 2019 (Syst. Biol). Gdi is used to compute the strength of divergence between pair of taxa for species delimitation analyses performed in bpp (Yang and Rannala 2010, PNAS) or ibpp (Solis-Lemus et al. 2015, Evolution).

Usage: 

wd =  working directory, defualt current directory

mcmc = file name containing mcmc output of bpp/ibpp

burnin = burnin for mcmc, 0.1 as default

out.name = output name, default gdi.output.txt


Example:

get.gdi(wd=getwd(), mcmc="bpp_mcmc.out", burnin=0.1, out.name="gdi.output.txt")

