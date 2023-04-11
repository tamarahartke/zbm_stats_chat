## stat chat 14 April 2023
## playing with the data from Dianzinga et al 2020: https://onlinelibrary.wiley.com/doi/10.1111/jbi.13959
## data from Dryad: https://datadryad.org/stash/dataset/doi:10.5061%2Fdryad.18931zcvk
## data pre-processing: saved as csv

# changelog ---------------------------------------------------------------
# framework 11 April 2023, TRH

# analysis plan -----------------------------------------------------------
## the original paper
  # alpha diversity: Hill numbers q = 0 (species richness) and q = 2 (Simpson index; evenness)
    iNEXT::iNEXT()
  # beta diversity: LCBD (local contribution to beta diversity) index
    adespatial::
  # linear mixed models (Gaussian errors) 
      # response variables: species richness, Simpson index, LCBD, total abundance (what about Hill numbers 0H and 2H?)
          # abundance, richness, Simpsons transformed: log[response + 0.5]
      # explanatory variables: plant identity, elevation, precipitation, fragmentation, habitat diversity, habitat amount, fragmentation:habitat diversity, fragmentation:habitat amount, habitat diversity:habitat amount
          # continuous expl variables scaled and centered (mean = 0, sd = 1)
      # random effect (spatial autocorrelation): site nested in transect, grid cell 
      # unit of replication: plant (did not nest plant within site; see methods for rationale)
    lme4::lmer()
    nloptr::lmerControl() # optimiser nloptwrap
    car::Anova()
    visreg::visreg()
      # model checking: visual inspection of residuals, predictmeans::CookD, spatial autocorrelation ape::Moran-I
  # NMDS: Hellinger-transformed data
    vegan::capscale()
  # multiscale analysis: not clear from the methods what these models looked like; check appendix
    MuMIn::r.squaredGLMM() # fdr correction

## what we want to do. some ideas...
    # glmmTMB, buildmer
    # DHARMa
    # emmeans, ggeffects

# housekeeping ------------------------------------------------------------
# install packages (in necessary) and load
  if (!require("Require")) install.packages("Require")
  my_packages <- c("tidyverse", "iNEXT", "adespatial", "lme4", "nloptr", "car", "visreg", "vegan", "MuMIn",
                   "glmmTMB", "buildmer", "DHARMa", "emmeans", "ggeffects")
  Require::Require(my_packages)
    
# set seed
  set.seed(42)

# data wrangling ----------------------------------------------------------
    
# read in the data
thrips <- read.csv2("DianzingaetalJBI2020.csv", stringsAsFactors = TRUE)
  # columns 1:41 are the community matrix (each column is one species of thrip)
  # columns 42:46 are date and location (transect, site, plant) data for each sample; they DO NOT combine to create unique identifiers!
  # column 47: elevation of the site in meters
  # columns 48:49: latitude and longitude. Achtung! Metadata from authors says these are in degrees, but they appear to actually be UTM coordinates 
  # column 50: grid cell for that sampling point (see Methods of paper)
  # columns 51:52: mean annual Temperature (°C) and rainfall (mm) for each site from MeteoFrance
  # columns 53:102: proportion of land uses sug (sugar cane), veg (vegetable crops), for (forest), plfor (forest plantation), gra (grassland), rock (bare rock), sava (savannah), orch (orchards), urb (urban) or habitat fragmentation (fragm; not sure how this was calculated) at 100 m, 300 m, 600 m, 1000 m, and 3000 m from the centre of the sampling site. 

# we will use the community matrix on its own a lot. Maybe we want to create a separate object for it?
  
# maybe we want to sum plants to create site-level observations?
  
# Hill numbers ------------------------------------------------------------

# observed species richness  
thrips$rich <- vegan::specnumber(thrips[ , 1:41])
  
# a first go at calculating Hill numbers
thrips_Hill <- iNEXT(t(thrips[ , 1:41]), q = c(0, 2)) # need to transpose the community matrix for iNEXT
  # vignette: http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_UserGuide.pdf
  # another walk-through: https://www.davidzeleny.net/anadat-r/doku.php/en:rarefaction_examples
  # even though we only asked for H0 and H2, it also gives us H1 in the output.
  # in the output object, we find estimates of each metric (qD) and sample coverage (SC), as well as their confidence intervals, based on rarefaction and extrapolation.
  thrips_Hill$DataInfo$S.obs == thrips$rich # this should match 
  thrips_Hill0 <- thrips_Hill$iNextEst$size_based[thrips_Hill$iNextEst$size_based$Method == "Rarefaction" & thrips_Hill$iNextEst$size_based$Order.q == 0, ] # the paper used rarefaction to estimate H0 and H2. Not clear how the authors chose the level, though. Minimum number of taxa and minimum number of individuals found in a samples is 1, maximum 7 and 107, respectively. iNext is giving us estimates for various m, up to (the total number of individuals observed in that samples - 1).
  summary(thrips$rich)
  summary(rowSums(thrips[, 1:41]))
# check out sample coverage
  # plot(thrips_Hill, type = 3) 
  summary(thrips_Hill0$SC) # we have some plots with really poor sample coverage!
 
# let's try a couple of different approaches and see if they are anything like what the paper shows
  # use estimateD because it is a bit faster, with easier-to-use output
  thrips_D_coverage <- estimateD(t(thrips[ , 1:41]), datatype = "abundance", base = "coverage") # compare based on coverage
  thrips_D_size <- estimateD(t(thrips[ , 1:41]), datatype = "abundance", base = "size") # compare based on sample size
  
  # I'm also working on a paper where we calculated Hill numbers like this:
  hill_0 <- iNEXT::ChaoRichness(x = t(thrips[ , 1:41]), datatype = "abundance") %>% 
    dplyr::select("Estimator") %>% 
    rowSums()
  hill_2 <- iNEXT::ChaoSimpson(x = t(thrips[ , 1:41]), datatype = "abundance") %>% 
    dplyr::select("Estimator") %>% 
    rowSums()
  
# discuss: how else could we do this? Do any of these approaches clearly match the results in the paper? Is it clear when they used observed species richness and when they used estimated richness? 