#########################################################################################################
####  Example ANC data extracted from Botswana 2014 Spectrum file. Data for Botswana Urban           ####
####  Spectrum files available on request from: http://apps.unaids.org/spectrum/                     ####
####  Code for reading Spectrum file available from: https://github.com/jeffeaton/read-epp-spectrum  ####
#########################################################################################################

## source("~/Documents/Code/R/read-epp-spectrum/read-epp-files.R")
## bw.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Botswana 2014/Botswana 2014_Nat 19_06_14-c"
## bw.eppd <- read.epp.data(paste(bw.path, ".xml", sep=""))

bw.urban.anc.prev <- rbind(c(NA, NA, NA, NA, NA, NA, 0.17, 0.149, 0.192, 0.278, 0.287, 0.313, 0.340, 0.391, 0.371, 0.362, 0.386, 0.382, 0.448, NA, 0.344, 0.353, 0.332, NA, 0.290, NA, 0.292, NA, NA, NA, NA, NA, NA, NA, NA, NA), 
                          c(NA, NA, NA, NA, NA, NA, 0.16, 0.237, 0.342, 0.297, 0.396, 0.431, 0.429, 0.430, 0.427, 0.444, 0.496, 0.402, 0.458, NA, 0.423, 0.367, 0.379, NA, 0.381, NA, 0.370, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                          c(NA, NA, NA, NA, NA, NA,   NA,    NA,    NA, 0.160,    NA, 0.217,    NA, 0.247,    NA, 0.407, 0.340, 0.331, 0.257, NA, 0.282, 0.250, 0.228, NA, 0.283, NA, 0.224, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                          c(NA, NA, NA, NA, NA, NA,   NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA, 0.321, 0.265, 0.279, NA, 0.321, 0.276, 0.267, NA, 0.288, NA, 0.158, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                          c(NA, NA, NA, NA, NA, NA,   NA,    NA,    NA,    NA,    NA, 0.239,    NA, 0.372,    NA, 0.304, 0.296, 0.292, 0.321, NA, 0.315, 0.317, 0.335, NA, 0.279, NA, 0.266, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                          c(NA, NA, NA, NA, NA, NA,   NA,    NA,    NA,    NA,    NA,    NA, 0.282,    NA, 0.320,    NA, 0.319, 0.398, 0.374, NA, 0.362, 0.299, 0.304, NA, 0.378, NA, 0.335, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                          c(NA, NA, NA, NA, NA, NA,   NA,    NA, 0.199,    NA, 0.299,    NA, 0.344,    NA, 0.418,    NA, 0.446, 0.367, 0.433, NA, 0.375, 0.368, 0.373, NA, 0.336, NA, 0.379, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                          c(NA, NA, NA, NA, NA, NA,   NA,    NA, 0.178,    NA, 0.389,    NA, 0.337,    NA, 0.313,    NA, 0.306, 0.346, 0.324, NA, 0.290, 0.284, 0.287, NA, 0.324, NA, 0.244, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                          c(NA, NA, NA, NA, NA, NA,   NA,    NA,    NA, 0.270,    NA, 0.378,    NA, 0.499,    NA, 0.503, 0.500, 0.481, 0.522, NA, 0.465, 0.411, 0.490, NA, 0.391, NA, 0.412, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                          c(NA, NA, NA, NA, NA, NA,   NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA,    NA, NA, 0.349, 0.200, 0.214, NA, 0.311, NA, 0.167, NA, NA, NA, NA, NA, NA, NA, NA, NA))

bw.urban.anc.n <- rbind(c(NA, NA, NA, NA, NA, NA,  58, 841, 881, 1205, 1307, 1232, 786, 1251, 480, 734, 576, 666, 555, NA, 722, 607, 692, NA, 678, NA, 572, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, NA, NA, NA, NA, NA, 128, 796, 803,  799,  626,  751, 802,  796, 576, 702, 494, 458, 593, NA, 659, 560, 702, NA, 527, NA, 547, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, NA, NA, NA, NA, NA,  NA,  NA,  NA,  508,   NA,  589,  NA,  300,  NA, 513, 347, 363, 359, NA, 313, 324, 319, NA, 316, NA, 302, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, NA, NA, NA, NA, NA,  NA,  NA,  NA,   NA,   NA,   NA,  NA,   NA,  NA,  NA, 192, 288, 306, NA, 368, 385, 441, NA, 421, NA, 242, NA, NA, NA, NA, NA, NA, NA, NA, NA),
c(NA, NA, NA, NA, NA, NA,  NA,  NA,  NA,   NA,   NA,  213,  NA,  534,  NA, 349, 540, 567, 522, NA, 551, 518, 487, NA, 486, NA, 262, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, NA, NA, NA, NA, NA,  NA,  NA,  NA,   NA,   NA,   NA, 393,   NA, 362,  NA, 336, 303, 326, NA, 379, 322, 403, NA, 409, NA, 463, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, NA, NA, NA, NA, NA,  NA,  NA, 267,   NA,  262,   NA, 289,   NA, 268,  NA, 341, 425, 389, NA, 534, 598, 667, NA, 405, NA, 458, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, NA, NA, NA, NA, NA,  NA,  NA, 258,   NA,  231,   NA, 266,   NA, 252,  NA, 232, 241, 246, NA, 243, 249, 242, NA, 274, NA, 169, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, NA, NA, NA, NA, NA,  NA,  NA,  NA,  307,   NA,  331,  NA,  471,  NA, 304, 351, 372, 302, NA, 300, 324, 310, NA, 263, NA, 282, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, NA, NA, NA, NA, NA,  NA,  NA,  NA,   NA,   NA,   NA,  NA,   NA,  NA,  NA,  NA,  NA,  NA, NA, 146, 143, 141, NA, 300, NA, 300, NA, NA, NA, NA, NA, NA, NA, NA, NA))

dimnames(bw.urban.anc.prev) <- dimnames(bw.urban.anc.n) <- list(c("Gaborone", "Francistown", "Southern", "South East", "Kweneng East", "Mahalapye", "Serowe Palapye", "Lobatse", "Selebi-Pikwe", "Jwaneng"), 1985:2020)

## M: time series of model mid-year prevalence outputs from 1970 through 2013.
M <- c(0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00111, 0.00123, 0.00136, 0.00150, 0.00166,
       0.00183, 0.00200, 0.00239, 0.00309, 0.00432, 0.00641, 0.00988, 0.01552, 0.02432, 0.03736,
       0.05551, 0.07899, 0.10703, 0.13798, 0.16956, 0.19942, 0.22556, 0.24661, 0.26194, 0.27150,
       0.27566, 0.27503, 0.27045, 0.26344, 0.25587, 0.24928, 0.24303, 0.23681, 0.23020, 0.22382,
       0.21841, 0.21386, 0.21015, 0.20715)

## ANC bias parameter (probit scale)
ancbias <- 0.2637549


#############################################
####  Example of likelihood calculation  ####
#############################################

## Load likelihood functions
source("anclik.R")

## Prepare data for likelihood calculation
bw.urban.anclikdat <- fnPrepareANCLikelihoodData(bw.urban.anc.prev, bw.urban.anc.n, anchor.year=1970L)

qM <- qnorm(M) # probit transformed prevalence

fnANClik(qM + ancbias, bw.urban.anclikdat)
fnANClik(qM + ancbias, bw.urban.anclikdat, VERSION="R")


## Compare performance of C versus R implementation
library(microbenchmark)

microbenchmark(fnANClik(qM + ancbias, bw.urban.anclikdat),
               fnANClik(qM + ancbias, bw.urban.anclikdat, VERSION="R"))



## Sample site-level random effects
b.site <- sample.b.site(qM + ancbias, bw.urban.anclikdat)

## Sample from clinic posterior predictive distribution
pred.site <- sample.pred.site(qM + ancbias, b.site, bw.urban.anclikdat)
