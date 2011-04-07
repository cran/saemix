####################################################################################
####				Tests for SaemixData				####
####################################################################################

# ECO TODO: verifier le nb de sujets ? au moins 2 sujets pour faire tourner saemix ? aussi warning si la reponse n'a que 2 modalites ?

############################
# Defining generic for new methods
setGeneric(name="read.saemixData",
  def=function(object){standardGeneric("read.saemixData")}
)

setGeneric(name="showall",
  def=function(object){standardGeneric("showall")}
)

setGeneric(name="psi",
  def=function(object,indiv.par){standardGeneric("psi")}
)

setGeneric(name="phi",
  def=function(object,indiv.par){standardGeneric("phi")}
)

setGeneric(name="eta",
  def=function(object,indiv.par){standardGeneric("eta")}
)

############################
# SaemixData
source("../R/SaemixData.R")

# SaemixRes
source("../R/SaemixRes.R")

# SaemixRes
source("../R/SaemixModel.R")

# SaemixObject
source("../R/SaemixObject.R")

####################################################################################
