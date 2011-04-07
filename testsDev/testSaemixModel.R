####################################################################################
####				Tests for SaemixModel				####
####################################################################################
cat("Beginning tests on SaemixModel class...\n")

setGeneric(name="showall",
  def=function(object){standardGeneric("showall")}
)

source("../R/SaemixModel.R")

model1cpt<-function(psi,id,xidep) { 
	  dose<-xidep[,1]
	  tim<-xidep[,2]  
	  ka<-psi[id,1]
	  V<-psi[id,2]
	  CL<-psi[id,3]
	  k<-CL/V
	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
	  return(ypred)
}

model1cpt.nodose<-function(psi,id,xidep) { 
	  tim<-xidep[,1] 
	  ka<-psi[id,1]
	  V<-psi[id,2]
	  CL<-psi[id,3]
	  k<-CL/V
	  ypred<-100*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
	  return(ypred)
}

wrong.model<-function(x) {
  x**2
}

# Only class
clmod<-new(Class="SaemixModel",model1cpt,psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))))

############ FAIL
# Missing model: not possible, should print a warning
y<-saemixModel()

# Non-existent model
y<-saemixModel(model="model2cpt")

y<-saemixModel(model=model2cpt)

# Wrong model  (only 1 argument)

y<-saemixModel(model=wrong.model)

# Missing psi0: not possible, should print a warning
y<-saemixModel(model=model1cpt)

#### Size mismatch
# Omega matrix of size 2
y<-saemixModel(model="model1cpt",psi0=matrix(c(1.,30,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))), covariance.model=matrix(c(1,0,0,1),ncol=2,byrow=TRUE))

# Omega not a square matrix
y<-saemixModel(model="model1cpt",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,
byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))), covariance.model=matrix(c(1,0,0,1,0,1),ncol=3,byrow=TRUE))

# Fixed estim not same size
y<-saemixModel(model="model1cpt",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),fixed.estim=c(1,0,1,1)) 

# Transform par not same size
y<-saemixModel(model="model1cpt",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,0))

############ SUCCESS
##### Minimal model (only name and psi0)
y<-saemixModel(model=model1cpt,psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))))

# same, with model name
y<-saemixModel(model="model1cpt",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,
byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))))

# only CI for fixed effects, given as a vector
y<-saemixModel(model="model1cpt",psi0=c(ka=1,V=20,CL=0.5))

# Only CI for fixed effects but covariate model given
y<-saemixModel(model="model1cpt",psi0=c(ka=1,V=20,CL=0.5), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE))

##### Giving all object specifications
y<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,
byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),
covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
error.model="combined")

##### Method show
y<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,
byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE),
fixed.estim=c(1,1,1),
covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
error.model="combined")

show(y)

showall(y)

##### Method plot

y<-saemixModel(model=model1cpt.nodose,description="One-compartment model with first-order absorption, dose=100", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

plot(y)

#####

cat("Model validity not checked (we don't check whether the model corresponds to Monolix format).\n")
cat("End of tests on SaemixModel class.\n")
####################################################################################
