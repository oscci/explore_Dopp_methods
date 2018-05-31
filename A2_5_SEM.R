#Max likelihood for laterality models
#by DVM Bishop, started 6th March 2018; updated 8th March 2018

#Prereg document at: https://osf.io/9uaw4/register/565fb3678c5e4a66b5582f67 

#Needs OpenMx, which you get with following command (not CRAN)
#source('http://openmx.psyc.virginia.edu/getOpenMx.R')
require(tidyverse)
require(OpenMx)
require(stargazer) #simple commands for nice tables
require(semTools)
library(DiagrammeR) #for the diagram
library('xlsx')
library(Hmisc) #added for correlation matrix

########################################################
# Select LI data type: 
#
# 1 = Peak L-R diff, using grand mean
# 2 = Mean L-R diff, using grand mean
# 3 = Peak L-R diff, median of all trials
# 4 = Mean L-R diff, median of all trials

#datatype <- as.numeric(readline("Which data type? 1=peak, 2=mean, 3=median peak, 4=median mean:   "))
datatype=2
#outdir <- paste0('C:/Users/zwoodhead/Dropbox/Project A2/A2_Data/SEM',datatype)

nsub <- 30 #change this when all data are collected


##########################################################################
# Read in data
# NB. Need to set working directory to location of data files - or else specify path
mydir<-''
#mydir <-'C:/Users/zwoodhead/Dropbox/Project A2/A2_Data/'
alltask <- read.csv(paste0('LI_data', datatype, '.csv'))

alltask <- cbind(alltask[1:nsub, 3:8], alltask[(nsub+1):(nsub*2), 3:8]) # Reshape into 12 columns

mylabels<-c('ListGen1','PhonDec1','SemDec1','SentGen1','SentComp1','Jabber1',
            'ListGen2','PhonDec2','SemDec2','SentGen2','SentComp2','Jabber2')

colnames(alltask)<-mylabels

# Exclude data with < 12 trials and outliers identified by H&I method
exclusions  <- read.csv(paste0("C:/Users/zwoodhead/Dropbox/Project A2/A2_Data/LI_exclusions", datatype, ".csv"))
exclusions  <- cbind(exclusions[1:nsub, 3:8], exclusions[(nsub+1):(nsub*2), 3:8]) # Reshape into 12 columns

for (t in 1:12){
  alltask[which(exclusions[ , t] > 0), t] = NA # Change excluded values to NA
}

#show means etc for tasks with time1 and time2 adjacent
stargazer(alltask[,c(1,7,2,8,3,9,4,10,5,11,6,12)],type='text')


#Check normality of data
## Have a look at the densities and do Shapiro-Wilks test and QQ plot
# Plots will be written to your working directory, called densities 1-6
#------------------------------------------------------------------

for (i in 1:6){
  png(filename=paste0(outdir,"/densities_",i,".png"))
  par(mfrow=c(2,2))
  for (j in 1:2){
    offset<- 6*(j-1)
  tempdata<-alltask[,(i+offset)]
  myshap<-shapiro.test(tempdata)
   plot(density(tempdata,na.rm=TRUE),main=mylabels[(i+offset)],xlim=c(-6,8),ylim=c(0,.3))
   text(-4,.1,paste('mean = \n',round(mean(tempdata,na.rm=TRUE),2)))
## Plot using a qqplot
   qqnorm(tempdata);qqline(tempdata, col = 2)
   text(.4,-1,paste('Shapiro\np = ',round(myshap$p.value,3)))
  }
  dev.off()
}


#------------------------------------------------------------------
#This bit of script was used at the start of OpenMx in case we wanted to try 'drop one' approach
#to test consistency of bifactor solution. It will drop cases specified in thisdrop.
#Could do this in a loop once we have an optimal approach.
thisdrop <- 0 #specify a number for participant to be dropped
dataRaw      <- mxData( observed=alltask, type="raw" )
if(thisdrop>0){
  dataRaw      <- mxData( observed=alltask[-thisdrop,], type="raw" )}


#------------------------------------------------------------------
#Structural equation models: start with models of means
#------------------------------------------------------------------

#Fully Saturated Model (from fig 4 in prereg document)
#This acts as baseline: it just models means and variance: no relation between variables, 
#and no consistency in LI over time
# We can then test how other models fit when we introduce constraints by
# equalising paths or by modeling covariances.
# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; free to vary by test and occasion
                        labels=c("e1","e2","e3","e4","e5","e6","e7","e8","e9","e10","e11","e12") )
# each has a different name, meaning they are estimated with different values
# means
means        <- mxPath( from="one", to=mylabels, arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =mylabels ) 
satModel <- mxModel("Saturated model", type="RAM",
                    manifestVars=mylabels,
                    dataRaw, resVars,  means)

mysatFit <- mxRun(satModel) #The mxRun command evaluates the model.
summarysat<-summary(mysatFit)
summarysat
#We'll create a Table, bigsummary, to hold fit stats for different models
myCFI <-round(fitMeasuresMx(mysatFit)[7],3)
myrmsea <-round(fitMeasuresMx(mysatFit)[21],3)
BICsat <- round(summarysat$BIC.Mx,1)
bigsummary <- data.frame(cbind(myCFI,myrmsea,BICsat,'Means/var only; no constraints'),row.names='satModel',stringsAsFactors=FALSE)
colnames(bigsummary)<-c('CFI','RMSEA','BIC','Comment')



# #----------------------------------------------------------------------
#Model 1: Tweak saturated model to check test-retest reliability of means
# Means and variances set to be the same for time1 and time2 for each measure
#This is achieved by giving the path the same name, e.g. meanA for A1 and A2
# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; same for each test on time 1 and time 2
                        labels=c("e1","e2","e3","e4","e5","e6","e1","e2","e3","e4","e5","e6") )
# means
means        <- mxPath( from="one", to=mylabels, arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("meanA","meanB","meanC",
                                  "meanD","meanE","meanF",
                                  "meanA","meanB","meanC",
                                  "meanD","meanE","meanF") ) 
myModel1 <- mxModel("Reliability model", type="RAM",
                    manifestVars=mylabels,
                    dataRaw, resVars,  means)

Model1Fit <- mxRun(myModel1) #The mxRun command evaluates the model.
summary1<-summary(Model1Fit)
summary1
myCFI <-round(fitMeasuresMx(Model1Fit)[7],3)
myrmsea <-round(fitMeasuresMx(Model1Fit)[21],3)
BIC1 <- round(summary1$BIC.Mx,1)
bigsummary[2,] <- c(myCFI,myrmsea,BIC1,'Means/var equated for time 1/2')
rownames(bigsummary)[2]<-'Reliability'


# -----------------------------------------------------------------------
# Compare Model1 with Fully Saturated model
bigsummary #show fit for both models
reliabtest<-mxCompare(mysatFit,Model1Fit)
# Make a message (pmessage) that tells us what the model comparison tells us. 
# Here we are looking for a nonsignificant p-value: that tells us that despite having fewer estimated parameters, fit does not suffer
# df is N observations (nsub*nvalues = 28*12) minus N estimated parameters (ep)
pmessage<-'Model 1 fit deteriorates relative to fully saturated! Means differ across test occasions' #default message
if(reliabtest$p[2]>.05){
  pmessage <- 'Model 1 fit does not deteriorate relative to fully saturated; ie means/vars equivalent across test occasions'}
reliabtest 
pmessage


# #---------------------------------------------------------------------------------------------------------------------------
# Now tweak reliability model to set all means and all vars to be the same
# This tests hypothesis that all measures are similarly lateralised
# Expect fit to worsen - assuming measures differ in extent of laterality
# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; same for each test on time 1 and time 2
                        labels=c("e1","e1","e1","e1","e1","e1","e1","e1","e1","e1","e1","e1") )
# means
means        <- mxPath( from="one", to=mylabels, arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("meanA","meanA","meanA",
                                  "meanA","meanA","meanA","meanA","meanA","meanA",
                                  "meanA","meanA","meanA") ) 
EquivmeansModel <- mxModel("Equivmeans model", type="RAM",
                           manifestVars=mylabels,
                           dataRaw, resVars,  means)

EquivmeansFit <- mxRun(EquivmeansModel) #The mxRun command evaluates the model.
summaryE<-summary(EquivmeansFit)

myCFI <-round(fitMeasuresMx(EquivmeansFit)[7],3)
myrmsea <-round(fitMeasuresMx(EquivmeansFit)[21],3)
BIC1a <- round(summaryE$BIC.Mx,1)
bigsummary[3,] <- c(myCFI,myrmsea,BIC1a,'Means/var equated all tests')
rownames(bigsummary)[3]<-'Equivmeans'

meanequivtest<-mxCompare(Model1Fit,EquivmeansFit)

#We predict that means differ, in which case p will be < .05
pmessage<-'Means do not differ between tasks' #default
if(meanequivtest$p[2]<.05){pmessage <- 'Means differ between tasks'}
meanequivtest
pmessage
bigsummary



# #---------------------------------------------------------------------------------------------------------------------------
#Now tweak meanequiv model to set means for taskAB, C and DEF (Model 2a in prereg document)
# This tests dorsal/ventral/mixed model
 
# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; same for each test on time 1 and time 2
                        labels=c("e1","e1","e2","e3","e3","e3","e1","e1","e2","e3","e3","e3") )
# means
means        <- mxPath( from="one", to=mylabels, arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("meanAB","meanAB","meanC",
                                  "meanDEF","meanDEF","meanDEF",
                                  "meanAB","meanAB","meanC",
                                  "meanDEF","meanDEF","meanDEF") ) 
myModel2a <- mxModel("model2a", type="RAM",
                     manifestVars=mylabels,
                     dataRaw, resVars,  means)

Model2aFit <- mxRun(myModel2a) #The mxRun command evaluates the model.
summary2a<-summary(Model2aFit)
myCFI <-round(fitMeasuresMx(Model2aFit)[7],3)
myrmsea <-round(fitMeasuresMx(Model2aFit)[21],3)
BIC2a <- round(summary2a$BIC.Mx,1)
bigsummary[4,] <- c(myCFI,myrmsea,BIC2a,'Means/var equated for AB and DEF')
rownames(bigsummary)[4]<-'DorsalVentral'

#compare with  reliability model (Model1) where task means free to vary
Model2atest<-mxCompare(Model1Fit,Model2aFit)
pmessage<-'Model2a fits as well as model with means free to vary any way'
if(Model2atest$p[2]<.05){
  pmessage <- paste0('Model2a (BIC=',BIC2a,') is worse fit than model with means not constrained (BIC=',BIC1,')')
}
#Add some output showing the mean estimates and testing if they fit the expected pattern
my2a<-summary(Model2aFit)$parameters[4:6,c(1,5,6)]
qmessage<-'Model predicts AB > DEF > C'
rmessage<-'Not confirmed'
if((my2a[1,2]>my2a[3,2])&&(my2a[3,2]>my2a[2,2])){rmessage<-'Confirmed'}
Model2atest
pmessage
qmessage
rmessage
my2a
bigsummary


# #---------------------------------------------------------------------------------------------------------------------------
#Now tweak meanequiv model to set means for taskBD, ACF and E (Model 2b)
# This tests lexical retrieval model

# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; same for each test on time 1 and time 2
                        labels=c("e1","e2","e1","e2","e3","e2","e1","e2","e1","e2","e3","e2") )
# means
means        <- mxPath( from="one", to=mylabels, arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("meanACF","meanBD","meanACF",
                                  "meanBD","meanE","meanACF",
                                  "meanACF","meanBD","meanACF",
                                  "meanBD","meanE","meanACF") ) 
myModel2b <- mxModel("model2b", type="RAM",
                     manifestVars=mylabels,
                     dataRaw, resVars,  means)

Model2bFit <- mxRun(myModel2b) #The mxRun command evaluates the model.
summary2b<-summary(Model2bFit)
myCFI <-round(fitMeasuresMx(Model2bFit)[7],3)
myrmsea <-round(fitMeasuresMx(Model2bFit)[21],3)
BIC2b <- round(summary2b$BIC.Mx,1)
bigsummary[5,] <- c(myCFI,myrmsea,BIC2b,'Means/var equated for ACF and BD')
rownames(bigsummary)[5]<-'LexicalRetrieval'

#compare with  reliability model (Model1) where task means free to vary
Model2btest<-mxCompare(Model1Fit,Model2bFit)
pmessage<-'Model2b fits as well as model with means free to vary any way'
if(Model2btest$p[2]<.05){
  pmessage <- paste0('Model2b (BIC=',BIC2b,') is worse fit than model with means not constrained (BIC=',BIC1,')')
}
#Add some output showing the mean estimates and testing if they fit the expected pattern
my2b<-summary(Model2bFit)$parameters[4:6,c(1,5,6)] #pull out just the rows/cols of interest
qmessage<-'Model predicts BD > ACF'
rmessage<-'Not confirmed'
if(my2b[2,2]>my2b[1,2]){rmessage<-'Confirmed'}
Model2btest
pmessage
qmessage
rmessage
my2b
bigsummary



#---------------------------------------------------------------------------------------------------------------------------
# The following models consider covariances
# Single factor model (means equalized for t1 and t2)
# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels=c("e1","e2","e3","e4","e5","e6","e11","e12","e13","e14","e15","e16") ) # variances can vary at each task/session

# latent variance - Factor1 is the single factor
latVar       <- mxPath( from="Factor1", arrows=2,
                        free=TRUE, values=1, labels ="varFactor1" )
# factor loadings
facLoads     <- mxPath( from="Factor1", to=mylabels, arrows=1,
                        free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), 
                        values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("l1","l2","l3","l4","l5","l6","l1","l2","l3","l4","l5","l6") )#same for each test on time 1 and 2
#The first path is fixed at one - others scaled relative to this

# means - one extra mean for the Factor, but this is set to NA
means        <- mxPath( from="one", to=c(mylabels,'Factor1'), arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T,FALSE), values=c(1,1,1,1,1,1,1,1,1,1,1,1,0),
                        labels =c("meanA","meanB","meanC",
                                  "meanD","meanE","meanF","meanA","meanB","meanC",
                                  "meanD","meanE","meanF",NA) ) #means constant from time 1 to time 2

oneFactorModel <- mxModel("Single Factor Model", type="RAM",
                          manifestVars=mylabels, latentVars="Factor1",
                          dataRaw, resVars, latVar, facLoads, means)

oneFactorFit<-mxRun(oneFactorModel)
#summary(oneFactorFit) #uncomment this to see results
summaryF1<-summary(oneFactorFit)
summaryF1
myCFI <-round(fitMeasuresMx(oneFactorFit)[7],3)
myrmsea <-round(fitMeasuresMx(oneFactorFit)[21],3)
BICF1 <- round(summaryF1$BIC.Mx,1)
bigsummary[6,] <- c(myCFI,myrmsea,BICF1,'= Reliability model + covariances')
rownames(bigsummary)[6]<-'SingleFactor'

#compare with  model where no correlation between measures
Model1Ftest<-mxCompare(oneFactorFit,Model1Fit)

pmessage<-'One factor model no better than model with no covariance between measures'
if(Model1Ftest$p[2]<.05){
  pmessage <- paste0('One factor model (BIC=',BICF1,') is better fit than model with no covariance between measures (BIC=',BIC1,')')
}
pmessage
bigsummary


#---------------------------------------------------------------------------------------------------------------------------
#Bifactor model 
# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels=c("e1","e2","e3","e4","e5","e6","e11","e12","e13","e14","e15","e16") ) # variances can vary at each task/session

# latent variances and covariance: NB assume totally independent, so covariance fixed at zero
latVars      <- mxPath( from=c("Factor1","Factor2"), arrows=2, connect="unique.pairs",
                        free=c(T,F,F), values=c(1,0,1), labels=c("varFactor1","cov","varFactor2") )
#changed the free statement from free =c(T,T,T) (that was error in prereg script: gives underidentified model)

# factor loadings for Factor1 #NB test A loading is fixed to one for this factor
facLoadsFactor1     <- mxPath( from="Factor1", to=mylabels, arrows=1,
                               free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), 
                               values=rep(1,12),
                               labels =c("k1","k2","k3","k4","k5","k6","k1","k2","k3","k4","k5","k6") )

# factor loadings for Factor2 #NB test A loading is fixed to zero for this factor
facLoadsFactor2     <- mxPath( from="Factor2", to=mylabels, arrows=1,
                               free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), 
                               values=c(0,rep(1,5),0,rep(1,5)),
                               labels =c("l1","l2","l3","l4","l5","l6","l1","l2","l3","l4","l5","l6") )

# means #estimated for all except the two factors
means        <- mxPath( from="one", to=c(mylabels,'Factor1','Factor2'), arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T,FALSE,FALSE), values=c(1,1,1,1,1,1,1,1,1,1,1,1,0,0),
                        labels =c("meanA","meanB","meanC",
                                  "meanD","meanE","meanF","meanA","meanB","meanC",
                                  "meanD","meanE","meanF",NA,NA) )

biFactorModel <- mxModel("BiFactor Model", type="RAM",
                         manifestVars=mylabels,
                         latentVars=c("Factor1","Factor2"),
                         dataRaw, resVars, latVars, facLoadsFactor1, facLoadsFactor2, means)

biFactorFit <- mxRun(biFactorModel)
summarybiF<-summary(biFactorFit)
summarybiF
myCFI <-round(fitMeasuresMx(biFactorFit)[7],3)
myrmsea <-round(fitMeasuresMx(biFactorFit)[21],3)
BICF2 <- round(summarybiF$BIC.Mx,1)
bigsummary[7,] <- c(myCFI,myrmsea,BICF2,'2 indepdt factors')
rownames(bigsummary)[7]<-'BiFactor'
mysummary<-summarybiF$parameters[1:10,c(1,3:6)]
mysummary$z<-mysummary$Estimate/mysummary$Std.Error
mcomp<-mxCompare(biFactorFit,oneFactorFit)
pmessage<-'Bi-factor model does not improve fit over one factor model'

if (mcomp$p[2]<.05){pmessage <- paste0('Bi-factor model (BIC=',BICF2,') is better fit than one factor model (BIC=',BICF1,')')}
pmessage
mcomp
bigsummary

#--------------------------------------------------------------------------------------------
# Explore reasons for poor fit

mymat.e<-mxGetExpected(biFactorFit,"covariance")
mymat.o <- cov(alltask,use="pairwise.complete.obs")
mymat.d<-abs(mymat.e-mymat.o) #absolute size of mismatch between obs and expected
mymat.d1<-mymat.e-mymat.o #size of mismatch with sign
cc<-rainbow(ncol(mymat.d1))
png(filename=paste0(outdir,"/covariances.png"))
heatmap(mymat.d1,keep.dendro=FALSE,Rowv=NA,Colv=NA,
        revC=TRUE,col = cm.colors(256),margins=c(5,5),
        main='Signed diff exp/obs cov')
dev.off()



#--------------------------------------------------------------------------------------------
# draw diagram of the bi-factor model, with NS paths omitted
# A shown in red as this has fixed paths to X1 (1) and X2 (0)
#--------------------------------------------------------------------------------------------
require(stringr)
# omxGraphviz(biFactorModel, dotFilename = "bifactor.dot")
# grViz("bifactor.dot") #this will generate a .dot file but it
# is messy, as it shows time 1 and time 2 measures, as well as means

# Script below shows time1/time2 combined and omits means for clarity
# The file for_graphviz is set up in advance and read in and  modified according to results
mybit<-read.csv('for_graphviz.csv',stringsAsFactors = FALSE,header=FALSE) #full list of all paths.
mybit2<-print.data.frame(mybit, 
                 quote=FALSE) #get rid of quotes

#we now want to a) remove rows that are NS and b) put in path coeffs for the rest
thisrow<-12 #NB: first row with path specification is col 13
thatrow<-0 #counter for the summary z scores: NB these exclude measure A! 
for (j in 1:2){#each *factor* (not each test occasion - these are collapsed in diagram)

for (i in 1:6){ #each task 
    thisrow<-thisrow+1

    if(i>1){ #measure A is fixed so not in the table
      thatrow<-thatrow+1
    if(mysummary$z[thatrow]<1.65)
    {mybit2[thisrow,]<-''} #delete this one
    else{
      pathlabel<-round(mysummary$Estimate[thatrow],2)
      bb<-mybit2[thisrow,]
      bb<-str_replace(bb,'xx',as.character(pathlabel))
      mybit2[thisrow,]<-bb
    }
    } #loop to here when i is 1: no action
  }
}
mybit2[19,]<-'' #delete path for A to X2: this one was fixed to 0
dotFilename<-'dottry.dot'
write.table(mybit2, dotFilename, append = FALSE,
            row.names = FALSE, col.names = FALSE,quote=FALSE)


grViz(dotFilename)

#NB If this figure is published, need to make following points:
# This is simplified path diagram. It shows just one measure per variable, 
# when in fact there were two, and it does not show means, though these were estimated.
# Also nonsignificant paths are omitted.

# Dorothy addition to calculate factors.
alltask$FacA.1 <- alltask$ListGen1+
                  2.07*alltask$PhonDec1+
                  1.95*alltask$SemDec1+
                  2.86*alltask$SentGen1+
                  2.07*alltask$SentComp1
alltask$FacA.2 <- alltask$ListGen2+
  2.07*alltask$PhonDec2+
  1.95*alltask$SemDec2+
  2.86*alltask$SentGen2+
  2.07*alltask$SentComp2

alltask$FacB.1 <-.73*alltask$SentComp1+.82*alltask$Jabber1
alltask$FacB.2 <-.73*alltask$SentComp2+.82*alltask$Jabber2
par(mfrow=c(2,2))
plot(alltask$FacA.1,alltask$FacA.2)
plot(alltask$FacB.1,alltask$FacB.2)
plot(alltask$FacA.1,alltask$FacB.1)
plot(alltask$FacA.2,alltask$FacB.2)

myncol<-ncol(alltask)
rcorr(as.matrix(alltask[,(myncol-3):myncol]))

alltask$meansess1 <- rowMeans(alltask[,1:6])
alltask$meansess2 <- rowMeans(alltask[,7:12])
cor(alltask$meansess1,alltask$meansess2)
plot(alltask$meansess1,alltask$meansess2)