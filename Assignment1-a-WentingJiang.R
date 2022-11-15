library(psych) #fa()
library(lavaan) #cfa()
library(paran) #procedure Horn

#function composite reliability
compositerel<-function(x){
  A<-(sum(x))^2
  B<-sum(1-x^2)
  return(A/(A+B))
}
setwd('/Users/wentingjiang/Desktop/Multi')
#load data
load("cosmetics.Rdata")

#compute centered data
cess2<-cosmetics[,0:9]
cess2<-scale(cess2,center=TRUE,scale=FALSE)

#####################################################################
#confirmatory factor analysis model with correlated latent variables
######################################################################
        '
#measurement model latent variables
cfa1<-' Att_organic=~1*Attitude_organic1+Attitude_organic2+Attitude_organic3         
        Att_packaging=~1*Attitude_packaging1+Attitude_packaging2 +Attitude_packaging3
        Att_crueltyfree=~1*Attitude_crueltyfree1+Attitude_crueltyfree2+Attitude_crueltyfree3
        Att_organic ~~Att_organic
        Att_packaging ~~Att_packaging
        Att_crueltyfree ~~Att_crueltyfree
         Att_organic ~~Att_packaging+Att_crueltyfree
         Att_crueltyfree~~Att_packaging
        '
#fit model on covariance matrix
fitcfa1<-cfa(cfa1,cess2)

#print fitmeasures
fitmeasures(fitcfa1,c("chisq","df", "pvalue", "cfi","tli","rmsea","srmr"))

#summary of results
summary(fitcfa1,fit.measures=TRUE)

#print standardized solution
standardizedSolution(fitcfa1)

#compute residual correlations (crossprod(t(fact3vm$loadings))+diag(fact3vm$uniquenesses))
cor_table <- residuals(fitcfa1, type = "cor")$cov

#compute number of residual correlations with absolute value larger than 0.05 below the diagonal
n<-sum(ifelse(abs(cor_table)>0.05,1,0))/2
print(n/(9*8/2))

#reliability factor scores
d<-standardizedSolution(fitcfa1)
#composite reliability Att_organic
compositerel(d[1:3,4])
#composite reliability Att_packaging
compositerel(d[4:6,4])
#composite reliability Att_crueltyfree
compositerel(d[7:9,4])

#overview table composite reliability
factorscore<-c("Att_organic","Att_packaging","Att_crueltyfree")
reliability<-round(c(compositerel(d[1:3,4]),compositerel(d[4:6,4]),compositerel(d[7:9,4])),3)
data.frame(factorscore,reliability)

#####################################################################
#confirmatory factor analysis model with correlated latent variables and correlated error terms 
######################################################################

#measurement model latent variables
cfa2<-' Att_organic=~1*Attitude_organic1+Attitude_organic2+Attitude_organic3         
        Att_packaging=~1*Attitude_packaging1+Attitude_packaging2 +Attitude_packaging3
        Att_crueltyfree=~1*Attitude_crueltyfree1+Attitude_crueltyfree2+Attitude_crueltyfree3
        Att_organic ~~Att_organic
        Att_packaging ~~Att_packaging
        Att_crueltyfree ~~Att_crueltyfree
        
         Attitude_organic1 ~~a*Attitude_packaging1
         Attitude_organic1 ~~a*Attitude_crueltyfree1
        Attitude_packaging1 ~~a*Attitude_crueltyfree1
         
         Attitude_organic2 ~~b*Attitude_packaging2
         Attitude_organic2~~b*Attitude_crueltyfree2
        Attitude_packaging2 ~~b*Attitude_crueltyfree2
        
        Attitude_organic3 ~~c*Attitude_packaging3
        Attitude_organic3 ~~c*Attitude_crueltyfree3
        Attitude_packaging3 ~~c*Attitude_crueltyfree3
        '
#fit model on covariance matrix
fitcfa1<-cfa(cfa2,cess2)

#print fitmeasures
fitmeasures(fitcfa1,c("chisq","df","pvalue" ,"cfi","tli","rmsea","srmr"))

#summary of results
summary(fitcfa1,fit.measures=TRUE)

#print standardized solution
standardizedSolution(fitcfa1)

#compute residual correlations
cor_table <- residuals(fitcfa1, type = "cor")$cov

#compute number of residual correlations with absolute value larger than 0.05 below the diagonal
n<-sum(ifelse(abs(cor_table)>0.05,1,0))/2
print(n/(9*8/2))
