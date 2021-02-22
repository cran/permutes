## ----ex1----------------------------------------------------------------------
library(permutes)
head(MMN)
nrow(MMN) #how many observations?
length(unique(MMN$Time)) #how many timepoints?

## ----permDry,eval=F-----------------------------------------------------------
#  perms <- permu.test(cbind(Fp1,AF3,F7,F3,FC1,FC5,C3,CP1,CP5,P7,P3,Pz,PO3,O1,Oz,O2,PO4,
#  	P4,P8,CP6,CP2,C4,FC6,FC2,F4,F8,AF4,Fp2,Fz,Cz) ~ Deviant * Session | Time,
#  	data=MMN)
#  ## (output not shown)

## ----perm,eval=F--------------------------------------------------------------
#  library(doParallel)
#  cl <- makeCluster(2) #or more
#  registerDoParallel(cl)
#  perms <- permu.test(cbind(Fp1,AF3,F7,F3,FC1,FC5,C3,CP1,CP5,P7,P3,Pz,PO3,O1,Oz,O2,PO4,
#  	P4,P8,CP6,CP2,C4,FC6,FC2,F4,F8,AF4,Fp2,Fz,Cz) ~ Deviant * Session | Time,data=MMN,
#  	parallel=TRUE)
#  ## (output not shown)

## ----interceptForDeterminism,cache=F,results='hide',echo=F--------------------
set.seed(10)
perms <- permu.test(cbind(Fp1,AF3,F7,F3,FC1,FC5,C3,CP1,CP5,P7,P3,Pz,PO3,O1,Oz,O2,PO4,
	P4,P8,CP6,CP2,C4,FC6,FC2,F4,F8,AF4,Fp2,Fz,Cz) ~ Deviant * Session | Time,
	data=MMN)

## ----permplot-----------------------------------------------------------------
plot(perms)

## ----es-----------------------------------------------------------------------
plot(perms,sig='p')

## ----agg----------------------------------------------------------------------
ROI <- c('Fp1','AF3','F7','F3','FC1','FC5','C3','CP1', 'CP5','CP6','CP2','C4',
	 'FC6','FC2','F4','F8','AF4','Fp2','Fz')
head(perms) #what do we have?
perms2 <- perms[perms$factor == 'Deviant' & perms$measure %in% ROI,]
perms2$sig <- perms2$p < .05
perms2 <- aggregate(sig ~ Time,perms2,sum)
plot(perms2$Time,perms2$sig) #look at all windows
print(perms2[perms2$sig > 10,'Time']) #our arbitrary criterion

## ----verif--------------------------------------------------------------------
print(unique(MMN$Time)[c(88,136)])

## ----lmer---------------------------------------------------------------------
data <- MMN[MMN$Time > 171 & MMN$Time < 265,]
data$amplitude <- rowMeans(data[,ROI])
data <- aggregate(amplitude ~ Deviant + Session + Subject,data,mean)
model <- perm.lmer(amplitude ~ Deviant * Session + (Deviant + Session | Subject),data)

## ----summary------------------------------------------------------------------
print(summary(model))

