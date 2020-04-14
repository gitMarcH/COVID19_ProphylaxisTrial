library(tidyverse)
library(pwr)
library(exact2x2)
library(doParallel)
library(purrr)
library(broom)

rm(list=ls())

cl <- makeCluster(12)
registerDoParallel(cl)

# global parameters
pow<-0.8
alpha<-0.05


# data simulation function
trialSim<-function(pPlcb=0.15,es=c(1,0.5,0.5),nTreatArmsInitial=1,nTreatArmsInterim=1,interimTimes=40,trialLength=100,followUp=14,nRecruitPerDay=10,dropTreat1External=NA,epiCurve=F,treatAllocRatios=c(1,1,1),pDropOut=0.1,stoppingRule=1,nExpPerTreatArm=298,alpha=0.05,pow=0.8){
  # pPlcb = attack rate over followUp (default=14) days in the placebo arm
  # es = vector of effect sizes of the various treatments in the trial (smaller number = larger effect: attack rate in treatment arms = pPlcb*es; first treatment = placebo and should be 1)
  # nTreatArmsInitial = initial number of treatment arms
  # nTreatArmsInterim = vector of arms added at each interim analysis timepoint
  # interimTimes = vector of times (in days since trial start) at which interim analyses are planned (note: interim analyses BEFORE trial participants for that day are recruited)
  # trialLength = time (in days) the trial is recruiting (+ 14 days until last patient visit)
  # followUp = follow-up time (in days) between recruitment and completing trial participation; developing sympoms anytime during those days, counts as disease endpoint
  # nRecruitPerDay = total number of participants that can be recruited per day
  # dropTreat1External = either NA or the day number at which treatment 1 is dropped due to external information
  # epiCurve = logical; if TRUE then will assume an epidemic curve which will modify the attack rates and the recruitment numbers per day
  # treatAllocRatios = vector of ratios (compared to placebo) to assign to a given treatment arm (first treatment = placebo)
  # pDropOut = probability of dropping out of study over the 14 days follow-up
  # stopping rule: one of 1 (as in Ventz et al with f=0.25, g=1.5), 2 (Ventz et al with f=0.25, g=2), 3 (Augu's MAMS trials; if one-sided p-value for no treatment benefit >0.2), 4 (Augu's MAMS trials; if one-sided p-value >0.4), 5 (posterior predictive probability of rejecting H0)
  # nExpPerTreat = number of participants desired to be randomised to each experimental arm (for treatment arms not dropped for futility)
  
  runParamsList<-list(
    pPlcb=pPlcb,
    es=es,
    nTreatArmsInitial=nTreatArmsInitial,
    nTreatArmsInterim=nTreatArmsInterim,
    interimTimes=interimTimes,
    trialLength=trialLength,
    followUp=followUp,
    nRecruitPerDay=nRecruitPerDay,
    dropTreat1External=dropTreat1External,
    epiCurve=epiCurve,
    treatAllocRatios=treatAllocRatios,
    pDropOut=pDropOut,
    stoppingRule=stoppingRule,
    nExpPerTreatArm=nExpPerTreatArm,
    alpha=alpha,
    pow=pow
  )
  
  nTot<-nRecruitPerDay*trialLength

  simDat<-data.frame(pid=1:nTot,dayRecruited=NA,treat=factor(rep(NA,nTot),levels=c("placebo",paste(sep="","treat",1:(nTreatArmsInitial+sum(nTreatArmsInterim))))),symptoms=rep(NA,nTot))

  curTreat<-1:nTreatArmsInitial
  curTreatMax<-nTreatArmsInitial
  treatsDroppedEarly<-integer(0)
  
  for(t in 1:trialLength){
    if(!is.na(dropTreat1External) & t==dropTreat1External){
      curTreat<-setdiff(curTreat,1)
    }
    
    if(is.element(el=t,set=interimTimes)){
      curInterim<-which(interimTimes==t)
      
      if(length(curTreat)>0){ 
        # do interim analysis (if current time is interim analysis time) as to whether to keep current treatments in the trial
        for(treatTmp in curTreat){
          idxComp<-which( (simDat$treat=="placebo" | simDat$treat==paste(sep="","treat",treatTmp)) & (t-simDat$dayRecruited+1)>=followUp)
          datTmp<-simDat[idxComp,]
          datTmp$treat<-factor(datTmp$treat)
          
          if(stoppingRule==1){
            p<-fisher.test(table(datTmp$treat,datTmp$symptoms),alternative="greater")$p.value
            tmp<-sum(datTmp$treat==paste(sep="","treat",treatTmp))/nExpPerTreatArm
            if(tmp>1){tmp<-1}
            if(p<0.5*(tmp^1.25)){
              curTreat<-setdiff(curTreat,treatTmp)
              treatsDroppedEarly<-c(treatsDroppedEarly,treatTmp)
              #cat(paste(sep="","Treatment ",treatTmp," dropped at time t = ",t," for futility with posterior probability of beneficial treatment effect p = ",p,".\n"))
            }
          }else if(stoppingRule==2){
            p<-fisher.test(table(datTmp$treat,datTmp$symptoms),alternative="greater")$p.value
            tmp<-sum(datTmp$treat==paste(sep="","treat",treatTmp))/nExpPerTreatArm
            if(tmp>1){tmp<-1}
            if(p<0.25*(tmp^1.5)){
              curTreat<-setdiff(curTreat,treatTmp)
              treatsDroppedEarly<-c(treatsDroppedEarly,treatTmp)
              #cat(paste(sep="","Treatment ",treatTmp," dropped at time t = ",t," for futility with posterior probability of beneficial treatment effect p = ",p,".\n"))
            }
          }else if(stoppingRule==3){
            p<-fisher.test(table(datTmp$treat,datTmp$symptoms),alternative="less")$p.value
            if(p>0.2){
              curTreat<-setdiff(curTreat,treatTmp)
              treatsDroppedEarly<-c(treatsDroppedEarly,treatTmp)
              #cat(paste(sep="","Treatment ",treatTmp," dropped at time t = ",t," for futility with current one-sided null hypothesis probability p = ",p,".\n"))
            }
          }else if(stoppingRule==4){
            p<-fisher.test(table(datTmp$treat,datTmp$symptoms),alternative="less")$p.value
            if(p>0.4){
              curTreat<-setdiff(curTreat,treatTmp)
              treatsDroppedEarly<-c(treatsDroppedEarly,treatTmp)
              #cat(paste(sep="","Treatment ",treatTmp," dropped at time t = ",t," for futility with current one-sided null hypothesis probability  p = ",p,".\n"))
            }
          }else if(stoppingRule==5){
            nNow<-sum(!is.na(datTmp$symptoms))
            nPlcbNow<-sum(!is.na(datTmp$symptoms) & datTmp$treat=="placebo")
            nPlcbCasesNow<-sum(!is.na(datTmp$symptoms) & datTmp$treat=="placebo" & datTmp$symptoms==1)
            nTreatNow<-sum(!is.na(datTmp$symptoms) & datTmp$treat==paste(sep="","treat",treatTmp))
            nTreatCasesNow<-sum(!is.na(datTmp$symptoms) & datTmp$treat==paste(sep="","treat",treatTmp) & datTmp$symptoms==1)
            
            # how many are likely to be recruited into the placebo arms and the current treatment arm in the remaining trial duration, assuming only placebo, current treatment and treatment(s) added at this point will be kept in the trial
            simDatTmp<-simDat[!is.na(simDat$dayRecruited),]
            nLeftToRecruit<-floor( (trialLength-t+1)*nRecruitPerDay*(1-pDropOut) )
            allocRatioTmp<-treatAllocRatios[c(1,treatTmp+1,(curTreatMax+1):(curTreatMax+nTreatArmsInterim[curInterim]))]
            allocRatioTmp<-allocRatioTmp/sum(allocRatioTmp)
            tmpAllocPlcb<-allocRatioTmp[1]
            tmpAllocCurTreat<-allocRatioTmp[2]
            nExpPlcb<-nPlcbNow+floor(sum((t-simDatTmp$dayRecruited+1)<followUp & simDatTmp$treat=="placebo")*(1-pDropOut))+round(nLeftToRecruit*tmpAllocPlcb)
            nExpTreat<-nTreatNow+floor(sum((t-simDatTmp$dayRecruited+1)<followUp & simDatTmp$treat==paste(sep="","treat",treatTmp))*(1-pDropOut))+round(nLeftToRecruit*tmpAllocCurTreat)
            nLeftPlcb<-nExpPlcb-nPlcbNow
            nLeftTreat<-nExpTreat-nTreatNow
            
            # needed maximum proportion of cases in treatment arm
            curPPlcb<-nPlcbCasesNow/nPlcbNow
            curPTreat<-nTreatCasesNow/nTreatNow
            #hReq<-pwr.2p2n.test(n1=nExpPlcb,n2=nExpTreat,sig.level=alpha,power=pow)$h
            #pReq<-(sin(asin(sqrt(curPPlcb))-hReq/2))^2
            # needed maximum proportion of cases in remaining subjects to be recruited
            #pReqLeft=(nExpTreat/(nExpTreat-nTreatNow))*(pReq-(nTreatNow/nExpTreat)*(nTreatCasesNow/nTreatNow))
            # maximum cases in outstanding participants to be recruited to this treatment arm
            #nCaseMax<-ceiling(pReqLeft*(nExpTreat-nTreatNow))
            # decision
            gr<-expand.grid(0:nLeftPlcb,0:nLeftTreat)
            fisherDf<-rbind(
              nPlcbNow+nExpPlcb-nPlcbCasesNow-gr[,1],
              nPlcbCasesNow+gr[,1],
              nTreatNow+nExpTreat-nTreatCasesNow-gr[,2],
              nTreatCasesNow+gr[,2]
            )
            fisherDf<-as.data.frame(fisherDf)
            jointBinomProba<-dbinom(gr[,1],size=nLeftPlcb,prob=curPPlcb)*dbinom(gr[,2],size=nLeftTreat,prob=curPTreat)
            fisherDf<-fisherDf[,jointBinomProba>1e-9]
            gr<-gr[jointBinomProba>1e-9,]
            jointBinomProba<-jointBinomProba[jointBinomProba>1e-9]
            postPredEachScenario<-function(x){
              fisher.test(matrix(byrow=T,ncol=2,x))
            }
            
            postPredP<-sum(
              jointBinomProba * 
                ((map_df(fisherDf, ~postPredEachScenario(.) %>%
                           broom::tidy() %>% 
                           .$p.value))[1,] < alpha & (nPlcbCasesNow+gr[,1])/(nExpPlcb) > (nTreatCasesNow+gr[,2])/(nExpTreat))
            )
            
            if(postPredP<0.1){
            #if(pbinom(round(nCaseMax),size=nExpTreat-nTreatNow,prob=nTreatCasesNow/nTreatNow)<0.01){
              curTreat<-setdiff(curTreat,treatTmp)
              treatsDroppedEarly<-c(treatsDroppedEarly,treatTmp)
              #cat(paste(sep="","Treatment ",treatTmp," dropped at time t = ",t," for futility with posterior predictive probability of rejecting the null at the final analysis p < ",0.05," (nNow=",nNow,", nPlcbNow=",nPlcbNow,", nPlcbCasesNow=",nPlcbCasesNow,", nTreatNow=",nTreatNow,", nTreatCasesNow=",nTreatCasesNow,").\n"))
            }
          }else{
            stop("stoppingRule parameter needs to be one of 1, 2, 3, 4 or 5.")
          }
        }
      }
      
      # add treatments if current time is interim analysis time
      maxTmp<-curTreatMax
      curTreatMax<-curTreatMax+nTreatArmsInterim[curInterim]
      curTreat<-c(curTreat,(maxTmp+1):curTreatMax)
    }
    
    # sample data
    if(length(curTreat)>0){
      treatProbs<-treatAllocRatios[c(1,1+curTreat)]/sum(treatAllocRatios[c(1,1+curTreat)])
      treatTmp<-sample(x=levels(simDat$treat)[c(1,1+curTreat)],size=nRecruitPerDay,replace=T,prob=treatProbs)
      treatProbsIndiv<-pPlcb*es[match(treatTmp,c("placebo",paste(sep="","treat",1:(nTreatArmsInitial+sum(nTreatArmsInterim)))))]
      #treatProbsIndiv<-ifelse(treatTmp=="placebo",pPlcb,pPlcb*es) # leftover when there was only a single effect size for all treatment arms
      simDat$dayRecruited[((t-1)*nRecruitPerDay+1):(t*nRecruitPerDay)]<-t
      simDat$treat[((t-1)*nRecruitPerDay+1):(t*nRecruitPerDay)]<-treatTmp
      simDat$symptoms[((t-1)*nRecruitPerDay+1):(t*nRecruitPerDay)]<-rbinom(nRecruitPerDay,size=1,prob=treatProbsIndiv)
      rm(treatTmp,treatProbs,treatProbsIndiv)
    }
  }
  
  # simulate drop-out
  simDat$symptoms[sample(1:nrow(simDat),size=round(nrow(simDat)*pDropOut))]<-NA
  
  return(list(simDat=simDat,runParamsList=runParamsList,treatsDroppedEarly=treatsDroppedEarly))
}


# analysis function
analyseTrial<-function(dat){
  treats<-setdiff(levels(dat$treat),"placebo")
  resMat<-data.frame(
    treat=treats,
    minTime=NA,
    maxTime=NA,
    nPlacebo=NA,
    nPlaceboNoMiss=NA,
    nPlaceboCase=NA,
    nTreat=NA,
    nTreatNoMiss=NA,
    nTreatCase=NA,
    pValueMiss=NA,
    attackRatePlcb=NA,
    attackRateTreat=NA,
    pValue=NA,
    stringsAsFactors=F
    )
  
  for(j in 1:nrow(resMat)){
    resMat$minTime[j]=min(dat$dayRecruited[dat$treat==resMat$treat[j]])
    resMat$maxTime[j]=max(dat$dayRecruited[dat$treat==resMat$treat[j]])
    
    resMat$nPlacebo[j]<-sum(dat$treat=="placebo" & dat$dayRecruited>=resMat$minTime[j] & dat$dayRecruited<=resMat$maxTime[j])
    resMat$nPlaceboNoMiss[j]<-sum(dat$treat=="placebo" & !is.na(dat$symptoms) & dat$dayRecruited>=resMat$minTime[j] & dat$dayRecruited<=resMat$maxTime[j])
    resMat$nPlaceboCase[j]<-sum(dat$treat=="placebo" & !is.na(dat$symptoms) & dat$symptoms==1 & dat$dayRecruited>=resMat$minTime[j] & dat$dayRecruited<=resMat$maxTime[j])
    
    resMat$nTreat[j]<-sum(dat$treat==resMat$treat[j] & dat$dayRecruited>=resMat$minTime[j] & dat$dayRecruited<=resMat$maxTime[j])
    resMat$nTreatNoMiss[j]<-sum(dat$treat==resMat$treat[j] & !is.na(dat$symptoms) & dat$dayRecruited>=resMat$minTime[j] & dat$dayRecruited<=resMat$maxTime[j])
    resMat$nTreatCase[j]<-sum(dat$treat==resMat$treat[j] & !is.na(dat$symptoms) & dat$symptoms==1 & dat$dayRecruited>=resMat$minTime[j] & dat$dayRecruited<=resMat$maxTime[j])
    
    resMat$attackRatePlcb[j]<-sum(dat$treat=="placebo" & !is.na(dat$symptoms) &dat$symptoms==1 & dat$dayRecruited>=resMat$minTime[j] & dat$dayRecruited<=resMat$maxTime[j])/resMat$nPlaceboNoMiss[j]
    resMat$attackRateTreat[j]<-sum(dat$treat==resMat$treat[j] & !is.na(dat$symptoms) &dat$symptoms==1 & dat$dayRecruited>=resMat$minTime[j] & dat$dayRecruited<=resMat$maxTime[j])/resMat$nTreatNoMiss[j]
    
    resMat$pValueMiss[j]<-fisher.test(matrix(byrow=T,ncol=2,c(resMat$nPlaceboNoMiss[j],resMat$nPlacebo[j]-resMat$nPlaceboNoMiss[j],resMat$nTreatNoMiss[j],resMat$nTreat[j]-resMat$nTreatNoMiss[j])))$p.value
    resMat$pValue[j]<-fisher.test(matrix(byrow=T,ncol=2,c(resMat$nPlaceboNoMiss[j]-resMat$nPlaceboCase[j],resMat$nPlaceboCase[j],resMat$nTreatNoMiss[j]-resMat$nTreatCase[j],resMat$nTreatCase[j])))$p.value
  }
  
  return(resMat)
}


# run several scenarios
Nsim<-1e3

gr<-expand.grid(
  c(0.1),
  c(0.1,0.125,0.15,0.2,0.3,0.4,0.5),
  c(5),
  c(0.5,0.67,0.75)
)

powMat<-data.frame(
  attackratePlcb=gr[,2],
  attackrateTreat=NA,
  effectSize=gr[,4],
  nPerArmParallel2ArmTrialNoMiss=NA,
  pDropOut=gr[,1],
  stoppingRule=gr[,3],
  nNoMissPlcbAvg=NA,
  nNoMissTreat1Avg=NA,
  nNoMissTreat2Avg=NA,
  nNoMissTreat3Avg=NA,
  propTreat1Dropped=NA,
  propTreat2Dropped=NA,
  powerTreat1=NA,
  powerTreat2=NA,
  powerTreat3=NA
)
#powMat<-rbind(powMat,c(0.6,NA,0.67,NA,0.1,1,rep(NA,5)),c(0.6,NA,0.67,NA,0.1,3,rep(NA,5)),c(0.6,NA,0.67,NA,0.1,5,rep(NA,5)))
powMat$attackrateTreat=round(digits=3,powMat$attackratePlcb*powMat$effectSize)

#res<-foreach(j = 1:nrow(powMat)) %dopar% {
for(j in 1:nrow(powMat)){
  #print(j)
  cat(file=paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/simRes_",Sys.Date(),"_Nsim",Nsim,".log"),append=T,paste(sep="",j,"\n"))
  #powMat$nPerArmParallel2ArmTrialNoMiss[j]<-ss2x2(p0=powMat$attackratePlcb[j],p1=powMat$attackrateTreat[j],sig.level=alpha,power=pow)$n0
  tmpSS<-ss2x2(p0=powMat$attackratePlcb[j],p1=powMat$attackrateTreat[j],sig.level=alpha,power=pow)$n0
  
  simResMat<-data.frame(
    run=1:Nsim,
    treat1Stopped=NA,
    treat2Stopped=NA,
    pTreat1=NA,
    pTreat2=NA,
    pTreat3=NA,
    effectDir1=NA,
    effectDir2=NA,
    effectDir3=NA,
    nNoMissPlcb=NA,
    nNoMissTreat1=NA,
    nNoMissTreat2=NA,
    nNoMissTreat3=NA
  )
  
  res<-foreach(i = 1:Nsim) %dopar% {
    #for(i in 1:Nsim){
    require(purrr)
    require(broom)
    datTmp<-trialSim(pPlcb=powMat$attackratePlcb[j],nTreatArmsInitial=2,treatAllocRatios=c(1,1,1,2),es=c(1,powMat$effectSize[j],1,powMat$effectSize[j]),interimTimes=50,trialLength=100,nRecruitPerDay=12,nExpPerTreatArm=tmpSS,stoppingRule=powMat$stoppingRule[j],pDropOut=powMat$pDropOut[j])
    # we assume 2 initial treatments, only one of which is effective and a third, effective, treatment is added at day 50
    resTmp<-analyseTrial(datTmp$simDat)
    
    tmp<-c(
      ifelse(length(datTmp$treatsDroppedEarly)==0,0,ifelse(is.element(el=datTmp$treatsDroppedEarly,set=1),1,0)),
      ifelse(length(datTmp$treatsDroppedEarly)==0,0,ifelse(is.element(el=datTmp$treatsDroppedEarly,set=2),1,0)),
      resTmp$pValue,
      ifelse(resTmp$attackRatePlcb>=resTmp$attackRateTreat,1,-1),
      sum(resTmp$nPlaceboNoMiss),
      resTmp$nTreatNoMiss[1],
      resTmp$nTreatNoMiss[2],
      resTmp$nTreatNoMiss[3]
    )
    return(tmp)
  }
  
  simResMat[,-1]<-matrix(unlist(res),byrow=T,ncol=12)
  
  powMat$propTreat1Dropped[j]<-sum(simResMat$treat1Stopped)/Nsim
  powMat$propTreat2Dropped[j]<-sum(simResMat$treat2Stopped)/Nsim
  powMat[j,c("powerTreat1","powerTreat2","powerTreat3")]<-colSums(simResMat[,c("pTreat1","pTreat2","pTreat3")]<0.05 & simResMat[,c("effectDir1","effectDir2","effectDir3")]==1)/Nsim
  powMat$nNoMissPlcbAvg[j]<-round(mean(simResMat$nNoMissPlcb))
  powMat$nNoMissTreat1Avg[j]<-round(mean(simResMat$nNoMissTreat1))
  powMat$nNoMissTreat2Avg[j]<-round(mean(simResMat$nNoMissTreat2))
  powMat$nNoMissTreat3Avg[j]<-round(mean(simResMat$nNoMissTreat3))
  
  cat(file=paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/simRes_",Sys.Date(),"_Nsim",Nsim,".log"),append=T,paste(sep="",paste(collapse=" ",powMat[j,]),"\n\n"))
}

tmp<-matrix(unlist(res),byrow=T,ncol=10)
powMat[,c("nPerArmParallel2ArmTrialNoMiss","propTreat1Dropped","propTreat2Dropped","powerTreat1","powerTreat2","powerTreat3","nNoMissPlcbAvg","nNoMissTreat1Avg","nNoMissTreat2Avg","nNoMissTreat3Avg")]<-tmp

write.csv(powMat,file=paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/simRes_",Sys.Date(),"_Nsim",Nsim,".csv"))
save(list=ls(),file=paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/simRes_",Sys.Date(),"_Nsim",Nsim,".RData"))

png(paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/simRes_",Sys.Date(),"_Nsim",Nsim,"_powerCurveTreat12_pDropOut0.1.png"),width=16,height=9,unit="in",res=600)
powMat %>%
  filter(pDropOut==0.1) %>%
  ggplot(mapping=aes(x=attackratePlcb,y=powerTreat1,col=as.factor(1-effectSize))) +
  geom_line(lwd=1.15) +
  geom_point(size=3) +
  theme(text=element_text(size=16)) +
  ggtitle("Statistical power to reject the null hypothesis of no effect for an effective initial treatment.") +
  ylab("Power") +
  xlab("Attack rate in control arm") +
  scale_color_manual(values=c("steelblue","orange","salmon"),name="Effect size") +
  #scale_linetype_manual(values=1:3,name="Stopping rule") +
  geom_abline(intercept=0.8,slope=0,lty=2,lwd=1.5,col="darkgrey")
dev.off()

png(paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/simRes_",Sys.Date(),"_Nsim",Nsim,"_powerCurveTreat3_pDropOut0.1.png"),width=16,height=9,unit="in",res=600)
powMat %>%
  filter(pDropOut==0.1) %>%
  ggplot(mapping=aes(x=attackratePlcb,y=powerTreat3,col=as.factor(1-effectSize))) +
  geom_line(lwd=1.15) +
  geom_point(size=3) +
  theme(text=element_text(size=16)) +
  ggtitle("Statistical power to reject the null hypothesis of no effect for an effective treatment added at day 50.") +
  ylab("Power") +
  xlab("Attack rate in control arm") +
  scale_color_manual(values=c("steelblue","orange","salmon"),name="Effect size") +
  #scale_linetype_manual(values=1:3,name="Stopping rule") +
  geom_abline(intercept=0.8,slope=0,lty=2,lwd=1.5,col="darkgrey")
dev.off()

png(paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/simRes_",Sys.Date(),"_Nsim",Nsim,"_powerCurveTreat12_pDropOut0.05.png"),width=16,height=9,unit="in",res=600)
powMat %>%
  filter(pDropOut==0.05) %>%
  ggplot(mapping=aes(x=attackratePlcb,y=powerTreat1,col=as.factor(1-effectSize),)) +
  geom_line(lwd=1.15) +
  geom_point(size=3) +
  theme(text=element_text(size=16)) +
  ggtitle("Statistical power to reject the null hypothesis of no effect for an effective initial treatment.") +
  ylab("Power") +
  xlab("Attack rate in control arm") +
  scale_color_manual(values=c("steelblue","orange","salmon"),name="Effect size") +
  #scale_linetype_manual(values=1:3,name="Stopping rule") +
  geom_abline(intercept=0.8,slope=0,lty=2,lwd=1.5,col="darkgrey")
dev.off()

png(paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/simRes_",Sys.Date(),"_Nsim",Nsim,"_powerCurveTreat3_pDropOut0.05.png"),width=16,height=9,unit="in",res=600)
powMat %>%
  filter(pDropOut==0.05) %>%
  ggplot(mapping=aes(x=attackratePlcb,y=powerTreat3,col=as.factor(1-effectSize))) +
  geom_line(lwd=1.15) +
  geom_point(size=3) +
  theme(text=element_text(size=16)) +
  ggtitle("Statistical power to reject the null hypothesis of no effect for an effective treatment added at day 50.") +
  ylab("Power") +
  xlab("Attack rate in control arm") +
  scale_color_manual(values=c("steelblue","orange","salmon"),name="Effect size") +
  #scale_linetype_manual(values=1:3,name="Stopping rule") +
  geom_abline(intercept=0.8,slope=0,lty=2,lwd=1.5,col="darkgrey")
dev.off()

png(paste(sep="","/Users/marc/work/MLW_LSTM/COVID19_Trial_Burke/sampleSize2ArmParallelTrial_NoDropOut_",Sys.Date(),".png"),width=16,height=9,unit="in",res=600)
powMat %>%
  ggplot(mapping=aes(x=attackratePlcb,y=nPerArmParallel2ArmTrialNoMiss,col=as.factor(1-effectSize),label=nPerArmParallel2ArmTrialNoMiss)) +
  geom_line(lwd=1.15) +
  geom_point(size=3) +
  geom_text(nudge_x=0.008,nudge_y=41) +
  theme(text=element_text(size=16)) +
  ggtitle("Sample size for 80% power for a conventional 2-parallel-arm trial (assuming no attrition).") +
  ylab("Sample size (per arm)") +
  xlab("Attack rate in control arm") +
  scale_color_manual(values=c("steelblue","orange","salmon"),name="Effect size") 
dev.off()

