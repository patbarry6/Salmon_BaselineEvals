####Genepop2Rubias####
library(tidyverse)


Genepop2rubias_mixture<-function(infile,outfile,digits){

dat<-readLines(infile)
PopIndex<-grep("pop|Pop|POP",dat)
if(1%in%PopIndex) PopIndex<-PopIndex[-1] #Remove the first index if Pop is in the file description
nloci<-PopIndex[1]-2
LociNames<-dat[2:(nloci+1)]


npops<-length(PopIndex)
dat<-dat[-c(1:(nloci+1),PopIndex)]
dat2<-matrix(unlist(str_split(dat,",")),nrow=length(dat),ncol=2,byrow=T)

#What are the names of the pops
pops<-dat%>%
  str_split(.,",")%>%
  lapply(.,`[[`,1)%>%
  {str_trim(unlist(.))}%>%
  str_split(.,"\\.")%>%
  lapply(.,`[[`,1)%>%
  unique()%>%
  unlist()

#Pull the pop designation for each individual. Should match rep(pops,each=popCounts)
Indpops<-dat%>%
  str_split(.,",")%>%
  lapply(.,`[[`,1)%>%
  {str_trim(unlist(.))}%>%
  str_split(.,"\\.")%>%
  lapply(.,`[[`,1)%>%
  unlist()

popCountIndex<-lapply(1:npops,function(x) which(Indpops==pops[x]))

Genos<-str_split(dat2[,2],pattern = " ")
Genos<-lapply(1:nrow(dat2), function(x) Genos[[x]][Genos[[x]]!=""])
GenoMat<-matrix(data=unlist(Genos),nrow=nrow(dat2),ncol=nloci,byrow=T)

HapLoc<-which(nchar(GenoMat[1,])==digits)
DipLoc<-which(nchar(GenoMat[1,])!=digits)

NewGenoMat<-matrix(data=NA,nrow=nrow(GenoMat),ncol=length(DipLoc)*2)

for (L in 1:length(DipLoc)){
  temp <- GenoMat[,L]
  temp2 <- lapply(1:(length(temp)),function(x) as.numeric(substring(temp[x],seq(1,digits*2,digits),seq(digits,digits*2,digits)))) #change width to 3 if uSats
  NewGenoMat[,((L*2-1):(L*2))]<-do.call(rbind, temp2)
}#over L loci

if(length(HapLoc)>0){
  HapList<-lapply(1:length(HapLoc), function(x) c(as.numeric(GenoMat[,HapLoc[x]]),rep(NA,length(GenoMat[,HapLoc[x]]))))
  HapMat<-matrix(data=unlist(HapList),nrow=nrow(NewGenoMat),ncol=length(HapLoc)*2)
  NewGenoMat<-cbind(NewGenoMat,HapMat)
}

NewGenoMat2<-as.data.frame(NewGenoMat)
NewGenoMat2[NewGenoMat2==0]<-"NA"

PopNumbs<-sapply(1:npops,function(x) length(popCountIndex[[x]]))
PopNames<-str_trim(unlist(lapply(1:npops,function(x) rep(pops[x],each=PopNumbs[x]))))
Ind_ID<-unlist(lapply(1:npops,function(x) seq(1,length(popCountIndex[[x]]),1)))
IND_POP_ID<-paste(PopNames,Ind_ID,sep="_")

Mixture<-cbind(rep("mixture",times=sum(PopNumbs)),rep("NA",times=sum(PopNumbs)),as.character(unlist(PopNames)),dat2[,1],NewGenoMat2)
colnames(Mixture)<-c("sample_type","repunit","collection","indiv",rep(unlist(LociNames),each=2))

write.csv(Mixture,outfile,quote = F, row.names=F)
cat(paste("Mixture file",infile,"successfully converted to", outfile,".",sep=" "))
}


##################################


Genepop2rubias_baseline<-function(infile,outfile,digits,ReportingGroupFile){
    if(file.exists(file.path(infile))){
        dat<-readLines(infile)
    }else{cat("Infle not found! Check to make sure the path to the file is correct.")
    }
    
    if(file.exists(file.path(ReportingGroupFile))){
        ReportingGroups<-read.csv(ReportingGroupFile,stringsAsFactors = F)
    } else{cat("ReportingGroupFile not found! Check to make sure the path to the file is correct.")
    }
    
    PopIndex<-grep("pop|Pop|POP",dat)
    if(1%in%PopIndex) PopIndex<-PopIndex[-1] #Remove the first index if Pop is in the file description
    nloci<-PopIndex[1]-2 # How many loci do we have?
    LociNames<-dat[2:(nloci+1)] #What are the loci names?
    
    
    npops<-length(PopIndex) #How many pops?
    dat<-dat[-c(1:(nloci+1))] #Lets get just the geno data
    IndIndex<-paste((grep("pop|Pop|POP",dat)+1),":",c((grep("pop|Pop|POP",dat)[-1])-1,length(dat)),sep="")
    
    PopIndex<-grep("pop|Pop|POP",dat)
    Genos<-dat[-PopIndex]
    
    dat2<-matrix(unlist(str_split(Genos,",")),nrow=length(Genos),ncol=2,byrow=T) #make the geno data a matrix
    
    
    #Pull the pop designation for each individual. Should match rep(pops,each=popCounts)
    Indpops<-lapply(1:nrow(ReportingGroups), function(x)
    rep(ReportingGroups$Pop[x],times=length(eval(parse(text=IndIndex[[x]])))))
    
    
    Genos<-str_split(dat2[,2],pattern = " ") #split the character vector of genotype for each individual
    Genos<-lapply(1:nrow(dat2), function(x) Genos[[x]][Genos[[x]]!=""])
    GenoMat<-matrix(data=unlist(Genos),nrow=nrow(dat2),ncol=nloci,byrow=T)
    
    #Let's figure out if there are any haploid loci in the bunch
    HapLoc<-which(nchar(GenoMat[1,])==digits)
    DipLoc<-which(nchar(GenoMat[1,])!=digits)
    
    NewGenoMat<-matrix(data=NA,nrow=nrow(GenoMat),ncol=length(DipLoc)*2)
    
    for (L in 1:length(DipLoc)){
        temp <- GenoMat[,L]
        temp2 <- lapply(1:(length(temp)),function(x) as.numeric(substring(temp[x],seq(1,digits*2,digits),seq(digits,digits*2,digits)))) #change width to 3 if uSats
        NewGenoMat[,((L*2-1):(L*2))]<-do.call(rbind, temp2)
    }#over L loci
    
    if(length(HapLoc)>0){
        HapList<-lapply(1:length(HapLoc), function(x) c(as.numeric(GenoMat[,HapLoc[x]]),rep(NA,length(GenoMat[,HapLoc[x]]))))
        HapMat<-matrix(data=unlist(HapList),nrow=nrow(NewGenoMat),ncol=length(HapLoc)*2)
        NewGenoMat<-cbind(NewGenoMat,HapMat)
    }
    
    NewGenoMat2<-as.data.frame(NewGenoMat)
    NewGenoMat2[NewGenoMat2==0]<-"NA"
    
    PopNumbs<-sapply(1:npops,function(x) length(Indpops[[x]]))
    PopNames<-str_trim(unlist(lapply(1:npops,function(x) rep(ReportingGroups$Pop[x],each=PopNumbs[x]))))
    IND_POP_ID<-lapply(1:npops,function(x) paste(rep(ReportingGroups$Pop[x],each=PopNumbs[x]),
    (1:length(rep(ReportingGroups$Pop[x],each=PopNumbs[x]))),sep="_"))%>%
    unlist()
    
    
    
    Baseline<-cbind(rep("reference",times=sum(PopNumbs)),as.character(PopNames),as.character(PopNames),str_trim(IND_POP_ID),NewGenoMat2)
    colnames(Baseline)<-c("sample_type","repunit","collection","indiv",rep(unlist(LociNames),each=2))
    
    #check to make sure to check that all the split individual pop identifiers are in the reporting group file
    if(all(Baseline$repunit%in%ReportingGroups[,1])){
        Baseline$repunit<-as.character(plyr::mapvalues(Baseline$repunit,ReportingGroups[,1],ReportingGroups[,2]))
    }else{cat(paste(paste(unique(Baseline$repunit[which((Baseline$repunit%in%ReportingGroups[,1])==F)]),collapse="/")," is/are missing from the ReportingGroupFile.",sep=""))
    }
    write.csv(Baseline,outfile,quote = F, row.names=F)
    cat(paste("Baseline file",infile,"successfully converted to", outfile,".",sep=" "))
}

