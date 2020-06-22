#remove specific locus from Genepop file#

LocusRM<-function(input="OncorSEAK.txt",locus2rm='One_lpp1-44',NewFileName='OncorSEAK_lpp1-44rm.txt',digits=3,PopNameAsPrefix=T,PopName=''){

  # input="Baseline36Loci.txt"
  # locus2rm=c("One_lpp1-44","One_MHC2_190")
  # NewFileName='Baseline34Loci.txt'
  # digits=3
  # PopNameAsPrefix=T
  
  
  
  Infile<-input
  dat<-readLines(Infile)
  #replace all the Pop with POP
  dat[grep(x=dat,pattern="POP|Pop")]<-"POP"
  nloci<-grep(x=dat,pattern="POP")[1]-2
  dat2<-str_split(dat,' ') #split lines on tab delimeter
  dat2<-lapply(1:length(dat2),function(x) dat2[[x]][(dat2[[x]]!='')]) #get rid of blank entries
  npops<-length(grep(x=dat,pattern="POP"))
  

  #If we use the prefix on the individual to identify pops
  if (PopNameAsPrefix==T){
    loci<-sapply(1:nloci+1,function(x) dat[[x]])
    AnalysisName<-dat2[[1]]
    
    dat2<-dat2[-(1:(nloci+1))]
    dat2[[length(dat2)+1]]<-"POP"
    
    
    dat2<-dat2[-which(dat2=="POP")]
    nInd<-length(dat2)
    
    
    DatMat<-matrix(data=unlist(dat2),nrow=nInd,ncol=nloci+2,byrow=T)
    
    
  PopNames<-unique(matrix(data=unlist(str_split(DatMat[,1],"\\.[0-9]+$")),nrow=nInd,ncol=2,byrow=T)[,1])
  
  PopIndex<-lapply(1:length(PopNames),function(P) which(sapply(str_split(DatMat[,1],"\\.[0-9]+$"),'[[',1)==PopNames[P]))
  DatMat<-DatMat[,-(which(loci%in%locus2rm)+2)]
  nloci<-nloci-length(locus2rm)
  
  #do we want 3 digit genepop format?
  DigScores<-unlist(lapply(1:nInd,function(i) paste(paste(DatMat[i,1:2],collapse=" "),paste(sapply((3:(nloci+2)),function(x) paste(formatC(as.numeric(substring(DatMat[i,x], seq(1, (digits + 1), digits),seq(digits, digits * 2, digits))),width=digits,format="d",flag="0"),collapse="")),collapse=" "),sep=" ")))
  
  #now to write all of this to a file!
  write.table(file=NewFileName,x=gsub(pattern='.txt',x=NewFileName,replacement=''),col.names=F,row.names=F,quote=F)
  write.table(file=NewFileName,x=loci[-which(loci%in%locus2rm)],append=T,col.names = F,row.names=F,quote=F)
  for (p in 1:npops){
    write.table(file=NewFileName,x="Pop",append=T,col.names = F,row.names=F,quote=F)
    write.table(file=NewFileName,x=DigScores[PopIndex[[p]]],append=T,col.names = F,row.names=F,quote=F)
  }
  
  }#if
  
  #if we specify pops in the function
  if(PopNameAsPrefix==F){
      loci<-sapply(1:nloci+1,function(x) dat[[x]])
      AnalysisName<-dat2[[1]]
      
      dat2<-dat2[-(1:(nloci+1))]
      dat2[[length(dat2)+1]]<-"POP"
      PopIndexLines<- sapply(1:npops,function(x) paste(which(dat2=="POP")[x],(which(dat2=="POP")[x+1]-2),sep=":"))
      IndivLineNums<-sapply(1:npops,function(x) paste(which(dat2=="POP")[x]+1,(which(dat2=="POP")[x+1])-1,sep=":"))
      PopIndex<-sapply(1:npops,function(x) paste(as.numeric(unlist(str_split(IndivLineNums[x],":")))-x,collapse =":"))
        
      dat2<-dat2[-which(dat2=="POP")]
      nInd<-length(dat2)
      
      
      DatMat<-matrix(data=unlist(dat2),nrow=nInd,ncol=nloci+2,byrow=T)
      
      
      PopNames<-PopName
    
    
   
    DatMat<-DatMat[,-(which(loci%in%locus2rm)+2)]
    nloci<-nloci-length(locus2rm)
    
    #do we want 3 digit genepop format?
    DigScores<-unlist(lapply(1:nInd,function(i) paste(paste(DatMat[i,1:2],collapse=" "),paste(sapply((3:(nloci+2)),function(x) paste(formatC(as.numeric(substring(DatMat[i,x], seq(1, (digits + 1), digits),seq(digits, digits * 2, digits))),width=digits,format="d",flag="0"),collapse="")),collapse=" "),sep=" ")))
    
    #now to write all of this to a file!
    write.table(file=NewFileName,x=gsub(pattern='.txt',x=NewFileName,replacement=''),col.names=F,row.names=F,quote=F)
    write.table(file=NewFileName,x=loci[-which(loci%in%locus2rm)],append=T,col.names = F,row.names=F,quote=F)
    for (p in 1:npops){
      write.table(file=NewFileName,x="Pop",append=T,col.names = F,row.names=F,quote=F)
      write.table(file=NewFileName,x=DigScores[eval(parse(text=PopIndex[p]))],append=T,col.names = F,row.names=F,quote=F)
    }
    
}#if
    }#function
  


