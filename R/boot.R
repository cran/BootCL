
BootCL.distribution<-function(chip,sampling.count,total.sampling.count=10000, windowsize=NULL)
{   seed<-sample(1,1:32767)
    if(chip=="HG.U133A")
    { data(HG.U133A); Chromosome.List<-HG.U133A;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="HG.U133B")
    { data(HG.U133B); Chromosome.List<-HG.U133B;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="MG.U74Av2")
    { data(MG.U74Av2); Chromosome.List<-MG.U74Av2;list.whole.count<-nrow(Chromosome.List)
    }  else if(chip=="MG.U74Bv2")
    { data(MG.U74Bv2); Chromosome.List<-MG.U74Bv2;list.whole.count<-nrow(Chromosome.List)
    }  else if(chip=="MG.U74Cv2")
    { data(MG.U74Cv2); Chromosome.List<-MG.U74Cv2;list.whole.count<-nrow(Chromosome.List)
    }  else if(chip=="RG.U34A")
    { data(RG.U34A); Chromosome.List<-RG.U34A;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="RG.U34B")
    { data(RG.U34B); Chromosome.List<-RG.U34B;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="RG.U34C")
    { data(RG.U34C); Chromosome.List<-RG.U34C;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="Hgfocus")
    { data(Hgfocus); Chromosome.List<-Hgfocus;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="Hgu133plus2")
    { data(Hgu133plus2); Chromosome.List<-Hgu133plus2;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="Hgu133A2")
    { data(Hgu133A2); Chromosome.List<-Hgu133A2;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="Mouse4302")
    { data(Mouse4302); Chromosome.List<-Mouse4302;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="Mouse430A2")
    { data(Mouse430A2); Chromosome.List<-Mouse430A2;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="hwgcod")
    { data(hwgcod); Chromosome.List<-hwgcod;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="mwgcod")
    { data(mwgcod); Chromosome.List<-mwgcod;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="rwgcod")
    { data(rwgcod); Chromosome.List<-rwgcod;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="hgug4112a")
    { data(hgug4112a); Chromosome.List<-hgug4112a;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="hgug4110b")
    { data(hgug4110b); Chromosome.List<-hgug4110b;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="mgug4122a")
    { data(mgug4122a); Chromosome.List<-mgug4122a;list.whole.count<-nrow(Chromosome.List)
    } else if(chip=="chicken")
    { data(chickenChip); Chromosome.List<-chickenChip;list.whole.count<-nrow(Chromosome.List)
    } else
    { print("Please enter VALID chip name!")
    } 
    data(total.count.ws)
    if(is.null(windowsize)) {  windowsize<-round(total.count.ws[chip,]$windowsize)}
    diffcount<-rep(0,total.sampling.count)
    RN1  = .C("distribution",
               as.integer(seed),
               as.integer(sampling.count), 
               as.integer(list.whole.count),
               as.integer(total.sampling.count),
               diff=as.integer(diffcount),
               as.integer(Chromosome.List$CHR.NAME),
               as.integer(Chromosome.List$START),
               as.integer(Chromosome.List$END),
               ws=as.integer(windowsize),
               PACKAGE="BootCL")

   return(list(dist=RN1$diff,windowsize=RN1$ws))
}

BootCL.Statistic<-function(chip,ID.data,windowsize=NULL)
{ 
  if(chip=="HG.U133A")
    { data(HG.U133A); Chromosome.List<-HG.U133A
    } else if(chip=="HG.U133B")
    { data(HG.U133B); Chromosome.List<-HG.U133B
    } else if(chip=="MG.U74Av2")
    { data(MG.U74Av2); Chromosome.List<-MG.U74Av2
    }  else if(chip=="MG.U74Bv2")
    { data(MG.U74Bv2); Chromosome.List<-MG.U74Bv2
    }  else if(chip=="MG.U74Cv2")
    { data(MG.U74Cv2); Chromosome.List<-MG.U74Cv2
    }  else if(chip=="RG.U34A")
    { data(RG.U34A); Chromosome.List<-RG.U34A
    } else if(chip=="RG.U34B")
    { data(RG.U34B); Chromosome.List<-RG.U34B
    } else if(chip=="RG.U34C")
    { data(RG.U34C); Chromosome.List<-RG.U34C
    } else if(chip=="Hgfocus")
    { data(Hgfocus); Chromosome.List<-Hgfocus
    } else if(chip=="Hgu133plus2")
    { data(Hgu133plus2); Chromosome.List<-Hgu133plus2
    } else if(chip=="Hgu133A2")
    { data(Hgu133A2); Chromosome.List<-Hgu133A2
    } else if(chip=="Mouse4302")
    { data(Mouse4302); Chromosome.List<-Mouse4302
    } else if(chip=="Mouse430A2")
    { data(Mouse430A2); Chromosome.List<-Mouse430A2
    } else if(chip=="hwgcod")
    { data(hwgcod); Chromosome.List<-hwgcod
    } else if(chip=="mwgcod")
    { data(mwgcod); Chromosome.List<-mwgcod
    } else if(chip=="rwgcod")
    { data(rwgcod); Chromosome.List<-rwgcod
    } else if(chip=="hgug4112a")
    { data(hgug4112a); Chromosome.List<-hgug4112a
    } else if(chip=="hgug4110b")
    { data(hgug4110b); Chromosome.List<-hgug4110b
    } else if(chip=="mgug4122a")
    { data(mgug4122a); Chromosome.List<-mgug4122a
    } else if(chip=="chicken")
    { data(chickenChip); Chromosome.List<-chickenChip
    } else
    { print("Please enter VALID chip name!")
    } 
  
    diff.count<-0;
    select.id<-NULL
    Affy.ID<-NULL
    for(i in 1:length(ID.data$Affy.list))
    { temp.id<-which(as.character(Chromosome.List$NAME)==as.character(ID.data$Access.list[i]))
      select.id<-c(select.id,temp.id)
      Affy.ID<-c(Affy.ID,rep(as.character(ID.data$Affy.list[i]),length(temp.id)))
    }
    select.list<-Chromosome.List[select.id,]
 
    chsort.id<-sort.list(select.list$CHR.NAME)
    select.list<-select.list[chsort.id,]
    Affy.ID<-Affy.ID[chsort.id]
    ch.list<-table(select.list$CHR.NAME)
    ch.name<-names(ch.list)
    sortCH.list<-NULL
    sortAffy.ID<-NULL
    for(i in 1:length(ch.list))
    {  temp.list<-select.list[select.list$CHR.NAME==ch.name[i],]
       temp.Affy<-Affy.ID[select.list$CHR.NAME==ch.name[i]]
       tempsort.id<-sort.list(temp.list$START)
       temp.list<-temp.list[tempsort.id,]
       sortCH.list<-rbind(sortCH.list,temp.list)
       sortAffy.ID<-c(sortAffy.ID,as.character(temp.Affy[tempsort.id]))
    }

    data(total.count.ws)
    if(is.null(windowsize)) {  windowsize<-round(total.count.ws[chip,]$windowsize)}

    n<-nrow(sortCH.list)
    conseq.state<-rep(0,n-1)
    for(i in 1:(n-1))
    {   if( sortCH.list$CHR.NAME[i+1] ==  sortCH.list$CHR.NAME[i])
        {   diff<- sortCH.list$START[i+1] -  sortCH.list$END[i]
#            if(diff>0 & diff<windowsize) 
             if(diff<windowsize) 
            {   diff.count<-diff.count+1
                conseq.state[i]<-1
            }
        }
    }
    return(list(Diff.count=diff.count,windowsize=windowsize,sampling.count=n,
conseq.state=conseq.state,Affy.ID=sortAffy.ID,Access.ID=sortCH.list$NAME))
}


BootCL.Pvalue<-function(Bstat,distribution)
{
    Pvalue<-sum(distribution$dist>=Bstat$Diff.count)/length(distribution$dist)
    return(list(Diff.count=Bstat$Diff.count, Pvalue=Pvalue,dist=distribution$dist,windowsize=distribution$windowsize,sampling.count=
     Bstat$sampling.count,conseq.state=Bstat$conseq.state,Affy.ID=Bstat$Affy.ID,Access.ID=Bstat$Access.ID))
}




BootCL.print<-function(BootP,cutoff=3,affyID.flag=TRUE)
{
   print.list<-NULL
   id<-1
   count<-1
   if(affyID.flag)
  { gene.name<-BootP$Affy.ID
  } else
  { gene.name<-BootP$Access.ID}

   for(i in 1:length(BootP$conseq.state))
   {  if(BootP$conseq.state[i]==1) 
      {   count<-count+1
      } else
      { 
         temp<-as.character(gene.name[(i-count+1):i])
         print.list<-c(print.list,list(temp))
         count<-1
       }
   }
   select.list<-NULL
   for(i in 1:length(print.list))
   { if(length(print.list[[i]])>=cutoff)
        select.list<-c(select.list,list(print.list[[i]]))
   }

   return(list(print.list=print.list,select.list=select.list))
}  
        
     


find.ID<-function(chip,sample,affyID.flag=TRUE)
{

  if(chip=="HG.U133A")
    { data(affy.hgu133a); xx<-affy.hgu133a;data(HG.U133A); Chromosome.List<-HG.U133A
    } else if(chip=="HG.U133B")
    { data(affy.hgu133b); xx<-affy.hgu133b;data(HG.U133B); Chromosome.List<-HG.U133B
    } else if(chip=="MG.U74Av2")
    { data(affy.mgu74av2); xx<-affy.mgu74av2;data(MG.U74Av2); Chromosome.List<-MG.U74Av2
    }  else if(chip=="MG.U74Bv2")
    { data(affy.mgu74bv2); xx<-affy.mgu74bv2;data(MG.U74Bv2); Chromosome.List<-MG.U74Bv2
    }  else if(chip=="MG.U74Cv2")
    { data(affy.mgu74cv2); xx<-affy.mgu74cv2;data(MG.U74Cv2); Chromosome.List<-MG.U74Cv2
    }  else if(chip=="RG.U34A")
    { data(affy.rgu34a); xx<-affy.rgu34a;data(RG.U34A); Chromosome.List<-RG.U34A 
    } else if(chip=="RG.U34B")
    { data(affy.rgu34b); xx<-affy.rgu34b;data(RG.U34B); Chromosome.List<-RG.U34B
    } else if(chip=="RG.U34C")
    { data(affy.rgu34c); xx<-affy.rgu34c;data(RG.U34C); Chromosome.List<-RG.U34C
    } else if(chip=="Hgfocus")
    { data(affy.hgfocus); xx<-affy.hgfocus;data(Hgfocus); Chromosome.List<-Hgfocus
    } else if(chip=="Hgu133plus2")
    { data(affy.hgu133plus2); xx<-affy.hgu133plus2;data(Hgu133plus2); Chromosome.List<-Hgu133plus2 
    } else if(chip=="Hgu133A2")
    { data(affy.hgu133a2); xx<-affy.hgu133a2;data(Hgu133A2); Chromosome.List<-Hgu133A2
    } else if(chip=="Mouse4302")
    { data(affy.mouse4302); xx<-affy.mouse4302;data(Mouse4302); Chromosome.List<-Mouse4302
    } else if(chip=="Mouse430A2")
    { data(affy.mouse430a2); xx<-affy.mouse430a2;data(Mouse430A2); Chromosome.List<-Mouse430A2
    } else if(chip=="hwgcod")
    { data(codelink.hwgcod); xx<-codelink.hwgcod;data(hwgcod); Chromosome.List<-hwgcod
      xx[,2]<-apply(as.matrix(xx[,2]),1,function(x){strsplit(x,split=".",fixed=TRUE)[[1]][1]})
    } else if(chip=="mwgcod")
    { data(codelink.mwgcod); xx<-codelink.mwgcod;data(mwgcod); Chromosome.List<-mwgcod
      xx[,2]<-apply(as.matrix(xx[,2]),1,function(x){strsplit(x,split=".",fixed=TRUE)[[1]][1]})
    } else if(chip=="rwgcod") 
    { data(codelink.rwgcod); xx<-codelink.rwgcod;data(rwgcod); Chromosome.List<-rwgcod
      xx[,2]<-apply(as.matrix(xx[,2]),1,function(x){strsplit(x,split=".",fixed=TRUE)[[1]][1]})
    } else if(chip=="hgug4112a")
    { data(agilent.hgug4112a); xx<-agilent.hgug4112a;data(hgug4112a); Chromosome.List<-hgug4112a
    } else if(chip=="hgug4110b")
    { data(agilent.hgug4110b); xx<-agilent.hgug4110b;data(hgug4110b); Chromosome.List<-hgug4110b
    } else if(chip=="mgug4122a")
    { data(agilent.mgug4122a); xx<-agilent.mgug4122a;data(mgug4122a); Chromosome.List<-mgug4122a
    }  else if(chip=="chicken")
    { data(affy.chickenChip); xx<-affy.chickenChip;data(chickenChip); Chromosome.List<-chickenChip
    } else
    { print("Please enter VALID chip name!")
    } 
    rownames(xx)<-as.character(xx[,1])
    if(affyID.flag==TRUE)
    {  if(is.null(dim(sample)))
       { Access.list<-xx[as.character(sample),2];n<-length(sample)
       } else
       { sample<-sample[,1]
         Access.list<-xx[as.character(sample),2];n<-length(sample)
       }
       if(sum(is.na(Access.list))==0)
       { probename.list<-sample
       } else
       { probename.list<-sample[-which(is.na(Access.list))]
       }
    }  else
    {   if(!is.null(dim(sample)))
       { sample<-sample[,1]
       }
       Access.list<-sample
       probename.list<-NULL
    }
   return(list(Access.list=Access.list,Affy.list=probename.list))
}


BootCL.plot<-function(BootP,xrange=NULL,freq.bootCL=NULL)
{  
   if(is.null(xrange)) {xrange<-c(0,max(BootP$dist)+5)}
   hist(BootP$dist,xlim=xrange,main="Bootstrapping distribution",freq=freq.bootCL,
#       sub=paste("Chromosomal spatial bias of ", chip ,"genes of #",BootP$sampling.count,sep=""),
       xlab=paste("P value : ",as.character(BootP$Pvalue)),ylab="",
       breaks=seq(0,BootP$sampling.count))
      abline(v=BootP$Diff.count,lty=1,col=2,lwd=2)
      text((BootP$Diff),-10,as.character(BootP$Diff),col="blue",cex=1)
}


