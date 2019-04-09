insetScale<-function(breaks,col,insetPos=c(.025,.015,.04,.25),main='',offset=1e-3,at=NULL,labels=NULL,cex=1){
  if(length(breaks)!=length(col)+1)stop('Number of breaks must be one more than colors')
  insetPos<-c(graphics::grconvertY(insetPos[1],'nfc','user'),graphics::grconvertX(insetPos[2],'nfc','user'),graphics::grconvertY(insetPos[3],'nfc','user'),graphics::grconvertX(insetPos[4],'nfc','user'))
  breakPos<-((breaks)-(min(breaks)))/max((breaks)-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  #add a bit of offset to avoid pdf viewers displaying breaks between exact rectangle border meeting
  offsetPos<-breakPos[-1]+c(rep(offset*diff(range(breakPos)),length(breakPos)-2),0)
  graphics::rect(breakPos[-length(breakPos)],insetPos[1],offsetPos,insetPos[3],col=col,xpd=NA,border=NA)
  graphics::rect(insetPos[2],insetPos[1],insetPos[4],insetPos[3],xpd=NA)
  if(is.null(at)){
    at<-pretty(breaks)
    at<-at[at<=max(breaks)&at>=min(breaks)]
  }
  if(is.null(labels))labels<-at
  convertPos<-(at-(min(breaks)))/((max(breaks))-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  graphics::segments(convertPos,insetPos[1],convertPos,insetPos[1]-diff(insetPos[c(1,3)])*.1,xpd=NA)
  graphics::text(convertPos,insetPos[1]-diff(insetPos[c(1,3)])*.175,labels,xpd=NA,adj=c(.5,1),cex=.85*cex)
  graphics::text(mean(insetPos[c(2,4)]),insetPos[3]+diff(insetPos[c(1,3)])*.45,main,xpd=NA,adj=c(.5,0),cex=cex)
  invisible(NULL)
}
slantAxis<-function(side,at,labels=at,srt=ifelse(side %in% c(1,4),-45,45),location=1.2,adj=ifelse(side==2,1,0),axisArgs=list(),...){
  do.call(graphics::axis,c(list(side,at,label=FALSE),axisArgs))
  if(side %in% c(1,3)){
    graphics::text(at, convertLineToUser(location,side), srt = srt, adj = adj, labels = labels, xpd = NA,...)
  }else{
    graphics::text(convertLineToUser(location,side),at, srt = srt, adj = adj, labels = labels, xpd = NA,...)
  }
  return(invisible(NULL))
}
conservativeBoundary<-function(boundaries,base=0){
  boundaries<-sort(boundaries)
  return(ifelse(all(boundaries>base),boundaries[1],ifelse(all(boundaries<base),boundaries[2],base)))
}
convertLineToUser<-function(line,axis=1){
  if(!(axis %in% 1:4))stop(simpleError('Undefined axis'))
  axisPair<-sort((c(axis-1,axis+1)%%4)+1)
  isHeight<-(axis%%2)==1
  isSecond<-axis>2
  thisMar<-graphics::par('mar')[axis]
  marWidth<-thisMar/sum(graphics::par('mar')[axisPair])*(graphics::par('fin')-graphics::par('pin'))[isHeight+1]
  widthPerLine<-marWidth/thisMar
  #find base line + add in if plot doesn't cover whole device e.g. graphics::par(mfrow=c(2,1))
  base<-ifelse(isSecond,graphics::par('fin')[isHeight+1]-widthPerLine*thisMar,widthPerLine*thisMar) + graphics::par('fig')[1+isHeight*2]*graphics::par('din')[isHeight+1]
  func<-if(isHeight)graphics::grconvertY else graphics::grconvertX
  out<-func(base+line*widthPerLine*ifelse(isSecond,1,-1),'inches','user')
  return(out)
}

vir<-read.csv('Viral Transmission MetadataTblS1.csv',stringsAsFactors=FALSE)
vir$RNA<-vir$NucAcid=='RNA'
allTransmissionVars <- c("Fecal.Oral","Arbovirus","Inhalation.Aerosols","Inhalation.Dust","Sexual","Eating","Oral.Bloodstream","Breastfeeding","Maternal.Fetal","Germ.line","Blood.Products","Contact.Skin.or.Eye")
transmissionVars<-allTransmissionVars[apply(vir[,allTransmissionVars],2,sum)>10]
attributeVars <- c("RNA","Spherical","Filamentous","Pleomorphic","Bullet.form","Lipid.Envelope")
vir[,c(attributeVars,transmissionVars)]<-apply(vir[,c(attributeVars,transmissionVars)],2,function(xx)as.logical(xx))

tree<-read.csv('ICTV Master Species List 2018a v1 - ICTV 2018 Master Species #33.csv',stringsAsFactors=FALSE)
taxa<-unique(tree[tree$Genus %in% vir$Virus.genus&tree$Genus!='',c('Phylum','Subphylum','Class','Order','Family','Genus')])
rownames(taxa)<-taxa$Genus
noNa<-vir[!apply(is.na(vir[,c(transmissionVars,attributeVars)]),1,any),]
rownames(noNa)<-sprintf('%s (%s)',noNa$Virus.name,noNa$Host.species)
speciesTaxa<-taxa[noNa$Virus.genus,]
speciesTaxa$Species<-noNa$Virus.name
speciesTaxa$Host<-paste(speciesTaxa$Species,noNa$Host.species)
taxaDist<-outer(1:nrow(speciesTaxa),1:nrow(speciesTaxa),function(xx,yy)mapply(function(xxx,yyy)sum(c(10,1,1,1)*(speciesTaxa[xxx,c('Family','Genus','Species','Host')]!=speciesTaxa[yyy,c('Family','Genus','Species','Host')])),xx,yy))
colnames(taxaDist)<-rownames(taxaDist)<-rownames(noNa)
taxaTree<-ape::nj(as.dist(taxaDist))

glmFits<-lapply(1:length(transmissionVars),function(ii){
  goodAtt<-sapply(attributeVars,function(xx)min(table(noNa[,xx],noNa[,transmissionVars[ii]])))>0
  if(sum(goodAtt)==0)return(NULL)
  phylolm::phyloglm(as.formula(sprintf('%s~%s',transmissionVars[ii],paste(attributeVars[goodAtt],collapse='+'))),noNa,taxaTree)
})
names(glmFits)<-transmissionVars
glmBounds<-t(do.call(rbind,lapply(glmFits,function(xx){
    if(is.null(xx))return(NA)
    zz<-as.data.frame(summary(xx)$coef)
    zz<-zz[rownames(zz)!='(Intercept)',]
    out<-structure(apply(cbind(zz$Estimate+2*zz$StdErr, zz$Estimate-2*zz$StdErr),1,conservativeBoundary),.Names=sub('TRUE','',rownames(zz)))
    out<-out[attributeVars]
})))
rownames(glmBounds)<-attributeVars
glmBounds<-glmBounds/log(2)


pagels<-parallel::mclapply(structure(attributeVars,.Names=attributeVars),function(att){
  lapply(structure(transmissionVars,.Names=transmissionVars),function(trans){
    phytools::fitPagel(taxaTree,structure(noNa[,att],.Names=rownames(noNa)),structure(noNa[,trans],.Names=rownames(noNa)))
  })
},mc.cores=length(attributeVars))
pagelP<-do.call(rbind,lapply(pagels,function(xx)sapply(xx,function(yy)yy$P)))
pagelP[,]<-p.adjust(pagelP,'fdr')
sigPs<-which(pagelP<.1,arr.ind=TRUE)

plotHeats<-function(glmBounds,pagelP){
  layout(matrix(c(0,1:2,0),ncol=4),width=c(.5,1,1,.07))
  par(mar=c(6,0,.1,1))
  maxBound<-max(abs(range(glmBounds,na.rm=TRUE)))
  breaks<-seq(-maxBound-.0001,maxBound+.0001,length.out=201)
  cols<-colorRampPalette(c('red','white','blue'))(200)
  image(1:nrow(glmBounds),1:ncol(glmBounds),glmBounds,breaks=breaks,col=cols,xaxt='n',yaxt='n',xlab='',ylab='')
  nas<-which(is.na(glmBounds),arr.ind=TRUE)
  rect(nas[,1]-.5,nas[,2]-.5,nas[,1]+.5,nas[,2]+.5,col='#DDDDDD',border=NA)
  abline(v=1:nrow(glmBounds)+.5,col='#00000022')
  abline(h=1:ncol(glmBounds)-.5,col='#00000077')
  axis(2,1:ncol(glmBounds),gsub('\\.',' ',colnames(glmBounds)),las=1)
  slantAxis(1,1:nrow(glmBounds),gsub('\\.',' ',rownames(glmBounds)),xpd=NA)
  box()
  insetScale(breaks,cols,insetPos=c(0.1,grconvertX(.02,'ndc','nfc'),.13,.03),at=-2:2,labels=2^(-2:2),main='Fold change in odds')
  cols<-structure(c(heat.colors(3),'white'),.Names=c('<.01','<.1','<.2','>.2'))
  image(1:nrow(pagelP),1:ncol(pagelP),pagelP,breaks=c(-.01,0.01,.1,.2,1.01),col=cols,xaxt='n',yaxt='n',xlab='',ylab='')
  abline(v=1:nrow(pagelP)+.5,h=1:ncol(pagelP)-.5,col='#00000033')
  slantAxis(1,1:nrow(pagelP),gsub('\\.',' ',rownames(pagelP)),xpd=NA)
  box()
  legend(grconvertX(0.01,'ndc','user'),grconvertY(-0.01,'ndc','user'),names(cols),fill=cols,xpd=NA,yjust=0,title='FDR adjusted p-value',ncol=4,bty='n',x.intersp=.1,y.intersp=.8)
}

pdf('figure/heat.pdf',height=4,width=8)
  plotHeats(glmBounds,pagelP)
dev.off()
png('figure/heat.png',height=1000,width=2000,res=250)
  plotHeats(glmBounds,pagelP)
dev.off()

