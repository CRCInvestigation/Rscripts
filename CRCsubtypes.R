# Short description: The following function can determine the molecular subtype of colon cancer samples, based on different classification algorithms. Work on Affymetrix HGU133A and HGU133Plus2.0 .CEL files
# author: Zsofia Sztupinszki; contact: sztup at hotmail.com
# Implemented algorithms
	#Budinska: Budinska, E. et al. Gene expression patterns unveil a new level of molecular heterogeneity in colorectal cancer. The Journal of pathology 231, 63-76, doi:10.1002/path.4212 (2013).
	#CCHS: Dekervel, J. et al. Hypoxia-driven gene expression is an independent prognostic factor in stage II and III colon cancer patients. Clinical cancer research : an official journal of the American Association for Cancer Research 20, 2159-2168, doi:10.1158/1078-0432.CCR-13-2958 (2014).
	#Chang95: Chang, W. et al. Gene expression profiling-derived immunohistochemistry signature with high prognostic value in colorectal carcinoma. Gut 63, 1457-1467, doi:10.1136/gutjnl-2013-305475 (2014).
	#CIN25: Carter, S. L., Eklund, A. C., Kohane, I. S., Harris, L. N. & Szallasi, Z. A signature of chromosomal instability inferred from gene expression profiles predicts clinical outcome in multiple human cancers. Nat Genet 38, 1043-1048, doi:ng1861 [pii] 10.1038/ng1861 (2006).
	#Cologuideex: Agesen, T. H. et al. ColoGuideEx: a robust gene classifier specific for stage II colorectal cancer prognosis. Gut 61, 1560-1567, doi:10.1136/gutjnl-2011-301179 (2012).
	#Cologuidepro: Sveen, A. et al. ColoGuidePro: a prognostic 7-gene expression signature for stage III colorectal cancer patients. Clin Cancer Res 18, 6001-6010, doi:10.1158/1078-0432.CCR-11-3302 (2012).
	#CRCassigner: Sadanandam, A. et al. A colorectal cancer classification system that associates cellular phenotype and responses to therapy. Nat Med 19, 619-625, doi:10.1038/nm.3175 (2013).
	#DeSousa: De Sousa, E. M. F. et al. Poor-prognosis colon cancer is defined by a molecularly distinct subtype and develops from serrated precursor lesions. Nat Med 19, 614-618, doi:10.1038/nm.3174 [pii] (2013).
	#Marisa: Marisa, L. et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. PLoS Med 10, e1001453, doi:10.1371/journal.pmed.1001453 (2013).
	#Merlos: Merlos-Suarez, A. et al. The intestinal stem cell signature identifies colorectal cancer stem cells and predicts disease relapse. Cell stem cell 8, 511-524, doi:10.1016/j.stem.2011.02.020 (2011).
	#MDA114: Oh, S. C. et al. Prognostic gene expression signature associated with two molecularly distinct subtypes of colorectal cancer. Gut 61, 1291-1298, doi:10.1136/gutjnl-2011-300812 (2012).
	#Popovici: Popovici, V. et al. Identification of a poor-prognosis BRAF-mutant-like population of patients with colon cancer. Journal of clinical oncology : official journal of the American Society of Clinical Oncology 30, 1288-1295, doi:10.1200/JCO.2011.39.5814 (2012).
	#Schetter: Schetter, A. J. et al. Association of inflammation-related and microRNA gene expression with cancer-specific mortality of colon adenocarcinoma. Clinical cancer research : an official journal of the American Association for Cancer Research 15, 5878-5887, doi:10.1158/1078-0432.CCR-09-0627 (2009).
	#Yuen3: Yuen, H. F. et al. TAZ expression as a prognostic indicator in colorectal cancer. PloS one 8, e54211, doi:10.1371/journal.pone.0054211 (2013).
	#V7RHS: Jiang, Y. et al. Development of a clinically feasible molecular assay to predict recurrence of stage II colon cancer. The Journal of molecular diagnostics : JMD 10, 346-354, doi:10.2353/jmoldx.2008.080011 (2008).
#input parameters:
	#foler: the working directory
	#method: one of the above mentioned classifier
	#normalized: optional, TRUE, if the data is already RMA normalized
		#rmafile: the normalied file (if the data is already normalized)
	#result_filename: optional, prefix for the result files
#example: 	CRCsubtypes("/samples",method="Yuen3",result_filename="newdataset")
#			CRCsubtypes("/samples",method="Budinska")
#License: GPL-2
CRCsubtypes<-function(folder, method,normalized=FALSE, rmafile=NULL, result_filename="colon"){
setwd(folder)
load("CRC.RData")
#install-packages
# source("https://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("simpleaffy")
# biocLite("hgu133plus2.db")
# biocLite("hgu133a.db")
# biocLite("impute")
# biocLite("DeSousa2013")

require(affy)
require(simpleaffy)
library(survplot)
require(tools)
library(citccmst)
require(MASS)
require(hgu133plus2.db)
require(hgu133a.db)
require(Hmisc)
library(impute)

##normalize data if necessary
if(normalized==T){
	myexprs<-read.table(rmafile, sep="\t", header=T, row.names=1)
	myexprs<-as.matrix(myexprs)
	mode(myexprs)<-"numeric"
	} else {
	celList<-list.celfiles(path=folder,full.names=F)
	myexprs<-justRMA(filenames = celList)
	myexprs <- exprs(myexprs)
}
###Budinska
if (method=="Budinska"){
path<-folder
select.probe=T
e<-myexprs
if(dim(myexprs)[1]>30000){
annotation<-annotation_plus2
} else {
annotation<-annotation_a
}

## 1. preparing EntrezID codes ##
EIDs <- annotation[match(rownames(e), annotation[,1]),]$EntrezID
# cleaning EIDs from // marks
symbol.end <- regexpr("//",as.character(EIDs))-2
symbol.length <- nchar(as.character(EIDs))
symbol.end[which(symbol.end==-3)] <- symbol.length[which(symbol.end==-3)]
names.temp <- substr(as.character(EIDs), 1, symbol.end)

## 2. preparing selection of probesets ##
#if(select.probe)
#{# computing most variable probesets per gene
sdrlm <- apply(e, 1, FUN=function(x) rlm(x~1)$s)
sdrlm_which_max_EID_affy <- tapply(sdrlm, names.temp , FUN=function(x) names(x[which.max(x)]))
#}else{load(paste(path,"/sdrlm_which_max_EID_affy.rdata", sep=""))}
# selecting most variable probesets and renaming rows to Entrez ID names
assign("temp", e[intersect(rownames(e),sdrlm_which_max_EID_affy),])
rownames(temp) <- names(sdrlm_which_max_EID_affy[match(intersect(rownames(e),sdrlm_which_max_EID_affy),sdrlm_which_max_EID_affy)])
assign(paste("e_uEID_",result_filename,sep=""), temp)
save(list=paste("e_uEID_",result_filename,sep=""), file=paste(path,"/e_uEID_",result_filename,".rdata",sep=""))

### 3. calculating metagenes ###
# loading genes for meta-gene construction as defined in Budinska et al., 2013
#load(file=paste(path,"/common.genes.40.rdata", sep=""))
assign(paste("e_uEID_short",result_filename,sep=""), temp[match(common.genes.40,rownames(temp)),])
save(list=paste("e_uEID_short",result_filename,sep=""), file=paste(path,"/e_uEID_short",result_filename,".rdata",sep="")) ###WTF
#loading gene modules
#load(file=paste(path, "/gene.modules.final.rdata", sep=""))
# calculating metagenes
#source(file=paste(path, "/MetaGenesFun.R", sep=""))
assign(paste("metamedian.",result_filename,sep=""),MetaGenesFun(get(paste("e_uEID_",result_filename,sep="")), groups=gene.modules.final[,c(1,2)], "median")) ###WTF
save(list=paste("metamedian.",result_filename,sep=""), file=paste(path,"/metamedian.",result_filename,".rdata",sep="")) ###WTF

# centering
dataset_medians <- apply(get(paste("metamedian.",result_filename,sep="")), 2, FUN=function(x) median(x, na.rm=T))
dataset_centered = c()
for (i in c(1:nrow(get(paste("metamedian.",result_filename,sep="")))))
{
dataset_centered <- cbind(dataset_centered, cbind(get(paste("metamedian.",result_filename,sep=""))[i,] )- cbind(dataset_medians))
}
save(dataset_centered, file=paste(path,"/dataset_centered.rdata",sep=""))

### 4. Applying LDA classifier for subtypes ###
#loading LDA classifier 
#load(file=paste(path,"/lda_res_centered.rdata", sep=""))
assign(paste("ldares.",result_filename,sep=""), predict(lda_res_centered, newdata=t(dataset_centered)))
save(list=paste("ldares.",result_filename,sep=""), file=paste(path,"/ldares.",result_filename,".rdata",sep=""))
subtypes_col = c("red","blue","black","orange","light grey")
subtypes_name = c("A","B","C","D","E")
subtypes_col <- subtypes_col[get(paste("ldares.",result_filename,sep=""))$class]
subtypes <- subtypes_name[get(paste("ldares.",result_filename,sep=""))$class]
names(subtypes)<-colnames(myexprs)
#names(subtypes)<-names(subtypes_col)<-rownames(paste("metamedian.",result_filename,sep=""))
write.table(subtypes,paste(result_filename,"_budinska.txt",sep=""),sep="\t")
}

###CCHS
else if (method=="CCHS"){
mf<-function(x){
a<-1.301+0.543*x["227896_at"]-0.416*x["221478_at"]+0.596*x["207574_s_at"]+0.538*x["209566_at"]-0.177*x["201746_at"]
ifelse(a<4.526,0,1)
}
mf2<-function(x){
a<-1.301+0.543*x["218264_at"]-0.416*x["221478_at"]+0.596*x["207574_s_at"]+0.538*x["209566_at"]-0.177*x["201746_at"]
ifelse(a<4.526,0,1)
}
if(dim(myexprs)[1]>30000){
myresults<-apply(2^myexprs,2,mf)
} else {
myresults<-apply(2^myexprs,2,mf2)}
write.table(myresults,paste(result_filename,"_CCHS.txt",sep=""), sep="\t")
}

###Chang95
else if (method=="Chang95"){
myprobes_chang<-c("210302_s_at","203695_s_at","209094_at","213909_at","203440_at","236335_at","1553808_a_at 226322_at","212531_at","203238_s_at","210512_s_at","204619_s_at","222484_s_at","209153_s_at","200098_s_at","205108_s_at","1555844_s_at 225540_at","202437_s_at","201261_x_at","226998_at","211964_at","229034_at","212446_s_at","221632_s_at","218272_at","213226_at","226908_at","238868_at","239952_at","212610_at","229055_at","205404_at","226492_at","205357_s_at","226930_at","202800_at","211980_at","200795_at","201141_at","223690_at","208760_at","204337_at","213075_at","223237_x_at","228968_at","202791_s_at","212233_at","220301_at","210809_s_at","228850_s_at","203083_at","201666_at","240289_at","213779_at","209747_at","205857_at","209082_s_at","226828_s_at","203382_s_at","203700_s_at","218353_at","205866_at","205695_at","202796_at","1554195_a_at 209875_s_at","205755_at","205422_s_at","223405_at","226691_at","236154_at","205258_at","201426_s_at","226705_at","203118_at","210004_at","205421_at","209687_at","205647_at","212464_s_at","213800_at","219054_at","202699_s_at","214170_x_at","210037_s_at","204416_x_at","202628_s_at","209646_x_at","1553243_at"," 206910_x_at","228224_at","228268_at","240189_at","209336_at")
me2<-myexprs[rownames(myexprs) %in% myprobes_chang,]
mode(me2)<-"numeric"
me2<-2^me2
mysum<-apply(me2,2, function(x) mean(x, na.rm=T))
myresults<-as.numeric(cut2(mysum, g=3)) -1
names(myresults)<-colnames(myexprs)
write.table(myresults,paste(result_filename,"_Chang95.txt",sep=""), sep="\t")
}

###CIN25
else if (method=="CIN25"){
cin70.score <- apply(myexprs[rownames(myexprs) %in% cin70$CIN70.probe[1:25],],2,mean)
dr<-range(cin70.score)
cin70.risk <- as.numeric(cut(x=cin70.score, breaks=c(dr[1]-1,median(cin70.score, na.rm=TRUE), dr[2]+1))) - 1
names(cin70.risk) <- names(cin70.score)
write.table(cin70.risk,paste(result_filename,"cin25risk.txt",sep=""), sep="\t",col.names=NA)
}


###Cologuidepro
else if (method=="Cologuidepro"){
myexprs<-2^myexprs
sveen7e1<-myexprs[rownames(myexprs) %in% sveen7[sveen7$level_associated_with_poor_prognosis=="High",]$probes,]
sveen7e2<-myexprs[rownames(myexprs) %in% sveen7[sveen7$level_associated_with_poor_prognosis=="Low",]$probes,]
sveen7eq1<-apply(sveen7e1,1,function(x) quantile(x,.8))
length(sveen7eq1)
sveen7eq2<-apply(sveen7e2,1,function(x) quantile(x,c(.2)))
sv1<-sveen7e1>sveen7eq1
sv2<-sveen7e2<sveen7eq2
sv3<-rbind(sv1,sv2)*1
mysum<-apply(sv3,2,sum)
mysum[mysum>2]<-"High"
mysum[mysum<3]<-"Low"
write.table(mysum,paste(result_filename,"_Cologuidepro.txt",sep=""), sep="\t", col.names=NA)
}

###Cologuideex
else if (method=="Cologuideex"){
myexprs<-2^myexprs
sveen7e1<-myexprs[rownames(myexprs) %in% agesen13[agesen13$level_associated_with_poor_prognosis=="High",]$probes,]
sveen7e2<-myexprs[rownames(myexprs) %in% agesen13[agesen13$level_associated_with_poor_prognosis=="Low",]$probes,]
sveen7eq1<-apply(sveen7e1,1,function(x) quantile(x,.8))
sveen7eq2<-apply(sveen7e2,1,function(x) quantile(x,c(.2)))
sv1<-sveen7e1>sveen7eq1
sv2<-sveen7e2<sveen7eq2
sv3<-rbind(sv1,sv2)*1
mysum<-apply(sv3,2,sum)
mysum[mysum>4]<-"High"
mysum[mysum<5]<-"Low"
write.table(mysum,paste(result_filename,"_Cologuideex.txt",sep=""), sep="\t", col.names=NA)
}


###CRCassigner
else if (method=="CRCassigner"){
#source("./screenExpr_DWD_normalize.R")
#source("./nmfconsensus_sztup.R")
if(dim(myexprs)[1]>30000){
x<-mget(rownames(myexprs),hgu133plus2SYMBOL)
} else {
x<-mget(rownames(myexprs),hgu133aSYMBOL)
}
x<-unlist(x)
rownames(myexprs) <- x
 m<-match(pamgenes2,rownames(myexprs))
 w<-which(!is.na(m))
 myexprs<-myexprs[m[w],]
filename1<-paste(result_filename,"_tmp.txt", sep="")
write.table(myexprs, filename1, sep="\t", col.names=NA)
screenExpr(filename1,sdCutoffs=0)
 fil1<-paste(result_filename,"_tmp_sd0.txt",sep="")  
 temp1<-read.delim(fil1,sep="\t",as.is=TRUE,row.names=1)
 dt1<-as.matrix(temp1)
 mode(dt1)<-"numeric"
 med<-apply(dt1,1,median)
 dataUnique12 <- t(apply(dt1,1,function(x){x-median(x)}))
 write.table(dataUnique12, paste(result_filename,"_ang_v2_",Sys.Date(),".txt",sep=""), sep="\t", col.names=NA)
nmfconsensus(input.ds=paste(result_filename,"_ang_v2_",Sys.Date(),".txt",sep=""),k.init=5,k.final=5,num.clusterings=20,maxniter=500,error.function="euclidean", non.interactive.run = T,directory=result_filename, myversion="v2",doc.string = "_")
}

###DeSousa
else if (method=="DeSousa"){
library(DeSousa2013)
library(affy)
library(frma)
library(sva)
library(rgl)
data(AMC, package="DeSousa2013")
data(dat, package="DeSousa2013")
data(uniGenes, package="DeSousa2013")
data(diffGenes.f, package="DeSousa2013")
data(classifier, package="DeSousa2013")
data(silh, package="DeSousa2013")
##medianra normalizalas##
 myexprs <- apply( myexprs,1,function(
 x){
    x-median(x)
 })
sdat<-t(myexprs)
datsel <- sdat[names(uniGenes), ]
rownames(datsel) <- uniGenes
datsel <- datsel[diffGenes.f, ]
pamcl <- pamClassify(datsel, signature, pam.rslt, thresh, postRth=1)
sdat.sig <- pamcl[["sdat.sig"]]
pred <- pamcl[["pred"]]
clu.pred <- pamcl[["clu.pred"]]
write.table(clu.pred, row.names=names(pamcl$clu.pred),file=paste(result_filename,"desousa.txt",sep=""), sep="\t")
}

###Merlos-Suarez
else if (method=="Merlos"){
myexprs<-t(apply(myexprs,1,function(x) (x-mean(x))/sd(x)))
#
myprobes<-as.character(merlos$Lgr5.ISC_affy)
myscore <- apply(myexprs[rownames(myexprs) %in% myprobes,],2,mean)
dr<-range(myscore)
mygroup11 <- as.numeric(cut(x=myscore, breaks=3)) - 1
myprobes<-as.character(merlos$EphB2.ISC_affy)
myscore <- apply(myexprs[rownames(myexprs) %in% myprobes,],2,mean)
dr<-range(myscore)
mygroup21 <- as.numeric(cut(x=myscore, breaks=3)) - 1
mysum<-cbind(mygroup11,mygroup21)
rownames(mysum) <- names(myscore)
colnames(mysum) <- c("Lgr5.ISC","EphB2.ISC")
write.table(mysum,paste(result_filename,"_merlos.txt",sep=""), sep="\t", col.names=NA)
}

###Marisa
else if (method=="Marisa") {
cithezannot<-myexprs[,1:2]
cithezannot<-cbind(cithezannot,rownames(cithezannot))
colnames(cithezannot)[3]<-"id"
citresult<-cit.assignCcmst(data=myexprs,
  data.annot=cithezannot,
  data.colId="id",
  data.colMap="id" ,
  citccmst.colMap="Probe.Set.ID",
 dist.method="dqda",
 plot=F
)
print(paste("Analyzing",i,"th sample"))
write.table(citresult,paste(result_filename,"_marisa.txt",sep=""), sep="\t")
}

###Popovici
else if (method=="Popovici") {
popovic<-function(myexprs){
ge1=mean(myexprs[names(myexprs) %in% as.character(pop$gene1_id)])
ge2=mean(myexprs[names(myexprs) %in% as.character(pop$gene2_id)])
BRAFscore<-ge2-ge1
if(ge1<ge2){
prediction<-"BRAFm"
}else {
prediction<-"WT"
}
}
popovic2<-function(myexprs){
ge1=mean(myexprs[names(myexprs) %in% as.character(pop$gene1_id)])
ge2=mean(myexprs[names(myexprs) %in% as.character(pop$gene2_id)])
BRAFscore<-ge2-ge1
if(BRAFscore>1.02){
prediction<-"BRAFm"
}else {
prediction<-"WT"
}
}
BRAFsig<-apply(myexprs,2,popovic)
BRAFscore<-apply(myexprs,2,popovic2)
BRAFsig<-cbind(BRAFsig,BRAFscore)
write.table(BRAFsig,paste(result_filename,"_popovic.txt",sep=""),sep="\t")
}

###Schetter
else if (method=="Schetter") {
myexprs<-2^myexprs
myexprs <- t(apply(myexprs,1,function(x){x-median(x,na.rm=T)}))
myexprs2<-myexprs[rownames(myexprs) %in% mygenes_schetter$affy[1:9],]
schetter_risk<-function(x){
names(x)<-mygenes_schetter$symbol[1:9]
risk1<-(0.855 * x["PRG1"]) + (0.720 * x["IL-10"]) + (0.458 * x["CD68"]) - (0.494 * x["IL-23a"]) - (0.635 * x["IL-12a"])
risk2<-(1.321 * x["PRG1"]) + (0.840 * x["ANXA1"]) + (0.123 * x["IL-23a"]) + (0.484 * x["IL-17a"]) + (0.367 * x["FOXP3"]) - (0.373 * x["HLA-DRA"])
risk<-cbind(risk1,risk2)
}
mr<-apply(myexprs2,2,function(x) schetter_risk(x))
mr.risk<-matrix(NA,2,dim(myexprs)[2])
for (j in 1:2){
dr<-range(as.numeric(mr[j,]))
mr.risk[j,] <- as.numeric(cut(x=as.numeric(mr[j,]), breaks=c(dr[1]-1,median(as.numeric(mr[j,]), na.rm=TRUE), dr[2]+1))) -1
}
mr.sum<-apply(mr.risk,2,sum)
mr.sum[mr.sum==2]<-"high risk"
mr.sum[mr.sum==1]<-"low risk"
mr.sum[mr.sum==0]<-"low risk"
names(mr.sum)<-colnames(myexprs2)
write.table(mr.sum,paste(result_filename,"_schetter.txt",sep=""), sep="\t", col.names=NA)
}

###MDA114
else if (method=="MDA114") {
me2<-myexprs[rownames(myexprs) %in% rownames(mda114),]
szorzo<-mda114[rownames(me2),]$weight
me3<-me2*szorzo
mysum<-apply(me3,2, sum)
if(dim(myexprs)[1]>30000){
cutoff<--125
} else {
cutoff<--65
}
mf<-function(x,cutoff){
if (x>(cutoff)){
myrisk<-"A"
}else{
myrisk<-"B"}
}
myresults<-sapply(mysum, function(x) mf(x,cutoff))
write.table(myresults,paste(result_filename,"_mda114.txt",sep=""), sep="\t")
}


###Yuen
else if (method=="Yuen3"){
me2<-myexprs[rownames(myexprs) %in% myprobes_yuen,]
mode(me2)<-"numeric"
mysum<-apply(me2,1, function(x)
as.numeric(cut(x=x, breaks=c(range(x)[1]-1,median(x), range(x)[2]+1))) -1)
myresults<-apply(mysum,1,sum)
names(myresults)<-colnames(myexprs)
write.table(myresults,paste(result_filename,"_yuen3.txt",sep=""), sep="\t")
}


###_V7RHS
else if (method=="V7RHS"){
jiang_rma<-function(x){
hk<-jianggenes[jianggenes$category=="HK",]$probes
hke<-x[as.character(hk)]
nemhk<-jianggenes[jianggenes$category=="nonHK",]$probes
nemhke<-x[as.character(nemhk)]
names(nemhke)<-jianggenes$Gene_symbol[1:7]
nemhk2<-6.64-(mean(hke)-nemhke)
RHS1 = -3.251*nemhk2["LILRB3"] - 3.156*nemhk2["YWHAH"] - 3.035*nemhk2["CHC1"] + 3.002*nemhk2["KLF5"] - 2.842*nemhk2["CAPG"] - 3.249*nemhk2["LAT"] - 2.835*nemhk2["EPM2A"]
nemhk2<-nemhk2*(-1)
RHS2 = -3.251*nemhk2["LILRB3"] - 3.156*nemhk2["YWHAH"] - 3.035*nemhk2["CHC1"] + 3.002*nemhk2["KLF5"] - 2.842*nemhk2["CAPG"] - 3.249*nemhk2["LAT"] - 2.835*nemhk2["EPM2A"]
nemhk2<-nemhke-mean(hke)
RHS3 = -3.251*nemhk2["LILRB3"] - 3.156*nemhk2["YWHAH"] - 3.035*nemhk2["CHC1"] + 3.002*nemhk2["KLF5"] - 2.842*nemhk2["CAPG"] - 3.249*nemhk2["LAT"] - 2.835*nemhk2["EPM2A"]
if(RHS1<0){
risk1<-"low"
} else {
risk1<-"high"
}
if(RHS2<0){
risk2<-"low"
} else {
risk2<-"high"
}
if(RHS3>0){
risk3<-"low"
} else {
risk3<-"high"
}
cbind(RHS1,RHS2,RHS3,risk1,risk2,risk3)
}
myexprs <- t(apply(myexprs,1,function(x){x-median(x,na.rm=T)}))
mr2<-apply(myexprs,2,jiang_rma)
mr.risk2<-matrix(NA,3,dim(myexprs)[2])
for (j in 1:3){
dr<-range(as.numeric(mr2[j,]))
mr.risk2[j,] <- as.numeric(cut(x=as.numeric(mr2[j,]), breaks=c(dr[1]-1,median(as.numeric(mr2[j,]), na.rm=TRUE), dr[2]+1))) -1
}
mr3<-rbind(mr2[6,],mr.risk2)
mr3<-rbind(mr2[6,],mr.risk2[1,])
mr3<-t(mr3)
#colnames(mr3)<-c("jiang1","jiang2", "jiang3","jiang4")
colnames(mr3)<-c("V7RHS_v1","V7RHS_v2")
mr3<-as.data.frame(mr3,stringsAsFactors=FALSE)
mr3$jiang2[mr3$jiang2==1]<-"low"
mr3$jiang2[mr3$jiang2==0]<-"high"
mr3$jiang3[mr3$jiang3==1]<-"high"
mr3$jiang3[mr3$jiang3==0]<-"low"
write.table(mr3,paste(result_filename,"_V7RHS.txt",sep=""),sep="\t", col.names=NA)
}
else {print("Method not found")}

}