system("cp SDF2PDF/* .")
#######################
system("rm -r *.temp")
system("rm fmcsR_cluster*")
system("rm *.smi")
#######################
#pdf(file="Test.pdf")
#######################
library(ChemmineR)#####
library(fmcsR)#########
library(plotrix)#######
#######################
#data(sdfsample)
#sdfset<-sdfsample
sdfset<-read.SDFset("all.sdf")
apset<-sdf2ap(sdfset)
#
sdfset<-sdfset[!sapply(as(apset,"list"),length)==1]
#
cid(sdfset)<-datablocktag(sdfset, tag="DRUGBANK_ID")
cluster<-cmp.cluster(apset, cutoff=c(0.7)) ##,0.4,0.6,0.7,0.8))
######################
#cid(sdfset)<-sdfid(sdfset)
###cid(sdfset)<-substring(gsub(" ","_",sdfid(sdfset)), 1, 20)
#cid(sdfset)<-gsub("\\=", "_", cid(sdfset))
#cid(sdfset)<-gsub("\\/", "_", cid(sdfset))
#cid(sdfset)<-gsub("\\?", "_", cid(sdfset))
#cid(sdfset)<-gsub(" ","_",cid(sdfset))
######################
#################
Work<-sort(unique(cluster$CLID_0.7))
##cluster$ids<-seq_along(cluster$ids)
#################
save.image("fmcsR_space.rda", compress=T)
#######################
DoTheWork<-function(f){
#######################
#load("fmcsR_space.rda")
#######################
#data(sdfsample)
#sdfset<-sdfsample
#
Numba<-length(sdfset[sort(cid(sdfset))%in%sort(subset(cluster,cluster$CLID_0.7==Work[f])$ids)])
if(Numba>2){
###
sdfset<-sdfset[sort(cid(sdfset))%in%sort(subset(cluster,cluster$CLID_0.7==Work[f])$ids)]
smiset<-sdf2smiles(sdfset)
######################
#fmcsR the sdfset
##################
Name<-paste("fmcsR_cluster_", Work[f], ".png", sep="") 
SDF<-paste("fmcsR_cluster_", Work[f], ".sdf", sep="") 
png(file=Name, width=2400, height=2400, units="px") 
table_name <-paste("fmcsR_cluster_", Work[f], "_SAR.txt", sep="") 
neighbors<-paste("fmcsR_cluster_", Work[f], ".rda", sep="")
html_name <- paste("fmcsR_cluster_", Work[f], "_SAR.html", sep="")
##################
write.SDF(sdfset, file=SDF, sig=T)
#sdfset<-read.SDFset(SDF)
##################
write.SMI(smiset, file=paste("fmcsR_cluster_", Work[f], ".sdf.smiset", sep=""), cid=T)
##################
WriteSmiOut<-function(a){
##################
write.SMI(smi=smiset[a], file=paste(cid(sdfset[a]), ".smi", sep=""), cid=T)
}
##################
a<-1:length(smiset) 
lapply(a,WriteSmiOut)
##################
#echo latex
######
system("for i in fmcsR_cluster_*.sdf
do
######
mkdir $i.temp
cp $i $i.temp/
cp $i.smiset $i.temp/
cp SDF2PDF/* $i.temp/
cp *.sh $i.temp/
cd $i.temp
babel $i $i.smi -m
cd ..
done")
######
system("for y in *.temp
do
cd $y
for s in *.smi
do
mol2chemfig -wo $s> $s.tex
done
./MakeLatex.sh
cd  ..
done")
######
##################
#echo homepage
##################
cid(sdfset)<-substring(sdfid(sdfset), 1, 50)
##################
library(fmcsR)####
##################
d <- sapply(cid(sdfset), function(x)
fmcsBatch(sdfset[x], sdfset, au=0, bu=0,
matching.mode="aromatic")[,"Overlap_Coefficient"])
##################
hc <- hclust(as.dist(1-d), method="complete")
hc[["labels"]] <- cid(sdfset)
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE, main="hclust of fmcsR tanimoto distances")
dev.off()
#################
#################
write.table(round(d[order(rev(hc$order)),]*100), file=table_name, sep=" & ", quote=FALSE)
Fix<-round(d[order(rev(hc$order)),]*100)
save(Fix, file=neighbors, compress=TRUE)
#################
#################
pngB<-paste("fmcsR_cluster_", Work[f], "_heatmap.png", sep="")
png(file=pngB, width=1200, height=1200, units="px") 
#################
library(plotrix)#
#################
#load("nnm.rda")
##
a<-Fix
##
n<-length(colnames(a))
##
#
color2D.matplot(a,c(0,1),c(0,0),c(0,0),show.legend=TRUE,show.values=FALSE,axes=FALSE)
#
axis(3,at=seq(1, n, 1)-0.50,labels=substring(colnames(d[order(rev(hc$order)),]), 1, 3))
axis(2,at=seq(1, n, 1)-0.50,labels=substring(rev(rownames(d[order(rev(hc$order)),])), 1, 3))
#
box()
##
dev.off()
##################
library(hwriter)##
##################
p=openPage(html_name)
##
#load("nnm.rda")
##
colnames(Fix)<-1:length(colnames(Fix))
##
#colors<=c('#000000','#FFFFFF')
##
hwrite(Fix, p, br=TRUE)
##
closePage(p)
##################

#################
#All DONE
#################
}
else
{}
#########################################
}
f<-1:length(Work)
lapply(f, DoTheWork)


#dev.off()
