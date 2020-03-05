wget https://www.cog-genomics.org/static/bin/plink/glist-hg19 -O glist-hg19
wget https://www.cog-genomics.org/static/bin/plink/glist-hg38
perl -p -i -e 's/ /\t/g' glist-hg19
perl -p -i -e 's/ /\t/g' glist-hg38

grep '\bIL4\b' glist-hg19 > IL4.fly

plink --bfile IL4 --pheno IL4.phen --mpheno 2 --maf 0.01 --assoc --allow-no-sex --make-set glist-hg19 --set-names "IL4" --make-set-border 5 --adjust

library("CMplot")
source("http://raw.githubusercontent.com/Shicheng-Guo/RobustSKAT/master/R/qqplotsource.R")
input<-read.table("plink.qassoc",head=T,as.is=T,check.names=F)
memo="IL4"
Chromosome<-input$CHR
Position<-input$BP
SNP<-input$SNP
P<-input$P
cminput<-data.frame(SNP,Chromosome,Position,P)
colnames(cminput)=c("SNP","Chromosome","Position","trait1")
CMplot(cminput,plot.type="b",memo="IL4",LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
write.csv(cminput,file="IL4.csv",quote=F,row.name=F,col.names=T)

awk '$1>0{print $1,$3-1,$3,$2,$5,$9}' OFS="\t" plink.qassoc | grep -v NA > plink.qassoc.hg19.bed

bedtools intersect -wao -a plink.qassoc.hg19.bed -b glist-hg19 | grep -v 

wget https://www.cog-genomics.org/static/bin/plink/glist-hg19 -O glist-hg19
wget https://www.cog-genomics.org/static/bin/plink/glist-hg38 -O glist-hg38
awk '{print $1,$2-5000,$3+5000,$4}' OFS="\t" glist-hg19 > glist-hg19.5k.bed
awk '{print $1,$2-5000,$3+5000,$4}' OFS="\t" glist-hg38 > glist-hg38.5k.bed
bedtools intersect -wao -a plink.qassoc.hg19.bed -b glist-hg19.5k.bed | grep -v '\-1' | awk '{print $1,$2,$3,$4,$5,$6,$10}' OFS="," > IL4.pvalue.csv