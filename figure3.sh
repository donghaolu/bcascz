##
## FIGURE 3 REGIONAL ASSOCIATION PLOT
##

## Created by Donghao Lu, donghao.lu@ki.se



## Prep input files for LocusZoom

R

rm(list=ls())
gc()


options(stringsAsFactors = F)


# read GWAS results
pgc <- read.table("scz_49EUR_MAF1_INFO6_rmindel_PRSice_input", header=T, sep="\t")
bcac <- read.table("bcac_EUR_MAF1_INFO6_rmindel_PRSice_input", header=T, sep="\t")


# merge
scz_bc <- merge(pgc, bcac, by=c("SNP","CHR","BP"))
names(scz_bc) <- c("SNP","CHR","BP","A1.scz","A2.scz","beta.scz","P.scz","SE.scz", "A1.bcac","A2.bcac","beta.bcac","P.bcac","SE.bcac")

# flip beta
scz_bc$beta.bcac <- ifelse(scz_bc$A1.scz==scz_bc$A1.bcac, scz_bc$beta.bcac, 0-scz_bc$beta.bcac)
scz_bc <- scz_bc[,-(9:10)]

# read top hits from Scz and BCa
scz <- read.table("SczEur95GenomicRiskLoci.txt",header=T,sep="\t")
#95 snps
bca <- read.table("BcaEur149GenomicRiskLoci.txt", header=T, sep="\t")
#149 snps


#Scz to BCa
datt <- merge(scz_bc, scz[,4:5], by.x=c("CHR","BP"), by.y=c("chr","pos"))
#95 obs

# Sign. after FDR adjustment
datt$padj <- p.adjust(datt$P.bcac, "fdr")
datt <- subset(datt, padj<0.05)
#7 obs


# output ref SNP +/- 1MB for Locuszoom
for (i in 1:nrow(datt)) {

  # position
  ch <- datt$CHR[i]
  ub <- datt$BP[i] + 1e6
  lb <- datt$BP[i] - 1e6
  dat <- subset(scz_bc, CHR==ch & BP>=lb & BP<=ub)
  # output scz
  write.table(dat[,c("SNP","P.scz")], paste("Scz2BCaChr",ch,datt$SNP[i],"Locus1mbBP",datt$BP[i],"scz.txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")
  # output bca
  write.table(dat[,c("SNP","P.bcac")], paste("Scz2BCaChr",ch,datt$SNP[i],"Locus1mbBP",datt$BP[i],"bca.txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")

}


#BCa to Scz
datt <- merge(scz_bc, bca[,4:5], by.x=c("CHR","BP"), by.y=c("chr","pos"))
#149 obs

# Sign. after FDR adjustment
datt$padj <- p.adjust(datt$P.scz, "fdr")
datt <- subset(datt, padj<0.05)
#10 obs


# output ref SNP +/- 1MB for Locuszoom
for (i in 1:nrow(datt)) {

  # position
  ch <- datt$CHR[i]
  ub <- datt$BP[i] + 1e6
  lb <- datt$BP[i] - 1e6
  dat <- subset(scz_bc, CHR==ch & BP>=lb & BP<=ub)
  # output scz
  write.table(dat[,c("SNP","P.scz")], paste("BCa2SczChr",ch,datt$SNP[i],"Locus1mbBP",datt$BP[i],"scz.txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")
  # output bca
  write.table(dat[,c("SNP","P.bcac")], paste("BCa2SczChr",ch,datt$SNP[i],"Locus1mbBP",datt$BP[i],"bca.txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")

}


q()
n




## LocusZoom Standalone on vector

# download package from web
wget https://statgen.sph.umich.edu/locuszoom/download/locuszoom_1.4.tgz
tar xvf locuszoom_1.4.tgz
#89G


# update database
locuszoom/bin/lzupdate.py --build hg19 --gencode 19 --gwas-cat


# plot locus chr19(p13.11)
# prep annote file
{
  echo -e 'snp\tstring\tcolor'
  echo -e 'rs2965183\tBCa\tblack'
  echo -e 'rs2905426\tScz\tblack'
} >chr19p13annot.txt

# BCa
# conditioning
./locuszoom --metal BCa2SczChr19rs2965183Locus1mbBP19545696bca.txt --markercol "SNP" --pvalcol "P.bcac" --refsnp "rs2965183" --add-refsnps "rs2905426" --flank 500KB --pop EUR --build hg19 --source 1000G_March2012 --plotonly --prefix "LocusChr19BCa" --denote-markers-file chr19p13annot.txt refsnpTextSize=.7 ymax=15 width=20 height=10 condLdColors="gray60,#FF00FF,#AFEEEE" showGenes=FALSE


# Scz
# conditioning
./locuszoom --metal BCa2SczChr19rs2965183Locus1mbBP19545696scz.txt --markercol "SNP" --pvalcol "P.scz" --refsnp "rs2965183" --add-refsnps "rs2905426" --flank 500KB --pop EUR --build hg19 --source 1000G_March2012 --plotonly --prefix "LocusChr19Scz" --denote-markers-file chr19p13annot.txt refsnpTextSize=.7 ymax=15 width=20 height=10 condLdColors="gray60,#FF00FF,#AFEEEE" showGenes=FALSE



##
## END OF FILE
##







