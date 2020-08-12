##
## FIGURE 2 PLOT TOP HITS OF SCZ ON BCA AND VICE VERSA
##


## Created by Donghao Lu, donghao.lu@ki.se



R


rm(list=ls())
gc()


setwd("/nfs/home/donlu/SczBca/")
options(stringsAsFactors = F)


# read GWAS results
dat <- read.table("cross_scz_bcac_GWAS_MAF1_INFO6_rmindel_original.txt.gz", header=T, sep="\t")
dat$BETA <- log(dat$OR)

# output files for PRS
pgc <- dat[, c("CHR","BP","SNP","A1","A2","BETA","P","SE")]
bcac <- dat[, c("CHR","BP","SNP","a1","a0","bcac_onco_icogs_gwas_beta","bcac_onco_icogs_gwas_P1df","bcac_onco_icogs_gwas_se")]
names(bcac) <- c("CHR","BP","SNP","A1","A2","BETA","P","SE")
bcac<- bcac[order(bcac$P),]
pgc <- pgc[order(pgc$P),]
write.table(pgc, "scz_49EUR_MAF1_INFO6_rmindel_PRSice_input", quote=F, row.names=F,col.names=T,sep="\t")
write.table(bcac, "bcac_EUR_MAF1_INFO6_rmindel_PRSice_input", quote=F, row.names=F,col.names=T,sep="\t")


# read top hits from Scz and BCa
scz <- read.table("SczEur95GenomicRiskLoci.txt",header=T,sep="\t")
#95 snps

bca <- read.table("BcaEur149GenomicRiskLoci.txt", header=T, sep="\t")
#149 snps


# keep top hits: 244 SNPs
snp <- unique(rbind(bca[,3:5], scz[,3:5]))
#244 snps
names(snp) <- c("rsID","CHR","BP")
pgc <- merge(pgc, snp[,2:3], by=c("CHR","BP"))
bcac <- merge(bcac, snp[,2:3], by=c("CHR","BP"))

# merge
dat <- merge(pgc, bcac, by=c("SNP","CHR","BP"))
names(dat) <- c("SNP","CHR","BP","scz.A1","scz.A2","scz.beta","scz.p","scz.se", "bca.A1","bca.A2","bca.beta","bca.p","bca.se")
#summary(dat$scz.A1==dat$bca.A1)

# flip beta
dat$bca.beta <- ifelse(dat$scz.A1==dat$bca.A1, dat$bca.beta, 0-dat$bca.beta)


# read top hits from Scz and BCa
scz <- read.table("SczEur95GenomicRiskLoci.txt",header=T,sep="\t")
#95 snps
bca <- read.table("BcaEur149GenomicRiskLoci.txt", header=T, sep="\t")
#149 snps

# exclude mhc regeion : chr6: 25,000,000-34,000,000 (hg 19)
ex <- dat$CHR==6 & (25000000<=dat$BP | dat$BP<=34000000)
#13
dat <- dat[!ex,]
#231 obs

#Scatter plot
library(gtx)
require(ggplot2)
library(gridExtra)
library(ggrepel)

#plot Scz to BCa
datt <- merge(dat, scz[,4:5], by.x=c("CHR","BP"), by.y=c("chr","pos"))
#91 obs

#FDR adjustment
datt$padj <- p.adjust(datt$bca.p, "fdr")
#summary(datt$padj<0.05)
#6 snps

#flip beta
datt$bca.beta <- ifelse(datt$scz.beta<0, 0-datt$bca.beta, datt$bca.beta)
datt$scz.beta <- ifelse(datt$scz.beta<0, 0-datt$scz.beta, datt$scz.beta)

#UB and LB
datt$bca.lb <- datt$bca.beta - 1.96*datt$bca.se
datt$bca.ub <- datt$bca.beta + 1.96*datt$bca.se

#overall effect estimate
tab <- grs.summary(datt$scz.beta, datt$bca.beta, datt$bca.se, 228951)
tab <- as.data.frame(tab)
tab <- tab[,c("m","ahat","aSE","pval","phet")]
tab$ahat <- round(tab$ahat,3)
tab$aSE <- round(tab$aSE,3)
tab$pval <- signif(tab$pval,3)
tab$phet <- signif(tab$phet,3)
names(tab) <- c("N of SNPs","Beta","SE","P value","Phet")

#code display
fig1 <- datt[,c("SNP","scz.beta","bca.beta","bca.ub","bca.lb","padj")]

fig1$SNP[fig1$padj>=0.05] <- NA
fig1$group <- ifelse(fig1$padj<0.05,1,0)

#plot ORs
fig1 <- fig1[order(fig1$group),]
fig<-ggplot(fig1, aes(x=scz.beta, y=bca.beta))+
  geom_point(data=subset(fig1, group==0), alpha=0.5, size=1, shape=1)+
  geom_errorbar(data=subset(fig1, group==0), aes(ymax=bca.ub,ymin=bca.lb),width=.001,alpha=0.2)+
  geom_abline(slope=tab[,2],col="black",linetype="dashed")+  
  scale_x_continuous(name="Effect on Schizophrenia (logOR)",breaks=seq(0,0.2,.02))+
  scale_y_continuous(name="Effect on Breast Cancer (logOR)",breaks=seq(-.04,.06,.02))+
  geom_point(data=subset(fig1, group==1), color="orange", size=2, shape=16) +
  geom_errorbar(data=subset(fig1, group==1), aes(ymax=bca.ub,ymin=bca.lb),width=.002,color="orange")+
  geom_label_repel(data=subset(fig1, group==1), aes(label=SNP), size=4, fill=NA, label.padding=.1) +
  theme_classic()+
  theme(axis.line.x=element_line(colour="grey50"),axis.line.y=element_line(colour="grey50"),
        legend.position = "none")
fig


#plot table
tb1 <- tableGrob(tab, rows="", theme=ttheme_minimal(base_size=12,padding = unit(c(10, 2), "mm")))
figtabA <- grid.arrange(fig, tb1,nrow=2, as.table=T, heights=c(5,1), top="A. Lead SNPs for Schizophrenia")
figtabA






#plot BCa to Scz
datt <- merge(dat, bca[,4:5], by.x=c("CHR","BP"), by.y=c("chr","pos"))
#140 obs

#FDR adjustment
datt$padj <- p.adjust(datt$scz.p, "fdr")
#summary(datt$padj<0.05)
#9 snps


#flip beta
datt$scz.beta <- ifelse(datt$bca.beta<0, 0-datt$scz.beta, datt$scz.beta)
datt$bca.beta <- ifelse(datt$bca.beta<0, 0-datt$bca.beta, datt$bca.beta)

#UB and LB
datt$scz.lb <- datt$scz.beta - 1.96*datt$scz.se
datt$scz.ub <- datt$scz.beta + 1.96*datt$scz.se

#overall effect estimate
tab <- grs.summary(datt$bca.beta, datt$scz.beta, datt$scz.se, 77096)
tab <- as.data.frame(tab)
tab <- tab[,c("m","ahat","aSE","pval","phet")]
tab$ahat <- round(tab$ahat,3)
tab$aSE <- round(tab$aSE,3)
tab$pval <- signif(tab$pval,3)
tab$phet <- signif(tab$phet,3)
names(tab) <- c("N of SNPs","Beta","SE","P value","Phet")

#code display
fig1 <- datt[,c("SNP","bca.beta","scz.beta","scz.ub","scz.lb","padj")]

fig1$SNP[fig1$padj>=0.05] <- NA
fig1$group <- ifelse(fig1$padj<0.05,1,0)


#plot ORs
fig1 <- fig1[order(fig1$group),]
fig<-ggplot(fig1, aes(x=bca.beta, y=scz.beta))+
  geom_point(data=subset(fig1, group==0), alpha=0.5, size=1, shape=1)+
  geom_errorbar(data=subset(fig1, group==0), aes(ymax=scz.ub,ymin=scz.lb),width=.001,alpha=0.2)+
  geom_abline(slope=tab[,2],col="black",linetype="dashed")+  
  scale_x_continuous(name="Effect on Breast Cancer (logOR)",breaks=seq(0,.3,.05))+
  scale_y_continuous(name="Effect on Schizophrenia (logOR)",breaks=seq(-.1,.1,.05))+
  geom_point(data=subset(fig1, group==1), color="orange", size=2, shape=16) +
  geom_errorbar(data=subset(fig1, group==1), aes(ymax=scz.ub,ymin=scz.lb),width=.002,color="orange")+
  geom_label_repel(data=subset(fig1, group==1), aes(label=SNP), size=4, fill=NA, label.padding=.1) +
  theme_classic()+
  theme(axis.line.x=element_line(colour="grey50"),axis.line.y=element_line(colour="grey50"),
        legend.position = "none")
fig


#plot table
tb1 <- tableGrob(tab, rows="", theme=ttheme_minimal(base_size=12,padding = unit(c(10, 2), "mm")))
figtabB <- grid.arrange(fig, tb1,nrow=2, as.table=T, heights=c(5,1), top="B. Lead SNPs for Breast Cancer")
figtabB


#panel figure
pfig <- grid.arrange(figtabA, figtabB, nrow=1, ncol=2)
ggsave('../output/pub_figure2.pdf',pfig, width=45, height=20, units="cm",scale=0.65) 




