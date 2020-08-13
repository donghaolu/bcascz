##################################################
### PRS analysis using PRSice V1.25 in terminal
##################################################

R -q --file=./PRSice_v1.25.XC.R --args \
plink ./plink_1.9_linux_160914 \
base /nfs/home/jieson/cross_trait_scz_bcar/PRSice/scz_49EUR_MAF1_INFO6_rmindel_PRSice_input \
target /nfs/home/jieson/cross_trait_scz_bcar/PRSice/bcac_EUR_MAF1_INFO6_rmindel_PRSice_input \
slower 0 supper 0.5 sinc 0.0001 \
barchart.levels 0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1,0.2,0.3,0.4,0.5 \
clump.snps T \
clump.ref /nfs/Stanley/GRS/files/1kg_p1v3_PLINK_cleaned/1kgEUR.noATCG.nomhc \
clump.p1 1 clump.p2 1 clump.kb 500 clump.r2 0.1 \
covary F quantiles T \
sumsum T size.targ 228951 \
figname SCZ_nomhc_BCarisk_500kb_test cleanup F 


R -q --file=./PRSice_v1.25.XC.R --args \
plink ./plink_1.9_linux_160914 \
base /nfs/home/jieson/cross_trait_scz_bcar/PRSice/bcac_EUR_MAF1_INFO6_rmindel_PRSice_input \
target /nfs/home/jieson/cross_trait_scz_bcar/PRSice/scz_49EUR_MAF1_INFO6_rmindel_PRSice_input \
slower 0 supper 0.5 sinc 0.0001 \
barchart.levels 0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1,0.2,0.3,0.4,0.5 \
clump.snps T \
clump.ref /nfs/Stanley/GRS/files/1kg_p1v3_PLINK_cleaned/1kgEUR.noATCG.nomhc \
clump.p1 1 clump.p2 1 clump.kb 500 clump.r2 0.1 \
covary F quantiles T \
sumsum T size.targ 77096 \
figname BCa_nomhc_SCZrisk_500kb cleanup T 

### output files "SCZ_nomhc_BCarisk_500kb_RAW_RESULTS_DATA.txt", "BCa_nomhc_SCZrisk_500kb_RAW_RESULTS_DATA.txt"

#################################
### DRAW PRS plot in R 
#################################

library(data.table)
library(foreign)
library(ggplot2)
require(gridExtra)
library(RColorBrewer)
colourCount = length(unique(R2.scz_prs_bca$r2))
getPalette = colorRampPalette(brewer.pal(9, "YlOrRd"))

### Import data in R
scz_prs_bca <- read.table("SCZ_nomhc_BCarisk_500kb_RAW_RESULTS_DATA.txt",header=T,stringsAsFactors=F)
bca_prs_scz <- read.table("BCa_nomhc_SCZrisk_500kb_RAW_RESULTS_DATA.txt",header=T,stringsAsFactors=F)

mytheme <-theme(legend.direction = "horizontal",
                legend.position="bottom",
                panel.grid.major.x = element_blank(), 
                axis.text.x=element_text(angle=25, hjust=1, size=12, face="bold"),
                axis.text.y=element_text(size=12, face="bold"), 
                axis.title=element_text(size=17, face="bold"), 
                strip.text = element_text(size=15, face="bold"),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12, face="bold"))

cthre <- c(5e-08, 5e-07, 5e-06, 5e-05, 1e-04, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5)
R2.scz_prs_bca <- scz_prs_bca[scz_prs_bca$thresh %in% cthre,]
R2.bca_prs_scz <- bca_prs_scz[bca_prs_scz$thresh %in% cthre,]

R2.scz_prs_bca$Variable <- c('scz_5e8', 'scz_5e7', 'scz_5e6', 'scz_5e5', 'scz_1e4', 'scz_001', 'scz_01', 'scz_05', 
				'scz_0_1','scz_0_2', 'scz_0_5')
R2.scz_prs_bca$SCZ_pT <- c("p<5e-08","p<5e-07","p<5e-06","p<5e-05","p<1e-04","p<0.001","p<0.01", 
							"p<0.05", "p<0.1", "p<0.2", "p<0.5")                            
# re-order the levels in the order of appearance in the data.frame
R2.scz_prs_bca$SCZ_pT1 <- factor(R2.scz_prs_bca$SCZ_pT, as.character(R2.scz_prs_bca$SCZ_pT))
R2.scz_prs_bca$P.value <- formatC(R2.scz_prs_bca$pval, format = "e", digits = 2)

R2.bca_prs_scz$Bca_pT <- c("p<5e-08","p<5e-07","p<5e-06","p<5e-05","p<1e-04","p<0.001","p<0.01", 
							"p<0.05", "p<0.1", "p<0.2", "p<0.5")
# re-order the levels in the order of appearance in the data.frame
R2.bca_prs_scz$Bca_pT <- factor(R2.scz_prs_bca$SCZ_pT, as.character(R2.scz_prs_bca$SCZ_pT))
R2.bca_prs_scz$P.value <- formatC(R2.bca_prs_scz$pval, format = "e", digits = 2)

### draw the plot
plot.R2.scz_prs_bca <- ggplot(R2.scz_prs_bca, aes(y=r2, x=SCZ_pT1)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_manual(values = getPalette(colourCount)) +
  mytheme + xlab("P-value threshold from schizophrenia GWAS") + ylim(0, max(R2.bca_prs_scz$r2)+0.0005) +
  geom_text(aes(label=P.value), vjust = -1.0, hjust = 0, angle = 45, cex = 3.8, parse=T, position=position_dodge(width = 0.9)) +
  ylab("Breast cancer variance explained")  
  
plot.R2.bca_prs_scz <- ggplot(R2.bca_prs_scz, aes(y=r2, x=Bca_pT)) +
  geom_bar(stat="identity",position="dodge") +
  #scale_fill_manual(values = getPalette(colourCount)) +
  mytheme + xlab("P-value threshold from breast cancer GWAS") + ylim(0, max(R2.bca_prs_scz$r2)+0.0005) +
  geom_text(aes(label=P.value), vjust = -1.0, hjust = 0, angle = 45, cex = 3.8, parse=T, position=position_dodge(width = 0.9)) +
  ylab("Schizophrenia variance explained")
   
### Save plots
pdf("Combined_PRS_rmmhc_500kb.pdf",width=20,height=10)
grid.arrange(plot.R2.scz_prs_bca, plot.R2.bca_prs_scz, ncol=2)
dev.off()  

png("Combined_PRS_rmmhc_500kb.png",width=1600,height=750,res=85)
grid.arrange(plot.R2.scz_prs_bca, plot.R2.bca_prs_scz, ncol=2)
dev.off()  

