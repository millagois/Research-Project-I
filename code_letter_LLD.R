####### Script used for generating the plots for the LLD letter 
#### Date: 04/12/2020 

#get ppackages
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggsignif)
library(mediation)
#files to use: 
LLD_gen <- read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-mgois/data_ready/LLD/LLD_genotype_final.txt", sep="\t", header=T)
LLD_dairy <- read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-mgois/data_ready/LLD/LLD_dairy.txt", sep="\t", header=T)
LLD_lac <- read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-mgois/data_ready/LLD/LLD_lactose_total.txt", sep="\t", header=T)
LLD_gut <- read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-mgois/data_ready/LLD/LLD_gut_complaints.txt", sep="\t", header=T)
LLD_bifido2 <- read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-mgois/data_ready/LLD/LLD_bifido_RD_corrected.txt", sep="\t", header = T)
LLD_intrinsic_factors<-  read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-mgois/data_ready/LLD/LLD_intrinsic_factors.txt", sep="\t", header = T)

#merge all dfs into one big one
LLD_df<- merge(LLD_gen, LLD_bifido2, by = 'LLDEEPID')
LLD_df <- merge(LLD_df, LLD_dairy, by= "LLDEEPID")
LLD_df <- merge(LLD_df, LLD_lac, by ="LLDEEPID")
LLD_df <- merge(LLD_df, LLD_gut, by= "LLDEEPID")
LLD_df <- merge(LLD_df, LLD_intrinsic_factors, by= "LLDEEPID")

names(LLD_df)[names(LLD_df) == "gdag"] <- "dairy"
names(LLD_df)[names(LLD_df)== "belching"] <- "burping"

#take out outliers that consume more than 1L of dairy products/day
LLD_df2 <- LLD_df[!(LLD_df$dairy >1000),]

#cleaning out the df
FDF <- LLD_df
FDF<- FDF[,c("LLDEEPID", "rs4988235_number", "rs4988235_var", "Read_Count",
              "Age", "Sex", "Bifido_genus_abundance", "dairy", "lactose", "pain",
             "discomfort","constipation", "diarrhea", "bloating", "flatulence", "burping",
             "combined", "nausea" , "antrop_BMI")] 

###transforming the data:             
invrank = function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))} 

FDF[,c("Bifido_genus_abundance", "dairy", "lactose", "combined", 
        "discomfort", "pain", "constipation","diarrhea","bloating",
        "flatulence", "burping", "nausea")] <- lapply(FDF[,c("Bifido_genus_abundance", "dairy", "lactose", "combined", 
                                                              "discomfort", "pain", "constipation","diarrhea",
                                                              "bloating", "flatulence", "burping", "nausea")], invrank)
###separating per genotype: 
FDF2 <- FDF

FDF2$LI_GG <- FDF2$rs4988235_number ==1
DF_LI <- FDF2[!(FDF2$LI_GG=='FALSE'),]
DF_LI <- DF_LI[!(is.na(DF_LI$combined)),]

DF_NOT_LI <- FDF2[!(FDF2$LI_GG=='TRUE'),]
DF_NOT_LI<-DF_NOT_LI[complete.cases(DF_NOT_LI), ]


FDF2$LNI_GT <- FDF2$rs4988235_number ==2
DF_LNI_GT <- FDF2[!(FDF2$LNI_GT=='FALSE'),]
DF_LNI_GT<-DF_LNI_GT[complete.cases(DF_LNI_GT), ]

FDF2$LT_TT <- FDF2$rs4988235_number ==3
DF_LT_TT <- FDF2[!(FDF2$LT_TT=='FALSE'),]
DF_LT_TT<-DF_LT_TT[complete.cases(DF_LT_TT), ]

##SNP-bifido (whole population)
wilcox.test(FDF2$Bifido_genus_abundance ~ FDF2$LI_GG)

kruskal.test(FDF2$Bifido_genus_abundance ~ FDF2$rs4988235_number)

png("letter_LLD_Bifido.png")
boxplot (Bifido_genus_abundance ~ rs4988235_var, data = FDF, col = c("#CC0000", "#333333", "#333333"), 
         main="LLD",
         xlab = "Genetic Variants (SNP rs4988235)", 
         ylab = "Bifidobacterium genus abundance") 

dev.off()


##LI individuals
###Dairy intake - Gut complaints 
png("LLD_LI_dairy_GC_combined.png")
scatter.smooth(DF_LI$dairy, DF_LI$combined, ylab= "Gastrointestinal Symptoms (score)", xlab = "Dairy Intake (g/day)", main= "LLD - LI Recessive Genotype")
dev.off()

cor.test(DF_LI$dairy, DF_LI$combined)
summary(lm(DF_LI$dairy ~ DF_LI$combined+DF_LI$antrop_BMI))
##Dairy intake - all gut complaints: 
cor_dairy_LI <- cor(DF_LI$dairy, DF_LI[,c("pain", "discomfort", "constipation", "diarrhea", 
                                          "bloating", "flatulence", "burping", "nausea", "combined")], use = 'pairwise.complete.obs')

pval_dairy_LI <- apply(DF_LI[,c("pain", "discomfort", "constipation", "diarrhea", 
                                "bloating", "flatulence", "burping", "nausea", "combined")], 2, function(x){ 
                                        cor.test(DF_LI$dairy, x, method =  "pearson", use = "pairwise")
                                })
p.vals <- sapply( pval_dairy_LI, "[[", "p.value")
p.vals

LI_dairy_cor <- as.data.frame(cor_dairy_LI)
LI_dairy_pv <- as.data.frame (p.vals)
LI_dairy_cor <- t(LI_dairy_cor)
LI_dairy_cor2 <- merge(LI_dairy_pv, LI_dairy_cor, by= "row.names")
LI_dairy_cor2$Row.names <- paste("LI-Dairy_", LI_dairy_cor2$Row.names)

###Lactose intake - Gut complaints 
png("LLD_LI_lac_GC_combined.png")
scatter.smooth(DF_LI$lactose, DF_LI$combined, ylab= "Gastrointestinal Symptoms (score)", xlab = "Lactose Intake (g/day)", main= "LLD - LI Recessive Genotype")
dev.off()

summary(lm(DF_LI$lactose ~ DF_LI$combined+DF_LI$antrop_BMI))

##Lactose intake - all gut complaints: 
cor_lac_LI <- cor(DF_LI$lactose, DF_LI[,c("pain", "discomfort", "constipation", "diarrhea", 
                                          "bloating", "flatulence", "burping", "nausea", "combined")], use = 'pairwise.complete.obs')

pval_lac_LI <- apply(DF_LI[,c("pain", "discomfort", "constipation", "diarrhea", 
                              "bloating", "flatulence", "burping", "nausea", "combined")], 2, function(x){ 
                                      cor.test(DF_LI$lactose, x, method =  "pearson", use = "pairwise")
                              })
p.vals <- sapply( pval_lac_LI, "[[", "p.value")
p.vals

LI_lac_cor <- as.data.frame(cor_lac_LI)
LI_lac_pv <- as.data.frame (p.vals)
LI_lac_cor <- t(LI_lac_cor)
LI_lactose_cor <- merge(LI_lac_pv, LI_lac_cor, by= "row.names")
LI_lactose_cor$Row.names <- paste("LI-Lactose_", LI_lactose_cor$Row.names)

###Bifidobacterium - Gut complaints 
png("LLD_LI_bifido_GC_combined.png")
scatter.smooth(DF_LI$Bifido_genus_abundance, DF_LI$combined, ylab= "Gastrointestinal Symptoms (score)", xlab = "Bifidobacterium abundance", main= "LLD - LI Recessive Genotype")
dev.off()

cor.test(DF_LI$Bifido_genus_abundance, DF_LI$combined)

##Bifido - bloating 
png("LLD_LI_bifido_GC_bloating.png")
scatter.smooth(DF_LI$Bifido_genus_abundance, DF_LI$bloating, ylab= "Bloating (rank transformation)", xlab = "Bifidobacterium abundance (rank transformation)", main= "LLD - LI Recessive Genotype")
dev.off()

##Bifido - burping 
png("LLD_LI_bifido_GC_burping.png")
scatter.smooth(DF_LI$Bifido_genus_abundance, DF_LI$burping, ylab= "Burping (rank transformation)", xlab = "Bifidobacterium abundance (rank transformation)", main= "LLD - LI Recessive Genotype")
dev.off()

##Bifido - discomfort 
png("LLD_LI_bifido_GC_discomfort.png")
scatter.smooth(DF_LI$Bifido_genus_abundance, DF_LI$discomfort, ylab= "Discomfort (rank transformation)", xlab = "Bifidobacterium abundance (rank transformation)", main= "LLD - LI Recessive Genotype")
dev.off()

##Bifido - Flatulence 
png("LLD_LI_bifido_GC_flatulence.png")
scatter.smooth(DF_LI$Bifido_genus_abundance, DF_LI$flatulence, ylab= "Flatulence (rank transformation)", xlab = "Bifidobacterium abundance (rank transformation)", main= "LLD - LI Recessive Genotype")
dev.off()

##Bifido - nausea 
png("LLD_LI_bifido_GC_nausea.png")
scatter.smooth(DF_LI$Bifido_genus_abundance, DF_LI$nausea, ylab= "Nausea (rank transformation)", xlab = "Bifidobacterium abundance (rank transformation)", main= "LLD - LI Recessive Genotype")
dev.off()

##Bifido - bloating 
png("LLD_LI_bifido_GC_pain.png")
scatter.smooth(DF_LI$Bifido_genus_abundance, DF_LI$pain, ylab= "Pain (rank transformation)", xlab = "Bifidobacterium abundance (rank transformation)", main= "LLD - LI Recessive Genotype")
dev.off()



###Bifidobacterium - all gut complaints: 
cor_bifido_LI <- cor(DF_LI$Bifido_genus_abundance, DF_LI[,c("pain", "discomfort", "constipation", "diarrhea", 
                                                            "bloating", "flatulence", "burping", "nausea", "combined")], use = 'pairwise.complete.obs')

pval_bifido_LI <- apply(DF_LI[,c("pain", "discomfort", "constipation", "diarrhea", 
                                 "bloating", "flatulence", "burping", "nausea", "combined")], 2, function(x){ 
                                         cor.test(DF_LI$Bifido_genus_abundance, x, method =  "pearson", use = "pairwise")
                                 })
p.vals <- sapply( pval_bifido_LI, "[[", "p.value")
p.vals

LI_bifido_cor <- as.data.frame(cor_bifido_LI)
LI_bifido_pv <- as.data.frame (p.vals)
LI_bifido_cor <- t(LI_bifido_cor)
LI_bifidobacterium_cor <- merge(LI_bifido_pv, LI_bifido_cor, by= "row.names")
LI_bifidobacterium_cor$Row.names <- paste("LI-Bifidobacterium Abundance_", LI_bifidobacterium_cor$Row.names)

###Dairy - Bifidobacterium abundance 
png("LLD_LI_dairy_bifido.png")
scatter.smooth(DF_LI$dairy, DF_LI$Bifido_genus_abundance, ylab="Bifidobacterium abundance", xlab = "Dairy intake (g/day)", main = "LLD - LI Recessive Genotype")
dev.off()

cor.test(DF_LI$Bifido_genus_abundance, DF_LI$dairy)

###Lactose - Bifidobacterium abundance 
png("LLD_LI_lac_bifido.png")
scatter.smooth(DF_LI$lactose, DF_LI$Bifido_genus_abundance, ylab="Bifidobacterium abundance", xlab = "Lactose intake (g/day)", main = "LLD - LI Recessive Genotype")
dev.off()

cor.test(DF_LI$Bifido_genus_abundance, DF_LI$lactose)

##merge all the correlation dfs
LI_cor_final <- rbind(LI_dairy_cor2, LI_lactose_cor, LI_bifidobacterium_cor)
write.table(LI_cor_final, file ="LLD_LI_pvals.txt",  col.names =T, row.names =T, append =F)

####Final plots
#Dairy - Bifido
summary(lm(DF_LI$dairy ~DF_LI$Bifido_genus_abundance+DF_LI$antrop_BMI))
png("letter_LLD_LI_bif_dairy.png")
plot (x= DF_LI$Bifido_genus_abundance, y=DF_LI$dairy, col= "#CC0000", pch=16 , cex=2 , xlab = "Bifidobacterium genus abundance (rank transformation)", ylab= "Dairy Intake (rank transformation)", main ="LLD - LI")
abline(lm(dairy ~ Bifido_genus_abundance, DF_LI), col = "#CC0000", lwd =4)
text (2.5, 2.6, "P-val: 0.038")
dev.off()

###dairy - GC (combined)
cor.test(DF_LI$combined, DF_LI$dairy)

png("letter_LLD_LI_dairy_GC.png")
plot (x= DF_LI$dairy, y=DF_LI$combined, col= "#CC0000", pch=16 , cex=2 , xlab = "Dairy Intake (rank transformation)", ylab= "Gastrointestinal Complaints (rank transformation)", main ="LLD - LI")
abline(lm(combined ~ dairy, DF_LI), col = "#CC0000",  lwd =4)
text (2.2, 2.4, "P-val: 0.33")
dev.off()

###dairy - GC (pain)
cor.test(DF_LI$pain, DF_LI$dairy)

png("letter_LLD_LI_dairy_pain.png")
plot (x= DF_LI$dairy, y=DF_LI$pain, col= "#CC0000", pch=16 , cex=2 , xlab = "Dairy Intake (rank transformation)", ylab= "Gastrointestinal Pain (rank transformation)", main ="LLD - LI")
abline(lm(pain ~ dairy, DF_LI), col = "#CC0000",  lwd =4)
text (2.2, 2.4, "P-val: 0.878")
dev.off()

###dairy - GC (discomfort)
cor.test(DF_LI$discomfort, DF_LI$dairy)

png("letter_LLD_LI_dairy_discomfort.png")
plot (x= DF_LI$dairy, y=DF_LI$discomfort, col= "#CC0000", pch=16 , cex=2 , xlab = "Dairy Intake (rank transformation)", ylab= "Gastrointestinal Discomfort (rank transformation)", main ="LLD - LI")
abline(lm(discomfort ~ dairy, DF_LI), col = "#CC0000",  lwd =4)
text (2.2, 2.4, "P-val: 0.796")
dev.off()

###dairy - GC (bloating)
cor.test(DF_LI$bloating, DF_LI$dairy)

png("letter_LLD_LI_dairy_bloating.png")
plot (x= DF_LI$dairy, y=DF_LI$bloating, col= "#CC0000", pch=16 , cex=2 , xlab = "Dairy Intake (rank transformation)", ylab= "Gastrointestinal Bloating (rank transformation)", main ="LLD - LI")
abline(lm(bloating ~ dairy, DF_LI), col = "#CC0000",  lwd =4)
text (2.2, 2.4, "P-val: 0.744")
dev.off()


###bifido - GC (combined)
cor.test(DF_LI$combined, DF_LI$Bifido_genus_abundance)

png("letter_LLD_LI_bifido_GC.png")
plot (x= DF_LI$Bifido_genus_abundance, y=DF_LI$combined, col= "#CC0000", pch=16 , cex=2 , xlab = "Bifidobacterium genus abundance (rank transformation)", ylab= "Gastrointestinal Complaints (rank transformation)", main ="LLD - LI")
abline(lm(combined ~ Bifido_genus_abundance, DF_LI), col = "#CC0000",  lwd =4)
text (2.2, 2.4, "P-val: 0.003")
dev.off()

###bifido - GC (pain)
cor.test(DF_LI$pain, DF_LI$Bifido_genus_abundance)

png("letter_LLD_LI_bifido_pain.png")
plot (x= DF_LI$Bifido_genus_abundance, y=DF_LI$pain, col= "#CC0000", pch=16 , cex=2 , xlab = "Bifidobacterium genus abundance (rank transformation)", ylab= "Gastrointestinal Pain (rank transformation)", main ="LLD - LI")
abline(lm(pain ~ Bifido_genus_abundance, DF_LI), col = "#CC0000",  lwd =4)
text (2.2, 2.4, "P-val: 0.019")
dev.off()

###bifido - GC (discomfort)
cor.test(DF_LI$discomfort, DF_LI$Bifido_genus_abundance)

png("letter_LLD_LI_bifido_discomfort.png")
plot (x= DF_LI$Bifido_genus_abundance, y=DF_LI$discomfort, col= "#CC0000", pch=16 , cex=2 , xlab = "Bifidobacterium genus abundance (rank transformation)", ylab= "Gastrointestinal Discomfort (rank transformation)", main ="LLD - LI")
abline(lm(discomfort ~ Bifido_genus_abundance, DF_LI), col = "#CC0000",  lwd =4)
text (2.2, 2.4, "P-val: 0.00035")
dev.off()

###bifido - GC (bloating)
cor.test(DF_LI$bloating, DF_LI$Bifido_genus_abundance)

png("letter_LLD_LI_bifido_bloating.png")
plot (x= DF_LI$Bifido_genus_abundance, y=DF_LI$bloating, col= "#CC0000", pch=16 , cex=2 , xlab = "Bifidobacterium genus abundance (rank transformation)", ylab= "Gastrointestinal Bloating (rank transformation)", main ="LLD - LI")
abline(lm(bloating ~ Bifido_genus_abundance, DF_LI), col = "#CC0000",  lwd =4)
text (2.2, 2.4, "P-val: 0.0065")
dev.off()

set.seed(1)
#mediation dairy intake + combined GC
bifido_dairy <- lm(DF_LI$Bifido_genus_abundance ~ DF_LI$dairy)
summary(bifido_dairy)

gut_bifidoanddairy <- lm(DF_LI$combined ~ DF_LI$Bifido_genus_abundance + DF_LI$dairy)
summary(gut_bifidoanddairy)

mediation_dairy <- mediate(bifido_dairy, gut_bifidoanddairy, treat = 'DF_LI$dairy', mediator = 'DF_LI$Bifido_genus_abundance')
summary(mediation_dairy) 

#mediation dairy intake + pain

pain_bifidoanddairy <- lm(DF_LI$pain ~ DF_LI$Bifido_genus_abundance + DF_LI$dairy)
summary(pain_bifidoanddairy)

mediation_dairy_pain <- mediate(bifido_dairy, gut_bifidoanddairy, treat = 'DF_LI$dairy', mediator = 'DF_LI$Bifido_genus_abundance')
summary(mediation_dairy_pain) 

