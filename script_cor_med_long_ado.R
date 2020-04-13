###This is the script used for the MMIT Research Project. It uses the LLD.metaphlan raw data to get the bacterial abundance and correlates it to the LLD info on gut complaints, dairy products intake and genotype of participants. 
##The aim is to check if there is a mediative role of microbiome in the relation between the lactose intolerant genotype and the display of gut complaints in lactose intolerant people. 
#Date of final analysis: 13 of April of 2020
#in this analysis I will use B. longum and B. adolescentis spp as the bifidobacteria abundance, since they were shown to be more significant in Esteban's findings. 
install.packages("mediation")
install.packages("ggpubr")
library(mediation)
library(data.table)
library(ggplot2)
library(ggsignif)

LCT<-read.table("C:/Users/Public/R/RPI/Linear_regressions/general_table_lct_snp.txt", fill= TRUE, header =TRUE)
dim(LCT)
str(LCT)
#I want to use the updated values for bifidobacteria abundance: 
updated_bacteria<- read.table("C:/Users/Public/R/RPI/LLD.metaphlan.raw.tsv", fill = TRUE, header = TRUE)
rownames(updated_bacteria)
bifido_updated <- updated_bacteria [171:195, ]
bifido_updated <- bifido_updated [c(7,22),]
#Now i will sum it up and have the total abundance of the bifido class 
bifido_new <- as.data.frame(t(bifido_updated))
rowSums(bifido_new)
bifido_new$total_abundance <- with(bifido_new, rowSums(bifido_new))


#create new dataframe with only the id and the total abundance:
total_abundance <- bifido_new$total_abundance
rownames(bifido_new)
#This function will get the first column of LCT and copy-paste it as row names, because at first the id names were not as rownames in the data frame. 
rownames(LCT) = LCT[,1]
#intersect the rownames and this function will return me a character vector with all the sample IDs that are present in both of them.  
updated_samples <- intersect(rownames(LCT), rownames(bifido_new))
#now try to order that
LCT_updated <- LCT[updated_samples, ]
bifido_new_ordered <- bifido_new[updated_samples,]
#now try merging 
LCT_merged <- merge(LCT_updated, bifido_new_ordered, by = 'row.names')
LCT_merged_final <- subset(LCT_merged, select = c(2:8, 14) )

#Checking milk consumption in total: 
boxplot(LCT_merged_final$both)
plot(density(LCT_merged_final$both, na.rm = TRUE), main = "Milk comsumption - Both")
# limiting the x axis to get a more  normal distribution (take out outliers)
plot(density(LCT_merged_final$both, na.rm = TRUE), xlim= c(0,1000), main = "Milk comsumption - Both")
plot(density(LCT_merged_final$both, na.rm = TRUE), xlim= c(0,1500), main = "Milk comsumption - Both") 
#take out outliers 
LCT_mf_ol_rm <- LCT_merged_final[!(LCT_merged_final$both >900),]

#Check the SNPs count to see if the recessive trait has really the lower occurance: 
SNPs <- data.table(LCT_mf_ol_rm$sampleID, genotype = LCT_mf_ol_rm$rs4988235_1_2_3)
genotypes <- ggplot(SNPs, aes(x=SNPs$genotype)) +geom_bar()
print(genotypes)
#select the recessive trait, named here as 1 (also the genotype of interest)
LCT_mf_ol_rm$rs_is1 = LCT_mf_ol_rm$rs4988235_1_2_3 ==1
sum(LCT_mf_ol_rm$rs_is1, na.rm = TRUE)
#make a new data frame separating populations: 
LCT_rs_is1 <- LCT_mf_ol_rm [!(LCT_mf_ol_rm$rs_is1 == 'FALSE'),]
LCT_rs_is1 <- LCT_rs_is1[complete.cases(LCT_rs_is1),]

LCT_rs_isNOT1 <- LCT_mf_ol_rm[!(LCT_mf_ol_rm$rs_is1 =='TRUE'), ]
LCT_rs_isNOT1 <- LCT_rs_isNOT1[complete.cases(LCT_rs_isNOT1),]


#Cheking bifido abundance 
boxplot(LCT_mf_ol_rm$total_abundance)
plot(density(LCT_mf_ol_rm$total_abundance, na.rm=TRUE), main = "bifido abundance")

#Check for gut complaints 
boxplot(LCT_mf_ol_rm$aggreg_compl_no_diar_const)
plot(density(LCT_mf_ol_rm$aggreg_compl_no_diar_const, na.rm = TRUE), main = "gut complaints")
# in this case we also se a clear outlier, but since gut complaints are interesting in this case, should we still take them out? 

#now correlating the variables to check for significant correlation. Always remember that the wilcoxon test is a test of dependency, so the samples are dependent. whereas the cor.test tests independent samples. 
#correlating SNP - milk consumption. 
scatter.smooth(LCT_mf_ol_rm$rs4988235_1_2_3, LCT_mf_ol_rm$both)
boxplot(LCT_mf_ol_rm$both ~LCT_mf_ol_rm$rs4988235_1_2_3)
wilcox.test(LCT_mf_ol_rm$both ~ LCT_mf_ol_rm$rs_is1)

#SNP - Bifido abundance (SIGNIFICANT!)
scatter.smooth(LCT_mf_ol_rm$rs4988235_1_2_3, LCT_mf_ol_rm$total_abundance)
boxplot(LCT_mf_ol_rm$total_abundance ~LCT_mf_ol_rm$rs4988235_1_2_3, ylab = "Bifido Abundance", xlab = "SNP" )
wilcox.test(LCT_mf_ol_rm$total_abundance ~ LCT_mf_ol_rm$rs_is1)

#SNP - gut complaints 
scatter.smooth(LCT_mf_ol_rm$rs4988235_1_2_3, LCT_mf_ol_rm$aggreg_compl_no_diar_const)
boxplot(LCT_mf_ol_rm$aggreg_compl_no_diar_const ~LCT_mf_ol_rm$rs4988235_1_2_3)
wilcox.test(LCT_mf_ol_rm$aggreg_compl_no_diar_const ~LCT_mf_ol_rm$rs_is1)

#Milk intake - gut complaints 
scatter.smooth(LCT_mf_ol_rm$both, LCT_mf_ol_rm$aggreg_compl_no_diar_const, ylab = "Gut complaints", xlab = "Milk intake", main = "Gut complaints - Milk intake")
cor.test(LCT_mf_ol_rm$aggreg_compl_no_diar_const, LCT_mf_ol_rm$both)

#bifido abundance - gut complaints 
scatter.smooth(LCT_mf_ol_rm$total_abundance, LCT_mf_ol_rm$aggreg_compl_no_diar_const, ylab = "Gut complaints", xlab = "Bifido abundance", main = "Gut complaints - Bifido abundance")
cor.test(LCT_mf_ol_rm$aggreg_compl_no_diar_const, LCT_mf_ol_rm$total_abundance)

#milk intake - bifido abundance 
scatter.smooth(LCT_mf_ol_rm$both, LCT_mf_ol_rm$total_abundance, ylab = "Bifido abundance", xlab = "Milk intake", main = "Bifido abundance - Milk intake")
cor.test(LCT_mf_ol_rm$total_abundance, LCT_mf_ol_rm$both)


#Now, I will use the separated dataframe to explore these relations 

#milk intake - gut complaints 
scatter.smooth(LCT_rs_is1$both, LCT_rs_is1$aggreg_compl_no_diar_const, ylab = "Gut complaints", xlab = "Milk intake", main = "Gut complaints - Milk intake", sub = "Recessive Genotype")
cor.test(LCT_rs_is1$aggreg_compl_no_diar_const, LCT_rs_is1$both)

scatter.smooth(LCT_rs_isNOT1$both, LCT_rs_isNOT1$aggreg_compl_no_diar_const)
cor.test(LCT_rs_isNOT1$aggreg_compl_no_diar_const, LCT_rs_isNOT1$both)

#bifido abundance - gut complaints (SIGNIFICANT!!)
scatter.smooth(LCT_rs_is1$total_abundance, LCT_rs_is1$aggreg_compl_no_diar_const, ylab = "Cut complaints", xlab = "Bifido Abundance", main = "Cut complaints - Bifido abundance", sub = "Recessive Genotype")
cor.test(LCT_rs_is1$total_abundance, LCT_rs_is1$aggreg_compl_no_diar_const)

scatter.smooth(LCT_rs_isNOT1$total_abundance, LCT_rs_isNOT1$aggreg_compl_no_diar_const)
cor.test(LCT_rs_isNOT1$total_abundance, LCT_rs_isNOT1$aggreg_compl_no_diar_const)


#milk intake - bifido abundance (SIGNIFICANT!!)
scatter.smooth(LCT_rs_is1$both, LCT_rs_is1$total_abundance, ylab = "Bifido Abundance", xlab = "Milk intake", main= "Bifido Abundance - Milk Intake", sub = "Recessive Genotype")
cor.test(LCT_rs_is1$both, LCT_rs_is1$total_abundance)


scatter.smooth(LCT_rs_isNOT1$both, LCT_rs_isNOT1$total_abundance)
cor.test(LCT_rs_isNOT1$both, LCT_rs_isNOT1$total_abundance)

#sour milk intake - bifido abundance (not significant but almost)
scatter.smooth(LCT_rs_is1$sourmilk, LCT_rs_is1$total_abundance, sub= "p-val: 0.08454" )
cor.test(LCT_rs_is1$sourmilk, LCT_rs_is1$total_abundance)

scatter.smooth(LCT_rs_isNOT1$sourmilk, LCT_rs_isNOT1$total_abundance)
cor.test(LCT_rs_isNOT1$sourmilk, LCT_rs_isNOT1$total_abundance)

#not-sourmilk intake - bifido abundance (check p-value because it doesnt show significance, but the image interpretation leads otherwise)
scatter.smooth(LCT_rs_is1$milk.blue., LCT_rs_is1$total_abundance, sub= "p-val: 0.1941")
cor.test(LCT_rs_is1$milk.blue., LCT_rs_is1$total_abundance)

scatter.smooth(LCT_rs_isNOT1$milk.blue., LCT_rs_isNOT1$total_abundance)
cor.test(LCT_rs_isNOT1$milk.blue., LCT_rs_isNOT1$total_abundance)


#now trying the mediation analysis between milk intake, bifido abundance and gut complaints 

bifido_milk <- lm(LCT_rs_is1$total_abundance ~ LCT_rs_is1$both)
summary(bifido_milk)

gut_bifidoandmilk <- lm(LCT_rs_is1$aggreg_compl_no_diar_const ~ LCT_rs_is1$total_abundance + LCT_rs_is1$both)
summary(gut_bifidoandmilk)

mediation_1 <- mediate(bifido_milk, gut_bifidoandmilk, treat = 'LCT_rs_is1$both', mediator = 'LCT_rs_is1$total_abundance')
summary(mediation_1) 


#another direction (checking again)
model.m <- lm(LCT_rs_is1$total_abundance  ~ LCT_rs_is1$both)
summary(model.m)
model.y <- lm(LCT_rs_is1$aggreg_compl_no_diar_const ~ LCT_rs_is1$total_abundance + LCT_rs_is1$both)
summary(model.y)
mediation3 <- mediate(model.m, model.y, treat = "LCT_rs_is1$both", mediator ="LCT_rs_is1$total_abundance")
summary(mediation3)


#milk as mediator
model_milk1 <- lm(both ~ total_abundance, LCT_rs_is1)
model_milk2 <- lm( aggreg_compl_no_diar_const ~ both + total_abundance, LCT_rs_is1)

mediation4 <- mediate(model_milk1, model_milk2, treat = "total_abundance", mediator = "both")
summary(mediation4)



