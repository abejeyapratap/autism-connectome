#!/usr/bin/env Rscript
rm(list = ls())
#cat("\014")
graphics.off()

library("optparse")

option_list = list(
  make_option("--scoreFile", type="character", default="allScores.txt", help="file containing scores [default= %default]", metavar="character"),
  make_option("--summaryOutputFile", type="character", default="summaryResults.txt", help="output file name for storing summary of results of mixed models [default= %default]", metavar="character"),
  make_option("--interceptOutputFile", type="character", default="intercept_", help="output file name for storing the intercepts of the fitted lines [default= %default]", metavar="character"),
  make_option("--discardTimePoint", type="character", default="none", help="discard patients at a certain time point [default= %default]", metavar="character"),
  make_option("--precision", type="integer", default=4, help="precision of outputfile in number of digits to report p values, Adj.R2, and slopes [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$scoreFile) || is.null(opt$summaryOutputFile) || is.null(opt$interceptOutputFile)){
  print_help(opt_parser)
  stop("Arguments for --scoreFile,--summaryResultsFile, and --interceptOutputFile must be supplied", call.=FALSE)
}


scoreFile<-opt$scoreFile
summaryOutputFile<-opt$summaryOutputFile
interceptOutputFile<-opt$interceptOutputFile
discardTimePoint<-opt$discardTimePoint
precision<-opt$precision

# scoreFile <- '/home/yusuf/repo/projects/tbiStructureAndFunction/results/brainLes++/01_matching/results/processed/raw_gmvmivol/direct/connectomeLevel/mixedModel_healthy/allScores.txt'
# summaryOutputFile<- '/home/yusuf/repo/projects/tbiStructureAndFunction/results/brainLes++/01_matching/results/processed/raw_gmvmivol/direct/connectomeLevel/mixedModel_healthy/similarityTime_model_summary_noneDiscarded.txt'
# interceptOutputFile<- '/home/yusuf/repo/projects/tbiStructureAndFunction/results/brainLes++/01_matching/results/processed/Schaefer216_det_raw_nonnormalized/direct/connectomeLevel/mixedModel_healthy/cognitiveScore_intercept_s2Discarded_'
# discardTimePoint<- 'none'
# precision <- 4

#setwd('/home/yusuf/Desktop/lmer')

library(lme4)
library(nlme)
library(lmerTest)
library(pbkrtest)
library(MuMIn)
Data <- read.table(scoreFile, header=TRUE)
if(!is.null(discardTimePoint) && (discardTimePoint=="s1" || discardTimePoint=="s2" || discardTimePoint=="s3") )
  Data <- Data[!(Data$timePoint==discardTimePoint),]
Data$similarityScore2 <- Data$similarityScore*Data$similarityScore
Data$similarityScore3 <- Data$similarityScore*Data$similarityScore*Data$similarityScore

cat('===========================================================================\n')
cat('=================== similarityScore ~ daysSinceInjury =====================\n')
cat('===========================================================================\n')
### Predict executive as a function of similarity score, but keep the slope of lines fixed for each subject
#### NOTE: daysSinceInjury needs to be part of the model to distinguish each sample for a single subject.
lmer0_dsi_interaction_frm=paste("scale(similarityScore) ~ scale(daysSinceInjury) + scale(I(ageBaseline)) + scale(I(pta)) + scale(daysSinceInjury)*scale(I(ageBaseline))*scale(I(pta)) + (1|subjectId)")
lmer0_dsi_frm=paste("similarityScore ~ daysSinceInjury + I(ageBaseline) + I(pta) + gender + (1|subjectId)")
lmer1_dsi_frm=paste("similarityScore ~ daysSinceInjury + I(ageBaseline) + gender + (1|subjectId)")
lmer2_dsi_frm=paste("similarityScore ~ daysSinceInjury + I(pta) + gender + (1|subjectId)")
lmer_null_frm=paste("similarityScore ~ daysSinceInjury + gender + (1|subjectId)")


lmer0_dsi_interaction=lmer(lmer0_dsi_interaction_frm, data=Data, REML=FALSE) # lme(similarityScore ~ daysSinceInjury + pta + ageBaseline, random = ~ 1|subjectId, data=Data)
lmer0_dsi=lmer(lmer0_dsi_frm, data=Data, REML=FALSE) # lme(similarityScore ~ daysSinceInjury + pta + ageBaseline, random = ~ 1|subjectId, data=Data)
lmer1_dsi=lmer(lmer1_dsi_frm, data=Data, REML=FALSE) # lme(similarityScore ~ daysSinceInjury + ageBaseline, random=~1|subjectId, data=Data)
lmer2_dsi=lmer(lmer2_dsi_frm, data=Data, REML=FALSE) # lme(similarityScore ~ daysSinceInjury + pta, random=~1|subjectId, data=Data)
lmer_null=lmer(lmer_null_frm, data=Data, REML=FALSE) # lme(similarityScore ~ 1, random =~ 1|subjectId, data=Data)


model0_dsi_interaction<-paste("# model0_dsi => ",lmer0_dsi_interaction_frm,"\n",sep="")
model0_dsi<-paste("# model0_dsi => ",lmer0_dsi_frm,"\n",sep="")
model1_dsi<-paste("# model1_dsi => ",lmer1_dsi_frm,"\n",sep="")
model2_dsi<-paste("# model2_dsi => ",lmer2_dsi_frm,"\n",sep="")
model_null<-paste("# model_null => ",lmer_null_frm,"\n",sep="")


write("models are fit\n",stderr())

cat('\n----------------------------------------------- Goodness of fit -------------------------------------------------------------\n')
rsquared_lmer0_dsi_interaction<-r.squaredGLMM(lmer0_dsi_interaction)
rsquared_lmer0_dsi<-r.squaredGLMM(lmer0_dsi)
rsquared_lmer1_dsi<-r.squaredGLMM(lmer1_dsi)
rsquared_lmer2_dsi<-r.squaredGLMM(lmer2_dsi)
rsquared_lmer_null<-r.squaredGLMM(lmer_null)

print(model0_dsi_interaction)
print(rsquared_lmer0_dsi_interaction)
cat('------------------------------------------\n')
print(model0_dsi)
print(rsquared_lmer0_dsi)
cat('------------------------------------------\n')
print(model1_dsi)
print(rsquared_lmer1_dsi)
cat('------------------------------------------\n')
print(model2_dsi)
print(rsquared_lmer2_dsi)
cat('------------------------------------------\n')
print(model_null)
print(rsquared_lmer_null)
cat('------------------------------------------\n')


write("goodness of fit calculated\n",stderr())

num_sim<-100
# cat('----------------------------------------------- pbkr test -------------------------------------------------------------\n')
pbkr_lmer_0_1_dsi<-PBmodcomp(lmer0_dsi,lmer1_dsi,nsim=num_sim)
pbkr_lmer_0_2_dsi<-PBmodcomp(lmer0_dsi,lmer2_dsi,nsim=num_sim)
pbkr_lmer_0_dsi_null<-PBmodcomp(lmer0_dsi,lmer_null,nsim=num_sim)
pbkr_lmer_1_dsi_null<-PBmodcomp(lmer1_dsi,lmer_null,nsim=num_sim)
pbkr_lmer_2_dsi_null<-PBmodcomp(lmer2_dsi,lmer_null,nsim=num_sim)

print(pbkr_lmer_0_1_dsi)
cat('------------------------------------------\n')
print(pbkr_lmer_0_2_dsi)
cat('------------------------------------------\n')
print(pbkr_lmer_0_dsi_null)
cat('------------------------------------------\n')
print(pbkr_lmer_1_dsi_null)
cat('------------------------------------------\n')
print(pbkr_lmer_2_dsi_null)

write("pbkr test is done\n",stderr())

cat('\n----------------------------------------------- ANOVA -------------------------------------------------------------\n')
anova_lmer_0_1_dsi<-anova(lmer0_dsi,lmer1_dsi) #is PTA significant in the presence of DSI and ageBaseline TO predict similarity
anova_lmer_0_2_dsi<-anova(lmer0_dsi,lmer2_dsi) #is ageBaseline significant in the presence of DSI and ageBaseline TO predict similarity
anova_lmer_0_dsi_null<-anova(lmer0_dsi,lmer_null) #is PTA+Age significant in the presence of DSI 
anova_lmer_1_dsi_null<-anova(lmer1_dsi,lmer_null) #is Age significant in the presence of DSI 
anova_lmer_2_dsi_null<-anova(lmer2_dsi,lmer_null) #is PTA significant minimally in the presence of DSI 

print(anova_lmer_0_1_dsi)
cat('------------------------------------------\n')
print(anova_lmer_0_2_dsi)
cat('------------------------------------------\n')
print(anova_lmer_0_dsi_null)
cat('------------------------------------------\n')
print(anova_lmer_1_dsi_null)
cat('------------------------------------------\n')
print(anova_lmer_2_dsi_null)

write("anova is done\n",stderr())

cat('\n----------------------------------------------- MODEL SUMMARY -------------------------------------------------------------\n')
summary_lmer0_dsi_interaction<-summary(lmer0_dsi_interaction)
summary_lmer0_dsi<-summary(lmer0_dsi)
summary_lmer1_dsi<-summary(lmer1_dsi)
summary_lmer2_dsi<-summary(lmer2_dsi)
summary_lmer_null<-summary(lmer_null)

print(summary_lmer0_dsi_interaction)
cat('------------------------------------------\n')
print(summary_lmer0_dsi)
cat('------------------------------------------\n')
print(summary_lmer1_dsi)
cat('------------------------------------------\n')
print(summary_lmer2_dsi)
cat('------------------------------------------\n')
print(summary_lmer_null)
cat('------------------------------------------\n')

### get summary of p values and put them into a separate file
####results with PBKR testing
models<-paste(model0_dsi,model1_dsi,model2_dsi,"scoreName\tpValue_LRT\tpValue_PBtest\tpValueDSI\tintercept\tslopeDSI\tslopeAge\tslopePTA\tAdj.R2\n",sep="")
lmer_0_null<-paste("matchingAccuracy_0_dsi_null\t",anova_lmer_0_dsi_null$`Pr(>Chisq)`[2],"\t",pbkr_lmer_0_dsi_null$test$p.value[2],"\t",summary_lmer0_dsi$coefficients[2,5],"\t",summary_lmer0_dsi$coefficients[1,1],"\t",summary_lmer0_dsi$coefficients[2,1],"\t",summary_lmer0_dsi$coefficients[3,1],"\t",summary_lmer0_dsi$coefficients[4,1],"\t",round(rsquared_lmer0_dsi[1],precision),"\n")
lmer_1_null<-paste("matchingAccuracy_1_dsi_null\t",anova_lmer_1_dsi_null$`Pr(>Chisq)`[2],"\t",pbkr_lmer_1_dsi_null$test$p.value[2],"\t",summary_lmer1_dsi$coefficients[2,5],"\t",summary_lmer1_dsi$coefficients[1,1],"\t",summary_lmer1_dsi$coefficients[2,1],"\t",summary_lmer1_dsi$coefficients[3,1],"\t","0","\t",round(rsquared_lmer1_dsi[1],precision),"\n")
lmer_2_null<-paste("matchingAccuracy_2_dsi_null\t",anova_lmer_2_dsi_null$`Pr(>Chisq)`[2],"\t",pbkr_lmer_2_dsi_null$test$p.value[2],"\t",summary_lmer2_dsi$coefficients[2,5],"\t",summary_lmer2_dsi$coefficients[1,1],"\t",summary_lmer2_dsi$coefficients[2,1],"\t","0","\t",summary_lmer2_dsi$coefficients[3,1],"\t",round(rsquared_lmer2_dsi[1],precision),"\n")

cat("# Summary results for the following linear mixed effects models",paste(models,lmer_0_null,lmer_1_null,lmer_2_null,sep=""),file=summaryOutputFile,sep='\n',append=FALSE)

### write coefficients of lines into a file
write.table(coef(lmer0_dsi)$subjectId,file=paste(interceptOutputFile,"similarityScore0.txt",sep=''), append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(coef(lmer1_dsi)$subjectId,file=paste(interceptOutputFile,"similarityScore1.txt",sep=''), append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(coef(lmer2_dsi)$subjectId,file=paste(interceptOutputFile,"similarityScore2.txt",sep=''), append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(coef(lmer3_dsi)$subjectId,file=paste(interceptOutputFile,"similarityScore3.txt",sep=''), append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE, quote=FALSE)
