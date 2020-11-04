library(minfiData)
library(minfi)
library(shinyMethyl)
library(cgageR)
library(ggthemes)

options(stringsAsFactors = F)

setwd("~/Dropbox/Parmar/Karolina_UPPMAX_TF-2613/")



#
# Read data
#

baseDir<-"~/Dropbox/Parmar/Karolina_UPPMAX_TF-2613/TF-2613_200928_IDAT/"
targets <- read.csv("TF-2613_200928_IDAT/TF-2613_200928_SampleID_SentrixID.EDIT_PS.csv")


RGSet <- read.metharray.exp(base = baseDir,targets = targets,recursive = T,verbose = T)
pd<-pData(RGSet)

manifest <- getManifest(RGSet)
MSet <- preprocessRaw(RGSet) 
GRset <- mapToGenome(RGSet)
beta <- getBeta(RGSet)

#
# QC
#

MSet.norm <- preprocessIllumina(RGSet, bg.correct = TRUE,
                                normalize = "controls",
                                reference = 2)
densityPlot(RGSet, sampGroups = pd$Gender, main = "Beta", xlab = "Beta")
mdsPlot(MSet.norm, numPositions = 10000,sampGroups = pd$Gender,sampNames = pd$SampleID)
qc <- getQC(MSet)
par(mfrow=c(1,1))

png("qc3.png")
densityPlot(MSet, sampGroups = merge1$SampleID,pal=tableau_color_pal(palette="Tableau 20")(20))
dev.off()

png("qc2.png",width = 800,h=600)
controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")
dev.off()


predictedSex <- getSex(GRset, cutoff = -2)$predictedSex

merge2<-data.frame(getSex(GRset),merge1)

ggplot(merge2,aes(x=xMed,y=yMed,col=Gender))+geom_point()+ggrepel::geom_text_repel(label=merge1$SampleID)+cowplot::theme_cowplot()



#
# Estimate age
#

est.ages <- getAgeR(beta,epitoc=TRUE,horvath=TRUE,hannum=TRUE,drift=FALSE,showStatusHannum=TRUE,
                    keepcpgs.epitoc=TRUE,keepcpgs.hannum=TRUE,keepres=FALSE,chrage=NULL)


merge1<-cbind(targets,est.ages$HorvathClock.output$Horvath.Est)



results <- data.frame(colnames(GRset),predictedSex,pd[1:15,],GenderMatch=predictedSex==toupper(substr(pd[1:15,"Gender"],1,1)),est.ages$HorvathClock.output$Horvath.Est)

results$AGE<-as.integer(results$AGE)

write.csv(data.frame(colnames(GRset),predictedSex,pd[1:15,],GenderMatch=predictedSex==toupper(substr(pd[1:15,"Gender"],1,1)),est.ages$HorvathClock.output$Horvath.Est),file="sampletable.csv")

corAll<-cor(results$AGE,results$Horvath.Est)

ggplot(results,aes(x=AGE,y=Horvath.Est))+geom_point()+geom_smooth(method="lm")+ggtitle("Corr: 0.67")+ggrepel::geom_text_repel(label=results$SampleID,col=ifelse(results$GenderMatch,"black","red"))+cowplot::theme_cowplot()
ggsave("Horvath.pdf")
        

results.correct_gender <- results[results$GenderMatch==TRUE,]

corCorrect<-cor(results.correct_gender$AGE,results.correct_gender$Horvath.Est)

summary(lm(results.correct_gender$AGE~results.correct_gender$Horvath.Est))

ggplot(results.correct_gender,aes(x=AGE,y=Horvath.Est))+geom_point()+geom_smooth(method="lm")+ggtitle(paste0("Corr: ",round(corCorrect,digits = 2)))+
  ggrepel::geom_text_repel(label=results.correct_gender$SampleID,col=ifelse(results.correct_gender$GenderMatch,"black","red"))+cowplot::theme_cowplot()
ggsave("Horvath.correctGender.pdf")

model<-lm(results.correct_gender$AGE~results.correct_gender$Horvath.Est)

as.formula(
  paste0("y ~ ", round(coefficients(model)[1],2), " + ", 
         paste(sprintf("%.2f * %s", 
                       coefficients(model)[-1],  
                       names(coefficients(model)[-1])), 
               collapse=" + ")
  )
)

#
# shinyMethyl
#

GRSet.norm <- preprocessQuantile(RGSet)
summary <- shinySummarize(RGSet)
summary.norm <- shinySummarize(GRSet.norm)

runShinyMethyl(summary, summary.norm)


myShinyMethylSet <- shinySummarize(RGSet)

runShinyMethyl(myShinyMethylSet)


