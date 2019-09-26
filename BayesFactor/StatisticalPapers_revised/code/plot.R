library(ggplot2)

theD <- read.csv(paste0("results/",n,"_",XGen,"_",epsilonDis,"_",betabGen,".csv"))
names(theD) <- c("SNR", "NEW", "GT", "EP", "DL")
myData <- NULL
for(i in names(theD)[-1]){
    tmp <- theD[,c(names(theD)[1],i)]
    tmp$method <- i
    myData <- rbind(myData,as.matrix(tmp))
}
myData <- as.data.frame(myData,stringsAsFactors=FALSE)
names(myData)[2] <- "pvalue"
myData[,2] <- as.numeric(myData[,2])
myData$method <- factor(myData$method,levels=names(theD)[-1],labels=names(theD)[-1])
myPlot <- ggplot(myData,aes(SNR,pvalue,color=method))+
    geom_boxplot()+
    geom_hline(yintercept = 0.05,linetype="dashed",color="purple")+
    ylab("Empirical powers")+
    scale_y_continuous(breaks=c(0.05,0.25,0.5,0.75,1))+
    guides(colour=guide_legend(title=NULL),linetype=guide_legend(title=NULL))+
    theme_bw()+
    theme(
          axis.text.y = element_text(size=rel(1.7)),
          axis.text.x = element_text(size=rel(1.7)),
          axis.title.y = element_text(size=rel(1.8)),
          axis.title.x = element_text(size=rel(1.8)),
          legend.position=c(0.15,0.8),
          legend.text= element_text(size=rel(1.55)),
          legend.title = element_text(size=rel(1.8)),
          aspect.ratio =0.6,
          plot.margin=unit(c(0,0,0,0),"cm")
          )
    #facet_wrap(~z,nrow=2)+

ggsave(paste0("figure/",n,"_",XGen,"_",epsilonDis,"_",betabGen,".eps"),myPlot,scale=0.8,width = 20, height = 12,units="cm")