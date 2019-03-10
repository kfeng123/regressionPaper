library(ggplot2)

for(XGen in c("uniform","equalCor","factor"))
    for(n in c(100))
        for(epsilonDis in c("t","chi"))
            for(betabGen in c("dense","sparse")){
                theD <- read.csv(paste0("tmp/",n,"_",XGen,"_",epsilonDis,"_",betabGen,".csv"))
                names(theD) <- c("SNR", "NEW", "GT", "EP")
                myData <- NULL
                for(i in names(theD)[-1]){
                    tmp <- theD[,c(names(theD)[1],i)]
                    tmp$method <- i
                    myData <- rbind(myData,as.matrix(tmp))
                }
                myData <- as.data.frame(myData,stringsAsFactors=FALSE)
                names(myData)[2] <- "pvalue"
                myData[,2] <- as.numeric(myData[,2])
                myData$method <- factor(myData$method,levels=c("NEW","GT","EP"),labels=c("NEW","GT","EP"))
                myPlot <- ggplot(myData,aes(SNR,pvalue,color=method))+
                    geom_boxplot()+
                    geom_hline(yintercept = 0.05,linetype="dashed",color="purple")+
                    ylab("Empirical powers")+
                    scale_y_continuous(breaks=c(0.05,0.25,0.5,0.75,1))+
                    guides(colour=guide_legend(title=NULL),linetype=guide_legend(title=NULL))+
                    theme_bw()+
                    theme(legend.position=c(0.1,0.85),
                          legend.text= element_text(size=rel(1.6)),
                          legend.title = element_text(size=rel(1.8))
                          )
                    #facet_wrap(~z,nrow=2)+
                
                ggsave(paste0("tmp/",n,"_",XGen,"_",epsilonDis,"_",betabGen,".eps"),myPlot,scale=0.8)
            }