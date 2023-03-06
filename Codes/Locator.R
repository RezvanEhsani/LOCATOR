

#' Clustering cells of interest across of all samples 

#' @param FeatureData Dataframe of features returend by ExtractFeatures.R function
#' @param MainData    Original data, required columns: Sample ID, 
#'                                                     X position, 
#'                                                     Y position, 
#'                                                     CellType, 
#'                                                     Meta CellType (1 if cell is tumor and 0 otherwise), 
#'                                                     ExistingClasses (if there are any)
#' @param SurvivalData Survival Data, required columns: Sample ID, 
#'                                                      SurvivalTime, 
#'                                                      Censored
#' @Cutoff   A cut off to make binary classification of samples for each clustered cell of interest.
#'           values: "Mean" (default), "Median", any numerical value between 0 and 1.
                                                 


#' @return a list contains 9 elements respectively
#'      (and colud called with the name, see the @example):
#'                  @ a dataframe same as feature dataframe but a extra columns called
#'                    "Clusters" which show clustering names for each cell of interest (CellClusterData)
#'                  @ a dataframe that rows are samples and columns are Clustered cells, and values show 
#'                    proportional clusterd cells within the sample (SampleData)
#'                  @ a heatmap shows z-score (standardized values) of majority of each features in the Clusters (CellHeatmap).
#'                  @ Tsne plot of clusterd cells of interest (CellTsne).
#'                  @ a boxplot of number of different clustered cells of interest across the samples (CellBoxplot).
#'                  @ a barplot of row counts of clustered cells of interest (CellBarplot).
#'                  @ MapPlot that maps clustered cells within original Image (MapPlot).
#'                  @ Alluvium that could be used to compare binary classification of LOCATOR and existing classes (AlluviumPlot).
#'                  @ Survival plots (SurvivalPlots).

#'  @example
#'  OutPuts <- locator(FeatureData = KerenFeatures, 
#'                     MainData = Data,
#'                     SurvivalData = SurvivalData,
#'                     Cutoff = 'Mean')



# Set some color
Cell_Colors <- c('#042ff2', '#fa1203', '#b27002', '#0285a2',
                 '#fe8d0d', '#119411', '#6f00e5', '#1cf7fe', 
                 '#6c84f2', '#bdbd38', '#03c403', '#006ce5', 
                 '#cb9dfc', '#bf05e8', '#02a288', '#768d22')



alluvium_func <- function(Dat, ExistClassName, j){
  
  ExistClass <- colnames(Dat)[2]
  ExistClassName <- String(ExistClassName)
  ExistClassName <- ensym(ExistClassName)
  
  dat <- Dat %>%
    group_by(NC, !!ExistClassName) %>%
    summarise(freq = n()) %>%
    ungroup()
  
  k_sample = length(unique(Dat$NC))
  
  ggallu <- ggplot(data = dat,aes(axis1 = !!ExistClassName, axis2 = NC, y = freq)) +
    geom_alluvium(aes(fill = NC),width =  1/16,alpha=.7) + #, knot.pos = 0
    geom_stratum(width =  1/16 , fill = "gray50") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),color = "White",angle = 0) +
    scale_x_discrete(limits = c(ExistClass, "NC"),expand = c(0.15, 0.05)) +
    scale_fill_manual(values = c("#b8b3b4", Cell_Colors[j-1])) +
    theme_void()+
    coord_flip( )+
    theme(legend.position="non")
  
  return(ggallu)
}

#Main function
locator <- function(FeatureData,
                    MainData,
                    SurvivalData,
                    Cutoff = 'Mean'){
    
    NCellTypes <- length(unique(MainData$CellType))
    Nfeatues <- 6 + 2*NCellTypes + 12 # 6 circel information + celltypes + connectivity of cell types + 12 graph features
    NExistClass <- ncol(FeatureData) - Nfeatues

    L1 <- 7
    L2 <- ncol(FeatureData)-NExistClass
    FeatureData[L1:L2] <- sapply(FeatureData[L1:L2], as.numeric)
    FeatureData[L1:L2] <- FeatureData[L1:L2]/rowSums(FeatureData[L1:L2])
    FeatureData <- cbind(FeatureData[1:6],
                         FeatureData[L1:(L1+NCellTypes-1)][,order(colnames(FeatureData[L1:(L1+NCellTypes-1)]))],
                         FeatureData[(L1+NCellTypes):ncol(FeatureData)])

    F_dat <- FeatureData[,c(L1:L2)]
    F_dat <- data.table(F_dat)
	  res <- Spectre::run.flowsom(F_dat ,
                            use.cols = names(F_dat)[c(1:ncol(F_dat))])
	  res <- data.frame(res)
  	res$FlowSOM_metacluster <- paste("NC", res$FlowSOM_metacluster, sep="")
  	res$FlowSOM_metacluster <- as.factor(res$FlowSOM_metacluster)

  	FeatureData$Clusters <- res$FlowSOM_metacluster
  	ClusterNames <- as.character(unique(FeatureData$Clusters))
  	FeatureData$Clusters <- factor(FeatureData$Clusters, ClusterNames[order(nchar(ClusterNames), ClusterNames)])
  	ClusteredData <- FeatureData

  	#Heatmap
    Dat_Feature1 <- ClusteredData[L1:Nfeatues]
    Dat_Feature1$Clusters <- ClusteredData$Clusters
    Dat_Feature1$Clusters <- factor(ClusteredData$Clusters)
    Dat_heat <- data.table(Dat_Feature1)
    Dat_heat <- data.frame(Dat_heat[, sapply(.SD, function(x) list(mean(x))), by=Clusters])
    
    Dat_heat <- Dat_heat[match(ClusterNames[order(nchar(ClusterNames), ClusterNames)], Dat_heat$Clusters), ] 
    rownames(Dat_heat) <- Dat_heat$Clusters
    
    Dat_heat <- Dat_heat[,c(2:ncol(Dat_heat))]
    Dat_heat <- scale(Dat_heat)
    Dat_heat[Dat_heat > 2] =2
    Dat_heat[Dat_heat < -2] =-2

    # Heatmap
    pheat <- pheatmap(as.matrix(Dat_heat),
                      name = "z-score",
                      border_color = "grey30",
                      gaps_col=c(NCellTypes,2*NCellTypes),
                      fontsize = 7,
                      cellwidth = 11,
                      cellheight=11,
                      show_rownames = TRUE,
                      show_colnames = T,
                      fontsize_row = 7,
                      fontsize_col =7,
                      display_numbers = F,
                      cluster_rows=F,
                      cluster_cols=F
                      )
    #Tsne plot    
    tmp <- Dat_Feature1
    myClusters <- Dat_Feature1$Clusters 

    duplicates <- duplicated(tmp) | duplicated(tmp, fromLast = TRUE)
    tsne_out<-Rtsne(tmp[!duplicates,1:ncol(tmp)],max_iter=1000,seed=Nseed, perplexity=30)
    tsne_data <- tmp[!duplicates,]
    tsne_data$TSNE1 <- tsne_out$Y[,1]
    tsne_data$TSNE2 <- tsne_out$Y[,2]
    tsne_data$Clusters <- myClusters[!duplicates]

    tsne_plot <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = Clusters )) +
                 geom_point(size = 1.5,alpha=0.8) +
                 theme_bw()+
                 ylab('t-SNE2')+
                 xlab('t-SNE1')+
                 #scale_color_brewer(palette="Set2")+
                 scale_color_manual(values=Cell_Colors)+
                 theme(axis.text.y = element_text( size = 12 ),
                 axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 12),
                 legend.position = "bottom",
                 legend.text=element_text(size=12))+
                 guides(color=guide_legend(nrow=4,override.aes = list(size=4)))

    # Box plot
    Melted <- melt(table(ClusteredData$sample_id, ClusteredData$Clusters))
    colnames(Melted) <- c("sample_id", "Clusters", "Counts")

    Boxplot <- ggplot(Melted, aes(x=Clusters, y=Counts, fill=Clusters)) +
               geom_boxplot(width=0.5)+
               scale_fill_manual(values=Cell_Colors)+
               theme_classic()+
               labs(title="", y="")+
               theme(legend.position = "none")+
               coord_flip()
                                   
    Melted2 <- melt(table(ClusteredData$Clusters))
    colnames(Melted2) <- c("Clusters", "NumberNeighborhoods")

    # Bar plot
    Barplot <-  ggplot(Melted2, aes(x=Clusters, y=NumberNeighborhoods, fill=Clusters)) + 
                geom_bar(stat="identity",width=0.7)+
                scale_fill_manual(values=Cell_Colors)+
                labs(title="", y="")+
                theme_classic()+  
                theme(legend.position = "none")+
                coord_flip()
    
    #
    SampleData <- dcast(ClusteredData, sample_id  ~ Clusters)
    L <- ncol(SampleData)
    SampleData[,2:L] <- SampleData[,2:L]/rowSums(SampleData[,2:L]) 
    
    if(NExistClass!=0) {
      uniDat <- unique(ClusteredData[,c(1, (L2+1):(ncol(ClusteredData)-1))])
      SampleData <- merge(SampleData, uniDat, by='sample_id')
    }
    
    #Map plot
    colnames(ClusteredData)[3] <- "cell_id"
    Data1 <- MainData[(MainData$CellType != CellTypesOfInterest),]
    Data1 <- Data1[,c("sample_id", "cell_id", "MetaCellType", "Pos_X", "Pos_Y")]
    colnames(ClusteredData)[5] <- "Pos_X"
    colnames(ClusteredData)[6] <- "Pos_Y"
    DatCells_s <- ClusteredData[,c("sample_id", "cell_id","Clusters", "Pos_X", "Pos_Y")]
    
    colnames(DatCells_s)[3] <- "MetaCellType"

    Data1 <- Data1[(Data1$sample_id %in% ClusteredData$sample_id),]
    MainDat <- rbind(Data1, DatCells_s)
    
    MainDat$MetaCellType[MainDat$MetaCellType==0] <- "NonTumor"
    MainDat$MetaCellType[MainDat$MetaCellType==1] <- "Tumor"
    ClusterNames <- unique(as.character(ClusteredData$Clusters))
    OrderedClusterNames <- ClusterNames[order(nchar(ClusterNames), ClusterNames)]
    MainDat$MetaCellType <- factor(MainDat$MetaCellType, levels=c('Tumor', 'NonTumor', OrderedClusterNames))
    cols <- c('#C0C0C0','#F9BFBF',Cell_Colors[1:length(ClusterNames)])
    names(cols) <- c('Tumor', 'NonTumor', OrderedClusterNames)

    Clusts <- unique(ClusteredData$Clusters)
    ID <- unique(MainDat$sample_id)

    plot_all <- list()
    colnames(MainDat)[3] <- "Type"
    
    for(i in 1:length(ID)){
      Datai <- MainDat[(MainDat$sample_id ==ID[i]),]
      p <- ggplot(Datai, aes(Pos_X, Pos_Y, color=Type))+
           geom_point(size=0.2, alpha=0.8)+
           scale_colour_manual(values = cols, drop = FALSE) +
           guides(color = guide_legend(override.aes = list(size = 3) ) )+
           theme_bw()+
           labs(title=ID[i], x ="", y = "") + 
           theme(legend.position="bottom")+
           theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
              axis.title.x = element_blank(), axis.title.y = element_blank(),
              legend.key = element_blank(),
              plot.title = element_text(size=7),
              legend.text = element_text(size = 8)
           )
      plot_all[[i]] <- p
    }
    
    MapPlot <- ggarrange(plotlist = plot_all, nrow = 3,ncol = 4,common.legend = TRUE, legend="bottom")
    
    
    #Alluvium Plots
    if(NExistClass!=0) {
        AlluviumPlot <- list()
        Ncluts <- length(unique(ClusteredData$Clusters))+1
        ExistCalssNames <- colnames(uniDat)[2:ncol(uniDat)]
        for(i in 1:length(ExistCalssNames)){
          AlluviumPlot_i <- list()
          for(j in 2:Ncluts){
            Dat <- as.data.frame(SampleData)
            if(Cutoff == 'Mean') mij <- mean(Dat[,j])
            if(Cutoff == 'Median') mij <- median(Dat[,j])
            if(Cutoff != 'Mean' & Cutoff != 'Median') mij <- Cutoff
            Dat$NC <- ifelse(Dat[,j] >= mij, 
                       paste(colnames(SampleData)[j]," +", sep="") ,
                       paste(colnames(SampleData)[j]," -", sep=""))
            EistClassAi <- Dat[,c('NC', ExistCalssNames[i])]
            EistClassAi <- na.omit(EistClassAi)
            EistClassAi[,2] <- factor(EistClassAi[,2])
            EistClassAi[,1] <- factor(EistClassAi[,1], levels=c( paste(colnames(SampleData)[j]," -", sep="") ,
                                                           paste(colnames(SampleData)[j]," +", sep="")))
      
            pij <- alluvium_func(EistClassAi, ExistClassName=colnames(EistClassAi)[2], j)
            AlluviumPlot_i[[j-1]] <- pij
          }
      AllAlluviumPlot_i  <- ggarrange(plotlist = AlluviumPlot_i, nrow = 2, ncol = 2)
      AlluviumPlot[[i]] <- AllAlluviumPlot_i
      }
    }
    else AlluviumPlot <- NULL
    
    # Survival Plots
    if(!is.null(SurvivalData)) {
        SurvivalData <- data.frame(SurvivalData)
        SurvivalData <- SurvivalData[,c('sample_id', 'SurvivalTime', 'Censored')]
        SurvivalData <- merge(SampleData, SurvivalData, by='sample_id')
        SurvPlots <- list()
        for(i in 2:Ncluts){
            Dat <- as.data.frame(SurvivalData)
            if(Cutoff == 'Mean') mi <- mean(Dat[,i])
            if(Cutoff == 'Median') mi <- median(Dat[,i])
            if(Cutoff != 'Mean' & Cutoff != 'Median') mi <- Cutoff
            Dat$NC <- ifelse(Dat[,i] >= mi, paste(colnames(SampleData)[i]," +", sep="") ,
                              paste(colnames(SampleData)[i]," -", sep="")) 
            
            Dat <- unique(Dat[,c('sample_id', 'NC', 'SurvivalTime', 'Censored')])
            Dat$SurvivalTime <- Dat$SurvivalTime  + runif(nrow(Dat), min=0.1, max=0.5)
            
            model_fit <- survfit(Surv(SurvivalTime, Censored) ~ NC, data = Dat)
            
            names(model_fit$strata) <- gsub("NC=", "", names(model_fit$strata))

            n1 <- model_fit$n[1]
            n2 <- model_fit$n[2]
            lab1 <- paste0(n1,"/",length(model_fit$n.censor[1:n1][model_fit$n.censor[1:n1]==0]))
            lab2 <- paste0(n2,"/",length(model_fit$n.censor[(n1+1):(n1+n2)][model_fit$n.censor[(n1+1):(n1+n2)]==0]))

            surv1 <- ggsurvplot(model_fit,  size = 0.25,censor.size=3,
                                break.time.by = 4, 
                                palette = c(Cell_Colors[i-1], "#b8b3b4"),
                                risk.table = FALSE, risk.table.y.text.col = FALSE,
                                font.x = c(9),
                                font.y = c(9),
                                pval = TRUE,
                                data = Dat)

            DF <- surv1$data.survplot
            if(DF[n1,"surv"] > DF[n1+n2,"surv"]) {
              ly2 <- DF[(n1+n2),"surv"]-0.05
              ly1 <- DF[n1,"surv"]+0.05
            }
            else{
              ly2 <- DF[(n1+n2),"surv"]+0.05
              ly1 <- DF[n1,"surv"]-0.05
              
            }
            surv1$plot <- surv1$plot+
              theme(legend.title=element_blank()) +
              ggplot2::annotate("text", 
                                x = DF[n1,"time"], y = ly1, # x and y coordinates of the text
                                label = lab1, size = 3)+
              ggplot2::annotate("text", 
                                x = DF[(n1+n2),"time"], y = ly2, # x and y coordinates of the text
                                label = lab2, size = 3)
            SurvPlots[[i-1]] <- surv1
        }
        SurvivalPlots  <- arrange_ggsurvplots(SurvPlots, print = TRUE, nrow = 2, ncol = 2)
    }
    else SurvivalPlots <- NULL
            
    return(list(ClusteredData = ClusteredData,
                SampleData    = SampleData,
                Heatmap   = pheat,
                Tsne      = tsne_plot,
                Boxplot   = Boxplot,
                Barplot   = Barplot,
                MapPlot   = MapPlot,
                AlluviumPlot = AlluviumPlot,
                SurvivalPlots = SurvivalPlots
                ))
    }

