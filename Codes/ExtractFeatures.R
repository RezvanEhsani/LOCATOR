
#' Extract three categories of features for cells of interest

#' @param Data                  a dataframe of samples, cells, markers and spatial information
#' @param SampleID_col          Index of Sample ID or Name
#' @param CellID_col            Index of Cell ID
#' @param CellType_col          Index of Cell type
#' @param CellTypesOfInterest   A vector give cell types for study
#' @param MarkersOfInterest_col A vector give column index of markers for study
#' @param Zscore                If Markers of interest already normalized by z-score
#'                              should set by FALSE and otherwise TRUE
#' @param PositivityCutoff      Cutoff to say a marker is positive, default is 0.5
#' @param X_col and Y_col       Column indexs of X and Y position of cells
#' @param r                     Radius of circle to define the marker or circle microenvironment
#' @param minCount              Cutoff to filter samples with few cells of interest
#' @param ExistingClass_col     A vector give column index of existing classification of samples for example: 
#'                              PAM50 classes of breast cancer



#' @return A data frame:
#'                      @ first 6 columns always are information of circles: "sample_id", 
#'                                                                            "CircelID", 
#'                                                                            "InterestCellID", 
#'                                                                            "InterestCellType",
#'                                                                            "XCircle", 
#'                                                                            "YCircle".
#'                      @ followed by composition of celltypes.
#'                      @ followed by connectivity score of celltypes.
#'                      @ followed by 12 graph feaures:  "nEdges", 
#'                                                       "nVertices", 
#'                                                       "AverageDegree", 
#'                                                       "EccentricityMean", 
#'                                                       "EccentricityMin", 
#'                                                       "EccentricityMax",
#'                                                       "Transitivity", 
#'                                                       "Entropy", 
#'                                                       "Diameter", 
#'                                                       "Density", 
#'                                                       "AlgebraicConnectivity", 
#'                                                       "SpectralRadius".
#'                      @ followed by existing classification of samples.



#' @examples for Keren Data
#'      SampleID_col = 2
#'      CellID_col = 1
#'      CellType_col = 6
#'      CellTypesOfInterest = "Macrophages"
#'      MarkersOfInterest_col = NULL
#'      Zscore = FALSE
#'      PositivityCutoff = 0.5
#'      MetaCellType_col = 3
#'      X_col = 4
#'      Y_col = 5
#'      r = 100
#'      minCount = 10
#'      ExistingClass_col = 7

#' RUN MAIN FUNCTION
#'      KerenFeatures <- Extract_Features(Data)

                                ############################################################################
#################################       PLEASE DO NOT RENAME OR REARRANGE THE FEATURE MATRIX COLUMNS       ####################################
                                ############################################################################
#

mkCircleDat <- function(Dat, CellsOfInterest, r){

    DS <- list()
    SampleIDs <- unique(Dat$sample_id)
    cat('Making circel data:', '\n')
    for(i in 1:length(SampleIDs)){
        Dati <- Dat[(Dat$sample_id == SampleIDs[i]), ]
        DatInteresti <- CellsOfInterest[(CellsOfInterest$sample_id == SampleIDs[i]), ]

        xyDati <- Dati[,c("Pos_X", "Pos_Y")]
        xyDatInteresti <- DatInteresti[,c("Pos_X", "Pos_Y")]
        
        Dist_ <- Dist(rbind(xyDatInteresti, xyDati))
        L <- nrow(xyDatInteresti)
        Dist_DatInteresti <- Dist_[1:L, (L+1):nrow(Dist_)]
        colnames(Dist_DatInteresti) <- Dati$cell_id
        InterestiID <- DatInteresti$cell_id
        InterestiCellType <- DatInteresti$CellType
        InterestiX <- xyDatInteresti$Pos_X
        InterestiY <- xyDatInteresti$Pos_Y
        
        D <- apply(Dist_DatInteresti, 1, function(x) names(which(x <= r)))
        if(!is.null(dim(D))) D <- as.list(data.frame(D)) 
        
        s1 <- lengths(D)
        CircleID <- unlist(lapply(1:length(s1), function(i) rep(i,s1[i])))
        InterestiID <- unlist(lapply(1:length(s1), function(i) rep(InterestiID[i],s1[i])))
        InterestiType <- unlist(lapply(1:length(s1), function(i) rep(InterestiCellType[i],s1[i])))
        XCircle <- unlist(lapply(1:length(s1), function(i) rep(InterestiX[i],s1[i])))
        YCircle <- unlist(lapply(1:length(s1), function(i) rep(InterestiY[i],s1[i])))
        
        S <- data.frame(cell_id=unlist(D))
        S$CircleID <- CircleID
        S$InterestCellID <- InterestiID
        S$InterestCellType <- InterestiType
        S$XCircle <- XCircle
        S$YCircle <- YCircle

        S1  <- merge(S, Dati, by = "cell_id")
        DS[[i]] <- S1
        setTxtProgressBar(txtProgressBar(min = 0, max = length(SampleIDs), style = 3), i)
        }
    cat('\n')
    return(do.call("rbind", DS))
}



Voronoi_Interaction <- function(Dat){

    AllCCIs <- data.frame()
    SampleIDs <- unique(Dat$sample_id)
    cat('Making spatial (Voronoi) interactions:', '\n')
    for(i in 1:length(SampleIDs)){

        coords <- Dat[(Dat$sample_id== SampleIDs[i]), c("cell_id", "Pos_X", "Pos_Y")]
        Int_voronoi <- voronoi_adjacency(coords, cell_id~Pos_X+Pos_Y)$Adjacencies
   
        colnames(Int_voronoi) <- coords$cell_id
        rownames(Int_voronoi) <- coords$cell_id
        CCI <- as.data.frame(as.table(Int_voronoi))
        CCI <- CCI[(CCI[,3]==TRUE),c(1,2)]

        CCI <- unique(CCI)
        CCI <- CCI[!duplicated(data.frame(t(apply(CCI,1,sort)))),]
        CCI$sample_id <- rep(SampleIDs[i],nrow(CCI))
        colnames(CCI) <- c("P1","P2","sample_id")

        AllCCIs <- rbind(AllCCIs, CCI)
        setTxtProgressBar(txtProgressBar(min = 0, max = length(SampleIDs), style = 3), i)
        }
    cat('\n')
    return(AllCCIs)
    }


getGraphProperties <- function(dt){

        dt1 <- dt[,c("P1", "P2")]
	g <- graph_from_data_frame(dt1, directed = F)
        E(g)$weight <- dt$weight
	# n of edges
	nEdges <- ecount(g)
	# n of vertices
	nVertices <- vcount(g)
        #
        AverageDegree <- mean(degree(g))
        #
        Dat_MAP1 <- dt[,c("P1","CellType1")]
        Dat_MAP2 <- dt[,c("P2","CellType2")]
        map = setNames(c(Dat_MAP1$CellType1, Dat_MAP2$CellType2) , c(Dat_MAP1$P1, Dat_MAP2$P2))
        V(g)$label <- unname(map[V(g)$name])
        freqs <- table(V(g)$label)/length(V(g)$label)
        Entropy <- -sum(freqs * log2(freqs))

        #
        Eccentricityg <- eccentricity(g)
        EccentricityMean <- mean(Eccentricityg)
        EccentricityMin <- min(Eccentricityg)
        EccentricityMax <- max(Eccentricityg)
        #
        # probability that the adjacent vertices of a vertex are connected. 
        Transitivity <- transitivity(g, type = "global", isolates = 'zero')
        #
        eigenvalues <- eigen(laplacian_matrix(g))$values

        #second-smallest eigenvalue. how well the overall graph is connected.
        AlgebraicConnectivity <- sort(abs(eigenvalues), decreasing = TRUE)[2]
        SpectralRadius <- max(abs(eigenvalues))
        #
        # length of the longest path (in number of edges) between two nodes
	Diameter <- diameter(g) 
	# average number of edges between any two nodes in the network.
	Distance <- mean_distance(g)
	# ratio of the number of edges and the number of possible edges
	Density <- edge_density(g)

        #
        sample_id <- unique(dt$sample_id)
        CircleID <- unique(dt$CircleID)

        dtout <- data.frame(sample_id, CircleID, nEdges, nVertices,
                            AverageDegree, EccentricityMean, EccentricityMin,
                            EccentricityMax, Transitivity, Entropy, Diameter,
                            Density, AlgebraicConnectivity, SpectralRadius
                            )
	return(dtout)
}


getConnectivityScores <- function(CircleObject){

       ConnectivityScores <- CircleObject %>%
                             select(sample_id, CircleID, CellType1, CellType2, weight) %>%
                             melt(id.vars=c("sample_id", "CircleID","weight")) %>%
                             group_by(sample_id, CircleID, value) %>%
                             mutate(weight = ifelse(weight == 0.5, 0.75, weight)) %>%
                             mutate(weight = ifelse(weight == 0.1, 0.25, weight)) %>%
                             mutate(weight = ifelse(weight == 1, 0.5, weight)) %>%
                             select(-variable) %>%
                             mutate( ConnectivityScore = sum(weight)) %>%
                             select(-weight) %>%
                             unique() %>%
                             dcast(sample_id + CircleID  ~ value, value.var=c("ConnectivityScore")) %>%
                             replace(is.na(.), 0) 
       colnames(ConnectivityScores)[3:ncol(ConnectivityScores)] <- paste("CS_", colnames(ConnectivityScores)[3:ncol(ConnectivityScores)], sep="")

    return(ConnectivityScores)
    }


CircleObjects <- function(CircleDat, CCI, r){

    cat('Making Circle objects:', '\n')
    CircleDat <- unique(CircleDat[,c("sample_id", "cell_id", "CircleID","CellType", "MetaCellType")])
    CCI <- CCI[order(CCI$sample_id),]
    CircleDat <- CircleDat[order(CircleDat$sample_id),]

    CCIsplits <- CCI %>%
                 group_by(sample_id)
    CCIsplits <- group_split(CCIsplits)

    Circlesplits <- CircleDat %>%
                 group_by(sample_id)
    Circlesplits <- group_split(Circlesplits)

    CCIDat <- list()
    for(i in 1:length(CCIsplits)){

        Circle_i <- Circlesplits[[i]]
        CCI_i <- CCIsplits[[i]]

        colnames(CCI_i)[1] <- "cell_id"
        CCI_i1 <- Reduce(function(x,y) merge(x,y,by=c("sample_id", "cell_id"),all.y=TRUE) ,list(CCI_i, Circle_i))
        CCI_i1 <- na.omit(CCI_i1)

        colnames(CCI_i1) <- c("sample_id", "P1", "P2", "CircleID", "CellType1", "MetaCellType1")
        colnames(CCI_i1)[3] <- "cell_id"
        CCIij <- Reduce(function(x,y) merge(x,y,by=c("sample_id", "cell_id", "CircleID"),all.y=TRUE) ,list(CCI_i1, Circle_i))
        CCIij <- na.omit(CCIij)
        
        colnames(CCIij)[c(2,7,8)] <- c("P2", "CellType2", "MetaCellType2")
        CCIij <- data.frame(CCIij)
        CCIij$weight[CCIij$CellType1 == CCIij$CellType2] <- 1
        CCIij$weight[CCIij$CellType1 != CCIij$CellType2 &
                 CCIij$MetaCellType1 == CCIij$MetaCellType2] <- 0.5
        CCIij$weight[CCIij$CellType1 != CCIij$CellType2 &
                 CCIij$MetaCellType1 != CCIij$MetaCellType2] <- 0.1
        CCIDat[[i]] <- CCIij
        setTxtProgressBar(txtProgressBar(min = 0, max = length(CCIsplits), style = 3), i)
        }
    cat('\n')
    CCIDat1 <- rbindlist(CCIDat, use.names=TRUE, fill=FALSE, idcol=NULL)
   
    return(CCIDat1)
    }

# Main function
Extract_Features <- function(Data){

            #
            Data$CellType <- factor(Data$CellType, levels=unique(Data$CellType))
            #remove duplicated cellID in each sample
          
            Data <- Data %>% 
                    group_by(sample_id) %>% 
                    filter(!duplicated(cell_id)) %>%
                    data.frame()
            
            # Find Cells of interest based on MarkersOfInterest_col and CellTypesOfInterest
            if(!(is.null(CellTypesOfInterest)) & is.null(MarkersOfInterest_col)) CellsOfInterest <- Data[(Data$CellType %in% CellTypesOfInterest), ]
            if(!(is.null(MarkersOfInterest_col))){
                if(Zscore == FALSE) Data[MarkersOfInterest_col] <- scale(Data[MarkersOfInterest_col])
                Data[MarkersOfInterest_col] <- ifelse(Data[MarkersOfInterest_col] >= PositivityCutoff, 1, 0)
                if(is.null(CellTypesOfInterest)) CellsOfInterest <- Data[(Data[MarkersOfInterest_col] == 1), ]
                else CellsOfInterest <- Data[(Data[MarkersOfInterest_col] == 1 & Data$CellType %in% CellTypesOfInterest), ]
            }
            
            # Remove samples based on minCount cutoff
            Dat1 <- CellsOfInterest %>%
                    group_by(sample_id) %>%
                    filter( n() >= minCount )%>%
                    data.frame()
            Dat <- Data[(Data$sample_id %in% unique(Dat1$sample_id)), ]

            # making Circles around the cells of interst
            CircleDat <- mkCircleDat(Dat, CellsOfInterest, r)

            # remove circles with less than 10 cells
            CircleDat <- CircleDat %>%
                         group_by(sample_id, CircleID) %>%
                         filter( n() >= 10 ) %>%
                         data.frame()

            # get Phenotype features for Circles
            DCircleDat <- CircleDat[,c("sample_id", "CircleID", "CellType")]
            DCircleDat$CellType <- factor(DCircleDat$CellType, unique(DCircleDat$CellType))

            Phenotypes <- suppressMessages(dcast(DCircleDat, sample_id +  CircleID ~ CellType))

            # making spatial cell cell interactions(CCI)
            xyCircleDat <- unique(CircleDat[,c("sample_id", "cell_id", "Pos_X", "Pos_Y")])
            CCI <-  Voronoi_Interaction(xyCircleDat)

            # Combine Circle Data and spatial CCI
            CircleObject <- CircleObjects(CircleDat, CCI)

            # grouping CircleObjects by sample_id and CircleID
            GroupedCircleObject <- CircleObject %>%
                                   group_by(sample_id, CircleID)
            GroupedCircleObject <- group_split(GroupedCircleObject)

            # get graph features
            cat('Making Graph features:', '\n')
            GraphLevelDts <- pbmclapply(GroupedCircleObject, getGraphProperties,mc.cores=1,
                                        ignore.interactive = getOption("ignore.interactive", T))
            GraphLevelDts <- do.call("rbind", GraphLevelDts)

            # get connectivity features
            cat('\n')
            cat('Making Connectivity features:', '\n')

            ConnectivityLevelDts <- getConnectivityScores(CircleObject)

            #Collecte all tree types of featres and
            # Add cell of interest information to CircleObject: InterestCellID, ... 
            DDCircleDat <- unique(CircleDat[,c("sample_id", "CircleID", "InterestCellID",
                                       "InterestCellType","XCircle", "YCircle")])

            Dat_List <- list(DDCircleDat, Phenotypes, ConnectivityLevelDts, GraphLevelDts)
            AllFeatures <- Reduce(function(x,y) merge(x,y,by=c("sample_id", "CircleID"),all=TRUE) ,Dat_List)

            # Add existing classes of samples for example:
            #                   molecular subtype of breast cancer
            if(!is.null(ExistingClass_col)){
                DatClasses <- Dat[, c(SampleID_col, ExistingClass_col)]
                DatClasses <- unique(DatClasses)
                AllFeatures <- merge(AllFeatures, DatClasses, by=c("sample_id"))
                }
                
            return(AllFeatures)
        }
