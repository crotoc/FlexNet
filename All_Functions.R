### Functions called from RWR-M.R

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# GENERAL FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 2.- We check the input parameters
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
existParameterAndFile <- function(Network,config){
    file_list <- lapply(config[[Network]], function(x){
        if(is.null(config[[x]])){
            if(opt$verbose)
                message("STEP: ",paste("This loaded layer", x,"IS NOT in config files"))
        }else{
            if(file.exists(config[[x]]) && !dir.exists(config[[x]])){
                message("STEP: ",paste("LoadedLayers",x,"IS FILE at:",paste(config[[x]],collapse = ", ")))
                a <- config[[x]]
                names(a) <- paste0(Network,"::",x,"::",x)
                a
            }else if(file.exists(config[[x]]) && dir.exists(config[[x]])) {
                if(opt$verbose)
                    message("STEP: ",paste("LoadedLayers",x,"IS DIR at:",paste(config[[x]],collapse = ", "),"\n"))
                file_list <- list.files(config[[x]],full.name = T)
                file_list <- file_list[grep(".gr$",file_list,perl = T)]
                print(Network)
                print(x)
                if(length(file_list)>1){
                    names(file_list) <- paste0(Network,"::",x,"::",basename(file_list))
                }else{
                    names(file_list) <- paste0(Network,"::",x,"::",x)
                }
                file_list
            }else{
                if(opt$verbose)
                    message("STEP: ",paste("The file for the loaded layer", x," IS NOT at",config[[x]]),"\n")
            }
        }
    }) %>% unlist
    return(file_list)
}

getNetworkSpecificParams <- function(Network,layers_list,Prefix,config, param.max = 1,param.min = 0){
    Layers <- config[[Network]]
    params <- as.numeric(config[[paste(Prefix,Network,sep = ".")]])
    if(Prefix == "Restart.Probability" || length(grep("Interlayer.Jump",Prefix)>0 )){
        all_files <- layers_list[[Network]]
        if(length(all_files) == length(params)){
            params_fill <- lapply(1:length(params),function(i){
                a <- rep(params[i],length(all_files[[i]]))
                names(a) <- names(all_files[[i]])
                a
            })
        }else if(length(all_files) < length(params)){
            if(opt$verbose){
                message("STEP: ","There are More entries of", paste(Prefix,Network,sep = "."), " than ", Network, "\n")
                message("STEP: ","Matched entries is used \n")
            }
            params_fill <- lapply(1:length(all_files),function(i){
                a <- rep(params[i],length(all_files[[i]]))
                names(a) <- names(all_files[[i]])
                a
            })
        }else if(length(all_files) > length(params)){
            stop("There are More networks of ",Network,  " than the entries of", paste(Prefix,Network,sep = "."), ". Please fix it\n")
        }
        params <- params_fill
    }
    
    message("STEP: ",paste(Prefix, Network,paste(unlist(params),collapse = ","))," ",length(unlist(params)))
    if(length(params)==1 && !is.list(params)){
        if ((params > param.max || params < param.min)){ stop("Incorrect params[[",Network,"]] it must be between ",param.min,"and",param.max)}
    }
    params <- unlist(params)
    return(params)
}



readNetworkSpecificParams <- function(Network,Prefix,config, param.max = 1,param.min = 0){
    Layers <- config[[Network]]
    params <- as.numeric(config[[paste(Prefix,Network,sep = ".")]])
    if(Prefix == "Restart.Probability" || length(grep("Interlayer.Jump",Prefix)>0 )){
        all_files <- lapply(Layers,function(x){existParameterAndFile(x,config)})
        if(length(all_files) == length(params)){
            params_fill <- lapply(1:length(params),function(i){
                a <- rep(params[i],length(all_files[[i]]))
                names(a) <- names(all_files[[i]])
                a
            })
        }else if(length(all_files) < length(params)){
            if(opt$verbose){
                message("STEP: ","There are More entries of", paste(Prefix,Network,sep = "."), " than ", Network, "\n")
                message("STEP: ","Matched entries is used \n")
            }
            params_fill <- lapply(1:length(all_files),function(i){
                a <- rep(params[i],length(all_files[[i]]))
                names(a) <- names(all_files[[i]])
                a
            })
        }else if(length(all_files) > length(params)){
            stop("There are More networks of ",Network,  " than the entries of", paste(Prefix,Network,sep = "."), ". Please fix it\n")
        }
        params <- params_fill
    }
    
    message("STEP: ",paste(Prefix, Network,paste(unlist(params),collapse = ","))," ",length(unlist(params)))
    if(length(params)==1 && !is.list(params)){
        if ((params > param.max || params < param.min)){ stop("Incorrect params[[",Network,"]] it must be between ",param.min,"and",param.max)}
    }
    params <- unlist(params)
    return(params)
}

check.parameters.length <- function(config,Network,Prefix,l){
    ## message("STEP: ",paste0(Prefix,".",Network),"\n")
    if(length(config[[paste0(Prefix,".",Network)]]) != l){
        stop(paste0(Prefix,".",Network), "must be equal to ",l)
    }
}


check.parameters.from.yaml.general <- function(opt){
    config <- opt$config
    r <- as.numeric(config$Global.Restart.Probability)
    if(opt$verbose)
        message("STEP: ",paste("Global restart probability (r):", r),"\n")
    if ((r > 1 || r <= 0)){ stop("Incorrect r, it must be between 0 and 1")}
    config$r <- r

    ## lambda inter network jump probability
    lambda <- config$Inter.Network.Jump.Probability.Table
    if(opt$verbose)
        message("STEP: ",paste("Inter Network Jump Probability Table (lambda):", paste(lambda,collapse = " ")),"\n")
    if(length(lambda)/(length(config$LoadedNetworks) * length(config$LoadedNetworks)) != 1) {stop("Incorrect Inter.Network.Jump.Probability, the sum of each lambda's component divided by number of networks must be 1","\n")}
    lambda <- matrix(lambda,nrow = length(config$LoadedNetworks), ncol = length(config$LoadedNetworks),byrow = T)
    colnames(lambda) <- config$LoadedNetworks
    rownames(lambda) <- config$LoadedNetworks
    config$lambda <- lambda


    ## eta weight for each type of seeds
    eta <- as.numeric(config$Networks.Restart.Probability)
    if(opt$verbose)
        message("STEP: ",paste("Networks Restart Probabily (eta):", paste(eta,collapse = " ")),"\n")
    ## if ((eta > 1 || eta < 0)){ stop("Incorrect eta, it must be between 0 and 1")}
    names(eta) <- config$LoadedNetworks
    config$eta <- eta

    ## output control
    if(opt$verbose)
        message("STEP: Output top K:")
    k <- lapply(config$LoadedNetworks,function(Network){readNetworkSpecificParams(Network,"Output.Top.K",config,100,0)})
    names(k) <- config$LoadedNetworks
    config$k <- k
    if(opt$verbose)
        message("\n")


    ## multiplex_layers_list <- lapply(config$LoadedNetworks,function(Network){out <- lapply(config[[Network]],function(x){suppressMessages(existParameterAndFile(x,config))});unlist(out)})
    ## names(multiplex_layers_list) <- config$LoadedNetworks
    ## config$multiplex_layers_list <- multiplex_layers_list

    ## bipartite_layers_list <-  lapply(unlist(lapply(combn(config$LoadedNetworks,2,simplify =F),function(x){c(paste0(x,collapse = "."),paste0(rev(x),collapse = "."))})),function(x){suppressMessages(existParameterAndFile(paste0("Bipartite.",x),config))})
    ## names(bipartite_layers_list) <- unlist(lapply(combn(config$LoadedNetworks,2,simplify =F),function(x){c(paste0(x,collapse = "_"),paste0(rev(x),collapse = "_"))}))
    ## bipartite_layers_list[sapply(bipartite_layers_list,is.null)] <- NULL
    ## config$bipartite_layers_list <- bipartite_layers_list

    return(config)
}




#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## MULTIPLEX RELATED FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 4.- We read the different layers that integrate the multiplex network.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.layer.number <- function(Networkid,config){
    layers_list <- lapply(config[[NetworkID]],function(x){existParameterAndFile(x,config)}) %>% unlist
    length(Layers)
}

get.adj.multiplex.by.networkid <- function(NetworkID,config,returnGraph=F){
    ## We read the different Networks (Layers) of our Multiplex network. We also simplify the networks
    ## by removing possible self loops and possible multiple nodes.
    layers_list <- lapply(config[[NetworkID]],function(x){existParameterAndFile(x,config)}) %>% unlist
    Layers <- lapply(1:length(layers_list),function(i){
        if(opt$verbose)
            message("STEP: Read #",i," ",layers_list[i],"\n")
        tb <- fread(layers_list[i],header = F)
        if(ncol(tb)==3){
            names(tb)[3] <- "weight"
        }
        Network <- graph.data.frame(tb,directed=FALSE)
        Network <- igraph::simplify(Network, remove.multiple = TRUE, remove.loops = TRUE)
        Network
    })

    ## save delta
    ## DELTA Weight for keeping in same layer or jumper to another homogenous layers
    if(opt$verbose)
        message("\nSTEP: delta:")
    check.parameters.length(config,NetworkID,"Interlayer.Jump.Probability",length(config[[NetworkID]]))
    delta <- readNetworkSpecificParams(NetworkID,"Interlayer.Jump.Probability",config)

    ## Allow DELTA. Whether allow jump from other layers to the specified layer.
    if(opt$verbose)
        message("\nSTEP: delta_allow")
    check.parameters.length(config,NetworkID,"Interlayer.Jump.Allow",length(config[[NetworkID]]));
    delta_allow <- readNetworkSpecificParams(NetworkID,"Interlayer.Jump.Allow",config)

    ## tau Weight of seeds for each layers in homogenous layers
    if(opt$verbose)
        message("\nSTEP: tau")
    tau <-  readNetworkSpecificParams(NetworkID,"Restart.Probability",config)

    ## We get a pool of nodes (Nodes in any of the layers.) 
    ## Some nodes will belong to just some layers, while others will be on all of them.
    pool_nodes <- sort(pool.of.nodes(Layers))
    
    ## We add to each layer the missing nodes with no connections (no edges for them)
    Layers <- add.missing.nodes(Layers, pool_nodes)
    names(Layers) <- layers_list

    ## We have to check that all the layers have the same number of Vertex. 
    ## We save N_Gene as the number of Nodes in all the layers (After adding those that were missing.)
    N <- Get.Number.Nodes(Layers)
    L <- length(Layers)
    LayerName <- gsub(".mat.origin.*","",names(layers_list),perl =T) %>% gsub("\\.","::",.,perl =T)

    ## get adjacency matrix before applying delta
    AdjacencyMatrix <- lapply(1:L,function(i){
        Adjacency_Layer <- matrix()
        if(opt$weighted && length(grep("weight",edge_attr_names(Layers[[i]])))>0){
            if(opt$verbose){
                message("STEP: use weighted network\n")
            }
            Adjacency_Layer <-  as_adjacency_matrix(Layers[[i]],sparse = TRUE,attr="weight")
        }else{
            if(opt$verbose){
                message("STEP: use unweighted network\n")
            }
            Adjacency_Layer <-  as_adjacency_matrix(Layers[[i]],sparse = TRUE)
        }
        ## the original code has a error, Adjacentcy_Layer should be column normalized to indicate the ration of (staying in the layer/jumping to other layers) as (1-delta)/delta 
        ## To make sure the ratio(staying/jumping) = delta, Adjacency_Layers should be col normalized before. If it's not, the ratio depends on the degree of each vertex.
        message("STEP: Column normalize AdjacencyMatrix for: ", names(Layers)[i],"\n")
        Adjacency_Layer <- colNorm(Adjacency_Layer)
        if(opt$verbose){
            message("STEP: ","The col sum of AdjacencyMatrix is: \n")
            message(paste0(capture.output(colSums(Adjacency_Layer) %>% round(.,5) %>% table), collapse = "\n"))
        }
        
        ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
        Adjacency_Layer <- Adjacency_Layer[match(pool_nodes,rownames(Adjacency_Layer)),match(pool_nodes,colnames(Adjacency_Layer))]
        Adjacency_Layer
    })

    ## Adjacency matrix of the multiplex network. (Transition Matrix for the multiplex system.)
    if(opt$verbose)
        message("STEP: ","Generating Multiplex Adjacency Matrix...",NetworkID,"\n")

    if(returnGraph){
        return(list(AdjacencyMatrix=AdjacencyMatrix,N=N,L=L,pool_nodes = pool_nodes,LayerName=LayerName,graphs=Layers,delta=delta, delta_allow=delta_allow,tau=tau))
    }else{
        return(list(AdjacencyMatrix=AdjacencyMatrix,N=N,L=L,pool_nodes = pool_nodes,LayerName=LayerName,delta=delta, delta_allow=delta_allow,tau=tau))
    }
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 5.- We generate a pool of nodes. We merge the nodes present in every layer.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
pool.of.nodes <- function(Layers){
    
    ## We get the number of layers
    Nr_Layers <- length(Layers)
    
    ## We get the nodes of all the layers of the multiplex network. We save them into a vector.
    Node_Names_all <- character()
    for (i in 1:Nr_Layers) {
        Node_Names_Layer <- V(Layers[[i]])$name
        Node_Names_all <-c(Node_Names_all,Node_Names_Layer)
    }

    ## We remove duplicates.
    Node_Names_all <- unique(Node_Names_all)
    
    return(Node_Names_all)
} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 6.- From the pool of nodes we add the missing proteins to each layer as 
##     isolated nodes.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
add.missing.nodes <- function (Layers,NodeNames) {
        
    ## We get the number of layers
    Nr_Layers <- length(Layers)
    
    ## We generate a new list of layers.
    Layers_New <- vector("list", Nr_Layers)
    
    ## We add to each layer the missing nodes of the total set of nodes, of the pool of nodes.
    for (i in 1:Nr_Layers){
        Node_Names_Layer <- V(Layers[[i]])$name
        Missing_Nodes <- NodeNames[which(!NodeNames %in% Node_Names_Layer)]
        Layers_New[[i]] <- add_vertices(Layers[[i]] ,length(Missing_Nodes), name=Missing_Nodes)
    }
    return(Layers_New)
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
    ## 7.- We check the total number of nodes in every layer and we return it. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
Get.Number.Nodes <- function(Layers_Allnodes) {
        
    ## We get the number of layers
    Nr_Layers <- length(Layers_Allnodes)
    vector_check <- numeric(length = Nr_Layers)  
    
    for (i in 1:Nr_Layers){
        vector_check[i] <- vcount(Layers_Allnodes[[i]])  
    }
    
    if (all(vector_check == vector_check[1])){
        if(opt$verbose)
            message("STEP: ","Number of nodes in every layer updated...\n")
        return(vector_check[1])  
    } else {
        stop("Not correct number of nodes in each Layer...\n")
    }
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
    ## 8.- We generate the multiplex adjacency matrix (Transition matrix) 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

colNorm <- function(AdjacencyMatrix){
    AdjacencyMatrix <- as.csc.matrix(AdjacencyMatrix)
    AdjacencyMatrix@x <- AdjacencyMatrix@x / rep.int(ifelse(colSums(AdjacencyMatrix)==0,1,colSums(AdjacencyMatrix)), diff(AdjacencyMatrix@p))
    AdjacencyMatrix
}

rowNorm <- function(AdjacencyMatrix){
    AdjacencyMatrix <- t(AdjacencyMatrix)
    AdjacencyMatrix <- as.csc.matrix(AdjacencyMatrix)
    AdjacencyMatrix@x <- AdjacencyMatrix@x / rep.int(ifelse(colSums(AdjacencyMatrix)==0,1,colSums(AdjacencyMatrix)), diff(AdjacencyMatrix@p))
    t(AdjacencyMatrix)
}


get.supra.adj.multiplex.by.networkid <- function(NetworkID,AdjacencyMatrixList,config){

    ## 1. The diagnal of the SupraAdjacencyMatrix should be the AdjacencyMatrix for each layer of gene networks
    ## 2. The off diagnal should be the idem matrix
    ## 3. Delta defines the probability of jumping to another layer, so 1-delta defines the probability of staying in a layer.
    ## 4. If a layer only can jump to a subset set of layers (L'), the weight for each layer of subset should be delta/L';
    ## If a layer can jump to any other layers, L' = L-1; See the demonstration in the codes.
    ## delta_allow_mat is generated to describe which layer can jump to any other layers and which layer can only jump a subset of layers.
    
    Adjacency_Layer_list <- AdjacencyMatrixList[[paste0(NetworkID,"_",NetworkID)]]$AdjacencyMatrix
    Adjacency_Layer_list_name <- AdjacencyMatrixList[[paste0(NetworkID,"_",NetworkID)]]$LayerName
    delta <- AdjacencyMatrixList[[paste0(NetworkID,"_",NetworkID)]]$delta
    delta_allow <- AdjacencyMatrixList[[paste0(NetworkID,"_",NetworkID)]]$delta_allow
    N <- AdjacencyMatrixList[[paste0(NetworkID,"_",NetworkID)]]$N
    tau <- AdjacencyMatrixList[[paste0(NetworkID,"_",NetworkID)]]$tau
    
    ## IDEM_MATRIX.
    Idem_Matrix <- Diagonal(N, x = 1)
    L <- length(Adjacency_Layer_list)
    
    SupraAdjacencyMatrixList <- list()

    ## if there is only 1 layer, the prob staying in the same layer is 1, so the delta is 0
    if(L==1){
        delta_allow_mat <- as.matrix(1)
    }else{
        ## to calculate delta weight matrix. Make sure that layers not allowed to jump has 0 
        delta_allow[delta_allow==0]=2
        delta_allow_mat <- delta_allow %o% delta_allow
        delta_allow_mat[delta_allow_mat==2] = 1
        delta_allow_mat[delta_allow_mat==2^2] = 0
        ## e.g. 1 indicates it can jump and 0 means it can not jump
        ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,]    1    1    1    1    1    1    1    1
        ## [2,]    1    1    1    1    1    1    1    1
        ## [3,]    1    1    0    0    0    0    0    0
        ## [4,]    1    1    0    0    0    0    0    0
        ## [5,]    1    1    0    0    0    0    0    0
        ## [6,]    1    1    0    0    0    0    0    0
        ## [7,]    1    1    0    0    0    0    0    0
        ## [8,]    1    1    0    0    0    0    0    0

        ## The diagnal of the delta_allow_mat is defined in the file, which is 1-delta.
        ## To calculate the delta weight for other layers,first mask the diagnal to 0
        diag(delta_allow_mat) <- 0
        ## [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
            ## [1,]    0    1    1    1    1    1    1    1
            ## [2,]    1    0    1    1    1    1    1    1
            ## [3,]    1    1    0    0    0    0    0    0
            ## [4,]    1    1    0    0    0    0    0    0
            ## [5,]    1    1    0    0    0    0    0    0
            ## [6,]    1    1    0    0    0    0    0    0
            ## [7,]    1    1    0    0    0    0    0    0
            ## [8,]    1    1    0    0    0    0    0    0


        ## calculate the column sum and the weight for each layers
        ## Col normalized the matrix
        delta_allow_mat <- sweep(delta_allow_mat,2,colSums(delta_allow_mat),"/")
        ## [,1]      [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,] 0.0000000 0.1428571  0.5  0.5  0.5  0.5  0.5  0.5
        ## [2,] 0.1428571 0.0000000  0.5  0.5  0.5  0.5  0.5  0.5
        ## [3,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [4,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [5,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [6,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [7,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [8,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        delta_allow_mat[is.na(delta_allow_mat)] <- 0

        ## multiple the delta to get the true weight
        delta_allow_mat <- sweep(delta_allow_mat,2,delta,"*")
        ## [,1]       [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,] 0.00000000 0.07142857 0.25 0.25 0.25 0.25 0.25 0.25
        ## [2,] 0.07142857 0.00000000 0.25 0.25 0.25 0.25 0.25 0.25
        ## [3,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [4,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [5,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [6,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [7,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [8,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00


        ## fill the diagnal with weight of staying in the layer
        diag(delta_allow_mat) <- 1-delta
        ## [,1]       [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,] 0.50000000 0.07142857 0.25 0.25 0.25 0.25 0.25 0.25
        ## [2,] 0.07142857 0.50000000 0.25 0.25 0.25 0.25 0.25 0.25
        ## [3,] 0.07142857 0.07142857 0.50 0.00 0.00 0.00 0.00 0.00
        ## [4,] 0.07142857 0.07142857 0.00 0.50 0.00 0.00 0.00 0.00
        ## [5,] 0.07142857 0.07142857 0.00 0.00 0.50 0.00 0.00 0.00
        ## [6,] 0.07142857 0.07142857 0.00 0.00 0.00 0.50 0.00 0.00
        ## [7,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.50 0.00
        ## [8,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.50

        ## col nomailize again to make sure
        delta_allow_mat <- sweep(delta_allow_mat,2,colSums(delta_allow_mat),"/")
        ## validate the col sum is 1
        round(colSums(delta_allow_mat),5) %>% table
        }

        if(opt$verbose){
            message("STEP: ","nomalized delta_allow_matrix is :\n")
        view(delta_allow_mat)
        }

    
    for (i in 1:L){
        Adjacency_Layer <- Adjacency_Layer_list[[i]]
        SupraAdjacencyMatrixList[[i]] <- do.call("cbind",lapply(1:L,function(j){
            if(i==j){
                ## For the diagnal of the SupraAdjacencyMatrix, using adj for each layer and multiple 1-delta
                delta_allow_mat[i,j] * Adjacency_Layer
            }else{
                ## For the off-diagnal of the SupraAdjacencyMatrix, using idem and multiple delta
                delta_allow_mat[i,j] * Idem_Matrix
            }
        }))

        ## We fill the diagonal blocks with the adjacencies matrix of each layer.
        if(L>1){
            rownames(SupraAdjacencyMatrixList[[i]]) <- paste(rownames(Adjacency_Layer),Adjacency_Layer_list_name[i],sep="::")
        }else{
            rownames(SupraAdjacencyMatrixList[[i]]) <- rownames(Adjacency_Layer)
        }
    }
        SupraAdjacencyMatrix <- do.call("rbind",SupraAdjacencyMatrixList)
        colnames(SupraAdjacencyMatrix) <- rownames(SupraAdjacencyMatrix)

        ## col normalized, to make sure the lambda weight of multiplex network is correct
        ## Before this step, if a gene has no edge, the col sum is 0. if a gene has edges, the col sum is 1.
        ##
        if(opt$verbose){
            message("Before the step, the col sum of SupraAdjacencyMatrix is: \n")
            message(paste0(capture.output(colSums(SupraAdjacencyMatrix) %>% round(.,5) %>% table), collapse = "\n"))
        }
        message("STEP: Column normalize SupraAdjacencyMatrix\n")
        SupraAdjacencyMatrix <- colNorm(SupraAdjacencyMatrix)

        if(opt$verbose){
            message("The col sum of SupraAdjacencyMatrix is: \n")
            message(paste0(capture.output(colSums(SupraAdjacencyMatrix) %>% round(.,5) %>% table), collapse = "\n"))
        }

        ## The Matrix should have dimension of (N_Gene * L_layers) * (N_Gene * L_layers)
        SupraAdjacencyMatrix %>% dim
        SupraAdjacencyMatrix %>% view
        
        ## validate
        ## A symetric matrix should be obtained and diagnol should be a and all other cells are i
        ## rownames(delta_allow_mat) <- NULL
        ## colnames(delta_allow_mat) <- NULL
        ## tmp <- lapply(1:L,function(i){
        ##     lapply(1:L,function(k){
        ##         N <- ncol(Adjacency_Layer_list[[i]])
        ##         if(all(SupraAdjacencyMatrix[((i-1)*N+1):(i*N),((k-1)*N+1):(k*N)] == delta_allow_mat[i,k] * Idem_Matrix)){
        ##             return("i")
        ##         }
        ##         if(all(SupraAdjacencyMatrix[((i-1)*N+1):(i*N),((k-1)*N+1):(k*N)] == delta_allow_mat[i,k] * Adjacency_Layer_list[[i]])){
        ##             return("a")
        ##         }
        ##     })
        ## })
        ## The colSums matrix should be the same to delta_allow_mat
        ## matrix(lapply(tmp,function(x){lapply(x,function(y){if(length(y$delta)==2){y$delta[y$delta!=0]}else{y$delta}})})  %>% unlist,nrow=L,byrow=T)
        ## [,1]    [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,] 0.50000 0.07143 0.25 0.25 0.25 0.25 0.25 0.25
        ## [2,] 0.07143 0.50000 0.25 0.25 0.25 0.25 0.25 0.25
        ## [3,] 0.07143 0.07143 0.50 0.00 0.00 0.00 0.00 0.00
        ## [4,] 0.07143 0.07143 0.00 0.50 0.00 0.00 0.00 0.00
        ## [5,] 0.07143 0.07143 0.00 0.00 0.50 0.00 0.00 0.00
        ## [6,] 0.07143 0.07143 0.00 0.00 0.00 0.50 0.00 0.00
        ## [7,] 0.07143 0.07143 0.00 0.00 0.00 0.00 0.50 0.00
        ## [8,] 0.07143 0.07143 0.00 0.00 0.00 0.00 0.00 0.50

        ## A symetric matrix should be obtained and diagnol should be a and all other cells are i
        ## matrix(lapply(tmp,function(x){lapply(x,function(y){y$ind})})  %>% unlist,nrow=L)
        ## [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,] "a"  "i"  "i"  "i"  "i"  "i"  "i"  "i" 
        ## [2,] "i"  "a"  "i"  "i"  "i"  "i"  "i"  "i" 
        ## [3,] "i"  "i"  "a"  "i"  "i"  "i"  "i"  "i" 
        ## [4,] "i"  "i"  "i"  "a"  "i"  "i"  "i"  "i" 
        ## [5,] "i"  "i"  "i"  "i"  "a"  "i"  "i"  "i" 
        ## [6,] "i"  "i"  "i"  "i"  "i"  "a"  "i"  "i" 
        ## [7,] "i"  "i"  "i"  "i"  "i"  "i"  "a"  "i" 
        ## [8,] "i"  "i"  "i"  "i"  "i"  "i"  "i"  "a" 

        
    return(list(SupraAdjacencyMatrix=SupraAdjacencyMatrix,N=N,L=L,pool_nodes = rownames(Adjacency_Layer_list[[1]]),LayerName=AdjacencyMatrixList[[paste0(NetworkID,"_",NetworkID)]]$LayerName,delta=delta,delta_allow=delta_allow,tau=tau))
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## HETEROGENEOUS NETWORK RELATED FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 9.- Function that generates a bipartite graph from gene-disease relations.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

get.bipartite.graph.from.edgelist <- function(V_Type_1,V_Type_2,edgelist){
    pars <- as.list(match.call()[-1])

    ## get graph from edgelist
    Bipartite_Network <- graph.data.frame(edgelist)

    ## specify types of vertices
    V_type <- data.table(V(Bipartite_Network)$name)
    V_type[V_type[[1]] %in% edgelist[[1]],c("type","group"):=list(FALSE,as.character(pars$V_Type_1))]
    V_type[V_type[[1]] %in% edgelist[[2]],c("type","group"):=list(TRUE,as.character(pars$V_Type_2))]
    V(Bipartite_Network)$group <- V_type$group
    V(Bipartite_Network)$type <- V_type$type

    ## Add missing vertices for each type
    V_Type_1_missing <- V_Type_1[! V_Type_1 %in% V(Bipartite_Network)$name[V(Bipartite_Network)$group == as.character(pars$V_Type_1)]]
    Bipartite_Network <- add_vertices(Bipartite_Network,nv = length(V_Type_1_missing),attr=list(name = V_Type_1_missing,type = FALSE, group=as.character(pars$V_Type_1)))
    V_Type_2_missing <- V_Type_2[! V_Type_2 %in% V(Bipartite_Network)$name[V(Bipartite_Network)$group == as.character(pars$V_Type_2)]]
    Bipartite_Network <- add_vertices(Bipartite_Network,nv = length(V_Type_2_missing),attr=list(name = V_Type_2_missing,type=TRUE,group = as.character(pars$V_Type_2)))
    Bipartite_Network
}

get.bipartite.graph.from.emptylist <- function(V_Type_1,V_Type_2){
    pars <- as.list(match.call()[-1])
    Bipartite_Network <- graph.data.frame(data.frame(V1=character(),V2=character()))
    Bipartite_Network <- add_vertices(Bipartite_Network,nv = length(V_Type_1),attr=list(name = V_Type_1,type = FALSE, group=as.character(pars$V_Type_1)))
    Bipartite_Network <- add_vertices(Bipartite_Network,nv = length(V_Type_2),attr=list(name = V_Type_2,type=TRUE,group = as.character(pars$V_Type_2)))
    Bipartite_Network
}

get.adj.bipartite.by.networkid <- function(Bipartite_NetworkID,AdjacencyMatrixList,config,returnGraph=FALSE){
    Bipartite_NetworkID <- unlist(Bipartite_NetworkID)
    if(opt$verbose)
        message("STEP: ","Reading bipartite graph between ",Bipartite_NetworkID,"\n")
    id <- paste(c("Bipartite",Bipartite_NetworkID),collapse = ".")

    rev <- FALSE
    if(is.null(config[[id]])){
        id <- paste(c("Bipartite",rev(Bipartite_NetworkID)),collapse = ".")
        if(!is.null(config[[id]])){
            Bipartite_NetworkID <- rev(Bipartite_NetworkID)
            rev <- TRUE
        }
    }
    
    ## get nodes names for row and column and Layer number 
    pool_nodes_Network_1 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[1],"_",Bipartite_NetworkID[1])]]$pool_nodes
    pool_nodes_Network_2 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[2],"_",Bipartite_NetworkID[2])]]$pool_nodes
    ## get layer number of network 1 and network 2
    L_Network_1 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[1],"_",Bipartite_NetworkID[1])]]$L
    L_Network_2 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[2],"_",Bipartite_NetworkID[2])]]$L

    if(is.null(config[[id]])){
        ## if id doesn't exist, get emtpy graph
        if(opt$verbose)
            message("STEP: ","No bipartite graph between",Bipartite_NetworkID,"\n")
        Bipartite_Network <- get.bipartite.graph.from.emptylist(pool_nodes_Network_1,pool_nodes_Network_2)
    }else{
        if(opt$verbose)
            message("STEP: ","bipartite graph between",Bipartite_NetworkID,": ", config[[id]], "\n")
        relation <- read.table(config[[id]], sep="\t", header=TRUE,stringsAsFactors = FALSE)
        relation <- relation[relation[,1] %in% pool_nodes_Network_1 & relation[,2] %in% pool_nodes_Network_2, ]
        relation %>% dim
        Bipartite_Network <- get.bipartite.graph.from.edgelist(pool_nodes_Network_1,pool_nodes_Network_2,relation)
    }


    AdjacencyMatrix <- as(as_incidence_matrix(Bipartite_Network, sparse =TRUE),"CsparseMatrix")
    AdjacencyMatrix <- AdjacencyMatrix[match(pool_nodes_Network_1,row.names(AdjacencyMatrix)),match(pool_nodes_Network_2,colnames(AdjacencyMatrix))]
    if(rev){
        AdjacencyMatrix <- t(AdjacencyMatrix)
    }
    ## To make sure the ratio(staying/jumping) = delta, Adjacency_Layers should be col normalized before. If it's not, the ratio depends on the degree of each vertex.
    AdjacencyMatrix <- colNorm(AdjacencyMatrix)
    if(opt$verbose){
        message("The col sum of bipartite AdjacencyMatrix is: \n")
        message(paste0(capture.output(colSums(AdjacencyMatrix) %>% round(.,5) %>% table), collapse = "\n"))
    }

    ## The Matrix should have dimension of (N_Gene * L_Gene) * (N_Disease * L_Disease)
    AdjacencyMatrix %>% dim
    AdjacencyMatrix %>% view

    ## save to graph list
    ## GraphList[[paste(Bipartite_NetworkID,collapse = "_")]] <- Bipartite_Network
    if(returnGraph){
        return(list(AdjacencyMatrix=AdjacencyMatrix,graphs = Bipartite_network))
    }else{
        return(list(AdjacencyMatrix=AdjacencyMatrix))
    }
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 10.- We expand the bipartite graph to fit the dimensions of our multilpex system.
##      From every layer, we can jump to the other subnetwork 
##      (disease similarity network in this case)
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

## universal version of expand bipartite
expand.bipartite.graph.universal <- function(Number_Type_1,Number_Layers_Type_1,Number_Type_2,Number_Layers_Type_2,LayerName_Type_1, LayerName_Type_2, Bipartite_matrix){
    ## The resulting matrix should be with dimension of  (Number_Type_1 * Number_Layers_Type_1) * (Number_Type_2 * Number_Layers_Type_2)

    ## Genrate matrix of  Number_Type_1 * (Number_Type_2 * Number_Layers_Type_2)
    SupraBipartiteMatrix <-  do.call("cbind",rep(list(Bipartite_matrix),Number_Layers_Type_2))

    ## Genrate matrix of  (Number_Type_1 * Number_Layers_Type_1) * (Number_Type_2 * Number_Layers_Type_2)
    SupraBipartiteMatrix <-  do.call("rbind",rep(list(SupraBipartiteMatrix),Number_Layers_Type_1))

    if(Number_Layers_Type_1>1){
        Row_Node_Names <- paste(rep(rownames(Bipartite_matrix),Number_Layers_Type_1),rep(LayerName_Type_1, each = Number_Type_1),sep="::")
    }else{
        Row_Node_Names <- rownames(Bipartite_matrix)
    }
    if(Number_Layers_Type_2>1){
        Col_Node_Names <- paste(rep(colnames(Bipartite_matrix),Number_Layers_Type_2),rep(LayerName_Type_2, each = Number_Type_2),sep="::")
    }else{
        Col_Node_Names <- colnames(Bipartite_matrix)
    }
    rownames(SupraBipartiteMatrix) <- Row_Node_Names
    colnames(SupraBipartiteMatrix) <- Col_Node_Names
    SupraBipartiteMatrix <- as(SupraBipartiteMatrix,"sparseMatrix")
    return(SupraBipartiteMatrix)
}

get.supra.adj.bipartite.by.networkid <- function(Bipartite_NetworkID,AdjacencyMatrixList,config){
    Bipartite_NetworkID <- unlist(Bipartite_NetworkID)

    ## get nodes names for row and column and Layer number 
    pool_nodes_Network_1 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[1],"_",Bipartite_NetworkID[1])]]$pool_nodes
    pool_nodes_Network_2 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[2],"_",Bipartite_NetworkID[2])]]$pool_nodes
    ## get layer number of network 1 and network 2
    L_Network_1 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[1],"_",Bipartite_NetworkID[1])]]$L
    L_Network_2 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[2],"_",Bipartite_NetworkID[2])]]$L

    ## layername of network 1 and 2
    LayerName_Network_1 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[1],"_",Bipartite_NetworkID[1])]]$LayerName
    LayerName_Network_2 <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[2],"_",Bipartite_NetworkID[2])]]$LayerName
    
    ## Sort row and col names
    Bipartite_matrix <- AdjacencyMatrixList[[paste0(Bipartite_NetworkID[1],"_",Bipartite_NetworkID[2])]]$AdjacencyMatrix
    message(paste0(capture.output(dim(Bipartite_matrix)), collapse = "\n"))

    ## We expand the biparite graph to fit the multiplex dimensions.
    ## The biparti matrix has now (N_Gene x N_disease)  dimensions. However we need it to have (N_Gene * L_Gene) x (N_Disease * L_Disease)
    ## The genes in all the layers have to point to the diseases
    N_Type_1 <-  length(pool_nodes_Network_1)
    N_Type_2 <-  length(pool_nodes_Network_2)

    SupraAdjacencyMatrix <- expand.bipartite.graph.universal(N_Type_1,L_Network_1,N_Type_2,L_Network_2,LayerName_Network_1,LayerName_Network_2,Bipartite_matrix)
    message(paste0(capture.output(dim(SupraAdjacencyMatrix)), collapse = "\n"))
    ## normalize to make sure weight of lamda is correct.
    if(opt$verbose){
        message("Before the step, the col sum of SupraAdjacencyMatrix is: \n")
        message(paste0(capture.output(colSums(SupraAdjacencyMatrix) %>% round(.,5) %>% table), collapse = "\n"))
    }
    message("STEP: Column normalize SupraAdjacencyMatrix\n")
    SupraAdjacencyMatrix <- colNorm(SupraAdjacencyMatrix)

    if(opt$verbose){
        message("The col sum of SupraAdjacencyMatrix is: \n")
        message(paste0(capture.output(colSums(SupraAdjacencyMatrix) %>% round(.,5) %>% table), collapse = "\n"))
    }
    return(list(SupraAdjacencyMatrix=SupraAdjacencyMatrix))
}





#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 11.- We compute the transition matrices for the multiplex-heterogeneous system. 
##      Those matrices will generate the final transition matrix.   
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

## 11.1.-Protein-Disease Transition Matrix.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.lambda <- function(config,AdjacencyMatrixList){
    Ls <- lapply(config$LoadedNetworks,function(x){AdjacencyMatrixList[[paste0(x,"_",x)]]$L}) %>% unlist
    names(Ls) <- config$LoadedNetworks

    lambda_mat_list <- list()
    for(i in 1:nrow(config$lambda)){
        for(j in 1:ncol(config$lambda)){
            lambda_mat <- matrix(1,nrow=Ls[rownames(config$lambda)[i]],ncol=Ls[rownames(config$lambda)[j]])
            rownames(lambda_mat) <- AdjacencyMatrixList[[paste0(rownames(config$lambda)[i],"_",rownames(config$lambda)[i])]]$LayerName
            colnames(lambda_mat) <- AdjacencyMatrixList[[paste0(rownames(config$lambda)[j],"_",rownames(config$lambda)[j])]]$LayerName
            lambda_mat <- sweep(lambda_mat,2,colSums(lambda_mat),FUN="/") * config$lambda[i,j]
            tmp <- list(lambda_mat)
            names(tmp) <- paste0(rownames(config$lambda)[i],"_",ncol=rownames(config$lambda)[j])
            lambda_mat_list <- c(lambda_mat_list,tmp)
            
        }
    }
    
    ## lambda_mat_list <- matrix(lambda_mat_list,nrow = nrow(config$lambda), ncol = ncol(config$lambda),byrow =T)
    lambda_mat_list
}

update.lambda <- function(lambda_mat_list,weights){
    allbipartite <- c(combn(names(weights),2,simplify =F),lapply(combn(names(weights),2,simplify =F),rev)) %>% do.call(rbind,.) %>% as.data.table
    colnames(allbipartite) <- c("To","From")
    allbipartite$id <- apply(allbipartite,1,function(x)paste0(x[1],"_",x[2]))
    ## only update bipartite

    lambda_mat_list_updated <- lapply(1:length(lambda_mat_list),function(i){
        ## print(names(lambda_mat_list)[i])
        if(length(grep(names(lambda_mat_list)[i],allbipartite$id))>0){
            lambda <- lambda_mat_list[[i]]
            item <- strsplit(names(lambda_mat_list)[i],"_",perl=T) %>% unlist
            ##
            if(nrow(lambda)==length(weights[[item[[1]]]])){
                if(opt$verbose){
                    message("STEP: add weights to row: ", item[[1]],"\n")
                }
                lambda <- sweep(lambda,1,weights[[item[1]]],"*")
            }
            if(ncol(lambda)==length(weights[[item[[2]]]])){
                if(opt$verbose){
                    message("STEP: add weights to col: ", item[[2]],"\n")
                }
                lambda <- sweep(lambda,2,weights[[item[2]]],"*")
            }
            lambda
        }else{
            lambda_mat_list[[i]]
        }
    })
    names(lambda_mat_list_updated) <- names(lambda_mat_list)
    lambda_mat_list_updated
}

get.transition <- function(delta_mat_list,lambda_mat_list,AdjacencyMatrixList,config){

    
}

## get transition for bipartite matrix
get.transition.universal <-  function(Type_To, Type_From,lambda,SupraAdjacencyMatrixList){
    ## 1. The transition matrix is constucted by this stacking mulplex adj and bipartite adj 
    ## 2. Note that Type_From is on the colnums of a table, so matrix id should be [Type_To]_[Type_From]
    TargetMatrixList <- list()
    if(opt$verbose)
        message("STEP: ","Extract target matrices for ",Type_From,"To",Type_To," \n",sep = "")
    TargetMatrixID <- paste(Type_To,"_",Type_From,sep = "")
    if(length(grep(TargetMatrixID,names(SupraAdjacencyMatrixList)))!=0){
        TargetMatrixList[[names(SupraAdjacencyMatrixList)[grep(TargetMatrixID,names(SupraAdjacencyMatrixList))]]] <- SupraAdjacencyMatrixList[[grep(TargetMatrixID,names(SupraAdjacencyMatrixList))]]$SupraAdjacencyMatrix
    }else{
        stop(TargetMatrixID,"matrix not found\n")
    }

    if(opt$verbose)
        message("STEP: ","Extract matrices starting from ",Type_From," \n",sep = "")
    for(i in which(grepl(paste("_",Type_From,sep = ""),names(SupraAdjacencyMatrixList),perl = T) & !grepl(TargetMatrixID,names(SupraAdjacencyMatrixList)))){
        TargetMatrixList[[names(SupraAdjacencyMatrixList)[i]]] <- SupraAdjacencyMatrixList[[i]]$SupraAdjacencyMatrix
        length(TargetMatrixList)
    }


    ## The column should be the same
    if(opt$verbose)
        message(paste0(capture.output(lapply(TargetMatrixList,dim)), collapse = "\n"))

    if(opt$verbose)
        message(paste0(capture.output(lapply(TargetMatrixList,sum)), collapse = "\n"))

    ## To compute the correct weight for target matrix
    ## Calculate whether there are edges for each gene in the three type of networks
    ## if there is no edge for a gene in a type of network, the node won't jump to that layer, the sum lambda of rest layers should be 1.
    Col_Sum_TargetMatrixList <- lapply(TargetMatrixList,function(x){colSums(x, na.rm = FALSE, dims = 1,sparseResult = TRUE)}) %>%  lapply(., as, "sparseMatrix") %>% do.call("cbind",.)
    colnames(Col_Sum_TargetMatrixList) <- names(TargetMatrixList)
    rownames(Col_Sum_TargetMatrixList) <- colnames(TargetMatrixList[[1]])
    ## make sure no value is great than 1
    Col_Sum_TargetMatrixList[Col_Sum_TargetMatrixList!=0]=1
    ## message(paste0(capture.output(head(Col_Sum_TargetMatrixList)), collapse = "\n"))
    ## diff(Col_Sum_TargetMatrixList@p)

    ## Extract the row of lambda whose row name is equal to Type_to and col name is equal to Type_from
    ##                     LoadedGeneLayers LoadedDiseaseLayers LoadedDrugLayers
    ## LoadedGeneLayers                   1                   1                1
    ## LoadedDiseaseLayers                1                   1                1
    ## LoadedDrugLayers                   1                   1                1

    if(opt$verbose)
        message("STEP: ","Extract jumping probability from",Type_From,"\n")
    lambda_str <- strsplit(colnames(Col_Sum_TargetMatrixList),"_")
    lambda_chosen <- lapply(lambda_str,FUN = function(x){lambda[grep(x[1],rownames(lambda)),grep(x[2],colnames(lambda))]}) %>% unlist
    names(lambda_chosen) <- lapply(lambda_str,function(x)paste(x,collapse = "_"))
    if(opt$verbose)
        message(paste0(capture.output(lambda_chosen), collapse = "\n"))

    ## If the a starting veterx has no edges to anohter vertex in a graph, then the weight of the graph should be zero in the lambda_chosen_corrected. For example if Gene can not jump to Disease but can jump to Gene and Drug, that means the weight of Gene to Disease should be 0.
    lambda_chosen_corrected <- sweep(Col_Sum_TargetMatrixList, 2, lambda_chosen,FUN = "*")
    lambda_chosen_corrected_norm <- sweep(lambda_chosen_corrected,1, rowSums(lambda_chosen_corrected) ,FUN = "/")
    lambda_chosen_corrected_norm[is.na(lambda_chosen_corrected_norm)] <- 0
    ## message(paste0(capture.output(head(lambda_chosen_corrected_norm)), collapse = "\n"))

    ## For vertices with > 0 edges, the row sum is 1; For the vertices with 0 edge, the row sum is 0
    if(opt$verbose)
        message(paste0(capture.output(table(rowSums(lambda_chosen_corrected_norm))), collapse = "\n"))
    
    Transition_Weight <- lambda_chosen_corrected_norm[,grep(TargetMatrixID,colnames(lambda_chosen_corrected_norm))]
    Transition_Weight <- Matrix::Diagonal(x=(Transition_Weight))
    Transition_Network <- TargetMatrixList[[TargetMatrixID]] %*% Transition_Weight
    dimnames(Transition_Network) <- list(row.names(TargetMatrixList[[1]]),colnames(TargetMatrixList[[1]]))
    ## Transition_Network %>% colSums %>% round(.,digits = 10) %>% table

    ## Validation: Column sums of Transition Networks should be equal to normalized lambda
    ## message(paste0(capture.output(all((Transition_Network %>% colSums %>% round(.,digits =5))==round(lambda_chosen_corrected_norm[,1],5))), collapse = "\n"))
    return(Transition_Network)
}

get.transition.multiplex.heterogenous <- function(config,SupraAdjacencyMatrixList){
    AllPossibleNetwork <- as.list(as.data.frame(rbind(rep(config$LoadedNetworks, each = length(config$LoadedNetworks)), config$LoadedNetworks)))
    TransitionMatrixList <- list()
    for(TransitionNetwork in AllPossibleNetwork){
        TransitionNetwork <- unlist(TransitionNetwork)
        TransitionMatrixList[[paste(TransitionNetwork,collapse = "_")]] <- get.transition.universal(TransitionNetwork[1],TransitionNetwork[2],config$lambda,SupraAdjacencyMatrixList)
    }

    lapply(TransitionMatrixList,dim)

    
    ## We generate the global transiction matrix and we return it.

    TransitionMultiplexHeterogeneousMatrix <- Matrix()
    rowMatrixOrder <- lapply(config$LoadedNetworks,function(x){paste(x,config$LoadedNetworks,sep = "_")})
    names(rowMatrixOrder) <- config$LoadedNetworks
    rowMatrix <- lapply(rowMatrixOrder,function(x){do.call(cbind,TransitionMatrixList[unlist(x)])})
    lapply(rowMatrix,dim)
    TransitionMultiplexHeterogeneousMatrix <- do.call(rbind,rowMatrix)

    TransitionMultiplexHeterogeneousMatrix %>% colSums %>% round(.,5) %>% table

    ## correct degrees of each nodes to reduce the effects. fix problem of a disease has too many genes. Also fix the hub genes
    ## row_degree <- rowSums(TransitionMultiplexHeterogeneousMatrix != 0)
    ## row_degree[grep("OMIM:103780",names(row_degree))]
    ## row_degree[grep("OMIM:114500",names(row_degree))]
    ## row_degree[grep("OMIM:181500",names(row_degree))]
    ## sort(row_degree,decreasing=T)[1:100] ## indicate the top diseases are all because of degree

    ##
    ## row_degree_correct <- log(row_degree)/row_degree
    ## sort(row_degree_correct,decreasing=T)[1:100] ## indicate the top diseases are all because of degree
    ## row_degree_correct[grep("OMIM:103780",names(row_degree_correct))]
    ## row_degree_correct[grep("OMIM:181500",names(row_degree_correct))]
    ## row_degree_correct[grep("OMIM:114500",names(row_degree_correct))]
    ## sort(row_degree_correct[grep("OMIM",names(row_degree_correct))],decreasing =T) %>% head
    ## sort(row_degree_correct[grep("OMIM",names(row_degree_correct))],decreasing =T) %>% tail
    ## row_degree[names(sort(row_degree_correct[grep("OMIM",names(row_degree_correct))],decreasing =T))]

    ## row_degree_correct

    ## initialize some parameter need for seed scoring initializaiton
    PoolNodesList <- lapply(SupraAdjacencyMatrixList,function(x){x$pool_nodes})
    PoolNodesList[sapply(PoolNodesList, is.null)] <- NULL
    names(PoolNodesList) <- sapply(names(PoolNodesList),function(x){unique(strsplit(x,"_",perl=T) %>% unlist)})

    LayerName <- lapply(config$LoadedNetworks,function(x){SupraAdjacencyMatrixList[[paste0(x,"_",x)]]$LayerName})
    names(LayerName) <- config$LoadedNetworks

    delta <- lapply(config$LoadedNetworks,function(x){SupraAdjacencyMatrixList[[paste0(x,"_",x)]]$delta})
    names(delta) <- config$LoadedNetworks

    delta_allow <- lapply(config$LoadedNetworks,function(x){SupraAdjacencyMatrixList[[paste0(x,"_",x)]]$delta_allow})
    names(delta_allow) <- config$LoadedNetworks

    tau <- lapply(config$LoadedNetworks,function(x){opt$SupraAdjacencyMatrixList[[paste0(x,"_",x)]]$tau})
    names(tau) <- config$LoadedNetworks

    return(list(TransitionMatrix = TransitionMultiplexHeterogeneousMatrix, PoolNodesList = PoolNodesList,LayerName=LayerName, delta=delta, delta_allow = delta_allow, tau = tau))
}



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 12.- We get the scores for all the seeds involved in the RWR-MH
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 3.- We check if the seed nodes are in our network.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

check.seeds.general <- function(Seeds, PoolNodesList){
    PoolNodesList %>% names

    All_seeds_ok <- lapply(PoolNodesList,function(x){data.frame(Seeds[which(Seeds[[1]] %in% x),])})
    
    All_seeds_ko <- Seeds[which(!Seeds[[1]] %in% do.call(rbind,All_seeds_ok)[,1]),]

    if(opt$verbose){
        message("Seeds OK: ")
        message(paste0(capture.output(All_seeds_ok), collapse = "\n"))
        message("Seeds KO: ")
        message(paste0(capture.output(All_seeds_ko), collapse = "\n"))
    }
    if (nrow(do.call(rbind,All_seeds_ok)) ==0){
        message("Seeds not found in our network")
        return(All_seeds_ok)
    } else {
        return(All_seeds_ok)
    }

}

get.seed.scores.weighted.universal <- function(SeedFileList,eta,tau) {
    Seeds_type_weight <- lapply(SeedFileList,function(x){ifelse(nrow(x)==0,0,1)}) %>% unlist
    eta_corrected <- Seeds_type_weight * eta
    eta_corrected_norm <- eta_corrected/sum(eta_corrected)
    tau_norm <- lapply(tau,function(x)x/sum(x))

    if(opt$verbose){
        message(paste0(capture.output(tau_norm), collapse = "\n"))
        message(paste0(capture.output(eta_corrected_norm), collapse = "\n"))
    }
    ## normalize weight for each gene ; if there is no weight, assign equal weight
    SeedFileList_weight <- lapply(SeedFileList,function(x){
        if(nrow(x)>0){
            if(ncol(x)==2){
                x[,2]=x[,2]/sum(x[,2])
            }else{
                x$weight=1;
                x$weight <- x$weight/sum(x$weight)
            }
            x
        }
    })

    
    ## lapply(SeedFileList,function(x){sum(x[,2])})
    SeedFileList_weight_param <- lapply(1:length(SeedFileList_weight),function(i){
        if(!is.null(SeedFileList_weight[[i]]) && nrow(SeedFileList_weight[[i]])>0){

            ## extract eta for the network type
            eta_corrected_norm_selected <- eta_corrected_norm[grep(names(SeedFileList_weight)[i],names(eta_corrected_norm))]
            SeedFile <- as.data.table(SeedFileList_weight[[i]])
            names(SeedFile)[1] <- "Seeds"

            SeedFile[, eta:=eta_corrected_norm_selected]

            ## extract tau for each layer of the network type
            tau_norm_selected <- tau_norm[grep(names(SeedFileList_weight)[i],names(tau_norm))]
            L <- length(tau_norm_selected[[1]])
            SeedFile_weight <- lapply(1:L,function(j){
                tmp <- copy(SeedFile)[,tau:=tau_norm_selected[[1]][j]]
                if(L>1){
                    tmp[[1]] <- paste(tmp[[1]],"_",j,sep = "")
                }
                tmp
            }) %>% rbindlist
            SeedFile_weight$Score <- SeedFile_weight[[2]] * SeedFile_weight$eta * SeedFile_weight$tau
            SeedFile_weight}
    })
    

    Seeds_Score <- rbindlist(SeedFileList_weight_param)
    Seeds_Score <- as.data.frame(Seeds_Score)
    return(Seeds_Score)
    }


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
                                        # 13.- RANDOM WALK WITH RESTART
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

Random_Walk_Restart <- function(Network_Matrix, r,SeedGenes, Threshold=NULL,max_iter=1000){

### We define the threshold and the number maximum of iterations for the randon walker.
    if(is.null(Threshold)){
        Threshold <- 1e-16
    }
    NetworkSize <- ncol(Network_Matrix)
    
### We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1
    
#### We define the prox_vector(The vector we will move after the first RW iteration. We start from The seed. We have to take in account
#### that the walker with restart in some of the Seed genes, depending on the score we gave in that file).
    prox_vector <- Matrix(0,nrow = NetworkSize,ncol=1)
    rownames(prox_vector) <- rownames(Network_Matrix)
    cn <- str_split(colnames(opt$TransitionMatrix),"::",n=4) %>% do.call(rbind,.) %>% data.table %>% .$V4
    
    if(length(which(cn %in% SeedGenes[,1]))!=nrow(SeedGenes)){
        message("All seeds should be in transition matrix")
                                        #quit("no")
    }
    prox_vector[which(cn %in% SeedGenes[,1])] <- (SeedGenes$Score)
    
    prox_vector  <- prox_vector/sum(prox_vector)
    restart_vector <-  prox_vector

    if(opt$verbose)
        message("r is: ",r,"\n")
    ##    while(iter <= max_iter){
    while(residue >= Threshold && iter <= max_iter){
        ## message("STEP: ","iter is: ",iter,"\n")
        old_prox_vector <- prox_vector
        prox_vector <- (1-r)*(Network_Matrix %*% prox_vector) + r * restart_vector
        residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
        ## message(paste0(capture.output(residue), collapse = "\n"))
        iter <- iter + 1; 
    }
    if(opt$verbose){
        message("STEP: ","final iter is: ",iter-1,"\n")
        message("STEP: ","final residule is: ",residue,"\n")
    }
    return(prox_vector) 
    } 

                                        #a <- as.data.table(SupraAdjacencyMatrixList[["LoadedGeneLayers_LoadedGeneLayers"]]$SupraAdjacencyMatrix[823,],keep.rownames=T)
                                        #a <- cbind(a,strsplit(a$V1,"_") %>% do.call(rbind,.))
                                        #names(a) <- paste0("V",1:4)
                                        #a[V2!=0,]

    ## a <- as.data.table(prox_vector[which(prox_vector!=0,arr.ind=T)[,1],],keep.rownames=T)
    ## a[order(-V2),][!grepl("OMIM|APOE",V1,perl=T),]
    ## a <- cbind(a,strsplit(a$V1,"_") %>% do.call(rbind,.))
    ## names(a) <- paste0("V",1:4)
    ## a <- dcast(a,formula=V3~V4,value.var="V2")

    ## a1 <- as.data.table(prox_vector[which(prox_vector!=0,arr.ind=T)[,1],],keep.rownames=T)
    ## a1 <- cbind(a1,strsplit(a1$V1,"_") %>% do.call(rbind,.))
    ## names(a1) <- paste0("V",1:4)
    ## a1 <- dcast(a1,formula=V3~V4,value.var="V2")




#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 14.- RANKING OF PROTEINS AFTER RANDOM WALK. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

rank_native <- function(PoolNodesList,Results,Seeds){
    ## PoolNodesList <- opt$TransitionMatrixList$PoolNodesList
    ## Results <-  Random_Walk_Results
    ## Seeds <- All_Seeds
    x <- PoolNodesList[[1]]

    Results_id <- str_split(rownames(Results),"::",n=4) %>% do.call(rbind,.) %>% data.table
    Results_dat <- data.table(Results_id,pvalue=as.matrix(Results)[,1])

    Results_Geomatric_Mean <- lapply(PoolNodesList,function(x){
        head(x)
        ## get Layer number for each PoolNodesList
        tmp <- Results_dat[V4 %in% x,] %>% dcast(.,as.formula(paste("V4","~","V3")),value.var = "pvalue" ) %>% as.data.table
        ## sort layers to match original layer name
        ## tmp <- tmp[,c(names(tmp)[1],mixedsort(names(tmp)[-1])),with=F]
        ## format names of columns
        names(tmp)[1] <- "id"
        ## calculate geomMean
        tmp[,geomMean := apply(.SD,1,function(x)exp(mean(log(x[x>0])))),.SDcols = names(tmp)[2:ncol(tmp)]]
        tmp[is.na(geomMean),geomMean:=0]
        ## add seeds indicator
        tmp[get(names(tmp)[1]) %in% Seeds[[1]],seeds :="Yes"]
        tmp[! get(names(tmp)[1]) %in% Seeds[[1]],seeds :="No"]
        ## sort by geomMean
        tmp <- tmp[order(-geomMean),]
        ## add rank
        tmp <- tmp[,rank:=1:.N]
        ## sort columns
        tmp <- tmp[,c(which(grepl("^(id|geomMean|rank|seeds)$",names(tmp),perl=T)),which(!grepl("^(id|geomMean|rank|seeds)$",names(tmp)))),with =F]
        tmp
    })

    return(Results_Geomatric_Mean)
}



####################################################################
#######		New funtions for bootstrapping
#####################################################################

####################################################################
#######		get multiplex graphs
#####################################################################
get_multiplex_network <- function(config){
    layers_list <- lapply(config[["LoadedNetworks"]],function(x){existParameterAndFile(x,config)})
    names(layers_list) <- config[["LoadedNetworks"]]

    Multiplex <- lapply(layers_list,function(layers){
        Layers <- lapply(1:length(layers),function(i){
            if(opt$verbose)
                message("STEP: Read #",i," ",layers[i],"\n")
            tb <- fread(layers[i],header = F)
            if(ncol(tb)==3){
                names(tb)[3] <- "weight"
            }
            Network <- graph.data.frame(tb,directed=FALSE)
            Network <- igraph::simplify(Network, remove.multiple = TRUE, remove.loops = TRUE)
            Network
        })

        ## We get a pool of nodes (Nodes in any of the layers.) 
        ## Some nodes will belong to just some layers, while others will be on all of them.
        pool_nodes <- sort(pool.of.nodes(Layers))
        Layers <- add.missing.nodes(Layers, pool_nodes)

        ## Layers <- lapply(1:length(Layers),function(i){
        ##     V(Layers[[i]])$name <- paste0(V(Layers[[i]])$name,"_",i)
        ##     Layers[[i]]
        ## })

        ## validation adding isolated nodes doesn't affect networks
        ## Layers_full <- add.missing.nodes(Layers, pool_nodes)
        ## identical_graphs(Layers[[1]],Layers_full[[1]])
        ## equal after removeing insolated nodes
        ## identical_graphs(Layers[[1]],delete.vertices(simplify(Layers_full[[1]]), degree(Layers_full[[1]])==0))
        names(Layers) <- names(layers)
        
        ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
        return(Layers)
    })
    ## format names
    names(Multiplex) <- paste0(names(layers_list),"_",names(layers_list))
    return(Multiplex)
}



####################################################################
#######		get bipartite graphs
#####################################################################

get_bipartite_network <- function(config){
    AllPossibleBipartite <- as.list(as.data.frame(combn(config$LoadedNetworks,2,simplify =T)))
    AllPossibleBipartite <- c(AllPossibleBipartite,lapply(AllPossibleBipartite,rev))
    names(AllPossibleBipartite) <- lapply(AllPossibleBipartite,function(x){paste0(x,collapse = "_")})
    
    pool_nodes <- config$pool_nodes
    Bipartite <- lapply(AllPossibleBipartite,function(Bipartite_NetworkID){
        Bipartite_Network <- NULL
        Bipartite_NetworkID <- unlist(Bipartite_NetworkID)
        id <- paste(c("Bipartite",Bipartite_NetworkID),collapse = ".")
        
        if(!is.null(config[[id]])){
            if(opt$verbose)
                message("STEP: ","Reading bipartite graph between ",Bipartite_NetworkID,": ",config[[id]],"\n")

            pool_nodes_Network_1 <- pool_nodes[[paste0(Bipartite_NetworkID[1],"_",Bipartite_NetworkID[1])]]
            pool_nodes_Network_2 <- pool_nodes[[paste0(Bipartite_NetworkID[2],"_",Bipartite_NetworkID[2])]]

            relation <- read.table(config[[id]], sep="\t", header=TRUE,stringsAsFactors = FALSE)
            
            relation <- relation[relation[,1] %in% pool_nodes_Network_1 & relation[,2] %in% pool_nodes_Network_2, ]
                                                                                                    if(opt$verbose)
                message("STEP: ","After filtering, getting ", relation %>% dim %>% .[1]," edges\n")
            
            Bipartite_Network <- get.bipartite.graph.from.edgelist(pool_nodes_Network_1,pool_nodes_Network_2,relation)
        }
        
        ## reverse edge list if oppsite exists.
        Bipartite_NetworkID_rev <- rev(Bipartite_NetworkID)
        id <- paste(c("Bipartite",Bipartite_NetworkID_rev),collapse = ".")
        if(!is.null(config[[id]])){
            if(opt$verbose)
                message("STEP: ","Reading bipartite graph between ",Bipartite_NetworkID_rev,": ",config[[id]],"\n")

            pool_nodes_Network_1 <- pool_nodes[[paste0(Bipartite_NetworkID_rev[1],"_",Bipartite_NetworkID_rev[1])]]
            pool_nodes_Network_2 <- pool_nodes[[paste0(Bipartite_NetworkID_rev[2],"_",Bipartite_NetworkID_rev[2])]]

            relation <- read.table(config[[id]], sep="\t", header=TRUE,stringsAsFactors = FALSE)
            relation <- relation[relation[,1] %in% pool_nodes_Network_1 & relation[,2] %in% pool_nodes_Network_2, ]
            relation %>% dim
            Bipartite_Network <- reverse_edges(get.bipartite.graph.from.edgelist(pool_nodes_Network_1,pool_nodes_Network_2,relation))
            ## reverse type if edge is reversed to make sure the row and col of incidance matrix
            V(Bipartite_Network)$type <- !V(Bipartite_Network)$type 
            
        }
        
        if(is.null(Bipartite_Network)){
            pool_nodes_Network_1 <- pool_nodes[[paste0(Bipartite_NetworkID[1],"_",Bipartite_NetworkID[1])]]
            pool_nodes_Network_2 <- pool_nodes[[paste0(Bipartite_NetworkID[2],"_",Bipartite_NetworkID[2])]]
            
            Bipartite_Network <- get.bipartite.graph.from.emptylist(pool_nodes_Network_1,pool_nodes_Network_2)
        }
        return(Bipartite_Network)
    })
    return(Bipartite)
}


####################################################################
#######		get full network ARRAY matrix. The cells store adjcency matrices
#####################################################################

get_full_network_mat <- function(config,GraphList){
    layers_list <- lapply(config[["LoadedNetworks"]],function(x){names(GraphList[[paste0(x,"_",x)]]);}) %>% unlist

    mat <- matrix(rep(list(),length(unlist(layers_list)) *length(unlist(layers_list))),nrow=length(unlist(layers_list)),ncol=length(unlist(layers_list)))
    colnames(mat) <- unlist(layers_list)
    rownames(mat) <- unlist(layers_list)

    

    AllPossibleComb <- lapply(config$LoadedNetworks,function(x){lapply(config$LoadedNetworks,function(y){c(x,y)}) %>% do.call(rbind,.)}) %>% do.call(rbind,.) %>% split(.,1:nrow(.))

    for(Network in AllPossibleComb){
        id <- paste0(Network,collapse = "_")
        ## fill in multplex network matrix to mat
        if(Network[1]==Network[2]){
            ## test Network <- AllPossibleComb[[1]]
            AllPossibleCombNetwork <- lapply(names(GraphList[[id]]),function(x){lapply(names(GraphList[[id]]),function(y){c(x,y)}) %>% do.call(rbind,.)}) %>% do.call(rbind,.) %>% split(.,1:nrow(.))
            for(x in AllPossibleCombNetwork){
                if(x[1]==x[2]){
                    if(opt$weighted && length(grep("weight",edge_attr_names(GraphList[[id]][[x[1]]])))>0){
                        if(opt$verbose){
                            message("STEP: use weighted network\n")
                        }
                        mat[[x[1],x[2]]] <- GraphList[[id]][[x[1]]] %>% as_adjacency_matrix(.,sparse = TRUE,attr="weight") %>% .[sort(rownames(.),index.return = TRUE)$ix,sort(colnames(.),index.return = TRUE)$ix]
                    }else{
                        if(opt$verbose){
                            message("STEP: use unweighted network\n")
                        }
                        mat[[x[1],x[2]]] <- GraphList[[id]][[x[1]]] %>% as_adjacency_matrix(.,sparse = TRUE) %>% .[sort(rownames(.),index.return=TRUE)$ix,sort(colnames(.),index.return=TRUE)$ix]
                    }
                }else{
                    m <- V(GraphList[[id]][[x[1]]]) 
                    n <- V(GraphList[[id]][[x[2]]]) 
                    Idem_Matrix <- Matrix(0,nrow=length(m),ncol=length(n))
                    diag(Idem_Matrix) <- 1
                    rownames(Idem_Matrix) <- m$name %>% sort
                    colnames(Idem_Matrix) <- n$name %>% sort
                    mat[[x[1],x[2]]] <-  Idem_Matrix 
                }
            }
        }else{
            ## fill in bipartite network
            ## test Network <- AllPossibleComb[[2]]
            AllPossibleCombNetwork <- lapply(names(GraphList[[paste0(Network[[1]],"_",Network[[1]])]]),function(x){lapply(names(GraphList[[paste0(Network[[2]],"_",Network[[2]])]]),function(y){c(x,y)}) %>% do.call(rbind,.)}) %>% do.call(rbind,.) %>% split(.,1:nrow(.))
            for(x in AllPossibleCombNetwork){
                ## test x <- AllPossibleCombNetwork[[1]]
                id <- paste0(Network,collapse = "_")
                bipartite_temp <- GraphList[[id]]
                if(opt$weighted && length(grep("weight",edge_attr_names(bipartite_temp)))>0){
                    if(opt$verbose){
                        message("STEP: use weighted network\n")
                    }
                    ## notice the row and col name of incidence matrix
                    mat[[x[1],x[2]]] <- bipartite_temp %>% as_incidence_matrix(.,sparse = TRUE,attr="weight") %>% .[sort(rownames(.),index.return = TRUE)$ix,sort(colnames(.),index.return=TRUE)$ix]
                }else{
                    if(opt$verbose){
                        message("STEP: use unweighted network\n")
                    }
                    mat[[x[1],x[2]]] <- bipartite_temp %>% as_incidence_matrix(.,sparse = TRUE) %>% .[sort(rownames(.),index.return =TRUE)$ix,sort(colnames(.),index.return=TRUE)$ix]
                }
            }
        }
    }
    return(mat=mat)
}


####################################################################
#######		get delta matrix. Same dimensions as network array matrix
#####################################################################

get_full_delta_mat <- function(config,layers_list){
    delta_mat <- matrix(0,nrow=length(unlist(layers_list)),ncol=length(unlist(layers_list)))
    colnames(delta_mat) <- unlist(layers_list)
    rownames(delta_mat) <- unlist(layers_list)

    AllPossibleComb <- lapply(config$LoadedNetworks,function(x){lapply(config$LoadedNetworks,function(y){c(x,y)}) %>% do.call(rbind,.)}) %>% do.call(rbind,.) %>% split(.,1:nrow(.))

    for(Network in AllPossibleComb){
        ## test code
        ## Network <- AllPossibleComb[[1]]
        ## fill in multplex network matrix to mat
        if(Network[1]==Network[2]){
            delta_multiplex <- get.delta.matrix.by.networkid(Network[1],layers_list,config)
            delta_mat[rownames(delta_multiplex),colnames(delta_multiplex)] <- delta_multiplex
        }else{
            ## test code
            ## Network <- AllPossibleComb[[2]]
            ## fill in bipartite network

            delta_mat[layers_list[[Network[1]]] %>% unlist,layers_list[[Network[2]]] %>% unlist] <- 1/length(layers_list[[Network[1]]] %>% unlist) 
        }
    }
    
    return(delta_mat=delta_mat)
}    

get.delta.matrix.by.networkid <- function(NetworkID,layers_list,config){

    ## 1. The diagnal of the SupraAdjacencyMatrix should be the AdjacencyMatrix for each layer of gene networks
    ## 2. The off diagnal should be the idem matrix
    ## 3. Delta defines the probability of jumping to another layer, so 1-delta defines the probability of staying in a layer.
    ## 4. If a layer only can jump to a subset set of layers (L'), the weight for each layer of subset should be delta/L';
    ## If a layer can jump to any other layers, L' = L-1; See the demonstration in the codes.
    ## delta_allow_mat is generated to describe which layer can jump to any other layers and which layer can only jump a subset of layers.

    L <- lapply(layers_list[[NetworkID]],function(x)length(x)) %>% unlist
    delta <- rep(config[[paste0("Interlayer.Jump.Probability.",NetworkID)]],times=L)
    delta_allow <- rep(config[[paste0("Interlayer.Jump.Allow.",NetworkID)]],times=L)


    ## if there is only 1 layer, the prob staying in the same layer is 1, so the delta is 0
    if(length(L)==1 && L==1){
        delta_allow_mat <- as.matrix(1)
        rownames(delta_allow_mat) <- layers_list[[NetworkID]] %>% unlist
        colnames(delta_allow_mat) <- layers_list[[NetworkID]] %>% unlist
    }else{
        ## to calculate delta weight matrix. Make sure that layers not allowed to jump has 0 
        delta_allow[delta_allow==0]=2
        delta_allow_mat <- delta_allow %o% delta_allow
        delta_allow_mat[delta_allow_mat==2] = 1
        delta_allow_mat[delta_allow_mat==2^2] = 0
        ## e.g. 1 indicates it can jump and 0 means it can not jump
        ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,]    1    1    1    1    1    1    1    1
        ## [2,]    1    1    1    1    1    1    1    1
        ## [3,]    1    1    0    0    0    0    0    0
        ## [4,]    1    1    0    0    0    0    0    0
        ## [5,]    1    1    0    0    0    0    0    0
        ## [6,]    1    1    0    0    0    0    0    0
        ## [7,]    1    1    0    0    0    0    0    0
        ## [8,]    1    1    0    0    0    0    0    0

        ## The diagnal of the delta_allow_mat is defined in the file, which is 1-delta.
        ## To calculate the delta weight for other layers,first mask the diagnal to 0
        diag(delta_allow_mat) <- 0
        ## [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,]    0    1    1    1    1    1    1    1
        ## [2,]    1    0    1    1    1    1    1    1
        ## [3,]    1    1    0    0    0    0    0    0
        ## [4,]    1    1    0    0    0    0    0    0
        ## [5,]    1    1    0    0    0    0    0    0
        ## [6,]    1    1    0    0    0    0    0    0
        ## [7,]    1    1    0    0    0    0    0    0
        ## [8,]    1    1    0    0    0    0    0    0


        ## calculate the column sum and the weight for each layers
        ## Col normalized the matrix
        delta_allow_mat <- sweep(delta_allow_mat,2,colSums(delta_allow_mat),"/")
        ## [,1]      [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,] 0.0000000 0.1428571  0.5  0.5  0.5  0.5  0.5  0.5
        ## [2,] 0.1428571 0.0000000  0.5  0.5  0.5  0.5  0.5  0.5
        ## [3,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [4,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [5,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [6,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [7,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        ## [8,] 0.1428571 0.1428571  0.0  0.0  0.0  0.0  0.0  0.0
        delta_allow_mat[is.na(delta_allow_mat)] <- 0

        ## multiple the delta to get the true weight
        delta_allow_mat <- sweep(delta_allow_mat,2,delta,"*")
        ## [,1]       [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,] 0.00000000 0.07142857 0.25 0.25 0.25 0.25 0.25 0.25
        ## [2,] 0.07142857 0.00000000 0.25 0.25 0.25 0.25 0.25 0.25
        ## [3,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [4,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [5,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [6,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [7,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00
        ## [8,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.00


        ## fill the diagnal with weight of staying in the layer
        diag(delta_allow_mat) <- 1-delta
        ## [,1]       [,2] [,3] [,4] [,5] [,6] [,7] [,8]
        ## [1,] 0.50000000 0.07142857 0.25 0.25 0.25 0.25 0.25 0.25
        ## [2,] 0.07142857 0.50000000 0.25 0.25 0.25 0.25 0.25 0.25
        ## [3,] 0.07142857 0.07142857 0.50 0.00 0.00 0.00 0.00 0.00
        ## [4,] 0.07142857 0.07142857 0.00 0.50 0.00 0.00 0.00 0.00
        ## [5,] 0.07142857 0.07142857 0.00 0.00 0.50 0.00 0.00 0.00
        ## [6,] 0.07142857 0.07142857 0.00 0.00 0.00 0.50 0.00 0.00
        ## [7,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.50 0.00
        ## [8,] 0.07142857 0.07142857 0.00 0.00 0.00 0.00 0.00 0.50

        ## col nomailize again to make sure
        delta_allow_mat <- sweep(delta_allow_mat,2,colSums(delta_allow_mat),"/")
        ## validate the col sum is 1
        if(opt$verbose)
            message(paste0(capture.output(round(colSums(delta_allow_mat),5) %>% table), collapse = "\n"))
    }
    
    if(opt$verbose){
        message("STEP: ","nomalized delta_allow_matrix is :\n")
        message(paste0(capture.output(view(delta_allow_mat,15)), collapse = "\n"))
    }
    rownames(delta_allow_mat) <- lapply(layers_list[[NetworkID]],function(x)x)  %>% unlist
    colnames(delta_allow_mat) <- lapply(layers_list[[NetworkID]],function(x)x)  %>% unlist
    return(delta_allow_mat)
}



####################################################################
#######		get lambda matrix. Same dimensions as network array matrix
#####################################################################

get_full_lambda_mat <- function(config, layers_list){
    lambda_mat <- matrix(0,nrow=length(unlist(layers_list)),ncol=length(unlist(layers_list)))
    colnames(lambda_mat) <- unlist(layers_list)
    rownames(lambda_mat) <- unlist(layers_list)
    
    AllPossibleComb <- lapply(config$LoadedNetworks,function(x){lapply(config$LoadedNetworks,function(y){c(x,y)}) %>% do.call(rbind,.)}) %>% do.call(rbind,.) %>% split(.,1:nrow(.))

    for(Network in AllPossibleComb){
        ## test code
        ## Network <- AllPossibleComb[[1]]
        ## fill in multplex network matrix to mat
        lambda_mat[layers_list[[Network[1]]] %>% unlist,layers_list[[Network[2]]] %>% unlist] <- config$lambda[Network[1],Network[2]]
    }

    return(lambda_mat=lambda_mat)
}



####################################################################
#######		get transition matrix from network array matrix, delta matrix and lambda matrix
#####################################################################

get_all_transition_mat <- function(network_mat,delta_mat,lambda_mat){
    mat_weight <- delta_mat * lambda_mat
    ## viewNoName(mat_weight,100)
    ## col norm adjacency matrix
    mat <- matrix(lapply(1:length(network_mat),function(i){colNorm(network_mat[[i]])}) %>% unlist,nrow=nrow(network_mat),ncol=ncol(network_mat))
    ## validation
    ## lapply(mat,function(x){colSums(x) %>% round(.,5) %>% table})

    ## multiple 
    mat <- matrix(lapply(1:length(mat),function(i){mat[[i]] * mat_weight[i]}) %>% unlist,nrow=nrow(mat),ncol=ncol(mat))
    colnames(mat) <- colnames(network_mat)
    rownames(mat) <- rownames(network_mat)

    ## validation
    ## lapply(mat,function(x){tmp <- colSums(x) %>% round(.,5) %>% unique;ifelse(length(tmp)==2,tmp[tmp!=0],tmp[tmp!=0])}) %>% unlist %>% matrix(.,nrow=8)
    ## should be the same to mat_weight
    ## [,1] [,2] [,3] [,4] [,5] [,6]    [,7]    [,8]
    ## [1,]  0.5  0.1  0.1  0.1  0.1  0.1 0.16667 0.16667
    ## [2,]  0.1  0.5  0.1  0.1  0.1  0.1 0.16667 0.16667
    ## [3,]  0.1  0.1  0.5  0.1  0.1  0.1 0.16667 0.16667
    ## [4,]  0.1  0.1  0.1  0.5  0.1  0.1 0.16667 0.16667
    ## [5,]  0.1  0.1  0.1  0.1  0.5  0.1 0.16667 0.16667
    ## [6,]  0.1  0.1  0.1  0.1  0.1  0.5 0.16667 0.16667
    ## [7,]  1.0  1.0  1.0  1.0  1.0  1.0 1.00000      NA
    ## [8,]  1.0  1.0  1.0  1.0  1.0  1.0      NA 1.00000



    mat <- lapply(1:nrow(mat),function(i){do.call(cbind,mat[i,])}) %>% do.call(rbind,.)
    mat <- colNorm(mat)
    ## validate;should be 0 and 1
    ## colSums(mat) %>% round(.,5) %>% table

    rownames(mat) <- lapply(1:nrow(network_mat),function(i){layername <- rownames(network_mat)[i] %>% gsub(".mat.*","",.,perl =T);paste0(layername,"::",rownames(network_mat[[i,1]]))}) %>% unlist
    colnames(mat) <- lapply(1:nrow(network_mat),function(i){layername <- rownames(network_mat)[i] %>% gsub(".mat.*","",.,perl =T);paste0(layername,"::",rownames(network_mat[[i,1]]))}) %>% unlist

    ## validation
    if(opt$verbose)
        message(capture.output(dim(mat),collapse=TRUE))
    return(mat)
}


####################################################################
#######		sampling disease-gene association
#####################################################################

sample_by_go <- function(gene_set,pro_p,N=20, topN=10){
    ## test: gene_set <- temp[mim_morbid=="OMIM:104300"]$hgnc_symbol
    ## if gene_set <=20, use all genes    
    if(length(gene_set)<=N){
        return(gene_set)
    }
    
    ## get overlaped and nonoverlaped genes with GO networks
    overlaped <- intersect(rownames(pro_p),gene_set)
    nonoverlaped <- gene_set[!gene_set %in% rownames(pro_p)]

    if(length(overlaped)==0){
        return(gene_set)
    }
    ## extract overlaped genes from GO networks
    gene_set_p <- pro_p[overlaped, overlaped] %>% as.data.table(., keep.rownames=T)

    ## make diagnol to be 0
    gene_set_p_mat <- as.matrix(gene_set_p[,-1])
    rownames(gene_set_p_mat) <- gene_set_p$rn
    diag(gene_set_p_mat) <- 0
    ## validation, should all be 0
    ## diag(gene_set_p_mat)

    ## remove gene rows with all zeros 
    gene_set_p_mat <- gene_set_p_mat[which(rowSums(gene_set_p_mat)!=0),which(colSums(gene_set_p_mat)!=0)]
    
    ## if gene_set <=20, use all genes    
    if(nrow(gene_set_p_mat)<=N){
        return(rownames(gene_set_p_mat))
    }

    gene_set_p <- as.data.table(gene_set_p_mat,keep.rownames =TRUE)

    ## filling rows when a gene is absent in go network
    ## if(length(nonoverlaped)>0){
    ##     row_filling <- lapply(nonoverlaped,function(x){y <- apply(gene_set_p[,-1],2,median)}) %>% do.call(rbind,.) %>% data.table(rn=nonoverlaped,.)
    ##     colnames(row_filling) <- colnames(gene_set_p)
    ##     gene_set_p <- rbind(gene_set_p,row_filling)
    ##     ## filling cols when a gene is absent in a go network
    ##     gene_set_p[,(nonoverlaped):=1/nrow(gene_set_p)]
    ## }
    ## validation, should be true
    ## all(colnames(gene_set_p)[-1] == gene_set_p$rn)


    ## extract melted table to plot the distribution
    ## gene_set_p_melt <- melt(gene_set_p,id.vars="rn")
    ## write.table(gene_set_p_melt,file="OMIM:104300",sep = "\t", quote=F,row.names=F)
    ## set.seed(1)

    ## sampling first gene
    gene_sampled <- sample(gene_set_p$rn,1)

    ## sampling the next one based on the sum probability of already sampled genes.

    for(i in 2:N){
        gene_sampled_p <- data.table(rn=gene_set_p$rn,rowsum=0,rank=0,gene_set_p[,gene_sampled,with=F])
        gene_sampled_p <- gene_sampled_p[!rn %in% gene_sampled,]
        if(nrow(gene_sampled_p)==0 || ncol(gene_sampled_p)==0){
            return(gene_sampled)
        }
        gene_sampled_p$rowsum <- rowSums(gene_sampled_p[,c(-1:-3)])
        gene_sampled_p[order(-rowsum),rank:=1:.N]
        ## gene_sampled_p[order(-rowsum),]
        gene_sampled_p_topN <- gene_sampled_p[rank<=topN,]
        gene_sampled_p_topN[,rowsum:=rowsum/sum(rowsum)]
        ## gene_sampled_p_topN[order(-rowsum),]
        ## temp <- sample(gene_sampled_p_topN$rn,10000000,prob=gene_sampled_p_topN$rowsum,replace=TRUE) %>% unlist %>% table
        ## temp/sum(temp)
        tryCatch({
            if(sum(gene_sampled_p_topN$rowsum)==0){
                next
            }
            gene_sampled <- c(gene_sampled,sample(gene_sampled_p_topN$rn,1,prob=gene_sampled_p_topN$rowsum,replace=TRUE))
        },  error = function(error_condition) {
            conditionMessage(error_condition)
            seed <-.Random.seed
            save(seed, gene_sampled_p, gene_sampled, gene_sampled_p_topN,gene_set,file = paste0(Sys.time(),".debug.rda"))
            break
        })        
    }
    
    return(gene_sampled)
}


substitute_bipartite_mat <- function(opt){
    network_mat <- copy(opt$SupraAdjacencyMatrix)
    ## extract disease gene bipartite matrix from the network_mat
    bipartite_disease_gene_mat <- network_mat[grep("LoadedDiseaseLayers",colnames(network_mat)),grep("LoadedGeneLayers",colnames(network_mat))]
    bipartite_disease_gene_edgelist <- lapply(bipartite_disease_gene_mat,function(x){data.table(rn = rownames(x)[x@i+1],cn = colnames(x)[x@j+1],weight=x@x)})
    ## validation, should return a single table
    ## identicalValue <- function(x,y) if (identical(x,y)) x else FALSE
    ## Reduce(f=identicalValue, bipartite_disease_gene_edgelist)
    ## Reduce(f=identicalValue, bipartite_disease_gene_edgelist)$rn %>% table %>% sort

    ## extract gene disease bipartite matrix from the network_mat
    ## bipartite_gene_disease_mat <- network_mat[grep("LoadedGeneLayers",names(colnames(network_mat))),grep("LoadedDiseaseLayers",names(colnames(network_mat)))]
    ## bipartite_gene_disease_edgelist <- lapply(bipartite_gene_disease_mat,function(x){data.table(cn = rownames(x)[x@i+1],rn = colnames(x)[x@j+1],weight=x@x)})
    ## ## validation, should return a single table
    ## identicalValue <- function(x,y) if (identical(x,y)) x else FALSE
    ## Reduce(f=identicalValue, bipartite_gene_disease_edgelist)
    ## Reduce(f=identicalValue, bipartite_gene_disease_edgelist)$rn %>% table %>% sort

    ## The disease-gene and gene-disease bipartite matrices should be the transpose of each other 
    ## validation: should be true
    ## all(Reduce(f=identicalValue, bipartite_gene_disease_edgelist)[,c("rn","cn","weight")] == Reduce(f=identicalValue, bipartite_disease_gene_edgelist)[,c("rn","cn","weight")])

    ## Sample 20 edges for each diseases from the edge list 
    bipartite_disease_gene_edgelist_sampled <- lapply(split(bipartite_disease_gene_edgelist[[1]],bipartite_disease_gene_edgelist[[1]]$rn),function(dis){cbind(unique(dis$rn),sample_by_go(dis$cn,opt$pro_p))}) %>% do.call(rbind,.) %>% data.table(.)
    ## validation: largest should be 20
    ## bipartite_disease_gene_edgelist_sampled  %>% .[[1]] %>% table %>% sort

    ## convert sampled bipartite edgelist to bipartite graph, use the full disease and gene lists.
    V_Type_1 <- bipartite_disease_gene_mat[[1]] %>% rownames
    V_Type_2 <- bipartite_disease_gene_mat[[1]] %>% colnames
    bipartite_disease_gene_sampled <- get.bipartite.graph.from.edgelist(V_Type_1,V_Type_2,bipartite_disease_gene_edgelist_sampled)

    ## get indicance matrix form sampled bipartite graph and sort the row and col names
    bipartite_disease_gene_sampled_mat <- bipartite_disease_gene_sampled %>% as_incidence_matrix(.,sparse = TRUE) %>% .[sort(rownames(.),index.return = TRUE)$ix,sort(colnames(.),index.return=TRUE)$ix]

    ## validation, largest value should be 20
    if(opt$verbose)
        message(paste0(capture.output(bipartite_disease_gene_sampled_mat  %>% rowSums %>% table %>% sort), collapse = "\n"))


    ## replace the disease-gene and gene-disease entries with sampled bipartite incidence matrix
    for( i in grep("LoadedDiseaseLayers",colnames(network_mat))){
        for( j in grep("LoadedGeneLayers",colnames(network_mat))){
            network_mat[[i,j]] <- bipartite_disease_gene_sampled_mat
        } 
    }

    for( j in grep("LoadedDiseaseLayers",colnames(network_mat))){
        for( i in grep("LoadedGeneLayers",colnames(network_mat))){
            network_mat[[i,j]] <- t(bipartite_disease_gene_sampled_mat)
        } 
    }
    return(network_mat)
}


substitute_bipartite_mat_randomSampling <- function(opt){
    network_mat <- copy(opt$SupraAdjacencyMatrix)
    ## extract disease gene bipartite matrix from the network_mat
    bipartite_disease_gene_mat <- network_mat[grep("LoadedDiseaseLayers",colnames(network_mat)),grep("LoadedGeneLayers",colnames(network_mat))]
    bipartite_disease_gene_edgelist <- lapply(bipartite_disease_gene_mat,function(x){data.table(rn = rownames(x)[x@i+1],cn = colnames(x)[x@j+1],weight=x@x)})
    ## validation, should return a single table
    ## identicalValue <- function(x,y) if (identical(x,y)) x else FALSE
    ## Reduce(f=identicalValue, bipartite_disease_gene_edgelist)
    ## Reduce(f=identicalValue, bipartite_disease_gene_edgelist)$rn %>% table %>% sort

    ## extract gene disease bipartite matrix from the network_mat
    ## bipartite_gene_disease_mat <- network_mat[grep("LoadedGeneLayers",names(colnames(network_mat))),grep("LoadedDiseaseLayers",names(colnames(network_mat)))]
    ## bipartite_gene_disease_edgelist <- lapply(bipartite_gene_disease_mat,function(x){data.table(cn = rownames(x)[x@i+1],rn = colnames(x)[x@j+1],weight=x@x)})
    ## ## validation, should return a single table
    ## identicalValue <- function(x,y) if (identical(x,y)) x else FALSE
    ## Reduce(f=identicalValue, bipartite_gene_disease_edgelist)
    ## Reduce(f=identicalValue, bipartite_gene_disease_edgelist)$rn %>% table %>% sort

    ## The disease-gene and gene-disease bipartite matrices should be the transpose of each other 
    ## validation: should be true
    ## all(Reduce(f=identicalValue, bipartite_gene_disease_edgelist)[,c("rn","cn","weight")] == Reduce(f=identicalValue, bipartite_disease_gene_edgelist)[,c("rn","cn","weight")])

    ## Sample 20 edges for each diseases from the edge list 
    bipartite_disease_gene_edgelist_sampled <- lapply(split(bipartite_disease_gene_edgelist[[1]],bipartite_disease_gene_edgelist[[1]]$rn),function(dis){cbind(unique(dis$rn),sample(dis$cn,size=min(length(dis$cn),20)))}) %>% do.call(rbind,.) %>% data.table(.)
    ## validation: largest should be 20
    ## bipartite_disease_gene_edgelist_sampled  %>% .[[1]] %>% table %>% sort

    ## convert sampled bipartite edgelist to bipartite graph, use the full disease and gene lists.
    V_Type_1 <- bipartite_disease_gene_mat[[1]] %>% rownames
    V_Type_2 <- bipartite_disease_gene_mat[[1]] %>% colnames
    bipartite_disease_gene_sampled <- get.bipartite.graph.from.edgelist(V_Type_1,V_Type_2,bipartite_disease_gene_edgelist_sampled)

    ## get indicance matrix form sampled bipartite graph and sort the row and col names
    bipartite_disease_gene_sampled_mat <- bipartite_disease_gene_sampled %>% as_incidence_matrix(.,sparse = TRUE) %>% .[sort(rownames(.),index.return = TRUE)$ix,sort(colnames(.),index.return=TRUE)$ix]

    ## validation, largest value should be 20
    if(opt$verbose)
        message(paste0(capture.output(bipartite_disease_gene_sampled_mat  %>% rowSums %>% table %>% sort), collapse = "\n"))


    ## replace the disease-gene and gene-disease entries with sampled bipartite incidence matrix
    for( i in grep("LoadedDiseaseLayers",colnames(network_mat))){
        for( j in grep("LoadedGeneLayers",colnames(network_mat))){
            network_mat[[i,j]] <- bipartite_disease_gene_sampled_mat
        } 
    }

    for( j in grep("LoadedDiseaseLayers",colnames(network_mat))){
        for( i in grep("LoadedGeneLayers",colnames(network_mat))){
            network_mat[[i,j]] <- t(bipartite_disease_gene_sampled_mat)
        } 
    }
    return(network_mat)
}
