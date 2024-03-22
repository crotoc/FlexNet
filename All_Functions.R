### Functions called from RWR-M.R

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# GENERAL FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## 2.- We check the input parameters
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
existParameterAndFile <- function(x,config){
    if(is.null(config[[x]])){
        if(opt$verbose)
            message("STEP: ",paste("This loaded layer", x,"IS NOT in config files"))
    }else{
        if(file.exists(config[[x]]) && !dir.exists(config[[x]])){
            message("STEP: ",paste("LoadedLayers",x,"IS FILE at:",paste(config[[x]],collapse = ", ")))
            a <- config[[x]]
            names(a) <- x
            a
        }else if(file.exists(config[[x]]) && dir.exists(config[[x]])) {
            if(opt$verbose)
                message("STEP: ",paste("LoadedLayers",x,"IS DIR at:",paste(config[[x]],collapse = ", "),"\n"))
            file_list <- list.files(config[[x]],full.name = T)
            file_list <- file_list[grep(".gr$",file_list,perl = T)]
            if(length(file_list)>1){
                names(file_list) <- basename(file_list)
            }else{
                names(file_list) <- x
            }
            file_list
        }else{
            if(opt$verbose)
                message("STEP: ",paste("The file for the loaded layer", x," IS NOT at",config[[x]]),"\n")
        }
    }
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
    AdjacencyMatrix@x <- AdjacencyMatrix@x / rep.int(ifelse(colSums(AdjacencyMatrix)==0,1,colSums(AdjacencyMatrix)), diff(AdjacencyMatrix@p))
    AdjacencyMatrix
}

rowNorm <- function(AdjacencyMatrix){
    AdjacencyMatrix <- t(AdjacencyMatrix)
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
        message("STEP: ","Reading bipartite graph between",Bipartite_NetworkID,"\n")
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
    if(length(which(colnames(Network_Matrix) %in% SeedGenes[,1]))!=nrow(SeedGenes)){
        message("All seeds should be in transition matrix")
        quit("no")
    }
    prox_vector[which(colnames(Network_Matrix) %in% SeedGenes[,1])] <- (SeedGenes$Score)
    
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
                                        # 14.- RANKING OF PROTEINS AFTER RANDOM WALK. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

rank_native <- function(PoolNodesList,Results,Seeds){
    ## PoolNodesList <- opt$TransitionMatrixList$PoolNodesList
    ## Results <-  Random_Walk_Results
    ## Seeds <- All_Seeds
    x <- PoolNodesList[[1]]
    
    Results_id <- gsub("(.*?)::(.+)","\\1\t\\2",rownames(Results),perl = T) %>% strsplit(.,"\t")
    Results_id <- lapply(Results_id,function(x){
        if(length(x)==1){x = c(x,1)}else{x}})
    Results_dat <- data.table(do.call(rbind,Results_id) ,pvalue=as.matrix(Results)[,1])

    Results_Geomatric_Mean <- lapply(PoolNodesList,function(x){
        head(x)
        ## get Layer number for each PoolNodesList
        tmp <- Results_dat[get(names(Results_dat)[1]) %in% x,] %>% dcast(.,as.formula(paste(names(Results_dat)[1],"~",names(Results_dat)[2])),value.var = "pvalue" ) %>% as.data.table
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



