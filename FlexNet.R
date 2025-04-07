suppressMessages(require(optparse))
if(!suppressWarnings(suppressMessages(require("myutils",character.only=TRUE)))==FALSE){
    suppressWarnings(suppressMessages(require(myutils)))
}else{
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/addGeom.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/formatPvalue.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/opendev.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/options.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/plotPLabelFormat.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/prettyPaginate.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/using.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/VAlignPlots.R" ,sep = ""))
    source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/mypkg/myutils/R/view.R" ,sep = ""))
}

cmdArgs <- commandArgs(TRUE)
suppressMessages(require(argparse))
parser <- ArgumentParser()
parser <- addArg(parser,c("General"))
## parser <- addArg(parser,c("Rmd"))
opt <- parser$parse_known_args()[[1]]
opt$cmdArgs <- cmdArgs

suppressMessages(require(ggplot2))
suppressMessages(require(dplyr))
suppressMessages(require(magrittr))
suppressMessages(require(reshape2))
suppressMessages(require(data.table))
suppressMessages(require(gtools))
suppressMessages(require(RColorBrewer))
suppressMessages(require(extrafont))
suppressMessages(require(gridExtra))
suppressMessages(require(grid))
suppressMessages(require(ggdendro))
suppressMessages(require(WGCNA))
suppressMessages(require(superheat))
suppressMessages(require(dendextend))
suppressMessages(require("ggpubr"))
suppressMessages(require("readr"))
suppressMessages(require("gtools"))
suppressMessages(using("doParallel"))
suppressMessages(using("parallel"))
suppressMessages(using("doRNG"))
require("igraph")
require("Matrix")
require("dnet")
require("supraHex")
require("yaml")
require("MatrixExtra")
require("stringr")
require("doFuture")
options(future.rng.onMisuse = "ignore")
suppressMessages(require("RhpcBLASctl"))
blas_set_num_threads(opt$threads)
omp_set_num_threads(opt$threads)

## loadfonts()
## theme_figure <- theme_classic() + theme(panel.background = element_rect(fill = 'transparent', colour = NA),legend.title=element_blank(),text=element_text(family=opt$fontfamily),axis.title.x = element_blank(),plot.title = element_text(hjust=0.5, face="bold",size=12),axis.line = element_blank())

if(!dir.exists(opt$dir_out)) dir.create(opt$dir_out,showWarnings=T,recursive=T)

print(opt$dir_out)


if(opt$verbose){
    cat("\n")
}

opt$output <- normalizeOutput(opt)
if(!dir.exists(dirname(opt$output))) dir.create(dirname(opt$output),showWarnings=T,recursive=T)


####################################################################
#######		script specific opt
#####################################################################
parser_script <- ArgumentParser()
parser_script$add_argument("--saveIntermediate", default="NULL",dest="saveIntermediate",help = "saveIntermediate is [default \"%(default)s\"]")
parser_script$add_argument("--weighted", action="store_true",dest="weighted", default=FALSE,help="weighted is [default %(default)s]")
parser_script$add_argument("--max_iter", type="integer", default=1000,help="max_iter [default %(default)s]", metavar="number")
parser_script$add_argument("--Seed.File.List", default="NULL",dest="Seed.File.List",help = "Seed.File.List is [default \"%(default)s\"]")
parser_script$add_argument("--Seed.File", default="NULL",dest="Seed.File",help = "Seed.File is [default \"%(default)s\"]")
parser_script$add_argument("--go", default="/fs0/chenr6/Database_fs0/Gibbs/network/go_propogation_probality_rp_0.3.RData",dest="go",help = "go is [default \"%(default)s\"]")
parser_script$add_argument("--specificity", default="/fs0/chenr6/project/netRepurpose/specificity/celltype.tissue.specificity",dest="specificity",help = "specificity is [default \"%(default)s\"]")
parser_script$add_argument("--cmd",nargs="+", default="1",dest="cmd",help = "cmd is [default \"%(default)s\"]")
parser_script$add_argument("--singleSamplingMode", action="store_true",dest="singleSamplingMode", default=FALSE,help="singleSamplingMode is [default %(default)s]")
parser_script$add_argument("--threadsSampling", type="integer", default=1, help="threadsSampling [default %(default)s]")

if(!is.null(opt$cmdArgs)){
    opt_script <- addArg(parser_script,c("General"))$parse_known_args(opt$cmdArgs)[[1]]
}else{
    opt_script <- addArg(parser_script,c("General"))$parse_known_args()[[1]]
}


for(opt_name in names(opt_script)){
    if(! opt_name %in% names(opt)){
        opt[opt_name] = opt_script[opt_name]
    }
}


if(opt$test){
    save.image(file = paste(opt$output,".rda",sep = ""));
    message(paste(opt$output,".rda",sep = ""))
    stop("test rda saved...")
}


####################################################################
#######		functions
#####################################################################

####################################################################
#######		generate AdjacencyMatrixList from .gr network file
#####################################################################
### explanation of parameters of model

## | Parameters | explain                                                |
## |------------|--------------------------------------------------------|
## | r          | is restarting probability of global                    |
## | delta      | possibility of jumping to a homolog layer in multiplex |
## | tau        | importance of each layer in multiplex                  |
## | lambda     | possibility of jumping to other types of networks      |
## | eta        | importance of each type of networks                    |

####################################################################
#######		Reading config file
#####################################################################

## opt$yamlx <- '/fs0/chenr6/project/netRepurpose/config.v2.yaml.test'
## opt$yaml <- '/fs0/chenr6/project/netRepurpose/script_test/config.v2.yaml.test'
## opt$yaml <- '/fs0/chenr6/project/netRepurpose/script_test/twoLayers.mix.withWeight.yaml'
## opt$yaml <- '/fs0/chenr6/project/netRepurpose/script_test/threeLayers.drugOnly.noWeight.yaml'
## opt$yaml <- '/fs0/chenr6/project/netRepurpose/perDrug/rwrmh_v2_20231204/dir_drug_1000/config.yaml'

## load("/panfs/accrepfs.vampire/nobackup/cgg/chenr6/project/netRepurpose/perDrug/dir_debugHeatmapProblem/dir_debug.oneSC/1.rda")
## opt$yaml <- '/fs0/chenr6/chenr6/myscript/Rscripts/rmd/netRepurpose/exampleYAML/debug/debug.FlexNet.yaml'



generate_GraphList <- function(opt,save=TRUE){
    config <- opt$config

    if(is.null(config$saveGraphList)){
        config$saveGraphList <- paste0(opt$saveIntermediate,".graph.","rda")
    }

    if(file.exists(config$saveGraphList)){
        if(opt$verbose)
            message("\nSTEP: ",config$saveGraphList, " exists. Skip generating and loading directly")
        if(length(grep("generate_SupraAdjacencyMatrix",opt$cmd))==0){
            return(opt)
        }else{
            load(config$saveGraphList)
        }
    }else{
        ## Add some network info to the save name.
        if(opt$verbose)
            message("\nSTEP: ",config$saveGraphList, " doesn't exist. Change name by adding network information. Generating...")

        ## Read the layers of our multiplex networks

        GraphList <- list()
        message("STEP: ","Reading Multiplex Layers...")
        GraphList <- c(GraphList,get_multiplex_network(config))

        ## Read the layers of our bipartite networks
        config$pool_nodes <- lapply(GraphList,function(x){pool.of.nodes(x) %>% sort})

        GraphList <- c(GraphList,get_bipartite_network(config))
        
        if(opt$verbose)
            message("STEP: ","save GraphList to ",config$saveGraphList, "...")
        if(save){
            save(GraphList,file = config$saveGraphList)
        }
    }
    opt$GraphList <- GraphList
    gc()
    if(opt$verbose)
        message("STEP: generating GraphList done.\n")
    opt
}


#####################################################################
#######	generate SupraAdjacencyMatrix if not exists
#####################################################################


generate_SupraAdjacencyMatrix <- function(opt,save=TRUE){
    config <- opt$config

    if(is.null(config$saveGraphList)){
        config$saveGraphList <- paste0(opt$saveIntermediate,".graph.rda")
    }


    if(is.null(opt$GraphList) && file.exists(config$saveGraphList)){
        load(config$saveGraphList)
        opt$GraphList <- GraphList
        rm(GraphList)
        gc()
    }

    
    if(is.null(config$saveSupraAdjacencyMatrix)){
        config$saveSupraAdjacencyMatrix <- paste0(opt$saveIntermediate,".supra.rda")
    }

    if(file.exists(config$saveSupraAdjacencyMatrix)){
        if(opt$verbose)
            message("\nSTEP: ",config$saveSupraAdjacencyMatrix, " exists. Skip generating and loading directly")
        if(length(grep("generate_TransitionMatrix",opt$cmd))==0){
            return(opt)
        }else{
            load(config$saveSupraAdjacencyMatrix)
        }
    }else{
        ## Add some network info to the save name.
        if(opt$verbose)
            message("\nSTEP: ",config$saveSupraAdjacencyMatrix, " doesn't exist. Change name by adding network information. Generating...")

        
        SupraAdjacencyMatrix <- get_full_network_mat(config,opt$GraphList)
        
        if(opt$verbose)
            message("STEP: ","save SupraAdjacencyMatrix to ",config$saveSupraAdjacencyMatrix, "...")
        if(save){
            save(SupraAdjacencyMatrix,file = config$saveSupraAdjacencyMatrix)
        }
    }
    opt$SupraAdjacencyMatrix <-  SupraAdjacencyMatrix
    if(opt$verbose)
        message("STEP: generating SupraAdjacencyMatrix done.\n")
    opt$GraphList <- NULL
    gc()
    opt
}


####################################################################
#######		generate transition matrix
#####################################################################
generate_TransitionMatrix  <- function(opt,save=TRUE){
    config <- opt$config

    if(is.null(config$saveSupraAdjacencyMatrix)){
        config$saveSupraAdjacencyMatrix <- paste0(opt$saveIntermediate,".supra.rda")
    }

    if(is.null(opt$SupraAdjacencyMatrix) && file.exists(config$saveSupraAdjacencyMatrix)){
        load(config$saveSupraAdjacencyMatrix)
        opt$SupraAdjacencyMatrix <- SupraAdjacencyMatrix
        rm(SupraAdjacencyMatrix)
        gc()
    }

    if(is.null(config$saveTransitionMatrix)){
        if(is.null(opt$i)){
            opt$i <-0
        }
        config$saveTransitionMatrix <- paste0(opt$saveIntermediate,".transition.",opt$i,".rda")
    }


    if(file.exists(config$saveTransitionMatrix)){
        if(opt$verbose)
            message("\nSTEP: ",config$saveTransitionMatrix, " exists. Skip generating and loading directly")
        if(length(grep("FlexNet",opt$cmd,perl=T))==0){
            return(opt)
        }else{
            load(config$saveTransitionMatrix)
        }
    }else{
        ## Add some network info to the save name.
        if(opt$verbose)
            message("\nSTEP: ",config$saveTransitionMatrix, " doesn't exist. Change name by adding network information. Generating...")

        if(opt$verbose)
            message("STEP: ","Transforming from adjacency matrices to transition matrix...\n")

        ## get layers from SupraAdjancency Matrix
        rn <- str_split(rownames(opt$SupraAdjacencyMatrix),"::",n=3) %>% do.call(rbind,.) %>% data.table %>% .[,V4:=rownames(opt$SupraAdjacencyMatrix)]        
        layers_list <-lapply(unique(rn$V1),function(x){tmp <- lapply(unique(rn[V1==x,]$V2),function(y){rn[V1==x & V2==y,]$V4});tmp;names(tmp) <- unique(rn[V1==x,]$V2);tmp})
        names(layers_list) <- unique(rn$V1)
        

        ## sort lists
        delta_mat <- get_full_delta_mat(config,layers_list)
        lambda_mat <- get_full_lambda_mat(config,layers_list)

        TransitionMatrix <-get_all_transition_mat(opt$SupraAdjacencyMatrix,delta_mat,lambda_mat)

        if(opt$verbose)
            message("STEP: ","save TransitionMatrix to ",config$saveTransitionMatrix, "...")
        if(save){
            save(TransitionMatrix,file = config$saveTransitionMatrix)
        }    
    }
    opt$TransitionMatrix <- TransitionMatrix

    if(opt$verbose)
        message("STEP: generating TransitionMatrix done.\n")

    opt$SupraAdjacencyMatrix <- NULL
    gc()
    opt

}

####################################################################
#######		integrate specifisity
#####################################################################
update_TransitionMatrixSpecificity <- function(opt){
    celltype_tissue_specificity <- fread(opt$specificity)

    transitionmatrix_row <-  str_split(rownames(opt$TransitionMatrix),"::",n=4) %>% do.call(rbind,.) %>% data.table(.,rn=rownames(opt$TransitionMatrix))
    transitionmatrix_row <- strsplit(transitionmatrix_row$V3,"\\.",perl =T) %>% lapply(.,function(x){if(length(x)==1){c(x,x)}else{x}}) %>% do.call(rbind,.) %>% data.table(transitionmatrix_row,.)
    dim(transitionmatrix_row)
    length(opt$TransitionMatrix %>% row.names)
    names(transitionmatrix_row) <- c("NetworkType","NetworkID_1","NetworkID_2","gene","rn","tissue","cell_type")
    transitionmatrix_row[,c("tissue","cell_type")] %>% unique

    transitionmatrix_specificity <- merge(transitionmatrix_row,celltype_tissue_specificity,by.x=c("gene","tissue","cell_type"),by.y=c("gene","variable","cell_type"),all.x=T,sort=F)
    transitionmatrix_specificity[,specificity:=cell_specificity * tissue_specificity]
    ## drug and disease has no specificity
    transitionmatrix_specificity[NetworkType %in% c("LoadedDiseaseLayers","LoadedDrugLayers"),specificity:=1]
    transitionmatrix_specificity[NetworkType %in% c("LoadedDiseaseLayers","LoadedDrugLayers"),]$specificity %>% table
    ## non scNetwork should have no specificity
    transitionmatrix_specificity[NetworkType %in% c("LoadedGeneLayers") & ! NetworkID_1 %in% "scNetwork",specificity:=1]
    transitionmatrix_specificity[NetworkType %in% c("LoadedGeneLayers") & ! NetworkID_1 %in% "scNetwork",]$specificity %>% table
    ## the rest NAs should be the genes don't have specificity value, assign to 0
    transitionmatrix_specificity[is.na(specificity),]$cell_type %>% unique
    transitionmatrix_specificity[is.na(specificity),]$tissue %>% unique
    transitionmatrix_specificity[is.na(specificity),specificity:=0]

    ## no NA
    transitionmatrix_specificity[is.na(specificity),]

    ## limit specificity value to quantile> 0.99
    transitionmatrix_specificity[specificity!=0]$specificity %>% quantile
    transitionmatrix_specificity[specificity>quantile(transitionmatrix_specificity$specificity[transitionmatrix_specificity$specificity>0],0.99),specificity:=quantile(transitionmatrix_specificity$specificity[transitionmatrix_specificity$specificity>0],0.99)]
    transitionmatrix_specificity[specificity!=0]$specificity %>% quantile

    opt$TransitionMatrix %>% dim
    transitionmatrix_specificity %>% dim
    all(rownames(opt$TransitionMatrix) ==  transitionmatrix_specificity$rn)

    opt$TransitionMatrix  <- Diagonal(length(transitionmatrix_specificity$specificity), x = transitionmatrix_specificity$specificity) %*% opt$TransitionMatrix
    ## temp  <- Diagonal(length(transitionmatrix_specificity$specificity), x = transitionmatrix_specificity$specificity) %*% opt$TransitionMatrix

    ## temp1 <- colNorm(temp)
    ## validation
    ## row.names(opt$TransitionMatrix)[grep("APOE",row.names(opt$TransitionMatrix),perl =T)]

    ## apoe <- grep("APOE",row.names(opt$TransitionMatrix),perl =T)
    ## tmpv <- opt$TransitionMatrix[,apoe] * transitionmatrix_specificity$specificity
    ## tmpv[tmpv!=0 &!is.na(tmpv)]
    ## apply(tmpv,2,function(x){x[x!=0 &!is.na(x)]})

    ## all(tmp[,apoe]==tmpv,na.rm=T)
     opt$TransitionMatrix <- colNorm(opt$TransitionMatrix)

     ##validate no missing value
     colSums(opt$TransitionMatrix) %>% round(.,5) %>% table

    rownames(opt$TransitionMatrix) <- transitionmatrix_specificity$rn
    opt
}

## ####################################################################
## ####### bootstrapping by sampling gene-disease association based on
## ####### GO networks
## #####################################################################
## load("/fs0/chenr6/Database_fs0/Gibbs/network/go_propogation_probality_rp_0.3.RData")
## opt <- generate_SupraAdjacencyMatrix(opt)
## opt$pro_p <- pro_p
## rm(pro_p)
## gc()
## set.seed(1)
## rm(.Random.seed, envir=globalenv())
sampleGeneDiseaseAssocation <- function(opt){
    load(opt$go)
    opt <- generate_SupraAdjacencyMatrix(opt)
    opt$pro_p <- pro_p
    rm(pro_p)
    gc()

    if(opt$verbose)
        message("STEP: Read GO RWR")
    
    if(opt$verbose){
        message("STEP: read SupraAdjancencyMatrix")
    }

    config <- opt$config

    if(is.null(config$saveSupraAdjacencyMatrix)){
        config$saveSupraAdjacencyMatrix <- paste0(opt$saveIntermediate,".supra.rda")
    }

    if(is.null(opt$SupraAdjacencyMatrix) && file.exists(config$saveSupraAdjacencyMatrix)){
        load(config$saveSupraAdjacencyMatrix)
        opt$SupraAdjacencyMatrix <- SupraAdjacencyMatrix
        rm(SupraAdjacencyMatrix)
        gc()
    }

    ## singleSamplingMode is used for generate a single transition rda using sampling
    ## FALSE indicates to generate a batch number. Use FALSE to speed up because it only load GO network and SupraAdjacencyMatrix once;
    if(opt$singleSamplingMode){
        Is = opt$number
    }else{
        Is = 1:opt$number
    }

    substitute <- function(opt,i){
        opt$i <- i
        opt$config$saveTransitionMatrix <- paste0(opt$saveIntermediate,".transition.",opt$i,".rda")        
        if(!file.exists(opt$config$saveTransitionMatrix)){
            if(opt$verbose){
                message("STEP: Sampling...",i, "")
            }
            opt$SupraAdjacencyMatrix <- substitute_bipartite_mat(opt)
        }
        if(opt$verbose){
            message("STEP: update transition")
        }
        opt <- generate_TransitionMatrix(opt)
        if(opt$verbose){
            message("STEP: FlexNet")
        }
        ## opt <- FlexNet(opt)
        rm(opt)
        gc()
    }
    
    options(future.globals.maxSize=220 * 1024*1024*1024) ##30G
    plan(multicore, workers = opt$threadsSampling)

    ## registerDoParallel(opt$threadsSampling)
    foreach(i=Is) %dofuture% { 
        substitute(opt, i)
    }
}


sampleRandomGeneDiseaseAssocation <- function(opt){
    if(opt$verbose){
        message("STEP: read SupraAdjancencyMatrix")
    }

    config <- opt$config

    if(is.null(config$saveSupraAdjacencyMatrix)){
        config$saveSupraAdjacencyMatrix <- paste0(opt$saveIntermediate,".supra.rda")
    }

    if(is.null(opt$SupraAdjacencyMatrix) && file.exists(config$saveSupraAdjacencyMatrix)){
        load(config$saveSupraAdjacencyMatrix)
        opt$SupraAdjacencyMatrix <- SupraAdjacencyMatrix
        rm(SupraAdjacencyMatrix)
        gc()
    }

    
    options(future.globals.maxSize=200 * 1024*1024*1024) ##30G
    plan(multicore, workers = opt$threadsSampling)

    ## singleSamplingMode is used for generate a single transition rda using sampling
    ## FALSE indicates to generate a batch number. Use FALSE to speed up because it only load GO network and SupraAdjacencyMatrix once;
    if(opt$singleSamplingMode){
        Is = opt$number
    }else{
        Is = 1:opt$number
    }

    substitute <- function(opt, i){
        opt$i <- i
        opt$config$saveTransitionMatrix <- paste0(opt$saveIntermediate,".transition.",opt$i,".rda")        
        if(!file.exists(opt$config$saveTransitionMatrix)){
            if(opt$verbose){
                message("STEP: Sampling... ",i, "")
            }
            opt$SupraAdjacencyMatrix <- substitute_bipartite_mat_randomSampling(opt)
        }
        if(opt$verbose){
            message("STEP: update transition")
        }
        opt <- generate_TransitionMatrix(opt)
        if(opt$verbose){
            message("STEP: FlexNet")
        }
        ## opt <- FlexNet(opt)
        rm(opt)
        gc()

    }
    registerDoParallel(opt$threadsSampling)
    foreach(i=Is) %dorng% {
        substitute(opt,i)
    }
}



## ####################################################################
## ####### WE PREPARE THE SEEDS SCORES AND WE PERFORM THE RWR-MH. WE WRITE
## ####### THE OUTPUT FILE CONTAINING HGNC GENE SYMBOL AND THEIR RESULTING SCORES. 		
## #####################################################################

### Now, that we have the diseases and the proteins of our system we can check the nodes.
FlexNet <- function(opt){
    config <- opt$config

    if(is.null(config$saveTransitionMatrix)){
        config$saveTransitionMatrix <- paste0(opt$saveIntermediate,".transition.",opt$i,".rda")
    }

    if(is.null(opt$TransitionMatrix) && file.exists(config$saveTransitionMatrix)){
        if(opt$verbose){
            message("STEP: ","load FlexNet...", config$saveTransitionMatrix)
        }
        load(config$saveTransitionMatrix)
        opt$TransitionMatrix <- TransitionMatrix
        rm(TransitionMatrix)
        gc()
    }

    if(opt$verbose){
        message("STEP: ","updating transition matrix using specificity")
    }

    opt <- update_TransitionMatrixSpecificity(opt)
    gc()

    if(opt$verbose){
        message("STEP: ","Start FlexNet...")
        message("STEP: ","Checking Seed nodes...")}



    if(is.null(config$Seed.File.List) && file.exists(config$Seed.File)){
        SeedFileList=config$Seed.File
    }else{
        if(opt$verbose)
            message("STEP: ","Run FlexNet on multiple seed files")
        SeedFileList <- fread(config$Seed.File.List,header=FALSE,sep="\t",stringsAsFactors = FALSE)$V1
    }

    SeedFileList <- lapply(SeedFileList,function(x){if(file.exists(x)){message("STEP: Seed file ",x, " is existing");x}else{message("STEP: Seed file ",x, " is not existing")}}) %>% unlist
    SeedFilelist <- SeedFileList[!is.null(SeedFileList)]


    ## need for RWR
    PoolNodesList <- split(opt$TransitionMatrix %>% rownames %>% strsplit(.,"::",perl =T) %>% lapply(., function(x){x[length(x)]}) %>% unlist, opt$TransitionMatrix %>% rownames %>% strsplit(.,"::",perl =T) %>% lapply(., function(x){x[1]}) %>% unlist)

    ## get layer information from transition matrix
    rn <- str_split(rownames(opt$TransitionMatrix),"::",n=4) %>% do.call(rbind,.) %>% data.table %>% .[,V5:=apply(.SD[,1:3],1,function(x){paste(x,collapse = "::")})]    %>% .[,V6:=rownames(opt$TransitionMatrix)]        
    layers_list <-lapply(unique(rn$V1),function(x){tmp <- lapply(unique(rn[V1==x,]$V2),function(y){unique(rn[V1==x & V2==y,]$V5)});tmp;names(tmp) <- unique(rn[V1==x,]$V2);tmp})
    names(layers_list) <- unique(rn$V1)

    ## get tau
    tau <- lapply(names(layers_list),function(x){getNetworkSpecificParams(x,layers_list,"Restart.Probability",config)})
    names(tau) <- names(layers_list)

    options(future.globals.maxSize=1024 * 1024*1024*1024) ##30G
    plan(multicore, workers = opt$threads)
    rlt <- foreach(Seed_File=SeedFileList) %dofuture% { 
        All_Seeds <- read.csv(Seed_File,header=FALSE,sep="\t",dec=".",stringsAsFactors = FALSE)

        Seed_File_List <- check.seeds.general(All_Seeds,PoolNodesList)
        
        ## validation
        ## all(Seeds == unlist(All_seeds_ok))

        ## eta is the importance of each type of networks
        ##We compute the restart probability of each seed, based on eta and tau.
        Seeds_Score <- get.seed.scores.weighted.universal(Seed_File_List,config$eta,tau)
        if(nrow(Seeds_Score)==0){
            return(NULL)
        }

        ## eta is  [1,2,3] for gene, disease,and drug
        ##validation
        ## gene eta =  1/6
        ## Seeds_Score[grepl("^[A-Z].*",Seeds_Score$Seeds,perl = T) & !grepl("OMIM",Seeds_Score$Seeds) ,]$Score %>% sum

        ##disease eta = 2/6
        ## Seeds_Score[grepl("^[A-Z].*",Seeds_Score$Seeds,perl = T) & grepl("OMIM",Seeds_Score$Seeds) ,]$Score %>% sum

        ##drug eta = 3/6
        ## Seeds_Score[!grepl("^[A-Z].*",Seeds_Score$Seeds,perl = T) & !grepl("OMIM",Seeds_Score$Seeds) ,]$Score %>% sum

        ## total score is 1 for all nodes
        ## Seeds_Score$Score %>% sum

        if(opt$verbose)
            print(Seeds_Score)

        if(opt$verbose)
            message("STEP: ","Performing Random Walk...")

        ## print(names(opt$TransitionMatrix))
        ## We perform the Walks on the multiplex-heterogeneous network.
        Random_Walk_Results <- Random_Walk_Restart(Network_Matrix = opt$TransitionMatrix,r = config$r,SeedGenes = Seeds_Score,max_iter = opt$max_iter)
        print("here")
        final_rank <- list()
### We remove the seed genes from the ranking, we sort by score and we write the results.
        final_rank <- rank_native(PoolNodesList,Random_Walk_Results,All_Seeds)
        tmp <- lapply(1:length(final_rank),function(i){
            if(opt$verbose)
                message("STEP: saving to", paste0(opt$dir_out,"/",basename(Seed_File),".",opt$i,".",names(final_rank)[i]))
            write.table(final_rank[[i]],file = paste0(opt$dir_out,"/",basename(Seed_File),".",opt$i,".",names(final_rank)[i]),sep="\t",row.names = FALSE, dec=".",quote=FALSE)
        })
        final_rank
    }
    ##names(rlt) <- SeedFileList
    if(opt$verbose)
        message("STEP: FlexNet done.\n")
    ##opt$final_rank <- rlt
    opt
}


FlexNetBootstrapping <- function(opt){
    ## singleSamplingMode is used for load a single transition rda. 
    ## FALSE indicates run Flexnet on input using a batch of transition rda.
    ## TRUE can speed up by using more cores
    if(opt$verbose){
        message("STEP: FlexNetBootstrapping...")
    }

    if(opt$singleSamplingMode){
        Is = opt$number
    }else{
        Is = 1:opt$number
    }
    options(future.globals.maxSize=200 * 1024*1024*1024) 
    plan(multicore, workers = opt$threadsSampling)
    ## registerDoParallel(opt$threadsSampling)
    foreach(i=Is) %dofuture% { 
        opt$i <- i
        opt$config$saveTransitionMatrix <- paste0(opt$saveIntermediate,".transition.",opt$i,".rda")
        
        if(file.exists(opt$config$saveTransitionMatrix)){
            if(opt$verbose){
                message("STEP: FlexNet...",opt$config$saveTransitionMatrix, "")
            }

            opt <- FlexNet(opt)
            rm(opt)
            gc()
        }else{
            if(opt$verbose)
                message("STEP: ",opt$config$saveTransitionMatrix, " not exist")
        }
    }
}


####################################################################
#######		Main function
#####################################################################

## source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/rmd/netRepurpose/Functions/All_Functions.v3.R",sep = ""))
## source("/fs0/chenr6/chenr6/myscript/Rscripts/rmd/netRepurpose/Functions/All_Functions.v2.R")
source(paste0(Sys.getenv("BASEDIR"),"/myscript/Rscripts/mypkg/FlexNet_bootstrap/All_Functions.R"))

opt$config <- read_yaml(file = opt$yaml)
if(opt$verbose)
    message("STEP: ","Read parameters from yaml",opt$yaml,"\n")
opt$config <- check.parameters.from.yaml.general(opt)
## overide opts in config
if(opt$saveIntermediate != "NULL"){
    opt$config$saveIntermediate <- opt$saveIntermediate
    if(!dir.exists(dirname(opt$saveIntermediate))) dir.create(dirname(opt$saveIntermediate),showWarnings=T,recursive=T)
}

## if seed.file.list is s
if(opt$Seed.File.List !="NULL"){
    opt$config$Seed.File.List <- opt$Seed.File.List
}

## If seed.file is specified, use it instead of the file specified in yaml
if(opt$Seed.File !="NULL"){
    opt$config$Seed.File <- opt$Seed.File
}

for(cmd in opt$cmd){
    if(grepl("generate_GraphList",cmd,perl =T,ignore.case = TRUE)){
        opt <- generate_GraphList(opt)
    }else if(grepl("generate_SupraAdjacencyMatrix",cmd,perl =T,ignore.case = TRUE)){
        opt <- generate_SupraAdjacencyMatrix(opt)
    }else if(grepl("generate_TransitionMatrix",cmd,perl =T,ignore.case = TRUE)){
        opt <- generate_TransitionMatrix(opt)
    }else if(grepl("FlexNet$",cmd,perl =T,ignore.case = TRUE)){
        opt <- FlexNet(opt)
    }else if(grepl("sampleGeneDiseaseAssocation",cmd,perl =T, ignore.case = TRUE)){
        ## set.seed(1)
        sampleGeneDiseaseAssocation(opt)
    }else if(grepl("sampleRandomGeneDiseaseAssocation",cmd,perl =T, ignore.case = TRUE)){
        ## set.seed(1)
        sampleRandomGeneDiseaseAssocation(opt)
    }else if(grepl("FlexNetBootstrapping",cmd,perl =T, ignore.case = TRUE)){
        FlexNetBootstrapping(opt)
    }
    
}

## plots <- heatmap.ggplot(opt=opt)

####################################################################
#######		save RData
#####################################################################
## if(opt$save_RData){
##     cat("save plot to RData\n")
##     save(x_plot,opt,x,plots,file="x_plot.RData")
## }

## if(is.na(opt$output)){
##     x11()
##     cat("Plotting QQ plot ...\n")
##     x_plot
##     message("Click the graph to exit..."); loc=locator(1);
##     loc = locator(1);
##     res = dev.off();

## } else {
##     opendev(opt)
##     lapply(plots,function(x){grid.newpage();grid.draw(x)})
##     res=dev.off()
## }
