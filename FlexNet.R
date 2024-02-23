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


## loadfonts()
## theme_figure <- theme_classic() + theme(panel.background = element_rect(fill = 'transparent', colour = NA),legend.title=element_blank(),text=element_text(family=opt$fontfamily),axis.title.x = element_blank(),plot.title = element_text(hjust=0.5, face="bold",size=12),axis.line = element_blank())

if(!dir.exists(opt$dir_out)) dir.create(opt$dir_out,showWarnings=T,recursive=T)

if(opt$verbose)
    message(opt$dir_out)


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
parser_script$add_argument("--Bipartite.LoadedGeneLayers.LoadedDiseaseLayers", default="NULL",dest="Bipartite.LoadedGeneLayers.LoadedDiseaseLayers",help = "Bipartite.LoadedGeneLayers.LoadedDiseaseLayers is [default \"%(default)s\"]")

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



generate_AdjacencyMatrixList <- function(opt,save=TRUE){
    config <- opt$config

    if(is.null(config$saveAdjacencyMatrixList)){
        config$saveAdjacencyMatrixList <- paste0(opt$saveIntermediate,".adj.","rda")
    }

    if(file.exists(config$saveAdjacencyMatrixList)){
        if(opt$verbose)
            message("\nSTEP: ",config$saveAdjacencyMatrixList, " exists. Skip generating and loading directly")
        load(config$saveAdjacencyMatrixList)
    }else{
        ## Add some network info to the save name.
        if(opt$verbose)
            message("\nSTEP: ",config$saveAdjacencyMatrixList, "doesn't exist. Change name by adding network information. Generating...")

        ##		Read the layers of our multiplex networks
        ## The structure of AdjacencyMatrixList is like following:
        ## For multiplex networks, there are four slots:

        ## | slot                 | description                                                       |
        ## | AdjacencyMatrix | The adjacency matrix for a specific type of multiplex network     |
        ## | L                    | Number of layers for this type of multiplex network               |
        ## | N                    | Number of nodes for this type of multiplex network                |
        ## | pool_nodes           | The sorted union of nodes name for this type of multiplex network |
        ## | graphs               | List of igraph objects for each layer in the multiplex network    |


        AdjacencyMatrixList <- list()

        if(opt$verbose)
            message("STEP: ","Reading Multiplex Layers...")
        for(NetworkID in config$LoadedNetworks){
            if(opt$verbose)
                message("STEP: ","Reading Multiplex Layers: ",NetworkID)
            AdjacencyMatrixList[[paste(NetworkID,NetworkID,sep = "_")]] <- get.adj.multiplex.by.networkid(NetworkID,config, returnGraph=F)
        }
        AdjacencyMatrixList %>% names

        ##		 Read the layers of our bipartite networks
        ## The structure of AdjacencyMatrixList is like following:

        ## For bipartite networks, there are two slots:
        ## | slot                 | description                                                   |
        ## | AdjacencyMatrix | The adjacency matrix for a specific type of bipartite network |
        ## | graphs               | igraph object for each layer in the bipartite network         |
        if(opt$verbose)
            message("STEP: ","Reading Bipartite Layers...")
        AllPossibleBipartite <- as.list(as.data.frame(combn(config$LoadedNetworks,2,simplify =T)))
        for(Bipartite_NetworkID in AllPossibleBipartite){
            if(opt$verbose)
                message("STEP: ","Reading Bipartite Layers: ",Bipartite_NetworkID)
            ## The starting vertex in on the column
            AdjacencyMatrixList[[paste(Bipartite_NetworkID,collapse = "_")]] <- get.adj.bipartite.by.networkid(Bipartite_NetworkID,AdjacencyMatrixList,config)
            ## save the transpose of  the bipartite matrix 
            AdjacencyMatrixList[[paste(rev(Bipartite_NetworkID),collapse = "_")]] <- get.adj.bipartite.by.networkid(rev(Bipartite_NetworkID),AdjacencyMatrixList,config)
        }
        AdjacencyMatrixList %>% names

        ## Save adjacency matrix

        if(opt$verbose)
            message("STEP: ","save AdjacencyMatrixList to ",config$saveAdjacencyMatrixList, "...")
        if(save){
            save(AdjacencyMatrixList,file = config$saveAdjacencyMatrixList)
        }
    }
    opt$AdjacencyMatrixList <- AdjacencyMatrixList
    if(opt$verbose)
        message("STEP: generating AdjacencyMatrixList done.\n")
    opt
}


#####################################################################
#######	generate SupraAdjacencyMatrixList if not exists
#####################################################################
generate_SupraAdjacencyMatrixList <- function(opt,save=TRUE){
    config <- opt$config

    if(is.null(config$saveAdjacencyMatrixList)){
        config$saveAdjacencyMatrixList <- paste0(opt$saveIntermediate,".adj.rda")
    }


    if(is.null(opt$AdjacencyMatrixList) && file.exists(config$saveAdjacencyMatrixList)){
        load(config$saveAdjacencyMatrixList)
        opt$AdjacencyMatrixList <- AdjacencyMatrixList
        rm(AdjacencyMatrixList)
        gc()
    }

    
    if(is.null(config$saveSupraAdjacencyMatrixList)){
        config$saveSupraAdjacencyMatrixList <- paste0(opt$saveIntermediate,".supra.rda")
    }

    if(file.exists(config$saveSupraAdjacencyMatrixList)){
        if(opt$verbose)
            message("\nSTEP: ",config$saveSupraAdjacencyMatrixList, " exists. Skip generating and loading directly")
        load(config$saveSupraAdjacencyMatrixList)
        
    }else{
        ## Add some network info to the save name.
        if(opt$verbose)
            message("\nSTEP: ",config$saveSupraAdjacencyMatrixList, "doesn't exist. Change name by adding network information. Generating...")

        SupraAdjacencyMatrixList <- list()
        if(opt$verbose)
            message("STEP: ","Reading Multiplex Layers...")
        for(NetworkID in config$LoadedNetworks){
            if(opt$verbose)
                message("STEP: ","Reading Multiplex Layers: ",NetworkID)
            SupraAdjacencyMatrixList[[paste(NetworkID,NetworkID,sep = "_")]] <- get.supra.adj.multiplex.by.networkid(NetworkID,opt$AdjacencyMatrixList,config)
        }
        SupraAdjacencyMatrixList %>% names

        if(opt$verbose)
            message("STEP: ","Reading Bipartite Layers...")
        AllPossibleBipartite <- as.list(as.data.frame(combn(config$LoadedNetworks,2,simplify =T)))
        if(length(AllPossibleBipartite)>0){
            for(Bipartite_NetworkID in AllPossibleBipartite){
                if(opt$verbose)
                    message("STEP: ","Reading Bipartite Layers: ",Bipartite_NetworkID)
                ## The starting vertex in on the column
                SupraAdjacencyMatrixList[[paste(Bipartite_NetworkID,collapse = "_")]] <- get.supra.adj.bipartite.by.networkid(Bipartite_NetworkID,opt$AdjacencyMatrixList,config)
                ## transpose the bipartite matrix and save it to reversed ID 
                SupraAdjacencyMatrixList[[paste(rev(Bipartite_NetworkID),collapse = "_")]] <- get.supra.adj.bipartite.by.networkid(rev(Bipartite_NetworkID),opt$AdjacencyMatrixList,config)
            }
        }
        SupraAdjacencyMatrixList %>% names
        lapply(SupraAdjacencyMatrixList,dim)
        
        if(opt$verbose)
            message("STEP: ","save SupraAdjacencyMatrixList to ",config$saveSupraAdjacencyMatrixList, "...")
        if(save){
            save(SupraAdjacencyMatrixList,file = config$saveSupraAdjacencyMatrixList)
        }
    }
    opt$SupraAdjacencyMatrixList <-  SupraAdjacencyMatrixList
    if(opt$verbose)
        message("STEP: generating SupraAdjacencyMatrixList done.\n")
    opt$AdjacencyMatrixList <- NULL
    gc()
    opt
}


####################################################################
#######		generate transition matrix
#####################################################################
generate_TransitionMatrixList  <- function(opt,save=TRUE){
    config <- opt$config

    if(is.null(config$saveSupraAdjacencyMatrixList)){
        config$saveSupraAdjacencyMatrixList <- paste0(opt$saveIntermediate,".supra.rda")
    }

    if(is.null(opt$SupraAdjacencyMatrixList) && file.exists(config$saveSupraAdjacencyMatrixList)){
        load(config$saveSupraAdjacencyMatrixList)
        opt$SupraAdjacencyMatrixList <- SupraAdjacencyMatrixList
        rm(SupraAdjacencyMatrixList)
        gc()
    }

    
    if(is.null(config$saveTransitionMatrixList)){
        config$saveTransitionMatrixList <- paste0(opt$saveIntermediate,".transition.rda")
    }

    if(file.exists(config$saveTransitionMatrixList)){
        if(opt$verbose)
            message("\nSTEP: ",config$saveTransitionMatrixList, " exists. Skip generating and loading directly")
        load(config$saveTransitionMatrixList)
        
    }else{
        ## Add some network info to the save name.
        if(opt$verbose)
            message("\nSTEP: ",config$saveTransitionMatrixList, "doesn't exist. Change name by adding network information. Generating...")

        if(opt$verbose)
            message("STEP: ","Transforming from adjacency matrices to transition matrix...\n")
        TransitionMatrixList <- get.transition.multiplex.heterogenous(config,opt$SupraAdjacencyMatrixList)

        if(opt$verbose)
            message("STEP: ","save TransitionMatrixList to ",config$saveTransitionMatrixList, "...")
        if(save){
            save(TransitionMatrixList,file = config$saveTransitionMatrixList)
        }    
    }
    opt$TransitionMatrixList <- TransitionMatrixList

    if(opt$verbose)
        message("STEP: generating TransitionMatrix done.\n")

    opt$SupraAdjacencyMatrixList <- NULL
    gc()
    opt

}

## ####################################################################
## ####### WE PREPARE THE SEEDS SCORES AND WE PERFORM THE RWR-MH. WE WRITE
## ####### THE OUTPUT FILE CONTAINING HGNC GENE SYMBOL AND THEIR RESULTING SCORES. 		
## #####################################################################

### Now, that we have the diseases and the proteins of our system we can check the nodes.
FlexNet <- function(opt){
    config <- opt$config

    if(is.null(config$saveTransitionMatrixList)){
        config$saveTransitionMatrixList <- paste0(opt$saveIntermediate,".transition.rda")
    }

    if(is.null(opt$TransitionMatrixList) && file.exists(config$saveTransitionMatrixList)){
        if(opt$verbose){
            message("STEP: ","load FlexNet...")
        }
        load(config$saveTransitionMatrixList)
        opt$TransitionMatrixList <- TransitionMatrixList
        rm(TransitionMatrixList)
        gc()
    }

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
    registerDoParallel(opt$threads)

    rlt <- foreach(Seed_File=SeedFileList) %dopar% { 
        All_Seeds <- read.csv(Seed_File,header=FALSE,sep="\t",dec=".",stringsAsFactors = FALSE)
        Seed_File_List <- check.seeds.general(All_Seeds,opt$TransitionMatrixList$PoolNodesList)
        
        ## validation
        ## all(Seeds == unlist(All_seeds_ok))

        ## eta is the importance of each type of networks
        ##We compute the restart probability of each seed, based on eta and tau.
        Seeds_Score <- get.seed.scores.weighted.universal(Seed_File_List,config$eta,opt$TransitionMatrixList$tau)
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
            message(paste0(capture.output(Seeds_Score), collapse = "\n"))

        if(opt$verbose)
            message("STEP: ","Performing Random Walk...")

        ## print(names(opt$TransitionMatrixList))
        ## We perform the Walks on the multiplex-heterogeneous network.
        Random_Walk_Results <- Random_Walk_Restart(Network_Matrix = opt$TransitionMatrixList$TransitionMatrix,r = config$r,SeedGenes = Seeds_Score,max_iter = opt$max_iter)

        final_rank <- list()
### We remove the seed genes from the ranking, we sort by score and we write the results.
        final_rank <- rank_native(opt$TransitionMatrixList$PoolNodesList,Random_Walk_Results,All_Seeds)
        names(final_rank[["LoadedGeneLayers"]])[-1:-4] <- opt$TransitionMatrixList$LayerName[["LoadedGeneLayers"]]
        tmp <- lapply(1:length(final_rank),function(i){
            write.table(final_rank[[i]],file = paste0(opt$dir_out,"/",basename(Seed_File),".",names(final_rank)[i]),sep="\t",row.names = FALSE, dec=".",quote=FALSE)
        })
        final_rank
    }
    names(rlt) <- SeedFileList
    if(opt$verbose)
        message("STEP: FlexNet done.\n")
    opt$final_rank <- rlt
    opt
}


####################################################################
#######		Main function
#####################################################################
require("igraph")
require("Matrix")
require("dnet")
require("supraHex")
require("yaml")
require("data.table")

## source(paste(Sys.getenv("BASEDIR"), "/myscript/Rscripts/rmd/netRepurpose/Functions/All_Functions.v3.R",sep = ""))
source("/fs0/chenr6/chenr6/myscript/Rscripts/mypkg/FlexNet/All_Functions.R")

opt$config <- read_yaml(file = opt$yaml)
if(opt$verbose)
    message("STEP: ","Read parameters from yaml",opt$yaml,"\n")
opt$config <- check.parameters.from.yaml.general(opt)
## overide opts in config
if(opt$saveIntermediate != "NULL"){
    opt$config$saveIntermediate <- opt$saveIntermediate
    if(!dir.exists(dirname(opt$saveIntermediate))) dir.create(dirname(opt$saveIntermediate),showWarnings=T,recursive=T)
}

if(opt$Seed.File.List !="NULL"){
    opt$config$Seed.File.List <- opt$Seed.File.List
}
if(opt$Bipartite.LoadedGeneLayers.LoadedDiseaseLayers!="NULL"){
    opt$config$Bipartite.LoadedGeneLayers.LoadedDiseaseLayers <- opt$Bipartite.LoadedGeneLayers.LoadedDiseaseLayers
}

for(cmd in opt$cmd){
    if(grepl("generate_AdjacencyMatrixList",cmd,perl =T,ignore.case = TRUE)){
        opt <- generate_AdjacencyMatrixList(opt)
    }else if(grepl("generate_SupraAdjacencyMatrixList",cmd,perl =T,ignore.case = TRUE)){
        opt <- generate_SupraAdjacencyMatrixList(opt)
    }else if(grepl("generate_TransitionMatrixList",cmd,perl =T,ignore.case = TRUE)){
        opt <- generate_TransitionMatrixList(opt)
    }else if(grepl("FlexNet",cmd,perl =T,ignore.case = TRUE)){
        opt <- FlexNet(opt)
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
