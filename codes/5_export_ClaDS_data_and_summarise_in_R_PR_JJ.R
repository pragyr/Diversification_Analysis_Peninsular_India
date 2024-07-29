library(parallel)
library(foreach)
library(doParallel)
library(doFuture)

# PLEASE CHECK THE FOLLOWING TWO LINES
# if you are using the "Rscript" command in linux/mac terminal
setwd(paste0(getwd(),"/../data/all_trees"))

# if you are using RStudio IDE
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/all_trees_or_super_trees"))

# OR DIRECTLY PUT THE PATH AS A STRING

filenames <- list.files()
lineage.names <- substr(filenames,1,(nchar(filenames)-4))

# making a list of node numbers corresponding to our clades of interest in each super-tree containing that clade
nodes.of.interest <- list(c("Acanthaceae"=16),
                          c("Ahaetullinae"=14),
                          c("Ethmostigmus"=89,"Rhysida_longipes_clade"=110,"Rhysida_immarginata_clade"=115,"Digitipes"=122),
                          c("Ceropegiae"=68),
                          c("Cnemaspis"=84),
                          c("Dravidogecko"=8),
                          c("Eutropis"=13),
                          c("Fan-Throated"=18),
                          c("Geckoella"=15),
                          c("Gegeneophis"=15),
                          c("Gerrhopilus"=5),
                          c("Hemidactylus"=49),
                          c("Hemiphyllodactylus"=13),
                          c("Heterometrinae"=24),
                          c("Memecylon"=39),
                          c("Microhyla1"=88,"Microhyla2"=96),
                          c("Montecincla"=69),
                          c("Indonaia"=85,"Lamellidens"=102,"Parreysia"=97),
                          c("Nyctibatrachidae"=59),
                          c("Ophisops"=37),
                          c("Piper1"=235,"Piper2"=193),
                          c("Pseudophilautus"=25),
                          c("Ranixalidae"=22),
                          c("Raorchestes"=20),
                          c("Testudines"=15),
                          c("Uperodon"=7),
                          c("Uropeltids"=45)
                          )

clads.rates.matrix <- matrix(ncol=4)

colnames(clads.rates.matrix) <- c("Lineage","Mean_Speciation_Rate","Mean_Extinction_Rate","Mean_Net-diversification_Rate")

# loop over all the ClaDS .Rdata result files
for(j in 1:length(filenames)){
  tree <- read.tree(filenames[j])
  # load the results of the ClaDS analysis into R
  # the result file should have the same name as that of the relevant super-tree
  load(paste0("../ClaDS_output_Julia/",lineage.names[j],".Rdata"))
  nodes <- unlist(nodes.of.interest[[j]])
  l.nodes <- length(nodes)
  
  # loop over all the relevant clades that belong to a super-tree
  for(k in 1:l.nodes){
    descendant.nodes <- unname(c(nodes[k],Descendants(tree,nodes[k],"all")))
    #descendant.nodes <- descendant.nodes[-which(descendant.nodes < nodes[k])]
    
    node.ages <- tree.age(tree)
    node.ages <- node.ages[which(rownames(node.ages) %in% as.character(descendant.nodes)),]
    cr.age <- max(node.ages$ages)
    
    time.bins <- seq(cr.age,0,-(cr.age/100))
    
    # if the node number corresponds to the root node of the super-tree, i.e., the entire super-tree is indeed our clade of interest,
    # then the rate in the oldest/first time bin would be the initial speciation rate, estimated in ClaDS
    # else extract it from the branch rates
    if(cr.age == max(tree.age(tree)$ages)){
      mean.sp.rate.in.each.time.bin <- c(CladsOutput$lambda0_map)
    }else{
      mean.sp.rate.in.each.time.bin <- c(CladsOutput$lambdai_map[which(tree$edge[,2] == 
                                                                         as.numeric(rownames(node.ages)[which(node.ages$ages == cr.age)]))])
    }
    
    # calculate rates through time within this clade
    for(t in 2:length(time.bins)){
      # t <- 101
      indices.older <- which(node.ages$ages >= abs(time.bins[t]))
      if(length(indices.older) == 1 & indices.older[1] == which(node.ages$ages == cr.age)){
        mean.sp.rate.in.each.time.bin <- append(mean.sp.rate.in.each.time.bin,mean.sp.rate.in.each.time.bin[1])
      }else{
        # search for branches that pertain to each time bin
        start.nodes.older.than.this.time <- as.numeric(rownames(node.ages)[indices.older])
        nodes.connected.to.those.nodes <- tree$edge[which(tree$edge[,1] %in% start.nodes.older.than.this.time),2]
        
        
        ages.of.end.nodes <- node.ages$ages[match(nodes.connected.to.those.nodes,as.numeric(rownames(node.ages)))]
        end.nodes.younger.than.previous.time <- nodes.connected.to.those.nodes[which(ages.of.end.nodes < time.bins[t-1])]
        
        indices.of.relevant.edges <- which(unlist(as.vector(as_tibble(tree$edge) %>% unite("combns",V1:V2))) %in% 
                unlist(as.vector(as_tibble(
          expand.grid(start.nodes.older.than.this.time,end.nodes.younger.than.previous.time)) %>% unite("combns",Var1:Var2))))
        
        
        rates <- CladsOutput$lambdai_map[indices.of.relevant.edges]
        mean.sp.rate.in.each.time.bin <- append(mean.sp.rate.in.each.time.bin,mean(rates))
      }
    }
    
    mean.sp.rate.in.each.time.bin[101] <- mean(CladsOutput$lambdatip_map)
    mean.ex.rate.in.each.time.bin <- mean.sp.rate.in.each.time.bin*CladsOutput$eps_map
    mean.nd.rate.in.each.time.bin <- mean.sp.rate.in.each.time.bin - mean.ex.rate.in.each.time.bin
    
    clads.rates.matrix <- rbind(clads.rates.matrix,
                               c(names(nodes)[k],mean(mean.sp.rate.in.each.time.bin),mean(mean.ex.rate.in.each.time.bin),mean(mean.nd.rate.in.each.time.bin)))
  }
  
}
clads.rates.matrix <- clads.rates.matrix[-1,]

clads.rates.matrix <- clads.rates.matrix[order(clads.rates.matrix[,"Lineage"]),]

dir.create("../Summary and Graphs")
write.csv(clads.rates.matrix, "../Summary and Graphs/clads_mean_rates.csv")
