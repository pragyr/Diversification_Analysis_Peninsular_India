# run the ClaDS analysis on your system
# Julia needs to be installed on your system beforehand
system2("codes/4_run_ClaDS_in_Julia_PR_JJ.jl")

library(phangorn)
library(ape)
library(dispRity)


setwd("data/trees_for_ClaDS_analyses")

filenames <- list.files() # store the filenames
lineage.names <- substr(filenames,1,(nchar(filenames)-4)) # remove ".tre" extension to retain the names of the lineages/clades

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
                          c("Pseudophilautus"=70),
                          c("Ranixalidae"=22),
                          c("Raorchestes"=75),
                          c("Testudines"=15),
                          c("Uperodon"=7),
                          c("Uropeltids"=45)
)

clads.mean.rates.matrix <- matrix(ncol=4) # make an empty matrix with 4 columns and 1 row to store the mean diversification rates for the clades of interest
clads.tip.rates <- matrix(ncol=4) # make an empty matrix with 4 columns and 1 row to store the tip rates of all the clades of interest

# set the column names for both the matrices
colnames(clads.mean.rates.matrix) <- c("Lineage","Mean_Speciation_Rate","Mean_Extinction_Rate","Mean_Net-diversification_Rate")
colnames(clads.tip.rates) <- c("Lineage","Speciation_Rate","Extinction_Rate","Net-diversification_Rate")

# loop over all the ClaDS .Rdata result files
for(j in 1:length(filenames)){
  #j<-1
  tree <- read.tree(filenames[j])
  # load the results of the ClaDS analysis into R
  # the result file should have the same name as that of the relevant super-tree
  load(paste0("../ClaDS_output_Julia/",lineage.names[j],".Rdata"))
  nodes <- unlist(nodes.of.interest[[j]]) # obtain the node number(s) of the clade(s) of interest within each of the input trees
  l.nodes <- length(nodes) # counting the number of clades that are of our interest within one input tree
  
  # loop over all the relevant clades that belong to the input tree
  for(k in 1:l.nodes){
    #k <- 1
    descendant.nodes <- unname(c(nodes[k],Descendants(tree,nodes[k],"all"))) # store the node numbers of nodes (internal and terminal) that are descendants of nodes[k] including the node itself
    descendant.only.tips <- unlist(Descendants(tree,nodes[k],"tips")) # store the node numbers of only the terminal nodes (tips) that are descendants of nodes[k]
    
    sp.rate.tips <- CladsOutput$lambdatip_map[descendant.only.tips] # obtain the tip speciation rates from the ClaDS output
    ex.rate.tips <- CladsOutput$eps_map*sp.rate.tips # calculate the extinction rate by multiplying the constant turnover rates (eps_map) with the tip speciation rates
    nd.rate.tips <- sp.rate.tips - ex.rate.tips # calculate the net-diversification rates by taking the difference between the two rates
    
    clads.tip.rates <- rbind(clads.tip.rates,cbind(rep(names(nodes)[k],length(descendant.only.tips)),
                                                   sp.rate.tips,ex.rate.tips,nd.rate.tips)) # add the three diversification rates to the matrix

    
    node.ages <- tree.age(tree)
    node.ages <- node.ages[which(rownames(node.ages) %in% as.character(descendant.nodes)),]
    cr.age <- max(node.ages$ages)
    
    time.bins <- seq(cr.age,0,-(cr.age/100)) # splitting the crown age into 100 timebins
    
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
    
    clads.mean.rates.matrix <- rbind(clads.mean.rates.matrix,
                                c(names(nodes)[k],mean(mean.sp.rate.in.each.time.bin),mean(mean.ex.rate.in.each.time.bin),mean(mean.nd.rate.in.each.time.bin)))
  }
  
}
clads.mean.rates.matrix <- clads.mean.rates.matrix[-1,]
clads.mean.rates.matrix <- clads.mean.rates.matrix[order(clads.mean.rates.matrix[,"Lineage"]),]
clads.tip.rates <- clads.tip.rates[-1,]
clads.tip.rates <- clads.tip.rates[order(clads.tip.rates[,"Lineage"]),]
View(clads.mean.rates.matrix)
write.csv(clads.mean.rates.matrix, "../clads_mean_rates.csv",row.names = F)
write.csv(clads.tip.rates,"../clads_tip_rates.csv",row.names = F)

# set the working directory back to the root of the R project
setwd("../../")