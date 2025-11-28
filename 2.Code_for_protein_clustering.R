######
###### Code for clustering of blood proteome for hub protein identification
#-----------------------------------------------------------------------------------------------------------------------

##### Co-regulation analysis for proteins from Group 1/2/3
setwd('/PATH')

Raw_protein_level=read.csv("./File.csv",header = T)

##-----------------------------------------------
##### Generate R value for pairwise correlations (proteins from Group 1/2/3)
##-----------------------------------------------
Correlation_matrix<-matrix(,ncol=length(Raw_protein_level[,1]), nrow = length(Raw_protein_level[,1]))
colnames(Correlation_matrix)<-rownames(Raw_protein_level)
rownames(Correlation_matrix)<-rownames(Raw_protein_level)

for (i in 1:length(Raw_protein_level[,1]))
{
  for (k in 1:length(Raw_protein_level[,1]))
  {
    Correlation_matrix[i,k]=as.numeric(cor(t(Raw_protein_level[i,]),t(Raw_protein_level[k,]),use="p"))
  }
}

Summary_correlation_matrix=data.frame(Correlation_matrix)

write.table(Summary_correlation_matrix, '/FILE_matrix.txt')

##------------------------------------------------------------------------
##### Generate pairwise correlation Heatmap plot among candidate proteins
##-------------------------------------------------

library(gplots)
set.seed(100)

breaks = seq(-1,max(Correlation_matrix),length.out=length(Raw_protein_level[,1])^2)
gradient1 = colorpanel( sum( breaks[-1]<=0 ), "dodgerblue4", "steelblue2","white" )
gradient2 = colorpanel( sum( breaks[-1]>0 ), "white", "tomato2","firebrick4" )
hm.colors = c(gradient1,gradient2)

heatmap.2(Correlation_matrix,scale="none",breaks=breaks,col=hm.colors,
          trace="none", 
          margin=c(5,20),lwid=c(1.5,5.0))

graph<-heatmap.2(Correlation_matrix,scale="none",breaks=breaks,col=hm.colors,
                 trace="none", 
                 margin=c(5,20),lwid=c(1.5,5.0))
cluster<-as.hclust(graph$rowDendrogram)
modules<-as.data.frame(cutree(cluster,h=4.0))
colnames(modules)<-c("Cluster_by_heatmap")
#heatmap(AD_Correlation_matrix,scale="none")
#p<-heatmap(AD_Correlation_matrix,scale="none")
#p
write.table(modules, '/FILE_model_clustering_based_on_Heatmap.txt')



##-----------------------------------------------
##### Network of each group of proteins
##-----------------------------------------------
## Input nodes and links info
library(dplyr)
library(tidyverse)
setwd('FIL_Protein_groups_clustering')

Nodes_info=read.csv('./Nodes_Cluster_for_Network.csv',header = T)
Links_info=read.csv('./Edges_Cluster_for_Network.csv',header = T)

Nodes=Nodes_info[c(1)]
Links=Links_info[c(1,2,4)]
##------------------------
# Rename nodes and edges

Nodes <- Nodes %>% rowid_to_column("id")
Edges <- Links %>% 
  left_join(Nodes, by = c("Proteins" = "Proteins")) %>% 
  rename(from = id)
Edges <- Edges %>% 
  left_join(Nodes, by = c("Hubs" = "Proteins")) %>% 
  rename(to = id)

Nodes$Size=Nodes_info$Size

colnames(Nodes)<-c("id","Protein_ID","value")

edges<-select(Edges,from,to,Thickness)
nodes<-select(Nodes,id,value)

##------------------------
# Modify the plot by visNetwork
#detach(package:network)
#rm(routes_network)

library(visNetwork)
#nodes$label=Nodes$Protein_ID
nodes$color=Nodes_info$Color
nodes$group=Nodes_info$Group
edges <- mutate(edges, width = log2(as.numeric(Thickness)+1)*10)
edges<- mutate(edges, length = 0.005)
edges$color<-"gray"

visNetwork(nodes, edges) %>% 
  #  visNodes(color=list(background="lightblue",border="darkblue",highlight="yellow")) %>% 
  #visIgraphLayout(layout = "layout_with_fr") %>% 
  # visEdges(arrows = "to") %>% 
  #visGroups(groupname = "Hub",color="#808080") %>%
  #visGroups(groupname = "Member",color="#FC6666")
  visLayout(randomSeed = 16)


##------------------------
# Plot network (quickly and briefly)
library(network)
routes_network <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE)
class(routes_network)
routes_network
plot(routes_network,vertex.cex=1,vertex.col="burlywood1")
