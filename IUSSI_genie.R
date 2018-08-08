library(igraph)
library(GENIE3)
library(data.table)

#Load data
tpm <- read.csv("workshop-iussi2018/Mariana/tpm.csv")

#Load gene descriptions, including drosophila melanogaster orthologs
aName <- fread("workshop-iussi2018/Mariana/gene_descriptions.tsv",sep="\t")

#Get genes as rownames
rownames(tpm) = tpm$gene_id
tpm = tpm[,-c(1)]

#Filter out genes with zero expression
tpm = tpm[rowSums(tpm) > 0,]

#hyperbolic sine transformation. This normalizes the data; is similar to log but defined at 0
tpm = log(tpm + sqrt(tpm ^ 2 + 1))

#Filter for genes in the ILP-2 module
mrjp = aName$gene_id[grepl("Mrjp",aName$gene_id)]
vg = aName$gene_id[grepl("itellogenin",aName$description)][c(1,3)]
jh = aName$gene_id[grepl("juvenile",aName$description)][c(1,2)]
insulin = aName$gene_id[grepl("insulin",aName$description)]

genes_mod <- c(mrjp,vg,jh,insulin)

sub = tpm[genes_mod,]

#Construct a network! 
network <- GENIE3(as.matrix(sub))

#Transform matrix into list of links, sorted by weight
netList <- getLinkList(network)

#Have to cut off links somewhere. Formally, could generate "null" distribution of connection weights.
#Here, just pick a cut-off (this is a number you can play with)
cut_off <- quantile(netList$weight,0.85)
netList_cut <- netList[netList$weight > cut_off,]

#Filter list for connections containing a specific gene...pick whatever you want!
focus_name <- "Vg"
g_out <- netList_cut[netList_cut$regulatory.gene == focus_name,]
g_in <- netList_cut[netList_cut$target.gene == focus_name,]

#Find the gene description of genes connected to ILP-2
g_out = merge(g_out,aName,by.x="target.gene",by.y="gene_id",sort = FALSE)
g_in = merge(g_in,aName,by.x="regulatory.gene",by.y="gene_id",sort = FALSE)

#Combine regulatory and target genes
g_all <- rbind(g_out,g_in)

#Get full list of genes connected to ILP-2
genes <- unique(c(as.character(g_all$target.gene),as.character(g_all$regulatory.gene)))

#output to igraph
nodes <- aName[aName$gene_id %in% genes,]
nodes <- nodes[!duplicated(nodes$description)]
#nodes <- nodes[,c(2,1)]
edges <- netList_cut[netList_cut$regulatory.gene %in% nodes$gene_id & 
                       netList_cut$target.gene %in% nodes$gene_id,]
colnames(edges) = c("from","to","weight")

#Construct network and plot it
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T) 
plot(net, edge.arrow.size=.4,vertex.label = nodes$description,margin=c(0,0,0,0))
plot(net, edge.arrow.size=.4,vertex.label = nodes$description,margin=c(0,0,0,0),
     layout=layout_in_circle)


#Other things to try with igraph: vary number of genes you keep, change edge weight,
#add color, change the layout if desired. 
#Check out the following link for examples:

#http://www.r-graph-gallery.com/248-igraph-plotting-parameters/

#especially section 4 and beyond:
#http://kateto.net/networks-r-igraph



