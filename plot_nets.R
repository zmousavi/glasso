
setwd("/work-zfs/abattle4/zeinab/prog/scnet/glasso")


library(ggplot2)
library(igraph)
library(argparser)


source('/work-zfs/abattle4/zeinab/prog/scnet/glasso/scale_free_test.R')
##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-r", help="rna_seq: sc or bulk", default="bulk")
args <- add_argument(args, "-c", help="cell type Dendritic, CD8", default="Dend")
args <- add_argument(args, "-p", help="significance value", default=0.5)
args <- add_argument(args, "-k", help="significance value", default=0)

argv = parse_args(args)
       
cell_type = argv$c
rna_seq = argv$r
significance = argv$p
k_try = argv$k

plot_dir = "/work-zfs/abattle4/zeinab/prog/scnet/glasso/Plots/PL_test"
dir_path_main = "/work-zfs/abattle4/zeinab/prog/scnet/glasso/results/"
#dir_path = paste0(dir_path, cell_type)
file_path = paste0(dir_path_main, rna_seq, '/', cell_type)
result_path = "/work-zfs/abattle4/zeinab/prog/scnet/glasso/results/PL_test"


file_pattern = sprintf("^%s.*_network.rds", cell_type)
filenames <- list.files(path=file_path, pattern=file_pattern, full.names=TRUE)
dataframe = data.frame(matrix(ncol = 5, nrow = 0))
col_names = c("alpha", "edge_count", "r_squared", "slope", "pl_dist")
colnames(dataframe) = col_names 
edge_count_list = c()
r_squared_list = c()
for (file in filenames){
  net =  readRDS(file)
  #prec =(abs(net))
  #cor = solve(prec)

  #dont really need to be computing this:
  adj = ifelse(net!=0 ,1,0)
  diag(adj)= 1
  G.degree.histogram = (as.data.frame(table(rowSums(as.matrix(adj)))))
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
  degree = G.degree.histogram[,1]
  count = G.degree.histogram[,2]
  linear_model = lm(log(count/sum(count)) ~ log(degree))
r_squared = summary(linear_model)$r.squared
slope = summary(linear_model)$coefficients[2]
    
  alpha = as.numeric(strsplit(basename(file), "_")[[1]][3])
  
  #end of unnecessary block
  
######PUT BACK
  ####deg_name = sprintf('%s_%s_%s_deg_JUL30_%s.rds', cell_type, rna_seq, alpha, k_try)
#####saveRDS(G.degree.histogram, file = file.path(result_path, deg_name))  


  graph = as.undirected(graph_from_adjacency_matrix(adj))
  edge_count = gsize(graph)

  pl_res = tryCatch(scalefree_test(graph, significance), error=function(e) NULL)
  pl_dist = pl_res[[1]]
  if (is.null(pl_dist)){
pl_dist = FALSE
} 
 
  df <- data.frame()
  df <- data.frame(alpha, edge_count, r_squared, slope, pl_dist)
  colnames(df) = col_names 
  dataframe = rbind(dataframe, df)
  
  graph_type = sprintf("%s_%s_%s", cell_type, rna_seq, alpha)
#####PUT BACK
########  print (graph_type)
########  plot_net(graph, graph_type, plot_dir, significance)
  
print (dataframe)
    #print("alpha")
    #print (alpha)
    #print (pl_dist)
}


print(dataframe) 
df_name = sprintf('%s_%s_SF_tests_JUL30_%s.rds', cell_type, rna_seq, k_try) 
saveRDS(dataframe, file = file.path(file_path, df_name))



result_df = dataframe
alpha = result_df$alpha_list
r_squared = result_df$r_squared_list 
edge_count = result_df$edge_count_list 

plot_r2 = file.path(plot_dir,sprintf('%s_%s_R2_JUL30_%s.pdf', cell_type, rna_seq, k_try))
title_r2 = sprintf("%s_Power_Law", rna_seq)
pdf(plot_r2)
ggplot() + 
  geom_point(result_df, mapping= aes(x=alpha, y=r_squared, color= pl_dist)) +
  geom_line() +
  geom_hline(yintercept=0.8, color = "red") +
  ggtitle(title_r2)
dev.off()
#print (result_df)


plot_e  = file.path(plot_dir,sprintf('%s_%s_edges_Jul30_%s.pdf', cell_type, rna_seq, k_try))
title_e =  sprintf("%s_Total_edges", rna_seq)
pdf(plot_e)
ggplot() + 
geom_point(result_df, mapping= aes(x=alpha, y=edge_count, color= pl_dist)) +
  geom_line() +
  ggtitle(title_e)
dev.off()
