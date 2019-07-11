setwd("/work-zfs/abattle4/zeinab/prog/scnet")

library(feather)
library(argparser)
library(ioutil)
library(genomicsutil)
library(miscutil)
library(ggplot2)
#library(ape)
#library(mappabilityutil)
#library(infotheo)
#library(arules)
library(igraph)
#library(lme4)
#library(reshape2)
#library(limma)
#library(RColorBrewer)

source('/work-zfs/abattle4/zeinab/prog/scnet/ashis/aggregated_mst.R')
source('/work-zfs/abattle4/zeinab/prog/scnet/ashis/coexpression_network.R')
source('/work-zfs/abattle4/zeinab/prog/scnet/ashis/network.R')

##### parse arguments ######
args <- arg_parser("program");

args_bulk <- add_argument(args, "-expr", help="expression data", default="/work-zfs/abattle4/ashis/progres/scproc/ye_data_processed_backup/CD8_T_cells/CD8_T_cells_pseudobulk_hvg500_normalized_pc_corrected.feather")
args_sc <- add_argument(args, "-expr", help="expression data", default="/work-zfs/abattle4/ashis/progres/scproc/ye_data_processed_backup/CD8_T_cells/CD8_T_cells_hvg500_normalized_pc_corrected.feather")


args = args_bulk# args_sc#    
args <- add_argument(args, "-o", help="output dir", default="results/")

argv = parse_args(args)
expr_fn = argv$expr
out_dir = argv$o

# data_sep = ' '
n_genes = 500 #500 is max?
# n.max.edges = 500
# string_score_threshold = 400
# ppi_enrichment_iter = 100
# 

string_dir = paste0(out_dir, '/string')
if(!dir.exists(out_dir))
  stop(sprintf('output directory does not exist: %s.', out_dir))
if(!dir.exists(string_dir))
  dir.create(string_dir)

### read and process data
if(endsWith(x = expr_fn, suffix = '.feather')){
  expr_df_raw = read_feather_df(fn = expr_fn, rownames.col = 1)
} else {
  expr_df_raw = read_df(expr_fn, sep = data_sep)
}
dim(expr_df_raw)


if(any(dim(expr_df_raw)<10))
  warning(sprintf('small number of features/samples in expression matrix: %s x %s', nrow(expr_df_raw), ncol(expr_df_raw)))


expr_df = filter_expr_by_coeff_of_variation(expr.df = expr_df_raw, raw.df = expr_df_raw, n = n_genes, min.mean = -Inf, min.var = 1e-6)


alpha_list = seq(0.001, 1, 0.1)

edge_count_list = c()
r_squared_list = c()
for (alpha in alpha_list){
  
  net <- get_glasso_net(expr_df, lambda = alpha)
  prec =(abs(net))
  cor = solve(prec)
  g = graph_from_adjacency_matrix(cor, weighted=T)
  g_bulk = g
  # List of degrees
  G.degrees <- degree(g)
  
  # Let's count the frequencies of each degree
  G.degree.histogram <- as.data.frame(table(G.degrees))
  
  # Need to convert the first column to numbers, otherwise
  # the log-log thing will not work (that's fair...)
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
  
  degree = G.degree.histogram[,1]
  count = G.degree.histogram[,2]
  r_squared = summary(lm(count ~ degree))$r.squared
  r_squared_list = c(r_squared_list, r_squared)
  
  edge_count = gsize(g)
  edge_count_list = c(edge_count_list, edge_count)
  
  if (r_squared > 0.8){
    break
  }
  
 
  
}
print(paste0("bulk_edge_count: ", edge_count))
print(paste0("bulk_r_squared: ", r_squared))
print(paste0("bulk_alpha: ", alpha))





###SC
#####
##### parse arguments ######
args =  args_sc#    
args <- add_argument(args, "-o", help="output dir", default="results/")

argv = parse_args(args)
expr_fn = argv$expr
out_dir = argv$o

# data_sep = ' '
n_genes = 500 #500 is max?
# n.max.edges = 500
# string_score_threshold = 400
# ppi_enrichment_iter = 100
# 

string_dir = paste0(out_dir, '/string')
if(!dir.exists(out_dir))
  stop(sprintf('output directory does not exist: %s.', out_dir))
if(!dir.exists(string_dir))
  dir.create(string_dir)

### read and process data
if(endsWith(x = expr_fn, suffix = '.feather')){
  expr_df_raw = read_feather_df(fn = expr_fn, rownames.col = 1)
} else {
  expr_df_raw = read_df(expr_fn, sep = data_sep)
}
dim(expr_df_raw)


if(any(dim(expr_df_raw)<10))
  warning(sprintf('small number of features/samples in expression matrix: %s x %s', nrow(expr_df_raw), ncol(expr_df_raw)))


expr_df = filter_expr_by_coeff_of_variation(expr.df = expr_df_raw, raw.df = expr_df_raw, n = n_genes, min.mean = -Inf, min.var = 1e-6)


alpha_list = seq(0.001, 1, 0.1)

edge_count_list = c()
r_squared_list = c()
for (alpha in alpha_list){
  
  net <- get_glasso_net(expr_df, lambda = alpha)
  prec =(abs(net))
  cor = solve(prec)
  g = graph_from_adjacency_matrix(cor, weighted=T)
  g_sc = g
  # List of degrees
  G.degrees <- degree(g)
  
  # Let's count the frequencies of each degree
  G.degree.histogram <- as.data.frame(table(G.degrees))
  
  # Need to convert the first column to numbers, otherwise
  # the log-log thing will not work (that's fair...)
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
  
  degree = G.degree.histogram[,1]
  count = G.degree.histogram[,2]
  r_squared = summary(lm(count ~ degree))$r.squared
  r_squared_list = c(r_squared_list, r_squared)
  
  edge_count = gsize(g)
  edge_count_list = c(edge_count_list, edge_count)
  
  if (r_squared > 0.8){
    break
  }

  
}

print(paste0("sc_edge_count: ", edge_count))
print(paste0("sc_r_squared: ", r_squared))
print(paste0("sc_alpha: ", alpha))

common_g <- graph.intersection(g_sc, g_bulk, keep.all.vertices = FALSE)
common_g_size = gsize(common_g)

print(paste0("intersect_edge_count: ", common_g_size))

common_g_nodes = length(V(common_g))
print(paste0("intersect_number_nodes: ", common_g_nodes))

pdf("./Plots/Intersection_dist.pdf")
plot(degree.distribution(common_g))
dev.off()





