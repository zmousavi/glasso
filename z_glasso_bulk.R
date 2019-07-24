setwd("/work-zfs/abattle4/zeinab/prog/scnet/glasso")

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

#args_bulk <- add_argument(args, "-expr", help="expression data", default="/work-zfs/abattle4/ashis/progres/scproc/ye_data_processed_backup/CD8_T_cells/CD8_T_cells_pseudobulk_hvg500_normalized_pc_corrected.feather")
#args_sc <- add_argument(args, "-expr", help="expression data", default="/work-zfs/abattle4/ashis/progres/scproc/ye_data_processed_backup/CD8_T_cells/CD8_T_cells_hvg500_normalized_pc_corrected.feather")

args_bulk <- add_argument(args, "-expr", help="expression data", default="/work-zfs/abattle4/ashis/progres/scproc/ye_data_processed_backup/CD8_T_cells/CD8_T_cells_pseudobulk_biovar_normalized_pc_corrected.feather")
args_sc <- add_argument(args, "-expr", help="expression data", default="/work-zfs/abattle4/ashis/progres/scproc/ye_data_processed_backup/CD8_T_cells/CD8_T_cells_biovar_normalized_pc_corrected.feather")


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

#expr_df = filter_expr_by_coeff_of_variation(expr.df = expr_df_raw, raw.df = expr_df_raw, n = n_genes, min.mean = -Inf, min.var = 1e-6)

expr_df = expr_df_raw


alpha_list = seq(0.001, 1, 0.01)

#alpha_list = c(0.3, 0.8)
r_squared_list = c()
edge_count_list = c()
for (alpha in alpha_list){

  net <- get_glasso_net(expr_df, lambda = alpha)
  prec =(abs(net))
  cor = solve(prec)
  g = graph_from_adjacency_matrix(cor, weighted=T)
  # List of degrees
  G.degrees <- degree(g)
  
  # Let's count the frequencies of each degree
  G.degree.histogram <- as.data.frame(table(G.degrees))
  
  # Need to convert the first column to numbers, otherwise
  # the log-log thing will not work (that's fair...)
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
  
  degree = G.degree.histogram[,1]
  count = G.degree.histogram[,2]
  r_squared = summary(lm(log(count) ~ log(degree)))$r.squared
  print(r_squared)
  edge_count = gsize(g)
  edge_count_list = c(edge_count_list, edge_count)
  r_squared_list = c(r_squared_list, r_squared)
  rm (g) 
  cat("alpha = ", alpha)
  cat("r_squared = ", r_squared)
  cat("edge_count = ", edge_count)
  cat( "\n")
}
r_squared_list
edge_count_list
result_df <- cbind(alpha_list, r_squared_list, edge_count_list)
result_df = as.data.frame(result_df)

pdf("./Plots/R2_bulk.pdf")
ggplot(result_df, aes(x=alpha_list, y=r_squared_list)) + geom_point() +
geom_line() +
geom_hline(yintercept=0.8, color = "red") +
ggtitle("Bulk Power Law")
dev.off()
#print (result_df)

pdf("./Plots/Edge_bulk.pdf")
ggplot(result_df, aes(x=alpha_list, y=edge_count_list)) + geom_point() +
geom_line() +
ggtitle("Bulk Edge Weight")
dev.off()




#alpha = 0/raw data
#G <- graph.data.frame(datExpr0)

# List of degrees
#G.degrees <- degree(G)

# Let's count the frequencies of each degree
#G.degree.histogram <- as.data.frame(table(G.degrees))

# Need to convert the first column to numbers, otherwise
# the log-log thing will not work (that's fair...)
#G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])

#degree = G.degree.histogram[,1]
#count = G.degree.histogram[,2]
#r_squared = summary(lm(count ~ degree))$r.squared



# dim(expr_df) 
# #500 x 60291 for sc, 500x for bulk = #500x140 for bulk
# #run glasso
# prec = run_glasso_net_0.3(expr_df)
# cor = solve(prec)
# 
# g = graph_from_adjacency_matrix(cor, weighted=T)
# gsize(g) 
# #166666 for sc, 220932 for bulk


