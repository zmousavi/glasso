
setwd("/work-zfs/abattle4/zeinab/prog/scnet/glasso")

library(feather)
library(argparser)
library(ioutil)
library(genomicsutil)
library(miscutil)
library(ggplot2)
library(igraph)


source('/work-zfs/abattle4/zeinab/prog/scnet/ashis/aggregated_mst.R')
source('/work-zfs/abattle4/zeinab/prog/scnet/ashis/coexpression_network.R')
source('/work-zfs/abattle4/zeinab/prog/scnet/ashis/network.R')

##### parse arguments ######
args <- arg_parser("program");

args <- add_argument(args, "-expr", help="/work-zfs/abattle4/ashis/progres/scproc/ye_data_processed_backup/CD8_T_cells/CD8_T_cells_pseudobulk_biovar_normalized_pc_corrected.feather")
args <- add_argument(args, "-o", help="output dir", default="results/")
args <- add_argument(args, "-f", help="output dataframe", default="df")
args <- add_argument(args, "-a", help="alpha_list", default=c(0.9))


argv = parse_args(args)
expr_fn = argv$expr
out_dir = argv$o
out_file = argv$f
alpha_list = argv$a

n_genes = 500 


string_dir = paste0(out_dir, '/string')
if(!dir.exists(out_dir))
  stop(sprintf('output directory does not exist: %s.', out_dir))
if(!dir.exists(string_dir))
  dir.create(string_dir)

### read and process data
if(endsWith(x = expr_fn, suffix = '.feather')){
  expr_df = read_feather_df(fn = expr_fn, rownames.col = 1)
} else {
  expr_df = read_df(expr_fn, sep = data_sep)
}


if(any(dim(expr_df_raw)<10))
  warning(sprintf('small number of features/samples in expression matrix: %s x %s', nrow(expr_df), ncol(expr_df)))



r_squared_list = c()
edge_count_list = c()
time_list = c()
for (alpha in alpha_list){

  start_time <- Sys.time()
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
  r_squared = summary(lm(count ~ degree))$r.squared
  edge_count = gsize(g)
  edge_count_list = c(edge_count_list, edge_count)
  r_squared_list = c(r_squared_list, r_squared)
  end_time <- Sys.time()
  run_time = end_time - start_time
  time_list = c(time_list, run_time)
  
  rm(g)
  cat("alpha = ", alpha)
  cat("r_squared = ", r_squared)
  cat("edge_count = ", edge_count)
  cat( "\n")
}
r_squared_list
edge_count_list

result_df <- cbind(alpha_list, r_squared_list, edge_count_list, time_list)
result_df = as.data.frame(result_df)

file_path = file.path(string_dir, paste(out_file, ".rds", sep=""))
saveRDS(result_df, file = file_path)


# pdf("./Plots/R2_sc.pdf")
# ggplot(result_df, aes(x=alpha_list, y=r_squared_list)) + geom_point() +
# geom_line() +
# geom_hline(yintercept=0.8, color = "red") + 
# ggtitle("Single Cell Power Law")
# dev.off()
# 
# pdf("./Plots/Edge_sc.pdf")
# ggplot(result_df, aes(x=alpha_list, y=edge_count_list)) + geom_point() +
# geom_line() +
# ggtitle("Single Cell Edge Weight")
# dev.off()
