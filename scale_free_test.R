


#power-law vs alternate distributions
pl_alt_fit <- function(observed_data, dist_test, significance=0.05) {
  library("poweRlaw")
  data_ref <- displ$new(observed_data)
  est <- estimate_xmin(data_ref)
  data_ref$xmin <- est
  
  data_alt_object <-  get(dist_test)
  data_alt = data_alt_object$new(observed_data)
  data_alt$xmin <- est$xmin
  data_alt$pars <- estimate_pars(data_alt)
  
  comp <- compare_distributions(data_ref, data_alt)
  res = list("score" = comp$test_statistic, "sig" = comp$p_one_sided)
  return(res)
  
  plot(data_alt, ylab="CDF")
  lines(data_ref)
  lines(data_alt, col=2, lty=2)
}


scalefree_test <- function(graph, significance){
  require(poweRlaw)
  #get degree sequence
  degree_table = degree(graph)
  degree_array = c()
  for (vertex in 1:length(V(graph))){
    d_v = degree_table[[vertex]]
    degree_array = c(degree_array, d_v)
  }
  
  #construct power law object for data
  #observed_data = degree_array+1
  observed_data = degree_array

  data_pl = displ$new(observed_data)
 print ("DAATA_PL") 
 print (data_pl)
  #estimate xmin
  est <- estimate_xmin(data_pl)
	print("EST")
print(est)
  data_pl$xmin <- est
  print(est$gof)
 pl_lnorm <- "N/A"
pl_exp <- "N/A"
pl_pois <- "N/A"
 
 #how well does power law distribution fit the data 
  if (est$gof > significance){
    power_law = as.logical(FALSE)
    
  } else{
    #how likely is the data drawn from power law?
    print ("bootstrap")
    print("data_pl")
    print (data_pl)

    bs <- bootstrap_p(data_pl)
    bs_p = bs$p
    print ("bootstrap finished")
    print(bs_p)
    if (bs_p < significance) {
      power_law = as.logical(FALSE)
      print("here")
    } else{
      print ("pl_lnorm start")
      pl_lnorm = pl_alt_fit(observed_data, "dislnorm")
      print("pl_lnorm finished")
      pl_exp = pl_alt_fit(observed_data, "disexp")
      pl_pois = pl_alt_fit(observed_data, "dispois")
      
      if ((pl_lnorm$score < 0 & pl_lnorm$sig < significance) | (pl_exp$score < 0 & pl_exp$sig < significance ) |  (pl_pois$score < 0 & pl_pois$sig < significance)){
        power_law = as.logical(FALSE)
      } else {
        power_law = as.logical(TRUE)
      }
    }
    
  }
  
  res = list(power_dist = power_law, lnorm_dist = pl_lnorm, exp_dist = pl_exp, pois_dist = pl_pois)
 
   
  return(res)
  }



plot_net = function(graph, graph_type, out_dir, significance){
  require(igraph)
  number_nodes = length(V(graph))
degree_table = degree(graph)
degree_array = c()
for (vertex in 1:length(V(graph))){
  d_v = degree_table[[vertex]]
  degree_array = c(degree_array, d_v)
}

#compute R2 & slope
G.degree.histogram = (as.data.frame(table(degree(graph))))
G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
degree = G.degree.histogram[,1]
count = G.degree.histogram[,2]
linear_model = lm(log(count/sum(count)) ~ log(degree))
r_squared = summary(linear_model)$r.squared
slope = summary(linear_model)$coefficients[2]
  
dist_name = sprintf('%s_net_%s.pdf' , graph_type, number_nodes)
pdf(file.path(out_dir, dist_name))
V(graph)$color = "lightblue"
plot.igraph(graph, vertex.label=NA, vertex.size = 8)
title(bquote(atop(.(graph_type),
                  atop("Number of Vertices =" ~ .(number_nodes), 
                       atop(R^2  ==  .(r_squared), "slope =" ~.(slope) ) ))))
dev.off()

net_name = sprintf('%s_distribution_%s.pdf' , graph_type, number_nodes)
pdf(file.path(out_dir, net_name))
hist(degree_array, freq=TRUE, col="lightblue", main="", xlab="Vertex Degrees", ylab="Frequency")
title(bquote(atop("?" ,
                  atop("Number of Vertices =" ~ .(number_nodes), 
                       atop(R^2  ==  .(r_squared), "slope =" ~.(slope) ) ))))
dev.off()


}
