#important juniper functions 
#first created 2025 mar 31
#last updated 2025 apr 28
#starting new file which reverts to using HPD to compute CI and p 


# install.packages("cli") # had to manually delete cli.dll and install with install.packages to update to adequate version
# install.packages("devtools")
# install.packages("lmPerm")
# install.packages("tibble")
library(devtools)
# devtools::install_github('xavierdidelot/TransPhylo')
# devtools::install_github("broadinstitute/juniper0", dependencies = FALSE) #because t kept failing on TransPhylo install
library(coda)
library(juniper0)
library(stringr)
library(ggplot2)
library(cowplot)
library(plyr)
library(data.table)
library(dplyr)
library(lmPerm) #maybe for permutation test
library(tidyr)
library(tibble)

my_theme = theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        text = element_text(size =12,colour = "black"),
        axis.text.x = element_text(size=12,colour = "black"),
        axis.text.y = element_text(size=12,colour = "black"),
        legend.text=element_text(size=12),
        axis.line = element_line(),
        plot.margin = margin(1,1,1,1, "cm")) 

Rcpp::sourceCpp("cpp_subroutines.cpp")


independent_downsample = function(res, burnin_pct = 0.2){
  loglik = res[[1]][(1+length(res[[1]])*burnin_pct):length(res[[1]])]
  ess = floor(unname(coda::effectiveSize(loglik)))
  sample_times = floor(seq(1+length(res[[1]])*burnin_pct, length(res[[1]]), length.out=ess))
  my_res = list()
  my_res[[1]] = res[[1]][sample_times]
  my_res[[2]] = res[[2]][sample_times]
  my_res[[3]] = res[[3]]
  my_res[[4]] = res[[4]]
  my_res[[5]] = res[[5]]
  
  return(my_res)
}

get_n_kids = function(res){
  ess = length(res[[1]])
  names = res[[3]][-1]
  rooted = res[[4]]
  s_max = res[[5]]
  n_obs = length(names)
  
  ts = c()
  n_kids = c()
  Rs = c()
  pis = c()
  
  for (s in 1:ess){
    ts = c(ts, list(sapply(res[[2]][[s]]$seq[2:(n_obs+1)],function(v){v[1]})))
    Rs = c(Rs, res[[2]][[s]]$R)
    pis = c(pis, res[[2]][[s]]$pi)
    
    iter_n_kids = c()
    for(host in 2:(n_obs+1)){
      iter_n_kids = c(iter_n_kids, sum(res[[2]][[s]]$h == host, na.rm=T))
    }
    n_kids = c(n_kids, list(iter_n_kids))
  }
  average_n_kids = data.frame(case = names, 
                              n_kids = Reduce(`+`, n_kids)/length(n_kids))
  return(average_n_kids)
}

probability_of_transmission = function(res, metadata, variable, name_var="case", list=FALSE){
  if(!list){
    distributions = transmission_by_variable_distributions(res,metadata,variable,name_var)[[1]]
    results = list()
    for(level in 1:length(distributions)){
      results[[level]] = 1-sapply(distributions[[level]], function(x){x[1]})
    }
    names(results) = names(distributions)
    results_df = stack(results)
    return(results_df)
  }
  if(list){
    all_object_results = list()
    pb = txtProgressBar(min = 0,max=length(res),initial=0,style=3) 
    for(object in 1:length(res)){
      setTxtProgressBar(pb,object)
      a_res = res[[object]]
      distributions = transmission_by_variable_distributions(a_res,metadata,variable,name_var)[[1]]
      results = list()
      for(level in 1:length(distributions)){
        results[[level]] = 1-sapply(distributions[[level]], function(x){x[1]})
      }
      names(results) = names(distributions)
      all_object_results[[object]] = results
    }
    close(pb)
    if(is.null(names(res))){
      names(all_object_results) = seq(1,length(res))
    }
    else{
      names(all_object_results) = names(res)
    }
    all_results_df_list = lapply(all_object_results, stack)
    all_results_df = bind_rows(all_results_df_list, .id="source")
    colnames(all_results_df) = c("object","prob","level")
    return(all_results_df)
  }
}

transmission_by_variable_distributions = function(res,metadata,variable,name_var = "case",drop=NULL){ 
  ess = length(res[[1]])
  names = res[[3]][-1]
  rooted = res[[4]]
  s_max = res[[5]]
  n_obs = length(names)
  
  metadata_variable_data = metadata[[variable]][match(names, metadata[[name_var]])] 
  if(length(which(metadata_variable_data==""))>0 | length(which(is.na(metadata_variable_data)))>0){
    warning("metadata contains missing values, refactoring missing data as new 'missing' category")
    metadata_variable_data[which(metadata_variable_data=="")]="missing"
    metadata_variable_data[which(is.na(metadata_variable_data))]="missing"
  }
  metadata_variable_levels = unique(metadata_variable_data)
  n_metadata_variable_levels = length(unique(metadata_variable_levels))
  
  ts = c()
  n_kids = c()
  Rs = c()
  pis = c()
  
  for (s in 1:ess){
    ts = c(ts, list(sapply(res[[2]][[s]]$seq[2:(n_obs+1)],function(v){v[1]})))
    Rs = c(Rs, res[[2]][[s]]$R)
    pis = c(pis, res[[2]][[s]]$pi)
    
    iter_n_kids = c()
    for(host in 2:(n_obs+1)){
      iter_n_kids = c(iter_n_kids, sum(res[[2]][[s]]$h == host, na.rm=T))
    }
    n_kids = c(n_kids, list(iter_n_kids))
  }
  
  mean_trans_by_variable = lapply(1:n_metadata_variable_levels, function(x) numeric(ess))
  names(mean_trans_by_variable) = metadata_variable_levels
  p_trans_by_variable = lapply(1:n_metadata_variable_levels, function(x) rep(list(rep(0, 7)), ess))
  names(p_trans_by_variable) = metadata_variable_levels
  
  for (i in 1:ess) {
    mean_trans <- rep(0, n_obs)
    min_t <- min(ts[[i]])
    wbar0 <- rev(wbar(min_t - 1, 0, Rs[i] * 0.5 / (1 - 0.5), 1 - 0.5, pis[i], 5, 1, 5, 1, 0.1))
    rho <- Rs[i]
    psi <- 0.5
    ws <- wbar0[round(-ts[[i]]/0.1)] 
    norms <- rep(0, n_obs)
    for (j in 1:n_obs) {
      norms[j] <- exp(alpha(n_kids[[i]][j], psi, rho, ws[j])) # Takes in wbar on log scale
    }
    for (k in 1:7) {
      ps <- choose(k-1, n_kids[[i]]) * dnbinom(k-1, rho, psi) * exp(ws)^(k - 1 - n_kids[[i]]) / norms
      for(level in metadata_variable_levels){
        level_indexes = which(metadata_variable_data == level)
        p_trans_by_variable[[level]][[i]][k] = mean(ps[level_indexes]) #need to get index of cases which match current level 
      }
    }
    means <- rep(0, n_obs)
    for (j in 1:n_obs) {
      means[j] <- exp(alpha2(n_kids[[i]][j], psi, rho, ws[j])) / norms[j] # Takes in wbar on log scale
    }
    for(level in metadata_variable_levels){
      level_indexes = which(metadata_variable_data == level)
      mean_trans_by_variable[[level]][i] = mean(means[level_indexes]) #need to get index of cases which match current level 
    }
  }
  p_trans_by_variable = p_trans_by_variable[!(names(p_trans_by_variable) %in% drop)]
  mean_trans_by_variable = mean_trans_by_variable[!(names(mean_trans_by_variable) %in% drop)]
  return(list(p_trans_by_variable, mean_trans_by_variable))
}


########################################

plot_transmission_by_variable = function(transmission_list, 
                                         title = "", 
                                         variable_name = "variable", 
                                         legend_name = NULL,
                                         legend_pos = NULL,
                                         show_means = TRUE, 
                                         type="p",
                                         colors = c("#4CB3A2", "#E3B861")) {

  # Use variable_name as default legend title if legend_name not provided
  #if(is.null(legend_name)) legend_name <- variable_name

  # default legend pos is right
  if(is.null(legend_pos)) legend_pos <- "right"

  # Prepare plot data for probabilities
  plot_data_p = data.frame()
  for(level in 1:length(transmission_list[[1]])){
    plot_data_p = rbind(plot_data_p, data.frame(
      variable = names(transmission_list[[1]])[level],
      transmissions = 1 - sapply(transmission_list[[1]][[level]], function(x) x[1])
    ))
  }

  probability_plot = ggplot(data = plot_data_p, mapping = aes(x = transmissions, fill = variable, color = variable)) +
    geom_density(color = NA, alpha = 0.6) +
    xlab("Mean Probability of Initiating at Least One Transmission") +
    ylab("Density of Number\nIterations (Downsampled)") +
    ggtitle(title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.position = legend_pos,
          axis.line = element_line(),
          plot.margin = margin(1,1,1,1, "cm")) +
    labs(color = legend_name, fill = legend_name) +
    scale_x_continuous(expand = c(0, 0), labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(name = legend_name, 
                       values = c("Yes" = colors[1], "No" = colors[2],
                                  #"yes" = colors[1], "no/unknown" = colors[2], "no.unknown" = colors[2],
                                  "No/Unknown" = colors[2], "No.Unknown" = colors[2],
                                  "Female" = colors[1], "Male" = colors[2])) +
    scale_fill_manual(name = legend_name, 
                      values = c("Yes" = colors[1], "No" = colors[2],
                                 #"yes" = colors[1], "no/unknown" = colors[2], "no.unknown" = colors[2],
                                  "No/Unknown" = colors[2], "No.Unknown" = colors[2],
                                 "Female" = colors[1], "Male" = colors[2]))

  if(show_means){
    means = ddply(plot_data_p, "variable", summarise, variable_mean=mean(transmissions))
    probability_plot = probability_plot + geom_vline(data=means, linewidth=1.3,
                                                     aes(xintercept=variable_mean, color=variable), linetype="dashed")
  }

  # Prepare plot data for number of transmissions
  mean_transmissions = as.data.frame(transmission_list[[2]])
  plot_data_n = gather(mean_transmissions, variable, transmissions, 1:dim(mean_transmissions)[2])

  transmission_plot = ggplot(data = plot_data_n, mapping = aes(x = transmissions, fill = variable, color = variable)) +
    geom_density(color = NA, alpha = 0.6) +
    xlab("Mean Number Transmissions Initiated") +
    ylab("Density of Number\nIterations (Downsampled)") +
    ggtitle(title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.position = legend_pos,
          axis.line = element_line(),
          plot.margin = margin(1,1,1,1, "cm")) +
    labs(color = legend_name, fill = legend_name) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(name = legend_name, 
                       values = c("Yes" = colors[1], "No" = colors[2],
                                  #"yes" = colors[1], "no/unknown" = colors[2], "no.unknown" = colors[2],
                                  "No/Unknown" = colors[2], "No.Unknown" = colors[2],
                                  "Female" = colors[1], "Male" = colors[2])) +
    scale_fill_manual(name = legend_name, 
                      values = c("Yes" = colors[1], "No" = colors[2],
                                 #"yes" = colors[1], "no/unknown" = colors[2], "no.unknown" = colors[2],
                                 "No/Unknown" = colors[2], "No.Unknown" = colors[2],
                                 "Female" = colors[1], "Male" = colors[2]))

  if(show_means){
    means = ddply(plot_data_n, "variable", summarise, variable_mean=mean(transmissions))
    transmission_plot = transmission_plot + geom_vline(data=means, linewidth=1.3,
                                                       aes(xintercept=variable_mean, color=variable), linetype="dashed")
  }

  if(type == "p") return(probability_plot)
  if(type == "n") return(transmission_plot)
}

# plot_transmission_by_variable = function(transmission_list, title = "", variable_name = "variable", legend_name=legend_name, show_means = TRUE, type="p"){
#   plot_data_p = data.frame()
#   for(level in 1:length(transmission_list[[1]])){
#     plot_data_p = rbind(plot_data_p, data.frame(variable = names(transmission_list[[1]])[level],
#                                                 transmissions = 1-sapply(transmission_list[[1]][[level]], function(x){return(x[1])})))
#   }
#   probability_plot = ggplot(data = plot_data_p, mapping = aes(x = transmissions, fill = variable, color = variable)) +
#     # geom_histogram(position = "identity", color = NA, alpha = 0.6) +
#     geom_density(color = NA, alpha = 0.6) +
#     xlab("Mean Probability of Initiating at Least One Transmission") +
#     # ylab("Number Iterations (Downsampled)") +
#     ylab("Density of Number\nIterations (Downsampled)") +
#     ggtitle(title) +
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           legend.position = "right",
#           axis.line = element_line(),
#           plot.margin = margin(1,1,1,1, "cm")) +
#     labs(color = variable_name, fill = variable_name) + 
#     # scale_x_continuous(expand = c(0, 0), labels = scales::percent, limits = c(0, 1)) +
#     scale_x_continuous(expand = c(0, 0), labels = scales::percent, limits = c(0, 1)) + # limits = c(0.25, 0.75)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     scale_color_manual(name = variable_name, values = c("Yes" = "#40B0A6", "No" = "#E1BE6A",
#                                                         "yes" = "#40B0A6", "no/unknown" = "#E1BE6A", "no.unknown" = "#E1BE6A",
#                                                         "Female" = "#40B0A6", "Male" = "#E1BE6A")) +
#     scale_fill_manual(name = variable_name, values = c("Yes" = "#40B0A6", "No" = "#E1BE6A",
#                                                        "yes" = "#40B0A6", "no/unknown" = "#E1BE6A", "no.unknown" = "#E1BE6A",
#                                                        "Female" = "#40B0A6", "Male" = "#E1BE6A"))
#   if(show_means){
#     means = ddply(plot_data_p, "variable", summarise, variable_mean=mean(transmissions))
#     probability_plot = probability_plot + geom_vline(data=means, linewidth=1.3, 
#                                                      aes(xintercept=variable_mean, color=variable), linetype="dashed")
#   }
#   mean_transmissions = as.data.frame(transmission_list[[2]])
#   plot_data_n = gather(mean_transmissions, variable, transmissions, 1:dim(mean_transmissions)[2])
  
#   transmission_plot = ggplot(data = plot_data_n, mapping = aes(x = transmissions, fill = variable, color = variable)) +
#     # geom_histogram(position = "identity", color = NA, alpha = 0.6) +
#     geom_density(color = NA, alpha = 0.6) +
#     xlab("Mean Number Transmissions Initiated") +
#     # ylab("Number Iterations (Downsampled)") +
#     ylab("Density of Number\nIterations (Downsampled)") +
#     ggtitle(title) +
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           legend.position = "right",
#           axis.line = element_line(),
#           plot.margin = margin(1,1,1,1, "cm")) +
#     labs(color = variable_name, fill = variable_name) + 
#     # scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
#     scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) + # limits = c(0.5, 1.5)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     scale_color_manual(name = variable_name, values = c("Yes" = "#40B0A6", "No" = "#E1BE6A",
#                                                         "yes" = "#40B0A6", "no/unknown" = "#E1BE6A", "no.unknown" = "#E1BE6A",
#                                                         "Female" = "#40B0A6", "Male" = "#E1BE6A")) +
#     scale_fill_manual(name = variable_name, values = c("Yes" = "#40B0A6", "No" = "#E1BE6A",
#                                                        "yes" = "#40B0A6", "no/unknown" = "#E1BE6A", "no.unknown" = "#E1BE6A",
#                                                        "Female" = "#40B0A6", "Male" = "#E1BE6A"))
  
#   if(show_means){
#     means = ddply(plot_data_n, "variable", summarise, variable_mean=mean(transmissions))
#     transmission_plot = transmission_plot + geom_vline(data=means, linewidth=1.3, 
#                                                        aes(xintercept=variable_mean, color=variable), linetype="dashed")
#   }
#   if(type == "p"){
#     return(probability_plot)
#   }
#   if(type == "n"){
#     return(transmission_plot)
#   }
# }

retrieve_means_for_transmission_by_variable = function(transmission_list, title = "", variable_name = "variable", type="p"){
  plot_data_p = data.frame()
  for(level in 1:length(transmission_list[[1]])){
    plot_data_p = rbind(plot_data_p, data.frame(variable = names(transmission_list[[1]])[level],
                                                transmissions = 1-sapply(transmission_list[[1]][[level]], function(x){return(x[1])})))
  }
  means_p = ddply(plot_data_p, "variable", summarise, variable_mean=mean(transmissions))
  
  mean_transmissions = as.data.frame(transmission_list[[2]])
  plot_data_n = gather(mean_transmissions, variable, transmissions, 1:dim(mean_transmissions)[2])
  means_n = ddply(plot_data_n, "variable", summarise, variable_mean=mean(transmissions))
  if(type == "p"){
    return(means_p)
  }
  if(type == "n"){
    return(means_n)
  }
}

# plot_transmission_by_variable = function(transmission_list, title = "none", variable_name = "variable", show_means = TRUE, type="p"){
#   plot_data_p = data.frame()
#   for(level in 1:length(transmission_list[[1]])){
#     plot_data_p = rbind(plot_data_p, data.frame(variable = names(transmission_list[[1]])[level],
#                                                 transmissions = 1-sapply(transmission_list[[1]][[level]], function(x){return(x[1])})))
#   }
#   probability_plot = ggplot(data = plot_data_p, mapping = aes(x = transmissions, fill = variable, color = variable)) +
#     geom_density(color = NA, alpha = 0.6) +
#     xlab("Probability of Transmission") +
#     ylab("Density") +
#     ggtitle(title) + 
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size =12,colour = "black"),
#           axis.text.x = element_text(size=12,colour = "black"),
#           axis.text.y = element_text(size=12,colour = "black"),
#           legend.text=element_text(size=12),
#           axis.line = element_line(),
#           plot.margin = margin(1,1,1,1, "cm")) +
#     labs(color = variable_name, fill = variable_name) + 
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0))
#   if(show_means){
#     means = ddply(plot_data_p, "variable", summarise, variable_mean=mean(transmissions))
#     probability_plot = probability_plot + geom_vline(data=means, linewidth=1.3, 
#                                                      aes(xintercept=variable_mean, color=variable), linetype="dashed")
#   }
#   mean_transmissions = as.data.frame(transmission_list[[2]])
#   plot_data_n = gather(mean_transmissions, variable, transmissions, 1:dim(mean_transmissions)[2])
#   transmission_plot = ggplot(data = plot_data_n, mapping = aes(x = transmissions, fill = variable, color = variable)) +
#     geom_density(color = NA, alpha = 0.6) +
#     xlab("Transmissions per Case") +
#     ylab("Density") +
#     ggtitle(title) + 
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           text = element_text(size =12,colour = "black"),
#           axis.text.x = element_text(size=12,colour = "black"),
#           axis.text.y = element_text(size=12,colour = "black"),
#           legend.text=element_text(size=12),
#           axis.line = element_line(),
#           plot.margin = margin(1,1,1,1, "cm")) +
#     labs(color = variable_name, fill = variable_name) + 
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0))
#   if(show_means){
#     means = ddply(plot_data_n, "variable", summarise, variable_mean=mean(transmissions))
#     transmission_plot = transmission_plot + geom_vline(data=means, linewidth=1.3, 
#                                                        aes(xintercept=variable_mean, color=variable), linetype="dashed")
#   }
#   if(type == "p"){
#     return(probability_plot)
#   }
#   if(type == "n"){
#     return(transmission_plot)
#   }
# }

transmission_by_variable_comparison = function(res,metadata,variable,name_var="case",
                                               method="b",comparator=NA,drop=NULL){
 
  list = TRUE
  if(length(res) >= 4 && length(res[[4]]) == 1){
    list = FALSE
  }
  if(method=="p"){
    method_ind_list = 1
  }
  if(method=="n"){
    method_ind_list = 2
  }
  if(method=="b"){
    method_ind_list = c(1,2)
  }
  if(list){
    num_results = length(res)
  }
  if(!list){
    num_results = 1
  }
  if(list & is.null(names(res))){
    names(res) = seq(1:length(res))
  }
  pb = txtProgressBar(min = 0,max=num_results,initial=0,style=3) 
  all_results =  vector("list", 2)
  combined_data = list(list(),list())
  for(r in 1:num_results){
    setTxtProgressBar(pb,r)
    if(list){
      a_res = res[[r]]
      object_name = names(res)[r]
    }
    else{
      a_res = res
      object_name = 1
    }
    transmisison_results = transmission_by_variable_distributions(a_res,metadata,variable,name_var,drop)
    for(method_ind in method_ind_list){
      results = data.frame()
      data = transmisison_results[[method_ind]]
      res_comparator = comparator
      if(!is.na(comparator) & !(comparator %in% names(data))){
        warning("comparator group not found, picking arbitrarily")
        res_comparator = NA 
      }
      if(is.na(comparator)){
        comparator = names(data)[1]
        res_comparator = names(data)[1] 
      }
      comparator_data = data[[res_comparator]]
      remaining_data = data[names(data) != res_comparator]
      if(length(remaining_data)==0){
        warning("all transmissions occured in one variable status, skipping this results object")
        next 
      }
      for(i in 1:length(remaining_data)){
        if(method_ind==1){ #probabilities
          level_name = paste(names(remaining_data)[i],"_vs_",res_comparator,sep="")
          obs_dif = (1-sapply(remaining_data[[i]], function(x){return(x[1])})) -
                           (1-sapply(comparator_data, function(x){return(x[1])}))
          results = rbind(results, data.frame(object = object_name,
                                              comparison = level_name, 
                                              lower = unname(quantile(obs_dif, 0.025)), 
                                              mean = mean(obs_dif), 
                                              upper = unname(quantile(obs_dif, 0.975)),
                                              pval = min(mean(obs_dif>0), mean(obs_dif<0))))
                                              # pval = format(min(mean(obs_dif > 0), mean(obs_dif < 0)), scientific = TRUE)))
          combined_data[[method_ind]][[level_name]] = c(combined_data[[method_ind]][[level_name]],obs_dif)
        }
        else{ #mean transmissions 
          level_name = paste(names(remaining_data)[i],"_vs_",res_comparator,sep="")
          obs_dif = remaining_data[[i]] - comparator_data
          results = rbind(results, data.frame(object = object_name,
                                              comparison = level_name, 
                                              lower = unname(quantile(obs_dif, 0.025)), 
                                              mean = mean(obs_dif), 
                                              upper = unname(quantile(obs_dif, 0.975)),
                                              pval = min(mean(obs_dif>0), mean(obs_dif<0))))
                                              # pval = format(min(mean(obs_dif > 0), mean(obs_dif < 0)), scientific = TRUE)))
          combined_data[[method_ind]][[level_name]] = c(combined_data[[method_ind]][[level_name]],obs_dif)
        }
      }
      all_results[[method_ind]] = rbind(all_results[[method_ind]], results)
    }
  }
  close(pb)
  if(list){
    for(method_ind in method_ind_list){
      for(l in 1:length(combined_data[[method_ind]])){
        all_results[[method_ind]] = rbind(all_results[[method_ind]], data.frame(object = "combined",
                                                                                comparison = names(combined_data[[method_ind]])[l],
                                                                                lower = unname(quantile(combined_data[[method_ind]][[l]],0.025)),
                                                                                mean = mean(combined_data[[method_ind]][[l]]),
                                                                                upper = unname(quantile(combined_data[[method_ind]][[l]],0.975)),
                                                                                pval = min(mean(obs_dif>0), mean(obs_dif<0))))
                                                                                # pval = format(min(mean(obs_dif > 0), mean(obs_dif < 0)), scientific = TRUE)))
      }
    }
  }
  if(method=="p"){
    all_results = all_results[[1]]
  }
  if(method=="n"){
    all_results = all_results[[2]]
  }
  if(method=="b"){
    names(all_results) = c("p","n")
  }
  return(all_results)
}

plot_transmission_by_variable_comparison = function(comparison_results, drop=NULL){
  n_data = comparison_results$n
  p_data = comparison_results$p
  if(!is.null(drop)){
    n_data = n_data[!Reduce(`|`, lapply(drop, function(x) grepl(x, n_data$comparison))), ]
    p_data = p_data[!Reduce(`|`, lapply(drop, function(x) grepl(x, p_data$comparison))), ]
    
  }
  n_plot = ggplot(n_data, aes(x = object, y = mean, colour = comparison)) + 
    geom_errorbar(aes(ymax = lower, ymin = upper), position = "dodge") + 
    # labs(y="mean difference in number of transmissions vs female") + 
    geom_point(position = position_dodge(0.9)) + my_theme + 
    ggtitle("number of transmissions") + 
    geom_hline(yintercept = 0, linetype="dashed") + coord_flip() +
    theme(legend.position = "top", legend.direction = "vertical")
  p_plot = ggplot(p_data, aes(x = object, y = mean, colour = comparison)) + 
    geom_errorbar(aes(ymax = lower, ymin = upper), position = "dodge") + 
    # labs(y="mean difference in number of transmissions vs female") + 
    geom_point(position = position_dodge(0.9)) + my_theme + 
    ggtitle("probability of any transmission") + 
    geom_hline(yintercept = 0, linetype="dashed") + coord_flip() +
    theme(legend.position = "top", legend.direction = "vertical")
  out_plot = plot_grid(n_plot, p_plot, ncol = 2)
  return(out_plot)
}

per_iteration_analysis = function(results, metadata,variable,name_var = "case"){
  n_reps <- length(results[[1]])
  names <- results[[3]]
  s_max <- results[[5]]
  n_obs <- length(names)
  names_df = data.table(case = names)
  included_meta = metadata[names_df, on = name_var, nomatch=NULL]
  all_iterations = data.table()
  for (i in 1:n_reps) {
    h <- results[[2]][[i]]$h
    n <- results[[2]][[i]]$n
    w <- sapply(results[[2]][[i]]$seq, length) - 1
    h_obs <- c()
    for (j in 2:n_obs) {
      h_obs[j] <- h[j]
      while (h_obs[j] > n_obs) {
        h_obs[j] <- h[h_obs[j]]
      }
    }
    trans <- cbind(h[2:n], 2:n)
    any_trans <- data.table(case=names[unique(trans[trans[, 1] <= n_obs, ][,1])])
    iteration_meta = metadata[any_trans, on= name_var, nomatch=NULL]
    iteration_results = iteration_meta[,.N,by = variable][included_meta[,.N,by = variable], on = variable]
    iteration_results$prop = iteration_results$N/iteration_results$i.N
    all_iterations = rbind(all_iterations, iteration_results[,c(1,4)])
  }
  return(all_iterations)
}

get_transmission_pairs = function(results){
  n_reps <- length(results[[1]])
  names <- results[[3]]
  s_max <- results[[5]]
  n_obs <- length(names)
  all_pairs = data.frame()
  for (i in 1:n_reps) {
    h <- results[[2]][[i]]$h
    n <- results[[2]][[i]]$n
    w <- sapply(results[[2]][[i]]$seq, length) - 1
    h_obs <- c()
    for (j in 2:n_obs) {
      h_obs[j] <- h[j]
      while (h_obs[j] > n_obs) {
        h_obs[j] <- h[h_obs[j]]
      }
    }
    trans <- cbind(h[2:n], 2:n)
    direct_trans <- trans[trans[, 1] <= n_obs & trans[, 2] <= 
                            n_obs & w[2:n] == 0, ]
    all_pairs = rbind(all_pairs, direct_trans)
  }
  counts <- aggregate(rep(1, nrow(all_pairs)), by = list(a = all_pairs$V1, b = all_pairs$V2), FUN = sum)
  pair_outs = data.frame(case1 = names[counts$a], case2 = names[counts$b], weight = counts$x/n_reps)
  return(pair_outs)
}



get_infection_dates = function(results){
  if(length(results) >= 4 && length(results[[4]]) == 1){
    n_reps <- length(results[[1]])
    s_max <- results[[5]]
    names <- results[[3]]
    n_obs = length(names)
    ts <- matrix(ncol = n_obs - 1, nrow = 0)
    for (i in 1:n_reps) {
      ts <- rbind(ts, sapply(results[[2]][[i]]$seq[2:n_obs], 
                             function(v) {
                               v[1]
                             }))}
    mean_ts <- colMeans(ts) + s_max
    lower_ts <- apply(ts, 2, function(v) {
      quantile(v, 0.025)
    }) + s_max
    upper_ts <- apply(ts, 2, function(v) {
      quantile(v, 0.975)
    }) + s_max
    out = data.frame(case = names[2:n_obs],
                     lower = lower_ts,
                     mean = mean_ts,
                     upper = upper_ts)
    
  }
  else{
    out = do.call(rbind, lapply(results, get_infection_dates))
  }
  return(out)
}



