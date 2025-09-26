##### loading data ==================================================================================================

#---------- phylogeny

caudata_tree <- ape::read.nexus("./caudata_tree_stewart_wiens_2024_arranged_lad.nexus")

#---------- species data

caudata_data <- read.csv(file = "./data/caudata_data/stewart_wiens_2024/caudata_data_arranged_lad.csv",
                         header = TRUE,
                         sep = ",")

##### run bamm ============================================================================================================

# Generate a control file with the priors.
generateControlFile(file = './bamm/divcontrol.txt', 
                    type = 'diversification', 
                    params = list(
                      treefile = 'caudata_tree_stewart_wiens_2024_arranged_lad.tre', # the tree file have to be in the same directory as bamm
                      globalSamplingFraction = '1',
                      numberOfGenerations = '100000000',
                      overwrite = '1',
                      lambdaInitPrior = "5.887243505",
                      lambdaShiftPrior = "0.006579862",
                      muInitPrior = "5.887243505",
                      expectedNumberOfShifts = '10',
                      lambdaIsTimeVariablePrior = "0"
                    ))

#---------- The bamm data object

edata <- getEventData(caudata_tree,
                      eventdata = "./run_bamm/event_data.txt",
                      burnin = 0.25)

#---------- Assessing MCMC convergence

# Plot the log-likelihood trace of the MCMC.
mcmcout <- read.csv("./run_bamm/mcmc_out.txt",
                    header = T)
plot(mcmcout$logLik ~ mcmcout$generation)

#' Discard some as burnin.
burnstart <- floor(0.25 * nrow(mcmcout)) 
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

save(postburn,
     file = "./run_bamm/postburn.R")

# Check the effective sample sizes of the log-likelihood and the number of shift events
# present in each sample.
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# In general, we want these to be at least 200 
# (and 200 is on the low side, but might be reasonable for very large datasets).

#---------- Assessing rate shifts

# Posterior probabilities of models sampled.
post_probs <- table(postburn$N_shifts) / nrow(postburn)

# See which models are part of the set that were sampled.
names(post_probs)
# In general, any model that is not included in names(post_probs) was so lousy 
# that it was never even sampled. Thus, if you fail to observe a model ‘0’ in this set, 
# this means that you have such overwhelming evidence for diversification rate heterogeneity 
# in your data that this model probability is effectively 0 (bear in mind that a model with name ‘0’ 
# is model M_0, or a model with no rate shifts). The probability of model ‘0’ is the posterior 
# probability of a model with just a single evolutionary rate dynamic (no rate shifts).
# IMPORTANT: We suggest that (usually) the overall best model from a BAMM analysis 
# is the model with the highest Bayes factor relative to the null model.

# Compute the posterior odds ratio for two models.
post_probs['1'] / post_probs['2']

# Summarize the posterior distribution of the number of shifts.
shift_probs <- summary(edata)
# Shift_probs is now a dataframe giving the posterior probabilities of 
# each rate shift count observed during simulation of the posterior.

#---------- Prior distribution

# Matrix of pairwise Bayes Factors comparing all models with a posterior/prior greater than zero.
bfmat <- computeBayesFactors(mcmcout,
                             expectedNumberOfShifts = 10,
                             burnin = 0.25)
# Bayes factors greater than 20 generally imply strong evidence for one model over another; 
# values greater than 50 are very strong evidence in favor of the numerator model.
# There is no definitive Bayes factor criterion for “significance”, but many researchers 
# consider values greater than 12 to be consistent with at least some effect.

# visualizing the prior and posterior simultaneously
plotPrior(mcmcout,
          expectedNumberOfShifts = 10)

#---------- Generate a mean phylorate plot

# Prepare color palette.
library(viridisLite)
v_colors =  viridis(1000, 
                    option = "H")

par(mar=c(0.1,0.1,0.1,0.1),
    lend = 2)

plot.bammdata(edata,
              # method = "polar",
              lwd = 2,
              pal = v_colors,
              legend = T)

# View a phylorate plot for any sample from the posterior.
# Here we will plot the 25th sample from the whale posterior:
index <- 25
e2 <- subsetEventData(edata, 
                      index = index)
plot.bammdata(e2, 
              lwd = 2)
addBAMMshifts(e2, 
              cex = 2)

#---------- Estimate the credible set of rate shifts

css <- credibleShiftSet(edata,
                        expectedNumberOfShifts = 10,
                        threshold = 5,
                        set.limit = 0.95)

# Number of distinct shift configurations in the data.
css$number.distinct
summary(css)

# Generate phylorate plots for each of the N shift configurations with the highest 
# posterior probabilities.
plot.credibleshiftset(css)
# The text above each phylorate plot gives the posterior probability of each shift configuration.

#---------- Best shift configuration

best <- getBestShiftConfiguration(edata,
                                  expectedNumberOfShifts = 10)

#---------- Marginal shift probabilities

# Compute the marginal shift probabilities for each branch, then plot a new phylogenetic tree 
# where the branch lengths are scaled by the probability that they contain a shift event.
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs)

#---------- Rate-through-time analysis

plotRateThroughTime(edata, 
                    ratetype = "extinction")

#---------- Macroevolutionary cohort analysis

# Macroevolutionary cohort analysis provides a way of summarizing the extent to which species share 
# correlated macroevolutionary dynamics. The basic idea is to visualize the pairwise probabilities 
# that any two species share a common macroevolutionary rate regime.
cmat <- getCohortMatrix(edata)
cohorts(cmat, edata)

#---------- Cumulative shift probabilities

# The cumulative shift probability tree shows the cumulative probability, on each branch, 
# that a shift occurred somewhere between the focal branch and the root of the tree.
# The occurrence of such a shift implies that evolutionary dynamics on the focal branch are 
# decoupled from the “background” diversification or trait evolutionary process at the root of 
# the tree.
cst <- cumulativeShiftProbsTree(edata)
plot.phylo(cst)

# extract the rates values 
rates_bamm <- getMarginalBranchRateMatrix(edata,
                                          verbose = TRUE)

# calculate mean speciation and extinction values
mean_speciation_rates_bamm <- apply(rates_bamm$lambda_branch_matrix, 1, mean)
b <- rowMeans(rates_bamm$lambda_branch_matrix)
mean_extinction_rates_bamm <- apply(rates_bamm$mu_branch_matrix, 1, mean)

mean_rates_bamm <- data.frame(
  mean_speciation = mean_speciation_rates_bamm,
  mean_extinction = mean_extinction_rates_bamm
)

# calculate net diversification
mean_rates_bamm$diversification <- mean_rates_bamm$mean_speciation - mean_rates_bamm$mean_extinction

# tip rates
tip_rates_bamm <- getTipRates(edata,
                              returnNetDiv = TRUE,
                              statistic = "mean")
bamm_tip_div_rates <- as.data.frame(tip_rates_bamm$netdiv.avg)
bamm_tip_div_rates$sp <- rownames(bamm_tip_div_rates)
colnames(bamm_tip_div_rates) <- c("net_div", "Species")

##### run clads ============================================================================================================================================

#---------- run the parameter inference in Julia

out = infer_ClaDS(caudata_tree, 10000, print_state = 1000)

#---------- extract rates 

# clads output
load("./run_clads/clads_run_caudata.rds")
clads_caudata <- CladsOutput

# extract rates at tips
lambda_sp <- clads_caudata$lambdatip_map
data_lambda_sp <- data.frame(
  sp = clads_caudata$tree$tip.label,
  lambda = lambda_sp
)

# calculate tip rates
epsilon <- clads_caudata$eps_map
tip_lambda <- clads_caudata$lambdatip_map
tip_div_rates = (1 - epsilon) * tip_lambda

tip_rates_clads_df <- data.frame(
  sp = clads_caudata$tree$tip.label,
  div_rates = tip_div_rates,
  lambda = tip_lambda,
  mu = c(tip_lambda - tip_div_rates))

# calculate branch rates
epsilon <- clads_caudata$eps_map
lambda_all_branches <- clads_caudata$lambdai_map

data_clads_stewart <- data.frame(
  n = 1:length(lambda_all_branches),
  lambda_all_branch = lambda_all_branches,
  mu_all_bhancges = lambda_all_branches * epsilon,
  div_rates = (1 - epsilon) * lambda_all_branches
)

##### test for correlation between bamm and clads rates =====================================================================

#---------- bamm and clads results

# bamm rates
mean_rates_bamm <- read.csv(file = "./test_corr_bamm_clads/mean_rates_bamm.csv",
                            header = TRUE,
                            sep = ",")

# branches rates from bamm and clads
all_models_rates <- read.csv2(file = "./test_corr_bamm_clads/all_models_rates.csv")

# tip rates from bamm and clads
all_models_tip_rates <- read.csv2(file = "./test_corr_bamm_clads/all_models_tip_rates.csv")

#---------- all branches rates

data <- all_models_rates[,2:7]

# initialise empty dataframe for results
results <- data.frame(
  Rate = character(),
  Test_Value = numeric(),
  R_squared = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# function to make correlation analyses
analyze_correlation <- function(column1, column2, rate_name) {
  correlation_test <- cor.test(column1, column2)
  r_squared <- correlation_test$estimate^2
  test_value <- correlation_test$statistic
  p_value <- correlation_test$p.value
  
  return(data.frame(
    Rate = rate_name,
    Test_Value = test_value,
    R_squared = r_squared,
    P_value = p_value
  ))
}

# correlation tests
results <- rbind(
  results,
  analyze_correlation(data$clads_spec, data$bamm_spec, "Speciation"),
  analyze_correlation(data$clads_ext, data$bamm_ext, "Extinction"),
  analyze_correlation(data$clads_div, data$bamm_div, "Diversification")
)

print(results)

#---------- tip rates

data <- all_models_tip_rates[,2:7]

# initialise empty dataframe for results
results <- data.frame(
  Rate = character(),
  Test_Value = numeric(),
  R_squared = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# correlation tests
results <- rbind(
  results,
  analyze_correlation(data$clads_spec, data$bamm_spec, "Speciation"),
  analyze_correlation(data$clads_ext, data$bamm_ext, "Extinction"),
  analyze_correlation(data$clads_net_div, data$bamm_net_div, "Diversification")
)

print(results)
##### pgls using bamm and clads rates - genus level =======================================================================

# caudata phylogeny - genus level
caudata_genus_tree <- ape::read.nexus("./pgls_genus_level/caudata_genus_tree.nexus")

#---------- bamm and clads results

# bamm and clads rates per genus
caudata_genus_data <- read.csv(file = "./pgls_genus_level/caudata_genus_data.csv")
rownames(caudata_genus_data) <- caudata_genus_data$genus

# dataset without groups with extreme values
caudata_genus_tree_ex <- drop.tip(caudata_genus_tree, c("Bolitoglossa", "Isthmura"))
caudata_genus_tree_ex <- drop.tip(caudata_genus_tree, "Bolitoglossa")
plot(caudata_genus_tree_ex)

caudata_genus_data_ex <- caudata_genus_data[!(caudata_genus_data$genus %in% c("Bolitoglossa", "Isthmura")), ]
caudata_genus_data_ex <- caudata_genus_data[!(caudata_genus_data$genus %in% "Bolitoglossa"), ]

# dataset without groups with extreme values and species richness < 2
genus_inf2sp <- caudata_genus_data_ex$genus[caudata_genus_data_ex$species_richness <= 2]
caudata_genus_tree_ex_sup2 <- drop.tip(caudata_genus_tree_ex, genus_inf2sp)
plot(caudata_genus_tree_ex_sup2)

caudata_genus_data_ex_sup2 <- caudata_genus_data_ex[caudata_genus_data_ex$species_richness > 2, ]

# prepare results table
pgls_phylolm_genus_rslts <- list()

library(phylolm)
phy <- caudata_genus_tree
data <- caudata_genus_data

# Initialize empty results for this iteration
pgls_rslts <- data.frame(formula = character(),
                         aicc = numeric(),
                         estimate = numeric(),
                         adj_r_squared = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE) # Ensure characters are not converted to factors

#---------- Models to test
models_to_test <- list(
  pgls_crown = "log10(species_richness) ~ crown_age",
  pgls_stem = "log10(species_richness) ~ stem_age",
  pgls_clads_div = "log10(species_richness) ~ log10(clads_div)",
  pgls_clads_spec = "log10(species_richness) ~ log10(clads_spec)",
  pgls_clads_ext = "log10(species_richness) ~ log10(clads_ext)",
  pgls_bamm_div = "log10(species_richness) ~ log10(bamm_div)",
  pgls_bamm_spec = "log10(species_richness) ~ log10(bamm_spec)",
  pgls_bamm_ext = "log10(species_richness) ~ log10(bamm_ext)",
  pgls_mean_div = "log10(species_richness) ~ log10(mean_div)",
  pgls_mean_spec = "log10(species_richness) ~ log10(mean_spec)",
  pgls_mean_ext = "log10(species_richness) ~ log10(mean_ext)"
)

# Loop over models
for (model_name in names(models_to_test)) {
  model_formula <- models_to_test[[model_name]]
  
  pgls_model <- phylolm(as.formula(model_formula), data = data, phy = phy)
  sum_model <- summary(pgls_model)
  
  # Extract formula and AICc
  aicc <- sum_model$aic
  
  # Extract F, adjusted R², and p-value
  f <- sum_model$coefficients[2, 1]
  r <- sum_model$adj.r.squared
  p <- sum_model$coefficients[2, 4]
  
  # Add results to the dataframe
  pgls_rslts <- rbind(pgls_rslts,
                      data.frame(
                        formula = model_formula,
                        aicc = aicc,
                        estimate = f,
                        adj_r_squared = r,
                        p_value = p,
                        stringsAsFactors = FALSE)) # Ensure characters are not converted to factors
  
  
}

# Store results for this iteration
pgls_phylolm_genus_rslts <- pgls_rslts
pgls_phylolm_genus_rslts

##### pgls using bamm and clads rates - family levels ========================================================================

# caudata phylogeny - family level
caudata_family_tree <- ape::read.nexus("./pgls_family_level/caudata_family_tree.nexus")

#---------- bamm and clads results

# bamm and clads rates per family
caudata_family_data <- read.csv(file = "./caudata_family_data.csv")
rownames(caudata_family_data) <- caudata_family_data$family

# datasets without groups with extreme values
caudata_family_tree_ex <- drop.tip(caudata_family_tree, "Plethodontidae")
plot(caudata_family_tree_ex)

caudata_family_data_ex <- caudata_family_data[!(caudata_family_data$family %in% "Plethodontidae"), ]

# prepare results table
pgls_phylolm_family_rslts <- list()

library(phylolm)
phy <- caudata_family_tree
data <- caudata_family_data

# Initialize empty results for this iteration
pgls_rslts <- data.frame(formula = character(),
                         aicc = numeric(),
                         estimate = numeric(),
                         adj_r_squared = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE) # Ensure characters are not converted to factors

#---------- Models to test
models_to_test <- list(
  pgls_crown = "log10(species_richness) ~ crown_age",
  pgls_stem = "log10(species_richness) ~ stem_age",
  pgls_clads_div = "log10(species_richness) ~ log10(clads_div)",
  pgls_clads_spec = "log10(species_richness) ~ log10(clads_spec)",
  pgls_clads_ext = "log10(species_richness) ~ log10(clads_ext)",
  pgls_bamm_div = "log10(species_richness) ~ log10(bamm_div)",
  pgls_bamm_spec = "log10(species_richness) ~ log10(bamm_spec)",
  pgls_bamm_ext = "log10(species_richness) ~ log10(bamm_ext)",
  pgls_mean_div = "log10(species_richness) ~ log10(mean_div)",
  pgls_mean_spec = "log10(species_richness) ~ log10(mean_spec)",
  pgls_mean_ext = "log10(species_richness) ~ log10(mean_ext)"
)

# Loop over models
for (model_name in names(models_to_test)) {
  model_formula <- models_to_test[[model_name]]
  
  pgls_model <- phylolm(as.formula(model_formula), data = data, phy = phy)
  sum_model <- summary(pgls_model)
  
  # Extract formula and AICc
  aicc <- sum_model$aic
  
  # Extract F, adjusted R², and p-value
  f <- sum_model$coefficients[2, 1]
  r <- sum_model$adj.r.squared
  p <- sum_model$coefficients[2, 4]
  
  # Add results to the dataframe
  pgls_rslts <- rbind(pgls_rslts,
                      data.frame(
                        formula = model_formula,
                        aicc = aicc,
                        estimate = f,
                        adj_r_squared = r,
                        p_value = p,
                        stringsAsFactors = FALSE)) # Ensure characters are not converted to factors
  
  
}

# Store results for this iteration
pgls_phylolm_family_rslts <- pgls_rslts
pgls_phylolm_family_rslts

##### bamm and clads randomisation procedure ================================================================================

#---------- random selection of nodes

num_repeats <- 1000
num_rand_nodes <- 10

# create a dataframe to store the results
results <- data.frame(
  node = integer(),
  species_richness = integer(),
  stem_Age = numeric(),
  crown_Age = numeric()
)

library(phangorn)
# function to select non-nested nodes
select_unique_nodes <- function(tree) {
  # identify intern nodes without the root node
  internal_nodes <- setdiff((Ntip(tree) + 2):max(tree$edge[, 1]), max(tree$edge[, 1]))
  
  repeat {
    selected_nodes <- integer(0)
    internal_nodes_temp <- internal_nodes
    
    # random selection of n nodes
    while (length(selected_nodes) < num_rand_nodes && length(internal_nodes_temp) > 0) {
      # random selection of a node
      node <- sample(internal_nodes_temp, 1)
      
      # identify the nodes included in this node
      clade_descendants <- allDescendants(tree)[[node]]
      
      # check if there is a conflict with other selected nodes
      if (all(!selected_nodes %in% clade_descendants)) {
        # add the node to the selection
        selected_nodes <- c(selected_nodes, node)
        
        # remove all nested nodes
        internal_nodes_temp <- setdiff(internal_nodes_temp, c(clade_descendants, node))
      } else {
        # if the node include already selected node, sample another one
        internal_nodes_temp <- setdiff(internal_nodes_temp, node)
      }
    }
    
    # if n nodes have been successfully selected, get out of the loop
    if (length(selected_nodes) == num_rand_nodes) {
      break
    }
  }
  
  return(selected_nodes)
}

# loop to repeat the operation n time
for (i in 1:num_repeats) {
  # call the function to select nodes
  selected_nodes <- select_unique_nodes(tree)
  
  # plot the selected nodes on the phylogeny
  plot(tree, show.tip.label = FALSE, show.node.label = TRUE)
  nodelabels(node = selected_nodes, pch = 19, col = "red", cex = 1.5)
  
  # extract data for each node
  for (node in selected_nodes) {
    # number of tips in the clade
    descendants <- Descendants(tree, node, "tips")[[1]]
    num_tips <- length(descendants)
    
    # Crown Age 
    crown_age <- node.depth.edgelength(tree)[node]
    
    # Stem Age 
    parent_edge <- tree$edge.length[which(tree$edge[, 2] == node)]
    stem_age <- crown_age + parent_edge
    
    # store data into the data frame
    results <- rbind(results, data.frame(
      node = node,
      species_richness = num_tips,
      stem_Age = stem_age,
      crown_Age = crown_age
    ))
  }
}

head(results)
write.csv2(results,
           file = "./selected_nodes.csv",
           row.names = FALSE)

#---------- extract rate data

# extract data from ClaDS results
clads_caudata <- CladsOutput

epsilon <- clads_caudata$eps_map
lambda_all_branches <- clads_caudata$lambdai_map

rates_clads_df <- data.frame(
  lambda_all_branch = lambda_all_branches,
  mu_all_bhancges = lambda_all_branches * epsilon,
  div_rates = (1 - epsilon) * lambda_all_branches
)

# assemble rate data
rates <- data.frame(
  branch = tree$edge[,2],
  clads_spec = rates_clads_df$lambda_all_branch,
  clads_ext = rates_clads_df$mu_all_bhancges,
  clads_div = rates_clads_df$div_rates,
  bamm_spec = mean_rates_bamm$mean_speciation,
  bamm_ext = mean_rates_bamm$mean_extinction,
  bamm_div = mean_rates_bamm$diversification
)

# function to get subtrees for each randomly selected node
library(geiger)
get_subtree_branches <- function(phylo_tree, node_number) {
  # extract all tips from a node
  tips_in_subtree <- tips(phylo_tree, node_number)
  
  # drop all tips that are not in the subtree
  subtree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, tips_in_subtree))
  
  # Return the branches of the subtree
  return(subtree$edge[, 2])
}

clads_mean_spec_rates <- c()
clads_mean_ext_rates <- c()
clads_mean_div_rates <- c()
bamm_mean_spec_rates <- c()
bamm_mean_ext_rates <- c()
bamm_mean_div_rates <- c()

nodes <- results$node

for (node in nodes) {
  # extract branched of the subtrees from node index
  branches <- get_subtree_branches(tree, node)
  
  # filter table rates to get the rates for the specific branches
  sub_tree_clads_speciation_rates <- rates$clads_spec[rates$branch %in% branches]
  sub_tree_clads_extinction_rates <- rates$clads_ext[rates$branch %in% branches]
  sub_tree_clads_diversification_rates <- rates$clads_div[rates$branch %in% branches]
  sub_tree_bamm_speciation_rates <- rates$bamm_spe[rates$branch %in% branches]
  sub_tree_bamm_extinction_rates <- rates$bamm_ext[rates$branch %in% branches]
  sub_tree_bamm_diversification_rates <- rates$bamm_div[rates$branch %in% branches]
  
  # calculate mean rate for this subtree
  mean_clads_spec_rate <- mean(sub_tree_clads_speciation_rates)
  mean_clads_ext_rate <- mean(sub_tree_clads_extinction_rates)
  mean_clads_div_rate <- mean(sub_tree_clads_diversification_rates)
  mean_bamm_spec_rate <- mean(sub_tree_bamm_speciation_rates)
  mean_bamm_ext_rate <- mean(sub_tree_bamm_extinction_rates)
  mean_bamm_div_rate <- mean(sub_tree_bamm_diversification_rates)
  
  # add result to the list
  clads_mean_spec_rates <- c(clads_mean_spec_rates, mean_clads_spec_rate)
  clads_mean_ext_rates <- c(clads_mean_ext_rates, mean_clads_ext_rate)
  clads_mean_div_rates <- c(clads_mean_div_rates, mean_clads_div_rate)
  bamm_mean_spec_rates <- c(bamm_mean_spec_rates, mean_bamm_spec_rate)
  bamm_mean_ext_rates <- c(bamm_mean_ext_rates, mean_bamm_ext_rate)
  bamm_mean_div_rates <- c(bamm_mean_div_rates, mean_bamm_div_rate)
}

#---------- assemble and arrange data for pgls analyses

data_per_clade <- data.frame(
  node_index = results$node, 
  species_richness = results$species_richness,
  crown_age = results$crown_Age,
  stem_age = results$stem_Age,
  clads_spe = clads_mean_spec_rates,
  clads_ext = clads_mean_ext_rates,
  clads_div = clads_mean_div_rates,
  bamm_spe = bamm_mean_spec_rates,
  bamm_ext = bamm_mean_ext_rates,
  bamm_div = bamm_mean_div_rates
)

data_per_clade$mean_spe <- rowMeans(data_per_clade[ , c(5, 8)])
data_per_clade$mean_ext <- rowMeans(data_per_clade[ , c(6, 9)])
data_per_clade$mean_div <- rowMeans(data_per_clade[ , c(7, 10)])

#---------- create subtrees for each dataset

# split dataset into the number of repetition and create subtrees 
total_rows <- nrow(data_per_clade)

library(castor)
for (i in 1:num_repeats) {
  # Déterminer l'index de début et de fin pour la tranche
  start_index <- (i - 1) * num_rand_nodes + 1
  end_index <- min(start_index + 9, total_rows)
  
  # Extraire la tranche
  rep <- data_per_clade[start_index:end_index, ]
  
  # get the picked nodes
  nodes = rep$node_index - Ntip(tree)
  
  # get the tips that included in the picked nodes
  subtrees_tip = get_subtrees_at_nodes(tree, nodes)$new2old_tips
  
  # keep the basal tip in each picked clade
  pickedtips <- data.frame()
  
  for (j in 1:num_rand_nodes) {
    pickedtips <- rbind(pickedtips, as.data.frame(subtrees_tip[[j]][1]))
  }
  
  tree2 <- keep.tip(tree, pickedtips$`subtrees_tip[[j]][1]`)
  
  # save the tree
  write.nexus(tree2, 
              file = paste0("./random_data/subtree_", i, ".nexus"))
  
  # identifier le nombre de tips dans l'arbre
  num_tips <- Ntip(tree2)
  
  # changer le nom des tips
  nametips <- data.frame()
  for (k in 1:num_rand_nodes) {
    nametips <- rbind(nametips, as.data.frame(tips(tree, subtrees_tip[[k]][1])))
  }
  
  colnames(nametips) <- "representative_species"
  
  # ajouter le nom du clade dans le data frame
  data <- cbind(nametips, rep)
  rownames(data) <- data$representative_species # necessary if we do pgls using phylolm
  
  # save data
  write.csv(data,
            file = paste0("./random_data/random_data_", i, ".csv"))
}

##### pgls using bamm and clads rates - randomisation procedure =================================================================

#---------- load data

path <- "./random_data/"

# data nodes
tables <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
tables_random_data <- tables_random_data <- lapply(tables, function(file) {
  read.csv(file, row.names = 1)  # Assumes the first column contains the row names
})

# subtrees
library(ape)
sub_trees <- list.files(path = path, pattern = "\\.nexus$", full.names = TRUE)
subtrees <- lapply(sub_trees, read.nexus)

#---------- realise pgls

results_pgls_random_clades <- list()

library(phylolm)
for(i in 1:length(subtrees)){
  
  phy <- subtrees[[i]]
  data <- tables_random_data[[i]]
  
  # Initialize empty results for this iteration
  pgls_rand_clades_rslts <- data.frame(formula = character(),
                                       aicc = numeric(),
                                       estimate = numeric(),
                                       adj_r_squared = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE) # Ensure characters are not converted to factors
  
  # Models to test
  models_to_test <- list(
    pgls_crown = "log10(species_richness) ~ crown_age",
    pgls_stem = "log10(species_richness) ~ stem_age",
    pgls_clads_div = "log10(species_richness) ~ log10(clads_div)",
    pgls_clads_spec = "log10(species_richness) ~ log10(clads_spe)",
    pgls_clads_ext = "log10(species_richness) ~ log10(clads_ext)",
    pgls_bamm_div = "log10(species_richness) ~ log10(bamm_div)",
    pgls_bamm_spec = "log10(species_richness) ~ log10(bamm_spe)",
    pgls_bamm_ext = "log10(species_richness) ~ log10(bamm_ext)",
    pgls_mean_div = "log10(species_richness) ~ log10(mean_div)",
    pgls_mean_spec = "log10(species_richness) ~ log10(mean_spe)",
    pgls_mean_ext = "log10(species_richness) ~ log10(mean_ext)"
  )
  
  # Loop over models
  for (model_name in names(models_to_test)) {
    model_formula <- models_to_test[[model_name]]
    
    pgls_model <- phylolm(as.formula(model_formula), data = data, phy = phy)
    sum_model <- summary(pgls_model)
    
    # Extract formula and AICc
    aicc <- sum_model$aic
    
    # Extract F, adjusted R², and p-value
    f <- sum_model$coefficients[2, 1]
    r <- sum_model$adj.r.squared
    p <- sum_model$coefficients[2, 4]
    
    # Add results to the dataframe
    pgls_rand_clades_rslts <- rbind(pgls_rand_clades_rslts, 
                                    data.frame(formula = model_formula,
                                               aicc = aicc,
                                               estimate = f,
                                               adj_r_squared = r,
                                               p_value = p,
                                               stringsAsFactors = FALSE)) # Ensure characters are not converted to factors
    
      
  }
  
  # Store results for this iteration
  results_pgls_random_clades[[i]] <- pgls_rand_clades_rslts
}

# check if there are NA in results
na_indices <- c()
# Loop through each element in the results list
for (i in seq_along(results_pgls_random_clades)) {
  # Check if the current data frame contains any NA values
  if (any(is.na(results_pgls_random_clades[[i]]))) {
    # If it does, add the index to the na_indices vector
    na_indices <- c(na_indices, i)
  }
}

#---------- assemble results

all_results <- list()

for (table in results_pgls_random_clades) {
  # extract data
  formulas <- table[[1]]
  aics <- table[[2]]
  estimates <- table[[3]]
  r_squared <- table[[4]]
  p_values <- table[[5]]
  
  # collect result for each formula
  unique_formulas <- unique(formulas)  # get all unique formulas
  for (formula in unique_formulas) {
    # filter lines corresponding to the current formula
    indices <- which(formulas == formula)
    
    # calculate means
    mean_aic <- mean(aics[indices], na.rm = TRUE)
    mean_r_squared <- mean(r_squared[indices], na.rm = TRUE)
    mean_p_value <- mean(p_values[indices], na.rm = TRUE)
    
    # count the number of pvalue < 0.05
    count_p_less_than_0_05 <- sum(p_values[indices] < 0.05, na.rm = TRUE)
    total_count <- length(indices)  # total values for this formula
    
    # calculate the percentage of pvalue < 0.05
    percent_p_value_less_than_0_05 <- ifelse(total_count > 0, (count_p_less_than_0_05 / total_count) * 100, NA)
    
    # add results to the table
    all_results[[length(all_results) + 1]] <- data.frame(
      formula = formula,
      mean_aic = mean_aic,
      mean_r_squared = mean_r_squared,
      mean_p_value = mean_p_value,
      percent_p_value_less_than_0_05 = percent_p_value_less_than_0_05,
      stringsAsFactors = FALSE
    )
  }
}

# combine all results
results_summary <- do.call(rbind, all_results)

# calculate the final means and percentages
final_results <- aggregate(. ~ formula, data = results_summary, FUN = mean, na.rm = TRUE)

# order table
order <- c(
  "log10(species_richness) ~ crown_age",
  "log10(species_richness) ~ stem_age",
  "log10(species_richness) ~ log10(clads_div)",
  "log10(species_richness) ~ log10(clads_spe)",
  "log10(species_richness) ~ log10(clads_ext)",
  "log10(species_richness) ~ log10(bamm_div)",
  "log10(species_richness) ~ log10(bamm_spe)",
  "log10(species_richness) ~ log10(bamm_ext)",
  "log10(species_richness) ~ log10(mean_div)",
  "log10(species_richness) ~ log10(mean_spe)",
  "log10(species_richness) ~ log10(mean_ext)"
)

final_results[[1]] <- factor(final_results[[1]], levels = order, ordered = TRUE)
final_results <- final_results[order(final_results[[1]]), ]
final_results

##### run environmental birth-death models - all caudata ==========================================================================

#---------- phylogeny
tree <- caudata_tree
caudata_tot_time <- max(TreeSim::getx(caudata_tree))
tot_time <- caudata_tot_time
f = 1

library(RPANDA)

#---------- paleo-environmental data

# temperature condamine 2022
data_temp_cond_22 <- read.csv2(file = "./data/environmental_data/stewart_wiens_2024/temperatures_condamine2022_arranged.csv")
dof_temp_cond_22 <- smooth.spline(data_temp_cond_22[, "time"], data_temp_cond_22[, "temperature"])$df

# fragmentation zaffos
data_frag <- read.csv(file = "./data/environmental_data/stewart_wiens_2024/fragmentation_zaffos_arranged.csv")
dof_frag <- smooth.spline(data_frag[, "time"], data_frag[, "fragmentation"])$df

# detrended temperature condamine 2022
temp_cond_22_detrended <- read.csv2(file = "./data/environmental_data/stewart_wiens_2024/temperature_condamine2022_detrended.csv")
dof_temp_cond_22_det <- smooth.spline(temp_cond_22_detrended[, "time"], temp_cond_22_detrended[, "detrended_temperature"])$df

# detrended zaffos
frag_detrended <- read.csv(file = "./data/environmental_data/stewart_wiens_2024/fragmentation_zaffos_detrended.csv")
dof_frag_det <- smooth.spline(frag_detrended[, "time"], frag_detrended[, "detrended_fragmentation"])$df

#---------- yule

f.lamb <- function(t, y) {y[1]}
f.mu <- function(t, y) {0}
lamb_par <- c(0.1)
mu_par <- c()

model <- RPANDA::fit_bd(
  phylo = tree,
  tot_time = tot_time,
  f.lamb = f.lamb,
  f.mu = f.mu,
  lamb_par = lamb_par,
  mu_par = mu_par,
  f = f,
  cst.lamb = TRUE,
  cst.mu = FALSE,
  expo.lamb = FALSE,
  expo.mu = FALSE,
  fix.mu = TRUE,
  cond = "crown",
  dt = 1e-3
)

save(model,
     file = paste0("./caudata_Bcst_noD.Rdata"))

#---------- Benv noD

f.lamb <- function(t, x, y) {y[1] * exp(y[2] * x)}
f.mu <- function(t, x, y) {0}
lamb_par_init <- c(0.1, 0.01)
mu_par_init <- c()

library(RPANDA)

#---------- Benv Dcst 

# exponential model with speciation rate depending on global average temperature through time and constant 
# extinction rate
f.lamb <- function(t, x, y) {y[1] * exp(y[2] * x)}
f.mu <- function(t, x, y) {y[1]}
lamb_par_init <- c(0.1, 0.01)
mu_par_init <- c(0.02)

Btempcond22_Dcst <- fit_env(phylo = tree,
                            env_data = data_temp_cond_22,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.mu = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = dof_temp_cond_22)

save(Btempcond22_Dcst,
     file = "./caudata_Btempcond22_Dcst.Rdata")

Btempcond22_Dcst <- fit_env(phylo = tree,
                            env_data = data_temp_cond_22,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.mu = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = 3)

save(Btempcond22_Dcst,
     file = "./caudata_Btempcond22_Dcst_df3.Rdata")

Btempcond22_Dcst_df025 <- fit_env(phylo = tree,
                                  env_data = data_temp_cond_22,
                                  tot_time = tot_time,
                                  f.lamb = f.lamb,
                                  f.mu = f.mu,
                                  lamb_par = lamb_par_init,
                                  mu_par = mu_par_init,
                                  f = f,
                                  cst.mu = TRUE,
                                  cond = "crown",
                                  dt = 1e-3,
                                  df = dof_temp_cond_22*0.25)

save(Btempcond22_Dcst_df025,
     file = "./caudata_Btempcond22_Dcst_df025.Rdata")

Btempcond22_Dcst_df05 <- fit_env(phylo = tree,
                                 env_data = data_temp_cond_22,
                                 tot_time = tot_time,
                                 f.lamb = f.lamb,
                                 f.mu = f.mu,
                                 lamb_par = lamb_par_init,
                                 mu_par = mu_par_init,
                                 f = f,
                                 cst.mu = TRUE,
                                 cond = "crown",
                                 dt = 1e-3,
                                 df = dof_temp_cond_22*0.5)

save(Btempcond22_Dcst_df05,
     file = "./caudata_Btempcond22_Dcst_df05.Rdata")

Btempcond22_Dcst_df075 <- fit_env(phylo = tree,
                                  env_data = data_temp_cond_22,
                                  tot_time = tot_time,
                                  f.lamb = f.lamb,
                                  f.mu = f.mu,
                                  lamb_par = lamb_par_init,
                                  mu_par = mu_par_init,
                                  f = f,
                                  cst.mu = TRUE,
                                  cond = "crown",
                                  dt = 1e-3,
                                  df = dof_temp_cond_22*0.75)

save(Btempcond22_Dcst_df075,
     file = "./caudata_Btempcond22_Dcst_df075.Rdata")

Btempcond22det_Dcst <- fit_env(phylo = tree,
                               env_data = temp_cond_22_detrended,
                               tot_time = tot_time,
                               f.lamb = f.lamb,
                               f.mu = f.mu,
                               lamb_par = lamb_par_init,
                               mu_par = mu_par_init,
                               f = f,
                               cst.mu = TRUE,
                               cond = "crown",
                               dt = 1e-3,
                               df = dof_temp_cond_22_det)
save(Btempcond22det_Dcst,
     file = "./caudata_Btempcond22det_Dcst.Rdata")

Btempcond22det_Dcst_df025 <- fit_env(phylo = tree,
                                     env_data = temp_cond_22_detrended,
                                     tot_time = tot_time,
                                     f.lamb = f.lamb,
                                     f.mu = f.mu,
                                     lamb_par = lamb_par_init,
                                     mu_par = mu_par_init,
                                     f = f,
                                     cst.mu = TRUE,
                                     cond = "crown",
                                     dt = 1e-3,
                                     df = dof_temp_cond_22_det*0.25)
save(Btempcond22det_Dcst_df025,
     file = "./caudata_Btempcond22det_Dcst_df025.Rdata")

Btempcond22det_Dcst_df05 <- fit_env(phylo = tree,
                                    env_data = temp_cond_22_detrended,
                                    tot_time = tot_time,
                                    f.lamb = f.lamb,
                                    f.mu = f.mu,
                                    lamb_par = lamb_par_init,
                                    mu_par = mu_par_init,
                                    f = f,
                                    cst.mu = TRUE,
                                    cond = "crown",
                                    dt = 1e-3,
                                    df = dof_temp_cond_22_det*0.5)
save(Btempcond22det_Dcst_df05,
     file = "./caudata_Btempcond22det_Dcst_df05.Rdata")

Btempcond22det_Dcst_df075 <- fit_env(phylo = tree,
                                     env_data = temp_cond_22_detrended,
                                     tot_time = tot_time,
                                     f.lamb = f.lamb,
                                     f.mu = f.mu,
                                     lamb_par = lamb_par_init,
                                     mu_par = mu_par_init,
                                     f = f,
                                     cst.mu = TRUE,
                                     cond = "crown",
                                     dt = 1e-3,
                                     df = dof_temp_cond_22_det*0.75)
save(Btempcond22det_Dcst_df075,
     file = "./caudata_Btempcond22det_Dcst_df075.Rdata")

Bfrag_Dcst <- fit_env(phylo = tree,
                      env_data = data_frag,
                      tot_time = tot_time,
                      f.lamb = f.lamb,
                      f.mu = f.mu,
                      lamb_par = lamb_par_init,
                      mu_par = mu_par_init,
                      f = f,
                      cst.mu = TRUE,
                      cond = "crown",
                      dt = 1e-3,
                      df = dof_frag)

save(Bfrag_Dcst,
     file = "./caudata_Bfrag_Dcst.Rdata")

Bfrag_Dcst <- fit_env(phylo = tree,
                      env_data = data_frag,
                      tot_time = tot_time,
                      f.lamb = f.lamb,
                      f.mu = f.mu,
                      lamb_par = lamb_par_init,
                      mu_par = mu_par_init,
                      f = f,
                      cst.mu = TRUE,
                      cond = "crown",
                      dt = 1e-3,
                      df = 3)

save(Bfrag_Dcst,
     file = "./caudata_Bfrag_Dcst_df3.Rdata")

Bfrag_Dcst_df025 <- fit_env(phylo = tree,
                            env_data = data_frag,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.mu = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = dof_frag*0.25)

save(Bfrag_Dcst_df025,
     file = "./caudata_Bfrag_Dcst_df025.Rdata")

Bfrag_Dcst_df05 <- fit_env(phylo = tree,
                           env_data = data_frag,
                           tot_time = tot_time,
                           f.lamb = f.lamb,
                           f.mu = f.mu,
                           lamb_par = lamb_par_init,
                           mu_par = mu_par_init,
                           f = f,
                           cst.mu = TRUE,
                           cond = "crown",
                           dt = 1e-3,
                           df = dof_frag*0.5)

save(Bfrag_Dcst_df05,
     file = "./caudata_Bfrag_Dcst_df05.Rdata")

Bfrag_Dcst_df075 <- fit_env(phylo = tree,
                            env_data = data_frag,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.mu = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = dof_frag*0.75)

save(Bfrag_Dcst_df075,
     file = "./caudata_Bfrag_Dcst_df075.Rdata")

Bfragdet_Dcst <- fit_env(phylo = tree,
                         env_data = frag_detrended,
                         tot_time = tot_time,
                         f.lamb = f.lamb,
                         f.mu = f.mu,
                         lamb_par = lamb_par_init,
                         mu_par = mu_par_init,
                         f = f,
                         cst.mu = TRUE,
                         cond = "crown",
                         dt = 1e-3,
                         df = dof_frag_det)
save(Bfragdet_Dcst,
     file = "./caudata_Bfragdet_Dcst.Rdata")

Bfragdet_Dcst_df025 <- fit_env(phylo = tree,
                               env_data = frag_detrended,
                               tot_time = tot_time,
                               f.lamb = f.lamb,
                               f.mu = f.mu,
                               lamb_par = lamb_par_init,
                               mu_par = mu_par_init,
                               f = f,
                               cst.mu = TRUE,
                               cond = "crown",
                               dt = 1e-3,
                               df = dof_frag_det*0.25)
save(Bfragdet_Dcst_df025,
     file = "./caudata_Bfragdet_Dcst_df025.Rdata")

Bfragdet_Dcst_df05 <- fit_env(phylo = tree,
                              env_data = frag_detrended,
                              tot_time = tot_time,
                              f.lamb = f.lamb,
                              f.mu = f.mu,
                              lamb_par = lamb_par_init,
                              mu_par = mu_par_init,
                              f = f,
                              cst.mu = TRUE,
                              cond = "crown",
                              dt = 1e-3,
                              df = dof_frag_det*0.5)
save(Bfragdet_Dcst_df05,
     file = "./caudata_Bfragdet_Dcst_df05.Rdata")

Bfragdet_Dcst_df075 <- fit_env(phylo = tree,
                               env_data = frag_detrended,
                               tot_time = tot_time,
                               f.lamb = f.lamb,
                               f.mu = f.mu,
                               lamb_par = lamb_par_init,
                               mu_par = mu_par_init,
                               f = f,
                               cst.mu = TRUE,
                               cond = "crown",
                               dt = 1e-3,
                               df = dof_frag_det*0.75)
save(Bfragdet_Dcst_df075,
     file = "./caudata_Bfragdet_Dcst_df075.Rdata")

#---------- Bcst Denv 

# exponential model with constant birth, death varying with environamental variable
f.lamb <- function(t, x, y) {y[1] }
f.mu <- function(t, x, y) {y[1] * exp(y[2] * x)}
lamb_par_init <- c(0.1)
mu_par_init <- c(0.02, 0.01)

Bcst_Dtempcond22 <- fit_env(phylo = tree,
                            env_data = data_temp_cond_22,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.lamb = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = dof_temp_cond_22)

save(Bcst_Dtempcond22,
     file = "./caudata_Bcst_Dtempcond22.Rdata")

Bcst_Dtempcond22 <- fit_env(phylo = tree,
                            env_data = data_temp_cond_22,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.lamb = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = 3)

save(Bcst_Dtempcond22,
     file = "./caudata_Bcst_Dtempcond22_df3.Rdata")

Bcst_Dtempcond22_df025 <- fit_env(phylo = tree,
                                  env_data = data_temp_cond_22,
                                  tot_time = tot_time,
                                  f.lamb = f.lamb,
                                  f.mu = f.mu,
                                  lamb_par = lamb_par_init,
                                  mu_par = mu_par_init,
                                  f = f,
                                  cst.lamb = TRUE,
                                  cond = "crown",
                                  dt = 1e-3,
                                  df = dof_temp_cond_22*0.25)

save(Bcst_Dtempcond22_df025,
     file = "./caudata_Bcst_Dtempcond22_df025.Rdata")

Bcst_Dtempcond22_df05 <- fit_env(phylo = tree,
                                 env_data = data_temp_cond_22,
                                 tot_time = tot_time,
                                 f.lamb = f.lamb,
                                 f.mu = f.mu,
                                 lamb_par = lamb_par_init,
                                 mu_par = mu_par_init,
                                 f = f,
                                 cst.lamb = TRUE,
                                 cond = "crown",
                                 dt = 1e-3,
                                 df = dof_temp_cond_22*0.5)

save(Bcst_Dtempcond22_df05,
     file = "./caudata_Bcst_Dtempcond22_df05.Rdata")

Bcst_Dtempcond22_df075 <- fit_env(phylo = tree,
                                  env_data = data_temp_cond_22,
                                  tot_time = tot_time,
                                  f.lamb = f.lamb,
                                  f.mu = f.mu,
                                  lamb_par = lamb_par_init,
                                  mu_par = mu_par_init,
                                  f = f,
                                  cst.lamb = TRUE,
                                  cond = "crown",
                                  dt = 1e-3,
                                  df = dof_temp_cond_22*0.75)

save(Bcst_Dtempcond22_df075,
     file = "./caudata_Bcst_Dtempcond22_df075.Rdata")

Bcst_Dtempcond22det <- fit_env(phylo = tree,
                               env_data = temp_cond_22_detrended,
                               tot_time = tot_time,
                               f.lamb = f.lamb,
                               f.mu = f.mu,
                               lamb_par = lamb_par_init,
                               mu_par = mu_par_init,
                               f = f,
                               cst.lamb = TRUE,
                               cond = "crown",
                               dt = 1e-3,
                               df = dof_temp_cond_22_det)
save(Bcst_Dtempcond22det,
     file = "./caudata_Bcst_Dtempcond22det.Rdata")

Bcst_Dtempcond22det_df025 <- fit_env(phylo = tree,
                                     env_data = temp_cond_22_detrended,
                                     tot_time = tot_time,
                                     f.lamb = f.lamb,
                                     f.mu = f.mu,
                                     lamb_par = lamb_par_init,
                                     mu_par = mu_par_init,
                                     f = f,
                                     cst.lamb = TRUE,
                                     cond = "crown",
                                     dt = 1e-3,
                                     df = dof_temp_cond_22_det*0.25)
save(Bcst_Dtempcond22det_df025,
     file = "./caudata_Bcst_Dtempcond22det_df025.Rdata")

Bcst_Dtempcond22det_df05 <- fit_env(phylo = tree,
                                    env_data = temp_cond_22_detrended,
                                    tot_time = tot_time,
                                    f.lamb = f.lamb,
                                    f.mu = f.mu,
                                    lamb_par = lamb_par_init,
                                    mu_par = mu_par_init,
                                    f = f,
                                    cst.lamb = TRUE,
                                    cond = "crown",
                                    dt = 1e-3,
                                    df = dof_temp_cond_22_det*0.5)
save(Bcst_Dtempcond22det_df05,
     file = "./caudata_Bcst_Dtempcond22det_df05.Rdata")

Bcst_Dtempcond22det_df075 <- fit_env(phylo = tree,
                                     env_data = temp_cond_22_detrended,
                                     tot_time = tot_time,
                                     f.lamb = f.lamb,
                                     f.mu = f.mu,
                                     lamb_par = lamb_par_init,
                                     mu_par = mu_par_init,
                                     f = f,
                                     cst.lamb = TRUE,
                                     cond = "crown",
                                     dt = 1e-3,
                                     df = dof_temp_cond_22_det*0.75)
save(Bcst_Dtempcond22det_df075,
     file = "./caudata_Bcst_Dtempcond22det_df075.Rdata")

Bcst_Dfrag <- fit_env(phylo = tree,
                      env_data = data_frag,
                      tot_time = tot_time,
                      f.lamb = f.lamb,
                      f.mu = f.mu,
                      lamb_par = lamb_par_init,
                      mu_par = mu_par_init,
                      f = f,
                      cst.lamb = TRUE,
                      cond = "crown",
                      dt = 1e-3,
                      df = dof_frag)

save(Bcst_Dfrag,
     file = "./caudata_Bcst_Dfrag.Rdata")

Bcst_Dfrag <- fit_env(phylo = tree,
                      env_data = data_frag,
                      tot_time = tot_time,
                      f.lamb = f.lamb,
                      f.mu = f.mu,
                      lamb_par = lamb_par_init,
                      mu_par = mu_par_init,
                      f = f,
                      cst.lamb = TRUE,
                      cond = "crown",
                      dt = 1e-3,
                      df = 3)

save(Bcst_Dfrag,
     file = "./caudata_Bcst_Dfrag_df3.Rdata")

Bcst_Dfrag_df025 <- fit_env(phylo = tree,
                            env_data = data_frag,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.lamb = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = dof_frag*0.25)

save(Bcst_Dfrag_df025,
     file = "./caudata_Bcst_Dfrag_df025.Rdata")

Bcst_Dfrag_df05 <- fit_env(phylo = tree,
                           env_data = data_frag,
                           tot_time = tot_time,
                           f.lamb = f.lamb,
                           f.mu = f.mu,
                           lamb_par = lamb_par_init,
                           mu_par = mu_par_init,
                           f = f,
                           cst.lamb = TRUE,
                           cond = "crown",
                           dt = 1e-3,
                           df = dof_frag*0.5)

save(Bcst_Dfrag_df05,
     file = "./caudata_Bcst_Dfrag_df05.Rdata")

Bcst_Dfrag_df075 <- fit_env(phylo = tree,
                            env_data = data_frag,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.lamb = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = dof_frag*0.75)

save(Bcst_Dfrag_df075,
     file = "./caudata_Bcst_Dfrag_df075.Rdata")

Bcst_Dfragdet <- fit_env(phylo = tree,
                         env_data = frag_detrended,
                         tot_time = tot_time,
                         f.lamb = f.lamb,
                         f.mu = f.mu,
                         lamb_par = lamb_par_init,
                         mu_par = mu_par_init,
                         f = f,
                         cst.lamb = TRUE,
                         cond = "crown",
                         dt = 1e-3)
save(Bcst_Dfragdet,
     file = "./caudata_Bcst_Dfragdet.Rdata")

Bcst_Dfragdet_df025 <- fit_env(phylo = tree,
                               env_data = frag_detrended,
                               tot_time = tot_time,
                               f.lamb = f.lamb,
                               f.mu = f.mu,
                               lamb_par = lamb_par_init,
                               mu_par = mu_par_init,
                               f = f,
                               cst.lamb = TRUE,
                               cond = "crown",
                               dt = 1e-3,
                               df = dof_frag_det*0.25)
save(Bcst_Dfragdet_df025,
     file = "./caudata_Bcst_Dfragdet_df025.Rdata")

Bcst_Dfragdet_df05 <- fit_env(phylo = tree,
                              env_data = frag_detrended,
                              tot_time = tot_time,
                              f.lamb = f.lamb,
                              f.mu = f.mu,
                              lamb_par = lamb_par_init,
                              mu_par = mu_par_init,
                              f = f,
                              cst.lamb = TRUE,
                              cond = "crown",
                              dt = 1e-3,
                              df = dof_frag_det*0.5)
save(Bcst_Dfragdet_df05,
     file = "./caudata_Bcst_Dfragdet_df05.Rdata")

Bcst_Dfragdet_df075 <- fit_env(phylo = tree,
                               env_data = frag_detrended,
                               tot_time = tot_time,
                               f.lamb = f.lamb,
                               f.mu = f.mu,
                               lamb_par = lamb_par_init,
                               mu_par = mu_par_init,
                               f = f,
                               cst.lamb = TRUE,
                               cond = "crown",
                               dt = 1e-3,
                               df = dof_frag_det*0.75)
save(Bcst_Dfragdet_df075,
     file = "./caudata_Bcst_Dfragdet_df075.Rdata")

#---------- Benv Dcst mtb 

# exponential model with speciation rate depending on global average temperature through time and constant 
# extinction rate
f.lamb <- function(t, x, y) {y[1] * exp(-y[2] / x)}
f.mu <- function(t, x, y) {y[1]}
lamb_par_init <- c(0.001, 0)
mu_par_init <- c(0)

Btempcond22_Dcst <- fit_env(phylo = tree,
                            env_data = data_temp_cond_22,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.mu = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = dof_temp_cond_22)

save(Btempcond22_Dcst,
     file = "./caudata_Btempcond22_Dcst_mtb.Rdata")

Btempcond22_Dcst <- fit_env(phylo = tree,
                            env_data = data_temp_cond_22,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.mu = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = 3)

save(Btempcond22_Dcst,
     file = "./caudata_Btempcond22_Dcst_df3_mtb.Rdata")

Btempcond22_Dcst_df025_mtb <- fit_env(phylo = tree,
                                      env_data = data_temp_cond_22,
                                      tot_time = tot_time,
                                      f.lamb = f.lamb,
                                      f.mu = f.mu,
                                      lamb_par = lamb_par_init,
                                      mu_par = mu_par_init,
                                      f = f,
                                      cst.mu = TRUE,
                                      cond = "crown",
                                      dt = 1e-3,
                                      df = dof_temp_cond_22*0.25)

save(Btempcond22_Dcst_df025_mtb,
     file = "./caudata_Btempcond22_Dcst_df025_mtb.Rdata")

Btempcond22_Dcst_df05_mtb <- fit_env(phylo = tree,
                                     env_data = data_temp_cond_22,
                                     tot_time = tot_time,
                                     f.lamb = f.lamb,
                                     f.mu = f.mu,
                                     lamb_par = lamb_par_init,
                                     mu_par = mu_par_init,
                                     f = f,
                                     cst.mu = TRUE,
                                     cond = "crown",
                                     dt = 1e-3,
                                     df = dof_temp_cond_22*0.5)

save(Btempcond22_Dcst_df05_mtb,
     file = "./caudata_Btempcond22_Dcst_df05_mtb.Rdata")

Btempcond22_Dcst_df075_mtb <- fit_env(phylo = tree,
                                      env_data = data_temp_cond_22,
                                      tot_time = tot_time,
                                      f.lamb = f.lamb,
                                      f.mu = f.mu,
                                      lamb_par = lamb_par_init,
                                      mu_par = mu_par_init,
                                      f = f,
                                      cst.mu = TRUE,
                                      cond = "crown",
                                      dt = 1e-3,
                                      df = dof_temp_cond_22*0.75)

save(Btempcond22_Dcst_df075_mtb,
     file = "./caudata_Btempcond22_Dcst_df075_mtb.Rdata")

Btempcond22det_Dcst <- fit_env(phylo = tree,
                               env_data = temp_cond_22_detrended,
                               tot_time = tot_time,
                               f.lamb = f.lamb,
                               f.mu = f.mu,
                               lamb_par = lamb_par_init,
                               mu_par = mu_par_init,
                               f = f,
                               cst.mu = TRUE,
                               cond = "crown",
                               dt = 1e-3,
                               df = dof_temp_cond_22_det)
save(Btempcond22det_Dcst,
     file = "./caudata_Btempcond22det_Dcst_mtb.Rdata")

Btempcond22det_Dcst_df025_mtb <- fit_env(phylo = tree,
                                         env_data = temp_cond_22_detrended,
                                         tot_time = tot_time,
                                         f.lamb = f.lamb,
                                         f.mu = f.mu,
                                         lamb_par = lamb_par_init,
                                         mu_par = mu_par_init,
                                         f = f,
                                         cst.mu = TRUE,
                                         cond = "crown",
                                         dt = 1e-3,
                                         df = dof_temp_cond_22_det*0.25)
save(Btempcond22det_Dcst_df025_mtb,
     file = "./caudata_Btempcond22det_Dcst_df025_mtb.Rdata")

Btempcond22det_Dcst_df05_mtb <- fit_env(phylo = tree,
                                        env_data = temp_cond_22_detrended,
                                        tot_time = tot_time,
                                        f.lamb = f.lamb,
                                        f.mu = f.mu,
                                        lamb_par = lamb_par_init,
                                        mu_par = mu_par_init,
                                        f = f,
                                        cst.mu = TRUE,
                                        cond = "crown",
                                        dt = 1e-3,
                                        df = dof_temp_cond_22_det*0.5)
save(Btempcond22det_Dcst_df05_mtb,
     file = "./caudata_Btempcond22det_Dcst_df05_mtb.Rdata")

Btempcond22det_Dcst_df075_mtb <- fit_env(phylo = tree,
                                         env_data = temp_cond_22_detrended,
                                         tot_time = tot_time,
                                         f.lamb = f.lamb,
                                         f.mu = f.mu,
                                         lamb_par = lamb_par_init,
                                         mu_par = mu_par_init,
                                         f = f,
                                         cst.mu = TRUE,
                                         cond = "crown",
                                         dt = 1e-3,
                                         df = dof_temp_cond_22_det*0.75)
save(Btempcond22det_Dcst_df075_mtb,
     file = "./caudata_Btempcond22det_Dcst_df075_mtb.Rdata")

#---------- Bcst Denv mtb 

# exponential model with constant birth, death varying with environamental variable
f.lamb <- function(t, x, y) {y[1] }
f.mu <- function(t, x, y) {y[1] * exp(-y[2] / x)}
lamb_par_init <- c(0.0001)
mu_par_init <- c(0.0000002, 0.000000001)

Bcst_Dtempcond22 <- fit_env(phylo = tree,
                            env_data = data_temp_cond_22,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.lamb = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = dof_temp_cond_22)

save(Bcst_Dtempcond22,
     file = "./caudata_Bcst_Dtempcond22_mtb.Rdata")

Bcst_Dtempcond22 <- fit_env(phylo = tree,
                            env_data = data_temp_cond_22,
                            tot_time = tot_time,
                            f.lamb = f.lamb,
                            f.mu = f.mu,
                            lamb_par = lamb_par_init,
                            mu_par = mu_par_init,
                            f = f,
                            cst.lamb = TRUE,
                            cond = "crown",
                            dt = 1e-3,
                            df = 3)

save(Bcst_Dtempcond22,
     file = "./caudata_Bcst_Dtempcond22_df3_mtb.Rdata")

Bcst_Dtempcond22_df025_mtb <- fit_env(phylo = tree,
                                      env_data = data_temp_cond_22,
                                      tot_time = tot_time,
                                      f.lamb = f.lamb,
                                      f.mu = f.mu,
                                      lamb_par = lamb_par_init,
                                      mu_par = mu_par_init,
                                      f = f,
                                      cst.lamb = TRUE,
                                      cond = "crown",
                                      dt = 1e-3,
                                      df = dof_temp_cond_22*0.25)

save(Bcst_Dtempcond22_df025_mtb,
     file = "./caudata_Bcst_Dtempcond22_df025_mtb.Rdata")

Bcst_Dtempcond22_df05_mtb <- fit_env(phylo = tree,
                                     env_data = data_temp_cond_22,
                                     tot_time = tot_time,
                                     f.lamb = f.lamb,
                                     f.mu = f.mu,
                                     lamb_par = lamb_par_init,
                                     mu_par = mu_par_init,
                                     f = f,
                                     cst.lamb = TRUE,
                                     cond = "crown",
                                     dt = 1e-3,
                                     df = dof_temp_cond_22*0.5)

save(Bcst_Dtempcond22_df05_mtb,
     file = "./caudata_Bcst_Dtempcond22_df05_mtb.Rdata")

Bcst_Dtempcond22_df075_mtb <- fit_env(phylo = tree,
                                      env_data = data_temp_cond_22,
                                      tot_time = tot_time,
                                      f.lamb = f.lamb,
                                      f.mu = f.mu,
                                      lamb_par = lamb_par_init,
                                      mu_par = mu_par_init,
                                      f = f,
                                      cst.lamb = TRUE,
                                      cond = "crown",
                                      dt = 1e-3,
                                      df = dof_temp_cond_22*0.75)

save(Bcst_Dtempcond22_df075_mtb,
     file = "./caudata_Bcst_Dtempcond22_df075_mtb.Rdata")

Bcst_Dtempcond22det <- fit_env(phylo = tree,
                               env_data = temp_cond_22_detrended,
                               tot_time = tot_time,
                               f.lamb = f.lamb,
                               f.mu = f.mu,
                               lamb_par = lamb_par_init,
                               mu_par = mu_par_init,
                               f = f,
                               cst.lamb = TRUE,
                               cond = "crown",
                               dt = 1e-3,
                               df = dof_temp_cond_22_det)
save(Bcst_Dtempcond22det,
     file = "./caudata_Bcst_Dtempcond22det_mtb.Rdata")

Bcst_Dtempcond22det_df025_mtb <- fit_env(phylo = tree,
                                         env_data = temp_cond_22_detrended,
                                         tot_time = tot_time,
                                         f.lamb = f.lamb,
                                         f.mu = f.mu,
                                         lamb_par = lamb_par_init,
                                         mu_par = mu_par_init,
                                         f = f,
                                         cst.lamb = TRUE,
                                         cond = "crown",
                                         dt = 1e-3,
                                         df = dof_temp_cond_22_det*0.25)
save(Bcst_Dtempcond22det_df025_mtb,
     file = "./caudata_Bcst_Dtempcond22det_df025_mtb.Rdata")

Bcst_Dtempcond22det_df05_mtb <- fit_env(phylo = tree,
                                        env_data = temp_cond_22_detrended,
                                        tot_time = tot_time,
                                        f.lamb = f.lamb,
                                        f.mu = f.mu,
                                        lamb_par = lamb_par_init,
                                        mu_par = mu_par_init,
                                        f = f,
                                        cst.lamb = TRUE,
                                        cond = "crown",
                                        dt = 1e-3,
                                        df = dof_temp_cond_22_det*0.5)
save(Bcst_Dtempcond22det_df05_mtb,
     file = "./caudata_Bcst_Dtempcond22det_df05_mtb.Rdata")

Bcst_Dtempcond22det_df075_mtb <- fit_env(phylo = tree,
                                         env_data = temp_cond_22_detrended,
                                         tot_time = tot_time,
                                         f.lamb = f.lamb,
                                         f.mu = f.mu,
                                         lamb_par = lamb_par_init,
                                         mu_par = mu_par_init,
                                         f = f,
                                         cst.lamb = TRUE,
                                         cond = "crown",
                                         dt = 1e-3,
                                         df = dof_temp_cond_22_det*0.75)
save(Bcst_Dtempcond22det_df075_mtb,
     file = "./caudata_Bcst_Dtempcond22det_df075_mtb.Rdata")

##### run environmental birth-death models - per family =========================================================================

#---------- phylogenies

# load sub-trees
library(ape)
files <- list.files("./bd_models_per_family/", pattern = "\\.nexus$", full.names = TRUE)

for (file in files) {
  file_name <- gsub(".nexus", "", basename(file))  # Extraire le nom de fichier sans l'extension
  assign(file_name, read.nexus(file))
}

#---------- data tables

# load taxo tables
files <- list.files("./bd_models_per_family/", pattern = "\\.csv$", full.names = TRUE)

for (file in files) {
  file_name <- gsub(".csv", "", basename(file))  # Extraire le nom de fichier sans l'extension
  assign(file_name, read.csv(file))
}

#---------- temperature condamine 2022

# hynobiidae
temp_hynobiidae <- read.csv2(file = "./bd_models_per_family/temp_hynobiidae.csv")
dof_temp_hynobiidae <- smooth.spline(temp_hynobiidae[, "time"], temp_hynobiidae[, "temperature"])$df

# salamandridae
temp_salamandridae <- read.csv2(file = "./bd_models_per_family/temp_salamandridae.csv")
dof_temp_salamandridae <- smooth.spline(temp_salamandridae[, "time"], temp_salamandridae[, "temperature"])$df

# plethodontidae
temp_plethodontidae <- read.csv2(file = "./bd_models_per_family/temp_plethodontidae.csv")
dof_temp_plethodontidae <- smooth.spline(temp_plethodontidae[, "time"], temp_plethodontidae[, "temperature"])$df

#---------- continental fragmentation

# hynobiidae
frag_hynobiidae <- read.csv2(file = "./bd_models_per_family/frag_hynobiidae.csv")
dof_frag_hynobiidae <- smooth.spline(frag_hynobiidae[, "time"], frag_hynobiidae[, "fragmentation"])$df

# salamandridae
frag_salamandridae <- read.csv2(file = "./bd_models_per_family/frag_salamandridae.csv")
dof_frag_salamandridae <- smooth.spline(frag_salamandridae[, "time"], frag_salamandridae[, "fragmentation"])$df

# plethodontidae
frag_plethodontidae <- read.csv2(file = "./bd_models_per_family/frag_plethodontidae.csv")
dof_frag_plethodontidae <- smooth.spline(frag_plethodontidae[, "time"], frag_plethodontidae[, "fragmentation"])$df

#---------- detrended temperature condamine 2022

# hynobiidae
temp_hynobiidae_det <- read.csv2(file = "./bd_models_per_family/temp_hynobiidae_det.csv")
dof_temp_hynobiidae_det <- smooth.spline(temp_hynobiidae_det[, "time"], temp_hynobiidae_det[, "detrended_temperature"])$df

# salamandridae
temp_salamandridae_det <- read.csv2(file = "./bd_models_per_family/temp_salamandridae_det.csv")
dof_temp_salamandridae_det <- smooth.spline(temp_salamandridae_det[, "time"], temp_salamandridae_det[, "detrended_temperature"])$df

# plethodontidae
temp_plethodontidae_det <- read.csv2(file = "./bd_models_per_family/temp_plethodontidae_det.csv")
dof_temp_plethodontidae_det <- smooth.spline(temp_plethodontidae_det[, "time"], temp_plethodontidae_det[, "detrended_temperature"])$df

#---------- detrended continental fragmentation

# hynobiidae
frag_hynobiidae_det <- read.csv2(file = "./bd_models_per_family/frag_hynobiidae_det.csv")
dof_frag_hynobiidae_det <- smooth.spline(frag_hynobiidae_det[, "time"], frag_hynobiidae_det[, "detrended_fragmentation"])$df

# salamandridae
frag_salamandridae_det <- read.csv2(file = "./bd_models_per_family/frag_salamandridae_det.csv")
dof_frag_salamandridae_det <- smooth.spline(frag_salamandridae_det[, "time"], frag_salamandridae_det[, "detrended_fragmentation"])$df

# plethodontidae
frag_plethodontidae_det <- read.csv2(file = "./bd_models_per_family/frag_plethodontidae_det.csv")
dof_frag_plethodontidae_det <- smooth.spline(frag_plethodontidae_det[, "time"], frag_plethodontidae_det[, "detrended_fragmentation"])$df

#---------- calculate total time for each family

hynobiidae_tot_time <- max(TreeSim::getx(Hynobiidae_tree))
plethodontidae_tot_time <- max(TreeSim::getx(Plethodontidae_tree))
salamandridae_tot_time <- max(TreeSim::getx(Salamandridae_tree))

#---------- Set the sampling fraction

f = 1

#---------- set lists of data

tree_list <- list(hynobiidae = ladderize(Hynobiidae_tree, right = TRUE),
                  plethodontidae = ladderize(Plethodontidae_tree, right = TRUE),
                  salamandridae = ladderize(Salamandridae_tree, right = TRUE)
                  )

#---------- BD models 

library(RPANDA)

names_list <- names(tree_list)

for (i in names(tree_list)){

  tree <- tree_list[[i]]
  
  tot_time <- get(paste0(i, "_tot_time"))
  
  data_temp <- get(paste0("temp_",i))
  dof_temp <- get(paste0("dof_temp_",i))
  
  data_frag <- get(paste0("frag_",i))
  dof_frag <- get(paste0("dof_frag_",i)) 
  
  data_temp_det <- get(paste0("temp_",i,"_det"))
  dof_temp_det <- get(paste0("dof_temp_",i,"_det"))
  
  data_frag_det <- get(paste0("frag_",i,"_det"))
  dof_frag_det <- get(paste0("dof_frag_",i,"_det")) 
  
  cond <- "crown"

  ################################################################################################################

  f.lamb <- function(t, y) {y[1]}
  f.mu <- function(t, y) {y[1]}
  lamb_par <- c(0.01)
  mu_par <- c(0.001)

  model <- RPANDA::fit_bd(
    phylo = tree,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par,
    mu_par = mu_par,
    f = f,
    cst.lamb = TRUE,
    cst.mu = TRUE,
    expo.lamb = FALSE,
    expo.mu = FALSE,
    fix.mu = FALSE,
    cond = cond,
    dt = 1e-3
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dcst.Rdata"))

  ################################################################################################################

  f.lamb <- function(t, y) {y[1] * exp(y[2] * t)}
  f.mu <- function(t, y) {y[1]}
  lamb_par <- c(0.01, 0.001)
  mu_par <- c(0.001)

  model <- RPANDA::fit_bd(
    phylo = tree,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par,
    mu_par = mu_par,
    f = f,
    cst.lamb = FALSE,
    cst.mu = TRUE,
    expo.lamb = TRUE,
    expo.mu = FALSE,
    fix.mu = FALSE,
    cond = cond,
    dt = 1e-3
  )

  save(model,
       file = paste0("./",i,"_Bvar_Dcst.Rdata"))

  ################################################################################################################

  f.lamb <- function(t, y) {y[1]}
  f.mu <- function(t, y) {y[1] * exp(y[2] * t)}
  lamb_par <- c(0.01)
  mu_par <- c(0.001, 0.0001)

  model <- RPANDA::fit_bd(
    phylo = tree,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par,
    mu_par = mu_par,
    f = f,
    cst.lamb = TRUE,
    cst.mu = FALSE,
    expo.lamb = FALSE,
    expo.mu = TRUE,
    fix.mu = FALSE,
    cond = cond,
    dt = 1e-3
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dvar.Rdata"))
  
  ################################################################################################################

  f.lamb <- function(t, x, y) {y[1] * exp(y[2] * x)}
  f.mu <- function(t, x, y) {y[1]}
  lamb_par_init <- c(0.01, 0.001)
  mu_par_init <- c(0.001)

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp
  )

  save(model,
       file = paste0("./",i,"_Btemp_Dcst.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag
  )

  save(model,
       file = paste0("./",i,"_Bfrag_Dcst.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det
  )

  save(model,
       file = paste0("./",i,"_Btempdet_Dcst.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_frag_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag_det
  )

  save(model,
       file = paste0("./",i,"_Bfragdet_Dcst.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = 3
  )

  save(model,
       file = paste0("./",i,"_Btemp_Dcst_df3.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = 3
  )

  save(model,
       file = paste0("./",i,"_Bfrag_Dcst_df3.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Btemp_Dcst_df025.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Btemp_Dcst_df05.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Btemp_Dcst_df075.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Bfrag_Dcst_df025.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Bfrag_Dcst_df05.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Bfrag_Dcst_df075.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Btempdet_Dcst_df025.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Btempdet_Dcst_df05.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Btempdet_Dcst_df075.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag_det*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Bfragdet_Dcst_df025.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag_det*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Bfragdet_Dcst_df05.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag_det*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Bfragdet_Dcst_df075.Rdata"))
  
  ################################################################################################################
  
  f.lamb <- function(t, x, y) {y[1]}
  f.mu <- function(t, x, y) {y[1] * exp(y[2] * x)}
  lamb_par_init <- c(0.01)
  mu_par_init <- c(0.001, 0.0001)

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dtemp.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dfrag.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dtempdet.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_frag_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag_det
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dfragdet.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = 3
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_df3.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp_shav,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = 3
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dtempshav_df3.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = 3
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dfrag_df3.Rdata"))

  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_df025.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_df05.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_df075.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dfrag_df025.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dfrag_df05.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dfrag_df075.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtempdet_df025.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtempdet_df05.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtempdet_df075.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag_det*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dfragdet_df025.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag_det*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dfragdet_df05.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_frag_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_frag_det*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dfragdet_df075.Rdata"))
  
  ################################################################################################################
  
  f.lamb <- function(t, x, y) {y[1] * exp(-y[2] / x)}
  f.mu <- function(t, x, y) {y[1]}
  lamb_par_init <- c(0.001, 0)
  mu_par_init <- c(0.001)

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp
  )

  save(model,
       file = paste0("./",i,"_Btemp_Dcst_mtb.Rdata"))


  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det
  )

  save(model,
       file = paste0("./",i,"_Btempdet_Dcst_mtb.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = 3
  )

  save(model,
       file = paste0("./",i,"_Btemp_Dcst_mtb_df3.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp_shav,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = 3
  )

  save(model,
       file = paste0("./",i,"_Btempshav_Dcst_mtb_df3.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Btemp_Dcst_df025_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Btemp_Dcst_df05_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Btemp_Dcst_df075_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Btempdet_Dcst_df025_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Btempdet_Dcst_df05_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.mu = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Btempdet_Dcst_df075_mtb.Rdata"))
  
  ################################################################################################################
  
  f.lamb <- function(t, x, y) {y[1]}
  f.mu <- function(t, x, y) {y[1] * exp(-y[2] / x)}
  lamb_par_init <- c(0.01)
  mu_par_init <- c(0.001, 0)

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_mtb.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dtempdet_mtb.Rdata"))

  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = 3
  )

  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_mtb_df3.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_df025_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_df05_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtemp_df075_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.25
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtempdet_df025_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.5
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtempdet_df05_mtb.Rdata"))
  
  model <- fit_env(
    phylo = tree,
    env_data = data_temp_det,
    tot_time = tot_time,
    f.lamb = f.lamb,
    f.mu = f.mu,
    lamb_par = lamb_par_init,
    mu_par = mu_par_init,
    f = f,
    cst.lamb = TRUE,
    cond = cond,
    dt = 1e-3,
    df = dof_temp_det*0.75
  )
  
  save(model,
       file = paste0("./",i,"_Bcst_Dtempdet_df075_mtb.Rdata"))
  
}

##### extract results from environmental birth-death models - all caudata ================================================================

#---------- load data

# create a function to load data
load_object <- function(file) {
  env <- new.env()
  load(file, envir = env)
  return(env[[ls(env)]])
}

# list of files to load
path_to_files <- "./"

files <- list.files(path_to_files,
                    pattern = "\\.Rdata$",
                    full.names = TRUE)

list_caudata <- list()

# load and arrange data
for (file in files) {
  
  object <- load_object(file)
  
  object_name <- gsub("\\.Rdata$", "", basename(file))
  
  list_caudata[[object_name]] <- object
}

#---------- extract results for all caudata

# create a function to extract results
get_values <- function(list, name) {
  values <- if (name %in% names(list)) {
    value <- list[[name]]
    if (length(value) < 2) {
      c(value, NA)
    } else {
      value[1:2]
    }
  } else {
    rep(NA, 2)
  }
  values
}

# create dataframes for all results
rslts <- list_caudata

matrix_results <- sapply(rslts, function(list) {
  c(list$model, 
    list$LH, 
    list$aicc, 
    get_values(list, "lamb_par"), 
    get_values(list, "mu_par"))
})

# transpose the matrix
matrix_results_t <- t(matrix_results)

# convert into dataframe
table_results <- as.data.frame(matrix_results_t,
                               stringsAsFactors = FALSE)
table_results <- cbind(Model = row.names(table_results),
                       table_results,
                       row.names = NULL)
colnames(table_results) <- c("Model", "Model_type", "LH", "AICC", "Lambda", "Alpha", "Mu", "Beta")
table_results$LH <- as.numeric(table_results$LH)
table_results$AICC <- as.numeric(table_results$AICC)

# add delta aicc and order data
library(dplyr)
table_results <- table_results %>% mutate(Delta_AICC = AICC - min(AICC))
table_results <- table_results[order(table_results$AICC), ]

##### extract results from environmental birth-death models - per family =============================================================

#---------- load data

# create a function to load data
load_object <- function(file) {
  env <- new.env()
  load(file, envir = env)
  return(env[[ls(env)]])
}

# list of files to load
path_to_files <- "./"

files <- list.files(path_to_files,
                    pattern = "\\.Rdata$",
                    full.names = TRUE)

# list_ambystomatidae <- list()
list_hynobiidae <- list()
list_plethodontidae <- list()
list_salamandridae <- list()

# load and arrange data
for (file in files) {

  object <- load_object(file)
  
  object_name <- gsub("\\.Rdata$", "", basename(file))
  
  if (startsWith(object_name, "ambystomatidae")) {
    list_ambystomatidae[[object_name]] <- object
  } else if (startsWith(object_name, "hynobiidae")) {
    list_hynobiidae[[object_name]] <- object
  } else if (startsWith(object_name, "plethodontidae")) {
    list_plethodontidae[[object_name]] <- object
  } else if (startsWith(object_name, "salamandridae")) {
    list_salamandridae[[object_name]] <- object
  }
}

#---------- extract results 

# create a function to extract results
get_values <- function(list, name) {
  values <- if (name %in% names(list)) {
    value <- list[[name]]
    if (length(value) < 2) {
      c(value, NA)
    } else {
      value[1:2]
    }
  } else {
    rep(NA, 2)
  }
  values
}

# extract results 
per_family_results <- list(
  # ambystomatidae = list_ambystomatidae,
  hynobiidae = list_hynobiidae,
  plethodontidae = list_plethodontidae,
  salamandridae = list_salamandridae
)

# create dataframes for all results
for (i in names(per_family_results)){
  
  rslts <- per_family_results[[i]]
  
  matrix_results <- sapply(rslts, function(list) {
    c(list$model, 
      list$LH, 
      list$aicc, 
      get_values(list, "lamb_par"), 
      get_values(list, "mu_par"))
  })
  
  # transpose the matrix
  matrix_results_t <- t(matrix_results)
  
  # convert into dataframe
  table_results <- as.data.frame(matrix_results_t,
                                 stringsAsFactors = FALSE)
  table_results <- cbind(Model = row.names(table_results),
                         table_results,
                         row.names = NULL)
  colnames(table_results) <- c("Model", "Model_type", "LH", "AICC", "Lambda", "Alpha", "Mu", "Beta")
  table_results$LH <- as.numeric(table_results$LH)
  table_results$AICC <- as.numeric(table_results$AICC)
  
  # add delta aicc and order data
  library(dplyr)
  table_results <- table_results %>% mutate(Delta_AICC = AICC - min(AICC))
  table_results <- table_results[order(table_results$AICC), ]
  
  # save results
  write.csv2(table_results,
             file = paste0("./",i,"_models_results.csv"))
}

##### ancestral state estimation - all caudata ==============================================================================

#---------- data preparation

caudata_data_anc_state_rec <- caudata_data

# change colnames
colnames(caudata_data_anc_state_rec)[4] <- "lc"

# change lc categories
unique(caudata_data_anc_state_rec$lc)
caudata_data_anc_state_rec$lc[caudata_data_anc_state_rec$lc %in% c("pd1", "pd1*", "pd2", "pd2*", "pd3", "pd3*", "pd4", "pd4*")] <- "pd"
caudata_data_anc_state_rec$lc[caudata_data_anc_state_rec$lc %in% "bi*"] <- "bi"
caudata_data_anc_state_rec$lc[caudata_data_anc_state_rec$lc %in% "f-bi*"] <- "f-bi"
caudata_data_anc_state_rec$lc[caudata_data_anc_state_rec$lc %in% "dd*"] <- "dd"
caudata_data_anc_state_rec$lc[caudata_data_anc_state_rec$lc %in% c("ovi", "ovi*", "f-vila", "f-vila*", "vipu", "vipu*")] <- "vi"

# 5 states life cycle for corHMM
mapping <- c("bi" = 0, 
             "dd" = 1,
             "f-bi" = 2, 
             "pd" = 3,
             "vi" = 4)
caudata_data_anc_state_rec$lc_5states_corHMM <- mapping[caudata_data_anc_state_rec$lc]

# extract characters as vectors
lc_5_corHMM <- caudata_data_anc_state_rec[, c(3, 5)]

#---------- ase runs

library(corHMM)
library(foreach)
library(doParallel)

#---------- prepare data

tree <- caudata_tree
x_cor <- lc_5_corHMM

#---------- create a function to run analyses

run_corHMM <- function(tree, data, rate.cat, model, rate.mat, timeout_duration = 600) {
  obj <- list()
  class(obj) <- "try-error"
  while (inherits(obj, "try-error") || is.null(obj)) {
    obj <- tryCatch({
      R.utils::withTimeout({
        corHMM::corHMM(
          phy = tree,
          data = data,
          rate.cat = rate.cat,
          model = model,
          rate.mat = rate.mat,
          node.states = "joint",
          get.tip.states = TRUE,
          nstarts = 1
        )
      }, timeout = timeout_duration, onTimeout = "silent")
    }, error = function(e) {
      message("Error encountered: ", e$message)
      NULL
    })
    
    if (is.null(obj)) {
      message("Retrying analysis due to timeout or error...")
    }
  }
  return(obj)
}

#---------- ER

niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(50, type = "PSOCK")
registerDoParallel(cl = mc)

all_er_cor_1r <- foreach(i = 1:niter, .packages = c("corHMM", "R.utils")) %dopar% {
  run_corHMM(tree = tree,
             data = x_cor,
             rate.cat = 1,
             model = "ER",
             rate.mat = NULL)
}
stopCluster(mc) 

save(all_er_cor_1r, 
     file = "./ase/all_er_cor_1r.Rdata")

# best fitting model
min_aicc <- Inf
min_aicc_index <- NULL

list_rslts <- all_er_cor_1r
for (i in 1:length(list_rslts)) {
  if (list_rslts[[i]]$AICc < min_aicc) {
    min_aicc <- list_rslts[[i]]$AICc
    min_aicc_index <- i
  }
}

cor_best_er_1r <- list_rslts[[min_aicc_index]]

save(cor_best_er_1r, 
     file = "./ase/cor_best_er_1r.Rdata")

niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(50, type = "PSOCK")
registerDoParallel(cl = mc)

all_er_cor_2r <- foreach(i = 1:niter, .packages = c("corHMM", "R.utils")) %dopar% {
  run_corHMM(tree = tree,
             data = x_cor,
             rate.cat = 2,
             model = "ER",
             rate.mat = NULL)
}
stopCluster(mc) 

save(all_er_cor_2r, 
     file = "./ase/all_er_cor_2r.Rdata")

# best fitting model
min_aicc <- Inf
min_aicc_index <- NULL

list_rslts <- all_er_cor_2r
for (i in 1:length(list_rslts)) {
  if (list_rslts[[i]]$AICc < min_aicc) {
    min_aicc <- list_rslts[[i]]$AICc
    min_aicc_index <- i
  }
}

best_er_cor_2r <- list_rslts[[min_aicc_index]]

save(best_er_cor_2r, 
     file = "./ase/best_er_cor_2r.Rdata")

#---------- SYM

niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(50, type = "PSOCK")
registerDoParallel(cl = mc)

all_sym_cor_1r <- foreach(i = 1:niter, .packages = c("corHMM", "R.utils")) %dopar% {
  run_corHMM(tree = tree,
             data = x_cor,
             rate.cat = 1,
             model = "SYM",
             rate.mat = NULL)
}
stopCluster(mc) 

save(all_sym_cor_1r, 
     file = "./ase/all_sym_cor_1r.Rdata")

# best fitting model
min_aicc <- Inf
min_aicc_index <- NULL

list_rslts <- all_sym_cor_1r
for (i in 1:length(list_rslts)) {
  if (list_rslts[[i]]$AICc < min_aicc) {
    min_aicc <- list_rslts[[i]]$AICc
    min_aicc_index <- i
  }
}

cor_best_sym_1r <- list_rslts[[min_aicc_index]]

save(cor_best_sym_1r, 
     file = "./ase/cor_best_sym_1r.Rdata")

niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(50, type = "PSOCK")
registerDoParallel(cl = mc)

all_sym_cor_2r <- foreach(i = 1:niter, .packages = c("corHMM", "R.utils")) %dopar% {
  run_corHMM(tree = tree,
             data = x_cor,
             rate.cat = 2,
             model = "SYM",
             rate.mat = NULL)
}
stopCluster(mc)

save(all_sym_cor_2r, 
     file = "./ase/all_sym_cor_2r.Rdata")

# best fitting model
min_aicc <- Inf
min_aicc_index <- NULL

list_rslts <- Filter(function(x) contains_null(x), all_sym_cor_2r)

for (i in 1:length(list_rslts)) {
  if (list_rslts[[i]]$AICc < min_aicc) {
    min_aicc <- list_rslts[[i]]$AICc
    min_aicc_index <- i
  }
}

best_sym_cor_2r <- list_rslts[[min_aicc_index]]

save(best_sym_cor_2r, 
     file = "./ase/best_sym_cor_2r.Rdata")

#---------- ARD

niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(50, type = "PSOCK")
registerDoParallel(cl = mc)

all_ard_cor_1r <- foreach(i = 1:niter, .packages = c("corHMM", "R.utils")) %dopar% {
  run_corHMM(tree = tree,
             data = x_cor,
             rate.cat = 1,
             model = "ARD",
             rate.mat = NULL)
}
stopCluster(mc) 

save(all_ard_cor_1r, 
     file = "./ase/all_ard_cor_1r.Rdata")

# best fitting model
min_aicc <- Inf
min_aicc_index <- NULL

list_rslts <- all_ard_cor_1r
for (i in 1:length(list_rslts)) {
  if (list_rslts[[i]]$AICc < min_aicc) {
    min_aicc <- list_rslts[[i]]$AICc
    min_aicc_index <- i
  }
}

cor_best_ard_1r <- list_rslts[[min_aicc_index]]

save(cor_best_ard_1r, 
     file = "./ase/cor_best_ard_1r.Rdata")

niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(50, type = "PSOCK")
registerDoParallel(cl = mc)

all_ard_cor_2r <- foreach(i = 1:niter, .packages = c("corHMM", "R.utils")) %dopar% {
  run_corHMM(tree = tree,
             data = x_cor,
             rate.cat = 2,
             model = "ARD",
             rate.mat = NULL)
}

save(all_ard_cor_2r, 
     file = "./ase/all_ard_cor_2r.Rdata")

# best fitting model
min_aicc <- Inf
min_aicc_index <- NULL

list_rslts <- all_ard_cor_2r
for (i in 1:length(list_rslts)) {
  if (list_rslts[[i]]$AICc < min_aicc) {
    min_aicc <- list_rslts[[i]]$AICc
    min_aicc_index <- i
  }
}

best_ard_cor_2r <- list_rslts[[min_aicc_index]]

save(best_ard_cor_2r, 
     file = "./ase/best_ard_cor_2r.Rdata")

##### secsse ==============================================================================================================

#---------- data preparation

caudata_data_secsse <- caudata_data

# change colnames
colnames(caudata_data_secsse)[3] <- "lc"

# change lc categories
unique(caudata_data_secsse$lc)
caudata_data_secsse$lc[caudata_data_secsse$lc %in% c("pd1", "pd1*", "pd2", "pd2*", "pd3", "pd3*", "pd4", "pd4*")] <- "pd"
caudata_data_secsse$lc[caudata_data_secsse$lc %in% "bi*"] <- "bi"
caudata_data_secsse$lc[caudata_data_secsse$lc %in% "f-bi*"] <- "f-bi"
caudata_data_secsse$lc[caudata_data_secsse$lc %in% "dd*"] <- "dd"
caudata_data_secsse$lc[caudata_data_secsse$lc %in% c("ovi", "ovi*", "f-vila", "f-vila*", "vipu", "vipu*")] <- "vi"

mapping <- c(
  "bi" = 0,
  "dd" = 1,
  "f-bi" = 2,
  "pd" = 3,
  "vi" = 4
)
caudata_data_secsse$lc_5states <- mapping[caudata_data_secsse$lc]

#---------- load the phylogeny

caudata_tree <- ape::read.nexus("./data/trees/trees_stewart_wiens_2024/caudata_tree_stewart_wiens_2024_arranged_lad.nexus")
tree <- caudata_tree

plethodon_tree <- ape::read.nexus("./data/trees/trees_stewart_wiens_2024/families_trees/plethodon_tree.nexus")
tree <- plethodon_tree

salamandridae_tree <- ape::read.nexus("./data/trees/trees_stewart_wiens_2024/families_trees/salamandridae_tree.nexus")
tree <- salamandridae_tree

#---------- load the data

caudata_lc_5states <- read.csv2(file = "./data/data_secsse/tree_stewart_wiens_2024/caudata_lc_5states.csv")
caudata_lc_5states <- caudata_lc_5states[match(caudata_tree$tip.label, caudata_lc_5states$Species), , drop = FALSE]

plet_lc_4states <- read.csv2(file = "./data/data_secsse/tree_stewart_wiens_2024/plet_lc_4states.csv")
plet_lc_4states <- plet_lc_4states[match(plethodon_tree$tip.label, plet_lc_4states$Species), , drop = FALSE]

sal_lc_3states <- read.csv2(file = "./data/data_secsse/tree_stewart_wiens_2024/sal_lc_3states.csv")
sal_lc_3states <- sal_lc_3states[match(salamandridae_tree$tip.label, sal_lc_3states$Species), , drop = FALSE]

#---------- load functions

source(file = "./secsse/functions_for_secsse.R")

#---------- prepare priors

library(DDD)
BD_model <- bd_ML(brts = ape::branching.times(tree))

init_lambda <- BD_model$lambda0
init_mu <- BD_model$mu0
init_q <- init_lambda / 3

#---------- sorting traits

library(secsse)

# caudata 5 states
traits <- sortingtraits(caudata_lc_5states, caudata_tree)

# plethodontidae 5 states
traits <- sortingtraits(plet_lc_4states, plethodon_tree)

# salamandridae 5 states
traits <- sortingtraits(sal_lc_3states, salamandridae_tree)

#---------- run non-timezone analyses

# define scenarios to test
dependencies <- c("ETD", "CTD", "CR")
inheritances <- c("dual_inheritance", "single_inheritance", "dual_symmetric_transition", "dual_asymmetric_transition")
transitions <- c("SYM")
initparsopts <- c("opt", "opt_x2", "opt_x0-5")
diff.conceal <- TRUE
parsfix <- c("none", "lambda", "mu")

# create a grid with all possible combinations
param_grid <- expand.grid(
  dependency = dependencies,
  inheritance = inheritances,
  transition = transitions,
  initparsopt = initparsopts,
  diff.conceal = diff.conceal,
  parsfix = parsfix,
  stringsAsFactors = FALSE
)

# detect/determine the number of core available
library(doParallel)
library(foreach)

n_cores <- parallel::detectCores() - 11
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# run the analysis
num_threads <- 1
sampling_fraction <- c(1, 1, 1, 1, 1)

results <- foreach(i = 1:nrow(param_grid), .packages = c("secsse")) %dopar% {
  p <- param_grid[i, ]

  cat("Running:", p$dependency, p$inheritance, p$transition, p$initparsopt, "\n")

  result <- secsse_run(
    tree = tree,
    traits = traits,
    num_concealed_states = num_concealed_states,
    dependency = p$dependency,
    inheritance = p$inheritance,
    transition = p$transition,
    diff.conceal = as.logical(p$diff.conceal),
    init_lambda = init_lambda,
    init_mu = init_mu,
    init_q = init_q,
    initparsopt = p$initparsopt,
    parsfix = p$parsfix,
    sampling_fraction = sampling_fraction,
    num_threads = num_threads,
    save_dir_idpars = "./idparslists/",
    save_dir_secsse = "./models/"
  )

  return(result)
}

# stop the cluster
stopCluster(cl)

#---------- run timezone analyses

idparslist_timezones <- secsse_idparslist_timezones(
  num_timezones = 2,
  tree = tree,
  traits = traits,
  num_concealed_states = 5,
  dependency = c("ETD", "CTD", "CR"),
  inheritance = c("dual_inheritance"),
  transition = c("ER"),
  diff.conceal = TRUE
)

# set up run parameters
initparsopts <- c("opt", "opt_x2", "opt_x0-5") # trois versions
critical_t <- c(90)
parsfix <- c("none", "lambda", "mu")
sampling_fraction <- rep(1, length(unique(traits)))
num_threads <- 1
save_dir <- "./models/"

# create a folder if necessary
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# define scenarios to test
param_grid <- expand.grid(
  idx = seq_along(idparslist_timezones),
  initparsopt = initparsopts,
  parsfix = parsfix,
  stringsAsFactors = FALSE
)

# detect/determine the number of core available
library(doParallel)
library(foreach)

n_cores <- parallel::detectCores() - 11
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# folder for log files
log_dir <- "./logs/"
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

# run parallel analyses
foreach(
  i = 1:nrow(param_grid),
  .packages = c("secsse")
) %dopar% {
  # parameters
  row <- param_grid[i, ]
  idx <- row$idx
  opt <- row$initparsopt
  fixpars <- row$parsfix

  log_file <- file.path(log_dir, sprintf("log_task_%03d.txt", i))

  # clean open connexions
  while (sink.number() > 0) sink(NULL)

  # connect logs
  con <- file(log_file, open = "wt")
  sink(con)
  sink(con, type = "message")

  cat("===== TÂCHE", i, "=====\n")
  cat("Paramètres : idparslist =", idx, ", initparsopt =", opt, ", parsfix =", fixpars, "\n")
  cat("Heure de début :", format(Sys.time()), "\n\n")

  mu_values <- c(
    init_mu, 0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
    1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4
  )
  lambda_values <- c(
    init_lambda, 0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
    1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4
  )

  success <- FALSE

  for (mu in mu_values) {
    for (lambda in lambda_values) {
      tryCatch(
        {
          cat("Essai avec init_mu =", mu, ", init_lambda =", lambda, "\n")
          res <- secsse_timezones_run(
            tree = tree,
            traits = traits,
            num_concealed_states = num_concealed_states,
            idparslist = idparslist_timezones[[idx]],
            critical_t = critical_t,
            init_lambda = lambda,
            init_mu = mu,
            init_q = init_q,
            initparsopt = opt,
            parsfix = fixpars,
            sampling_fraction = sampling_fraction,
            num_threads = 1
          )

          success <- TRUE
          cat("✅ SUCCÈS avec init_mu =", mu, ", init_lambda =", lambda, "\n")

          save_name <- paste0("secsse_timezones_idparslist_", idx, "_", fixpars, "_", opt, ".rds")
          save_path <- file.path(save_dir, save_name)
          saveRDS(res, file = save_path)

          break
        },
        error = function(e) {
          cat("❌ Échec avec init_mu =", mu, ", init_lambda =", lambda, "→", conditionMessage(e), "\n")
        }
      )
      if (success) break
    }
    if (success) break
  }

  if (!success) {
    cat("⚠️ Tous les essais ont échoué pour la tâche", i, "\n")
  }

  cat("\nHeure de fin :", format(Sys.time()), "\n")

  # Fermeture propre
  sink(type = "message")
  sink()
  close(con)
}

stopCluster(cl)
cat("All analyses (", nrow(param_grid), ") are finished and saved.\n")

#---------- extract non-timezone results

# path to files
rslts_folder <- "./"

# list .rds files in the folder
file_rds <- list.files(path = rslts_folder, pattern = "\\.rds$", full.names = TRUE)

# initiate a list to store all results
results <- lapply(file_rds, function(file) {
  # load rds object
  object <- readRDS(file)

  # extract data
  ml <- if (!is.null(object$model_output$ML)) object$model_output$ML else NA
  total_pars <- if (!is.null(object$total_pars)) object$total_pars else NA
  init_mu <- if (!is.null(object$init_mu)) object$init_mu else NA
  init_lambda <- if (!is.null(object$init_lambda)) object$init_lambda else NA

  # get model name
  file_name <- basename(file)
  model_name <- gsub("^secsse_|\\.rds$", "", file_name)

  # return list of results
  list(
    model = model_name,
    ll = ml,
    k = total_pars,
    init_mu = init_mu,
    init_lambda = init_lambda
  )
})

# convert into a data frame
df_results <- do.call(rbind, lapply(results, as.data.frame))
rownames(df_results) <- NULL

#---------- arrange results 

# identify parameters
trans_pattern <- "SYM"

parse_name <- function(model) {
  parts <- strsplit(model, "_")[[1]]

  sym_index <- which(parts == trans_pattern)
  dependency <- parts[1]
  inheritance <- paste(parts[2:(sym_index - 1)], collapse = "_")
  transitions <- trans_pattern
  fixed_parameter <- gsub("^fix", "", parts[sym_index + 1])
  initial_parameters <- paste(parts[(sym_index + 2):length(parts)], collapse = "_")

  return(data.frame(
    dependency = dependency,
    inheritance = inheritance,
    transitions = transitions,
    fixed_parameter = fixed_parameter,
    initial_parameters = initial_parameters
  ))
}

parsed <- do.call(rbind, lapply(df_results$model, parse_name))

# assemble and arrange the table
rslts <- cbind(parsed, df_results[, -1])

#---------- rank models

table_rslts <- rslts

# calculate AIC
table_rslts$AIC <- 2 * table_rslts$k - 2 * table_rslts$ll

# calculate  AICc
n <- length(traits)
table_rslts$AICc <- table_rslts$AIC + (2 * table_rslts$k * (table_rslts$k + 1)) / (n - table_rslts$k - 1)

# calculate delta AIC
library(dplyr)
table_rslts <- table_rslts %>% mutate(delta_AIC = AIC - min(AIC))
table_rslts <- table_rslts %>% mutate(delta_AICc = AICc - min(AICc))

# calculate AICw
table_rslts <- table_rslts %>% mutate(AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)))

# order
table_rslts <- table_rslts[order(table_rslts$AICc), ]

#---------- keep best AIC from each model

# without time zones
table_rslts_f <- table_rslts %>%
  group_by(dependency, inheritance, transitions, fixed_parameter) %>%
  slice_min(delta_AIC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(delta_AIC)

table_rslts_f_df <- as.data.frame(table_rslts_f)

#---------- extract non-timezone results

# path to files
rslts_folder <- "./"

# list .rds files in the folder
file_rds <- list.files(path = rslts_folder, pattern = "\\.rds$", full.names = TRUE)

# initiate a list to store all results
results <- lapply(file_rds, function(file) {
  # load rds object
  object <- readRDS(file)

  # extract data
  ml <- if (!is.null(object$model_output$ML)) object$model_output$ML else NA
  total_pars <- if (!is.null(object$total_pars)) object$total_pars else NA
  init_mu <- if (!is.null(object$init_mu)) object$init_mu else NA
  init_lambda <- if (!is.null(object$init_lambda)) object$init_lambda else NA

  # get model name
  file_name <- basename(file)
  model_name <- gsub("^secsse_|\\.rds$", "", file_name)

  # return list of results
  list(
    model = model_name,
    ll = ml,
    k = total_pars,
    init_mu = init_mu,
    init_lambda = init_lambda
  )
})

# convert into a data frame
df_results <- do.call(rbind, lapply(results, as.data.frame))
rownames(df_results) <- NULL

#---------- arrange results

# identify parameters
trans_pattern <- "ER"

# identify parameters
parse_name <- function(model) {
  parts <- strsplit(model, "_")[[1]]

  dependency <- parts[3]
  inheritance <- "dual_inheritance"
  transitions <- trans_pattern
  fixed_parameter <- gsub("^fix", "", parts[4])
  initial_parameters <- paste(parts[5:length(parts)], collapse = "_")

  return(data.frame(
    dependency = dependency,
    inheritance = inheritance,
    transitions = transitions,
    fixed_parameter = fixed_parameter,
    initial_parameters = initial_parameters
  ))
}

parsed <- do.call(rbind, lapply(df_results$model, parse_name))

# assemble and arrange the table
rslts_tz <- cbind(parsed, df_results[, -1])

#---------- rank models

table_rslts <- rslts_tz

# calculate AIC
table_rslts$AIC <- 2 * table_rslts$k - 2 * table_rslts$ll

# calculate  AICc
n <- length(traits)
table_rslts$AICc <- table_rslts$AIC + (2 * table_rslts$k * (table_rslts$k + 1)) / (n - table_rslts$k - 1)

# calculate delta AIC
library(dplyr)
table_rslts <- table_rslts %>% mutate(delta_AIC = AIC - min(AIC))
table_rslts <- table_rslts %>% mutate(delta_AICc = AICc - min(AICc))

# calculate AICw
table_rslts <- table_rslts %>% mutate(AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)))

# order
table_rslts <- table_rslts[order(table_rslts$AICc), ]

#---------- keep best AIC from each model

# with time zones
table_rslts_f <- table_rslts %>%
  group_by(dependency, inheritance, transitions, fixed_parameter) %>%
  mutate(dependency = recode(dependency,
    "1" = "CTD_ETD",
    "2" = "CR_ETD",
    "3" = "ETD_CTD",
    "4" = "CR_CTD",
    "5" = "ETD_CR",
    "6" = "CTD_CR"
  )) %>%
  slice_min(delta_AICc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(delta_AICc)

table_rslts_f_df <- as.data.frame(table_rslts_f)