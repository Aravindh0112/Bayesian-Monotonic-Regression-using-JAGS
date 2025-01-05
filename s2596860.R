# Aravindh Sankar Ravisankar : s2596860

# This .R file compares two JAGS models using experimental data
# on the effect of a chemical on genetically engineered cells.
# Performs Bayesian inference on the models, setting up plots
# for visualization, computes DIC scores to assess and
# compare model perfomance

# Loading the required library for using JAGS
library(rjags)

# Reading the data from the vin.txt file
vin_jags_data <- read.table("vin.txt", header = TRUE)
n_experiments <- length(unique(vin_jags_data$exper))
n_doses <- length(unique(vin_jags_data$conc))
effect <- vin_jags_data$effect

# Organizing the data into a jags_list
jags_data <- list (
  'n_experiments' = n_experiments,
  'n_doses' = n_doses,
  'effect' = matrix(effect , nrow = n_doses , ncol = n_experiments)
)

# Running the data from the first JAGS model
mod <- jags.model("M1s2596860.jags", data = jags_data)
no_of_iterations <- 20000
burn_in = 2000

# Perform the JAGS sampling through coda.sampling method
jags_samples <- coda.samples(mod, variable.names = c("mu","Mj","tau","tau_o"), n.iter = no_of_iterations, burn_in = burn_in)
str(jags_samples)


# Selecting the stochastic variable columns for effective
# sample size computation (variable Mj)
selected_cols <- c("Mj[1]","Mj[2]","Mj[3]","Mj[4]","Mj[5]")

# Computing the effective sample size
eff_sample_sizes<- effectiveSize(jags_samples[,selected_cols])

# Computing the maximum and minimum effective sample size
max_eff_param <- which.max(eff_sample_sizes)
jags_samples[,max_eff_param]
min_eff_param <- which.min(eff_sample_sizes)
jags_samples[,min_eff_param]

# Plotting the graph for the maximum effective sample parameter
plot(jags_samples[,max_eff_param])

#Plotting the graph for the minimum effective sample parameter
plot(jags_samples[,min_eff_param])


# Loading the ggplot2 library
library(ggplot2)

# Fetching the effective samples
eff_samples <- jags_samples[[1]]
eff_samples

# Initializing the array to store the values of expected_effect
expected_effect <- array(NA ,dim = c(n_doses,n_experiments,no_of_iterations))
mean_effect <- c() # Initializing for storing mean of the effect values
lower_quantile <- c() # Initializing for storing the lower_quantile
upper_quantile <- c() # Initializing for storing the upper_quantile

# Executing a for loop for computing expected effect,
# lower quantile and upper quantile values.
  for (j in 1:n_experiments) {
    for (i in 1:n_doses) {
       expected_effect[i,j,] <- eff_samples[,paste0('Mj[',j,']')] * eff_samples[,paste0('mu[',i,']')]
       lower_quantile[(j-1)*9 + i] <- quantile(expected_effect[i,j,],0.025)
       upper_quantile[(j-1)*9 + i] <- quantile(expected_effect[i,j,] , 0.975)
    }
  }

# Executing a for loop for fetching the mean of the
# expected effect values computed from above.
for (j in 1:n_experiments) {
  for(i in 1:n_doses) {
    mean_effect[(j-1)*9 + i] <- mean(expected_effect[i,j,])
  }
}

# creating a dataframe for plotting the graph
plot_data <- data.frame(
  dose_index = 1:9,
  effect =  c(vin_jags_data$effect),
  mean_effect = mean_effect,
  lower_quantile = lower_quantile,
  upper_quantile = upper_quantile,
  exp_grp = rep(c("Exp1","Exp2","Exp3","Exp4","Exp5"), each = 9)
)

colors <- c("blue","green","red","cyan","orange")
colScale <- scale_colour_manual(name = "exp_grp",values = colors)

#Plotting the graph with ggplot with first plot being observed effect,
# and overlaying on it is the computed mean_expected effect alongside
# the lower quantile and upper quantile calculated (95% credible regions)
plot_1 <- ggplot(data = plot_data, aes(x = dose_index, y = effect,color = exp_grp)) +
          geom_line() +
          geom_point(aes(y = mean_effect, color = exp_grp)) +
          geom_ribbon(data = plot_data , aes(x = dose_index , ymin = lower_quantile , ymax = upper_quantile, fill = exp_grp),alpha = 0.3)+
          scale_color_manual(name = "exp_grp", values = colors)+
          scale_fill_manual(name = "exp_grp", values = colors)+
          labs(title = "Effect vs Dose Index", x = "Dose Index", y = "Effect") + colScale +
          theme_minimal()

# Displaying the plot
print(plot_1)

# Running the second updated JAGS model
mod <- jags.model("M2s2596860.jags", data = jags_data)

# Sampling on the updated JAGS model
updated_jags_samples <- coda.samples(mod, variable.names = c("mu","Mj","tau","tau_o"), n.iter = no_of_iterations, burn_in = burn_in)
str(updated_jags_samples)

# Producing a plot just like the one which we produced for 
# the earlier model

# Samples are fetched from the updated_jags_samples
updated_eff_samples <- updated_jags_samples[[1]]
updated_eff_samples

# Initialization of expected effect, mean_effect, 
# upper quantile and lower quantile
upd_expected_effect <- array(NA ,dim = c(n_doses,n_experiments,no_of_iterations))
upd_mean_effect <- c()
upd_lower_quantile <- c()
upd_upper_quantile <- c()

# Computing expected effect, lower quantile and upper quantile
for (j in 1:n_experiments) {
  for (i in 1:n_doses) {
    upd_expected_effect[i,j,] <- updated_eff_samples[,paste0('Mj[',j,']')] * updated_eff_samples[,paste0('mu[',i,',',j,']')]
    upd_lower_quantile[(j-1)*9 + i] <- quantile(upd_expected_effect[i,j,],0.025)
    upd_upper_quantile[(j-1)*9 + i] <- quantile(upd_expected_effect[i,j,] , 0.975)
  }
}
# Computing the mean of expected effect values
for (j in 1:n_experiments) {
  for(i in 1:n_doses) {
    upd_mean_effect[(j-1)*9 + i] <- mean(upd_expected_effect[i,j,])
  }
}
# Initializing the data frame for plotting the graph
updated_plot_data <- data.frame(
  dose_index = 1:9,
  effect =  c(vin_jags_data$effect),
  mean_effect = upd_mean_effect,
  lower_quantile = upd_lower_quantile,
  upper_quantile = upd_upper_quantile,
  exp_grp = rep(c("Exp1","Exp2","Exp3","Exp4","Exp5"), each = 9)
)

colors <- c("blue","green","red","cyan","orange")
colScale <- scale_colour_manual(name = "exp_grp",values = colors)

updated_plot <- ggplot(data = updated_plot_data, aes(x = dose_index, y = effect,color = exp_grp)) +
  geom_line() +
  geom_point(aes(y = mean_effect, color = exp_grp)) +
  geom_ribbon(data = updated_plot_data , aes(x = dose_index , ymin = lower_quantile , ymax = upper_quantile, fill = exp_grp),alpha = 0.3)+
  scale_color_manual(name = "exp_grp", values = colors)+
  scale_fill_manual(name = "exp_grp", values = colors)+
  labs(title = "Effect vs Dose Index", x = "Dose Index", y = "Effect") + colScale +
  theme_minimal()

print(updated_plot)

# Validating the model on the DIC score measures
model_1 <- jags.model("M1s2596860.jags",data = jags_data,n.chains=2)
model_2 <- jags.model("M2s2596860.jags",data = jags_data,n.chains=2)

dic.samples(model_1,n.iter=no_of_iterations)
dic.samples(model_2,n.iter=no_of_iterations)

# model_1 (M1s2596860.jags) is preferred as compared
# to the model_2 (M2s2596860.jags) because of the lower
# penalized deviance and mean deviance values.


