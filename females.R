# read in the ones we want
expos_f <- read.csv("data/expos_f_old.csv")
deaths_f <- read.csv("data/deaths_f_old.csv")

dim(deaths_f)

N <- dim(expos_f)[1]
n_years <- dim(expos_f)[2]

age  <- 1:N
years <- 1:n_years
# standardise age and years to aid convergence
age_standardised <- (age - mean(age)) / sd(age)
years_standardised <- (years -mean(years))/ sd(years)

# created a named list with all the elements mentioned in the stan model

stan_input_data <-list(N=N,
                       n_years=n_years,
                       age_standardised=age_standardised,
                       years_standardised=years_standardised,
                       expos=expos_f, 
                       # note the name in the list (expos) is different from the variable name expos_m
                       # The list name must match those declared in the stan model data section.
                       deaths=round(deaths_f))



stan_fit <- stan(file = "stan/old_age_save.stan", data = stan_input_data, iter = 1000,
                 chains=3, cores=1)

print(stan_fit)

traceplot(stan_fit, nrow=3)

stan_rhat(stan_fit)
stan_ac(stan_fit)


mu_samples <- as.matrix(stan_fit, par="mu")
dim(mu_samples) 
n_samples <- dim(mu_samples)[1]
# arrange our samples so that we have a 3-dimensional array
#  the dimensions are                 samples, age, year
mu_array <- array(mu_samples, dim=c(n_samples,   N, n_years))

# extract the mean ---------
# apply the function 'mean' to the array `mu_array`,
# averaging over dimension 1, leaving dimensions 2 and 3 intact.
mu_mean <- apply(mu_array, MARGIN=c(2,3), FUN=mean)

ages <- as.numeric(row.names(expos_f))
year_index <- 1

par(bty = 'l')
plot(x=ages, y=mu_mean[,year_index], type="l", ylab="Rate", xlab="Ages", ylim=c(0,2.5))
# type="l" means draw a continuous line, not a point.

# we can compute the real rates and plot them also. 
real_rates <- deaths_f/expos_f
points(x=ages,y=real_rates[,year_index], pch=19, col="red")



stan_fit <- stan(file = "stan/old_age_gen_quan.stan", data = stan_input_data, iter = 1000,
                 chains=3, cores=1)

print(stan_fit)


sim_deaths <- as.matrix(stan_fit, par="simulated_deaths")

# arrange our samples so that we have a 3-dimensional array
#  the dimensions are                 samples, age, year
sim_deaths_array <- array(sim_deaths, dim=c(n_samples,   N, n_years))

expos_matrix <- data.matrix(expos_f)
# divide every sample array by observed exposures to get implied rates
sim_rates_array <- sapply(1:n_samples, 
                          function(i) sim_deaths_array[i,,] / expos_matrix, 
                          simplify = "array")

# permute this matrix to get back to the right dimension order
sim_rates_array <- aperm(sim_rates_array, c(3,1,2))


# compute quantiles  ---------
# apply the function 'quantiles' to the array,
# passing additional argument probs to quantile function
sim_rates_q <- apply(sim_rates_array, MARGIN=c(2,3), FUN=quantile, 
                     probs=c(0.05,0.5,0.95), na.rm=T)


ages <- as.numeric(row.names(expos_f))
year_index <- 53

par(bty = 'l')
plot(x=ages, y=sim_rates_q[2,,year_index], type="l", ylab="Rate", xlab="Ages", ylim=c(0,2.5))
points(x=ages, y=sim_rates_q[1,,year_index], type="l", lty=3) # lty=3 gives a dotted line
points(x=ages, y=sim_rates_q[3,,year_index], type="l", lty=3) # lty=3 gives a dotted line


real_rates <- deaths_f/expos_f
points(x=ages,y=real_rates[,year_index], pch=19, col="red")
title(paste("Predicted Mortality Rates for ",1960 + year_index),
      sub="90% interval")
