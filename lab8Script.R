library(tidyverse)
library(stringr)
library(nleqslv)
library(patchwork)
data1 = read_csv("wb_API/Metadata3.csv") |>
  select(c("Country Name", "Country Code", "Indicator Name", "Indicator Code", "2022")) |>
  mutate(`2022` = `2022`/1000)  |>
  filter(!is.na(`2022`))

##################
#Method of Moments
##################
mom.beta = function(data, par)
{
  alpha = par[1]
  beta = par[2]
  mom.mean = mean(data)
  mom.mean.squared = mean(data^2)
  eq1 = alpha/(alpha + beta)
  eq2 = ((alpha + 1)*alpha)/((alpha + beta + 1)*(alpha + beta))
  
  return(c(eq1 - mom.mean, eq2 - mom.mean.squared))
}

(measure.of.moments = nleqslv(x = c(1,1),
                             fn = mom.beta,
                             data = data1$`2022`))
(alpha = measure.of.moments$x[1])
(beta = measure.of.moments$x[2])

#####################
#Method of Likelihood
#####################
llmle.beta = function(data, par, neg)
{
  alpha = par[1]
  beta = par[2]
  ll <- sum(log(dbeta(x=data, shape1=alpha, shape2=beta)), na.rm = TRUE) #Why are we doing the sum?
  #Ask prof. and look at notes in recall section
  return(ifelse(neg, -ll, ll))
}
mle.val <- optim(fn = llmle.beta,
              par = c(1,1),
              data = data1$`2022`,
              neg=T)
alpha = mle.val$par[1]
beta = mle.val$par[2]
#################################
#Estimating beta and alpha values
#################################
mle.beta.list = numeric(0)
mle.alpha.list = numeric(0)
mom.beta.list = numeric(0)
mom.alpha.list = numeric(0)
for(i in 1:50)
{
  set.seed(7272 + i)
  dat = rbeta(266, shape1 = 8, shape2 = 950)
  #Appends alpha mle
  mle.alpha.list = append(mle.alpha.list, optim(fn = llmle.beta,
                    par = c(1,1),
                    data = dat,
                    neg=T)$par[1])
  #Appends beta mle
  mle.beta.list = append(mle.beta.list, optim(fn = llmle.beta,
                               par = c(1,1),
                               data = dat,
                               neg=T)$par[2])
  #Appends alpha mom
  mom.alpha.list = append(mom.alpha.list, nleqslv(fn = mom.beta, 
                                 x = c(1,1),
                                 data = dat)$x[1])
  #Appends beta mom
  mom.beta.list = append(mom.beta.list, nleqslv(fn = mom.beta, 
                                x = c(1,1),
                                data = dat)$x[2])
}
view(mom.beta.list)
##################
#Creating Plots
#################
#NOTE: MIGHT WANT TO REMOVE HISTOGRAMS
#Plots mle alpha
mle.alpha.plot = ggplot() +
  geom_histogram(aes(x = mle.alpha.list, y = after_stat(density))) + 
  geom_density(aes(x = mle.alpha.list)) + 
  xlab("Alpha (MLE)") + 
  ylab("Density")
#Plots mle beta
mle.beta.plot = ggplot() +
  geom_histogram(aes(x = mle.beta.list, y = after_stat(density))) + 
  geom_density(aes(x = mle.beta.list)) + 
  xlab("Beta (MLE)") + 
  ylab("Density")
#Plots mom alpha
mom.alpha.plot = ggplot() +
  geom_histogram(aes(x = mom.alpha.list, y = after_stat(density))) + 
  geom_density(aes(x = mle.alpha.list)) + 
  xlab("Alpha (MOM)") + 
  ylab("Density")
#Plots mom beta
mom.beta.plot = ggplot() +
  geom_histogram(aes(x = mom.beta.list, y = after_stat(density))) + 
  geom_density(aes(x = mom.beta.list)) + 
  xlab("Beta (MLE)") + 
  ylab("Density")
(mle.alpha.plot + mle.beta.plot) /
(mom.alpha.plot + mom.beta.plot)

#######################
#Calculating indicators
#######################

indicator.table = indicator.table |>
  summary()


