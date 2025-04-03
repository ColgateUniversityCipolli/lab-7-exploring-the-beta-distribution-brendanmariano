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
(mom.alpha = measure.of.moments$x[1])
(mom.beta.val = measure.of.moments$x[2])
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
mle.alpha = mle.val$par[1]
mle.beta = mle.val$par[2]
#################################
#Plotting distributions
#################################
#Histogram of data
#MOM and MLE Distribution
mom.seq = tibble(x = seq(0,.022,.0001)) |>
  mutate(mom.dat = dbeta(x, shape1 = mom.alpha, shape2 = mom.beta.val)) |>
  mutate(mle.dat = dbeta(x, shape1 = mle.alpha, shape2 = mle.beta))
mom.plot = ggplot() +
  geom_histogram(aes(x = data1$`2022`, y = after_stat(density)), 
                 breaks = seq(0,.022, .0017), fill = "black") + 
  geom_line(data = mom.seq, aes(x = x, y = mom.dat, color = "MOM")) +
  geom_line(data = mom.seq, aes(x = x, y = mle.dat, color = "MLE")) +
  xlab("Proportion of Deaths relative to Entire Population") + 
  ylab("Density") +
  ggtitle("Population vs Estimation of Death Distribution") + 
  theme_bw()
mom.plot

#################################
#Estimating beta and alpha values
#################################
mle.beta.list = numeric(0)
mle.alpha.list = numeric(0)
mom.beta.list = numeric(0)
mom.alpha.list = numeric(0)
#Need to change to 1000
for(i in 1:1000)
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
##################
#Creating Plots
#################
#Plots mle alpha
mle.alpha.plot = ggplot() +
  geom_density(aes(x = mle.alpha.list), fill = "red") + 
  xlab("Alpha (MLE)") + 
  ylab("Density") + 
  ggtitle("MLE Alpha Distribution") + 
  theme_bw()
#Plots mle beta
mle.beta.plot = ggplot() +
  geom_density(aes(x = mle.beta.list), fill = "black") + 
  xlab("Beta (MLE)") + 
  ylab("Density") +
  ggtitle("MLE Beta Distribution") +
  theme_bw()
#Plots mom alpha
mom.alpha.plot = ggplot() +
  geom_density(aes(x = mom.alpha.list),fill = "blue") + 
  xlab("Alpha (MOM)") + 
  ylab("Density") +
  ggtitle("MOM Alpha Distribution") +
  theme_bw()
#Plots mom beta
mom.beta.plot = ggplot() +
  geom_density(aes(x = mom.beta.list), fill = "green") + 
  xlab("Beta (MLE)") + 
  ylab("Density") +
  ggtitle("MOM Beta Distribution") +
  theme_bw()
(mle.alpha.plot + mle.beta.plot) /
(mom.alpha.plot + mom.beta.plot)

#######################
#Calculating indicators
#######################
#Could make below more efficient by combing vectors and using group by
#Review
indicator.table = tibble() |>
  summarize(element = "Alpha (MLE)", 
          bias = mean(mle.alpha.list) - 8, 
          precision = 1/var(mle.alpha.list), 
          mse = var(mle.alpha.list) + bias^2) %>%
  bind_rows(summarize(.,element = "Beta (MLE)", 
                   bias = mean(mle.beta.list) - 950, 
                   precision = 1/var(mle.beta.list), 
                   mse = var(mle.beta.list) + bias^2))  %>%
  bind_rows(summarize(., element = "Alpha (MOM)",
                   bias = mean(mom.alpha.list) - 8,
                   precision = 1/var(mom.alpha.list),
                   mse = var(mom.alpha.list) + bias^2)) %>%
  bind_rows(summarize(., element = "Beta (MLE)",
                   bias = mean(mom.beta.list) - 950,
                   precision = 1/var(mom.beta.list),
                   mse = var(mom.beta.list) + bias^2))

indicator.table


