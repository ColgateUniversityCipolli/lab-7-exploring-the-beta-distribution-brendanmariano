library(tidyverse)
library(cumstats)
library(patchwork)
calculate_distribution <- function(alpha, beta)
{
  label = paste("Beta(",alpha,",", beta,")", sep = "")
  # print(label)
  q1.fig.dat <- tibble(x = seq(-0.25, 1.25, length.out=1000)) |>   # generate a grid of points
    mutate(beta.pdf = dbeta(x, alpha, beta),
           norm.pdf = dnorm(x, mean = alpha/(alpha+beta), 
           sd = sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))))  # Gaussian distribution with same mean and variance
    plot = ggplot(data= q1.fig.dat)+                                              # specify data
      geom_line(aes(x=x, y=beta.pdf, color=label)) +                 # plot beta dist
      geom_line(aes(x=x, y=norm.pdf, color="Gaussian(0.2857, 0.0255)")) +  # plot guassian dist
      geom_hline(yintercept=0)+                                            # plot x axis
      theme_bw()+                                                          # change theme
      xlab("x")+                                                           # label x axis
      ylab("Density")+                                                     # label y axis
      scale_color_manual("", values = c("orange", "grey"))+                 # change colors
      theme(legend.position = "bottom") +
      ggtitle(plot.titles[i]) 

  values = tibble(mean = numeric(1), variance = numeric(1), 
                  skewness = numeric(1), kurtosis = numeric(1)) |>
    mutate(mean = alpha/(alpha + beta)) |>
    mutate(variance = (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))) |>
    mutate(skewness = 2*(beta-alpha)*(alpha+beta+1)^(1/2)/((alpha+beta+2)*(alpha*beta)^(1/2))) |>
    mutate(kurtosis = 6*((alpha-beta)^2*(alpha+beta+1)-alpha*beta*(alpha+beta+2))
           /(alpha*beta*(alpha+beta+2)*(alpha+beta+3)))
    return (list(plot = plot, values = values))
}
alpha.vals = c(2,5,5,.5)
beta.vals = c(5,5,2,.5)
plot.titles = c("PDF of Beta(2,5)", "PDF of Beta(5,5)", 
                "PDF of Beta(5,2)", "PDF of Beta(.5,.5)")
total.data = tibble(mean = numeric(0), variance = numeric(0),
                    skewness = numeric(0), kurtosis = numeric(0))
all.plots = list()
for(i in 1:4)
{
  return.val = calculate_distribution(alpha.vals[i], beta.vals[i])
  total.data = bind_rows(total.data, return.val$values)
  all.plots[[i]] <- return.val$plot
}
# view(total.data)
(all.plots[[1]] + all.plots[[2]]) / 
(all.plots[[3]] + all.plots[[4]])

#Step 2:
beta.moment <- function(alpha, beta, k, centered)
{
  mean.funct = function(x) x^1*dbeta(x,alpha, beta)
  mean.val = integrate(mean.funct,lower = 0,upper = 1)$value
  if(k == 1)
  {
    mean.val
  }
  else if(k == 2)
  {
    variance.funct = function(x) (x-mean.val)^k*dbeta(x,alpha,beta)
    variance.val = integrate(variance.funct, 0,1)
    variance.val$value
  }
  else if(k == 3)
  {
    skewness.funct.top = function(x) (x-mean.val)^k*dbeta(x,alpha,beta)
    skewness.val.top = integrate(skewness.funct.top, 0,1)$value
    skewness.funct.bott = function(x) (x-mean.val)^(k-1)*dbeta(x,alpha,beta)
    skewness.val.bott = integrate(skewness.funct.bott, 0,1)$value
    skewness.val = skewness.val.top/(skewness.val.bott^(3/2))
    skewness.val
  }
  else if(k == 4)
  {
    kurt.funct.top = function(x) (x-mean.val)^k*dbeta(x,alpha,beta)
    kurt.val.top = integrate(kurt.funct.top, 0, 1)$value
    kurt.funct.bott = function(x) (x-mean.val)^(k-2) *dbeta(x,alpha,beta)
    kurt.val.bott = integrate(kurt.funct.bott, 0, 1)$value
    kurt.val = kurt.val.top/((kurt.val.bott)^2) -3
    kurt.val
  }
    
}
population.data = tibble()
for(i in 1:4)
{
  row = tibble() |>
  summarize(type = "Population", alpha = alpha.vals[i], beta = beta.vals[i],
            mean = beta.moment(alpha.vals[i], beta.vals[i],1,F),
            variance = beta.moment(alpha.vals[i], beta.vals[i], 2, T),
            skewness = beta.moment(alpha.vals[i], beta.vals[i], 3,T),
            excess.kurtosis = beta.moment(alpha.vals[i], beta.vals[i], 4, T))
  population.data = bind_rows(population.data, row)
}
# print(beta.moment(2,5,3, F))
# print(calculate_distribution(2,5)$values)
# view(population.data)
#Creates histograms

#################################################################################
#Making Data Summaries
#################################################################################
set.seed(7272)
sample.size = 500
plot.list = list() #For storing the plots
summaries.list = list() #For storing summaries
plot.titles = c("Sample of Beta(2,5)", "Sample of Beta(5,5)", 
                "Sample of Beta(5,2)", "Sample of Beta(.5,.5)")
#Should call calculate distribution
for(i in 1:4)
{
  total.distrib = tibble(type = character(), alpha = numeric(), beta = numeric(), mean = numeric(), variance = numeric(), 
                       skewness = numeric(), excess.kurtosis = numeric())
  beta.sample = tibble(x = rbeta(n = sample.size, shape1 = alpha.vals[i], shape2 = beta.vals[i])) 
  #Stores actual probability density function
  population.fig = tibble(x = seq(0, 1, length.out = 1000)) |>
    mutate(beta.pdf = dbeta(x, shape1 = alpha.vals[i], shape2 = beta.vals[i]))
  #Creates plots
  plot.list[[i]] = ggplot() + 
    geom_histogram(data = beta.sample, aes(x = x, y = after_stat(density)), fill = "royalblue", 
                   breaks = seq(0,1,.08)) +
    geom_density(data = beta.sample, aes(x = x, y = after_stat(density), color = "Sample")) + 
    geom_line(data = population.fig, aes(x = x, y = beta.pdf, color = "Population")) + 
    ylab("Density") +
    xlab("x") + 
    ggtitle(plot.titles[i]) + 
    theme(plot.title = element_text(size = 10), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)) +
    labs(color = "Color")
    
  #Making data summary table
  summary = beta.sample %>%
    summarize(type = "Sample", mean = mean(x), variance = var(x),
              skewness = skewness(x), excess.kurtosis = kurtosis(x),
              alpha = alpha.vals[i], beta = beta.vals[i]) 
    total.distrib = bind_rows(total.distrib, summary)
    total.distrib = bind_rows(total.distrib, population.data |> slice(i))
    summaries.list[[i]] = total.distrib
}

(plot.list[[1]] + plot.list[[2]]) /
(plot.list[[3]] + plot.list[[4]])

################################################################################
#Step 4:
################################################################################
###################
#Adds cum mean
###################
mean.plot = ggplot() +
  ylab("Cumulative Mean") + 
  xlab("Observations") + 
  geom_hline(yintercept = beta.moment(2,5, 1, F)) +
  ggtitle("Cumulative Mean with 49 samples") + 
  theme(plot.title = element_text(size = 10))
###################
#Adds cum variance
###################
var.plot = ggplot() +
  ylab("Cumulative Variance") +
  xlab("Observations") + 
  geom_hline(yintercept = beta.moment(2,5,2,T)) + 
  ggtitle("Cumulative Variance with 49 samples") + 
  theme(plot.title = element_text(size = 10))
###################
#Adds cum skewness
###################
skew.plot = ggplot() +
  ylab("Cumulative Skewness") +
  xlab("Observations") + 
  geom_hline(yintercept = beta.moment(2,5,3,T)) +
  ggtitle("Cumulative Skewness with 49 samples") + 
  theme(plot.title = element_text(size = 10))
###################
#Adds cum kurtosis
###################
kurt.plot = ggplot() +
  ylab("Cumulative Kurtosis") +
  xlab("Observations") + 
  geom_hline(yintercept = beta.moment(2,5,4,T)) + 
  ggtitle("Cumulative Kurtosis with 49 samples") + 
  theme(plot.title = element_text(size = 10))

#Adds lines to each graph
for(i in 2:50)
{
observations = 1:500
  set.seed(7272 + i)
  #Hardcoded alpha and beta values (can change if needed)
  random.data = tibble(x = rbeta(n = 500, shape1 = 2, shape2 = 5)) |>   
    mutate(cummean.1 = cummean(x))
  mean.plot = mean.plot +
    geom_line(data = random.data, aes(x=observations, y = cummean.1), color = i)
  # print(random.data)
  #Variance
  random.data = tibble(x = rbeta(n = 500, shape1 = 2, shape2 = 5)) |>   
    mutate(cumvar.1 = cumvar(x))
  var.plot = var.plot +
    geom_line(data = random.data, aes(x=observations, y = cumvar.1), color = i)
  # print(random.data)
  #Skewness
  random.data = tibble(x = rbeta(n = 500, shape1 = 2, shape2 = 5)) |>   
    mutate(cumskew.1 = cumskew(x))
  skew.plot = skew.plot +
    geom_line(data = random.data, aes(x=observations, y = cumskew.1), color = i)
  #Kurtosis
  random.data = tibble(x = rbeta(n = 500, shape1 = 2, shape2 = 5)) |>   
    mutate(cumkurt.1 = cumkurt(x))
  kurt.plot = kurt.plot +
    geom_line(data = random.data, aes(x=observations, y = cumkurt.1-3), color = i)
  
}
#Utilize Patchwork for displaying plots
(mean.plot + skew.plot)/(var.plot + kurt.plot)
################################################################################
#Step 5:Modeling the variation
################################################################################
cs.total.sample = tibble(sample.mean = numeric(), sample.var = numeric(), 
                                sample.skew = numeric(), sample.kurt = numeric())
for(i in 1:1000)
{
  #Generates random sample
  random.data = tibble(x = rbeta(n= 500, shape1 = 2, shape2 = 5))
  #Stores current sample
  cs.curr.sample = tibble() |>
    summarize(sample.mean = mean(pull(random.data,x)), sample.var = var(pull(random.data,x)),
           sample.skew = skewness(pull(random.data,x)), sample.kurt = kurtosis(pull(random.data,x))-3)
  cs.total.sample = cs.total.sample |>
    bind_rows(cs.curr.sample)
}
######################
#Creating plots
######################
#Mean plot
mean.plot = ggplot(data = cs.total.sample, aes(x = sample.mean, y = after_stat(density))) +
  geom_histogram(fill = "navy", breaks = seq(.26,.31,.003)) +
  geom_density(color = "gray") + 
  theme_bw() + 
  xlab("Mean") + 
  ylab("Density") + 
  ggtitle("PDF of Mean")

#Variance plot
variance.plot = ggplot(data = cs.total.sample, aes(x = sample.var, y = after_stat(density))) +
  geom_histogram(fill = "navy",breaks = seq(.02,.03,.0008)) +
  geom_density(color = "gray") + 
  theme_bw() + 
  xlab("Variance") + 
  ylab("Density") + 
  ggtitle("PDF of Variance")
#Skewness plot
skewness.plot = ggplot(data = cs.total.sample, aes(x = sample.skew, y = after_stat(density))) +
  geom_histogram(fill = "navy", breaks = seq(.3,.9,.04)) +
  geom_density(color = "gray") + 
  theme_bw() + 
  xlab("Skewness") + 
  ylab("Density") + 
  ggtitle("PDF of Skewness")
#Kurtosis plot
kurt.plot = ggplot(data = cs.total.sample, aes(x = sample.kurt, y = after_stat(density))) +
  geom_histogram(fill = "navy", breaks = seq(-.6,.6,.095)) +
  geom_density(color = "gray") + 
  theme_bw() + 
  xlab("Kurtosis") + 
  ylab("Density") + 
  ggtitle("PDF of Kurtosis")

(mean.plot + variance.plot) / (skewness.plot + kurt.plot)