library(tidyverse)
library(stringr)
library(nleqslv)

data1 = read_csv("wb_API/Metadata3.csv") |>
  select(c("Country Name", "Country Code", "Indicator Name", "Indicator Code", "2022")) |>
  mutate(`2022` = `2022`/1000)  |>
  filter(!is.na(`2022`))

view(data1)
#################
#Method of Moments
#################
view(pull(data,`2022`))
mom.mean = mean(pull(data1,`2022`))
mom.mean.squared = mean(pull(data1,`2022`)^2)
mom.beta = function(data, par)
{
  alpha = par[1]
  beta = par[2]
    
  eq1 = alpha/(alpha + beta)
  eq2 = ((alpha + 1)*alpha)/((alpha + beta + 1)*(alpha + beta))
  
  return(c(eq1 - mom.mean, eq2 - mom.mean.squared))
}

(measure.of.moments = nleqslv(x = c(1,1),
                             fn = mom.beta,
                             data = data$`2022`))
(alpha = measure.of.moments$x[1])
(beta = measure.of.moments$x[2])

#################
#Method of Likelihood
#################
llmle.beta = function(data, par, neg)
{
  alpha = par[1]
  beta = par[2]
  
  ll <- sum(log(dbeta(x=data, shape1=alpha, shape2=beta)), na.rm = TRUE) #Why are we doing the sum?
  
  return(ifelse(neg, -ll, ll))
}

(mle.val <- optim(fn = llmle.beta,
              par = c(1,1),
              data = data$`2022`,
              neg=T))
(alpha = mle.val$par[1])
(beta = mle.val$par[2])
  