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
