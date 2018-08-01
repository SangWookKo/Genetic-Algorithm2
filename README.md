Genetic-Algorithm2

## Refernece https://cran.r-project.org/web/packages/GA/vignettes/GA.html#introduction by Luca Scrucca


## Introduction

Genetic algorithms (GAs) are stochastic search algorithms inspired by the basic principles of biological evolution and natural selection. GAs simulate the evolution of living organisms, where the fittest individuals dominate over the weaker ones, by mimicking the biological mechanisms of evolution, such as selection, crossover and mutation.

* The R package GA provides a collection of general purpose functions for optimization using genetic algorithms. 
* The package includes a flexible set of tools for implementing genetic algorithms search in both the continous and discrete case, whether constrained or not.
* Users can easily define their own objective function depending on the problem at hand.
* Several genetic operators are available and can be combined to explore the best settings for the current task.
* Furthermore, users can define new genetic operators and easily evaluate their performance.
* Local search using general-purpose optimisation algorithms can be applied stochastically to exploit interesting regions.
* GAs can be run sequentially or in parallel, using an explicit master-slave paralleisation or a coarse-grain islands approach.

This document gives a quick tour of GA(version 3.1.1) functionalities. It was originally written in R Markdown, using the <U>knitr</U> package for production. Further details are provided in the papers Scrucca(2013) and Scrucca(2017). 
See also help(package="GA") for a list of available functions and methods.

``` 
library(GA)

Type 'citiation("GA")' for citing this R package in publications.
```

## Function optimisation in one dimension

Consider the function *f(x) = (x^2 + x) cos(x)* defined over the range -10 <= *x* <= 10:

```
f <- function(x) (x^2 + x) * cos(x)
lbound <- -10; ubound <- 10
curve(f, from = lbound, to ubound, n = 1000)
```

```
  GA <- ga(type = "real-valued", fitness = f, lower = c(th = lbound), upper = ubound)
  summary(GA)
```

```
  curve(f, from = lbound, to = ubound, n = 1000)
  points(GA@solution, GA@fitnessValue, col = 2, pch = 19)
```

## Function Optimisation in two dimensions

Consider the Rastrigin function, a non-convex function often used as a test problem for optimisation algorithms because it is a difficult problem due to its large number of local minima. 

In two dimensions it is defined as 

*f(x1,x2) = 20 + x1^2 + x2^2 - 10(cos(2 pie x1) + cos(2 pie x2)), with -5.12 <= xi <= 5.12 for i = 1,2. 
It has a global minimum at (0,0) where f(0,0) = 0.

```
Rastrigin <- function(x1, x2)
{
  20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}


x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
f <- outer(x1, x2, Rastrigin)
persp3D(x1, x2, f, theta = 50, phi = 20, color.palette = bl2gr.colors)
```

```
  filled.contour(x1, x2, f, color.palette = bl2gr.colors)
```

A GA minimisation search is obtained as follows( note the minus sign used in the definition of the local fitness function):

``` 
    GA <- ga(type = "real-valued",
             fitness = function(x) - Rastrigin(x[1], x[2]),
             lower = c(-5.12, -5.12) , upper = c(5.12,5.12),
             popSize = 50, maxiter = 1000, run = 100)
    summary(GA)         
```

```
filled.contour(x1, x2, f, color.palette = bl2gr.colors, plot.axes = { axis(1); axis(2); points(GA@solution[,1], GA@solution[,2], pch = 3, cex = 2, col = "white", lwd= 2)} )
```

The GA search process can be visualised by defining a monitoring function as follows:

``` 
  monitor <- function(obj)
  {
    contour(x1, x2, f, drawlabels = FALSE, col = grey(0.5))
    title(paste("iteration = ", obj@iter), font.main = 1)
    points(obj@population, pch = 20, col = 2)
    Sys.sleep(0.2)  
  }
    ```

  
  
  
