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

GA <- ga(type = "real-valued",
         fitness = function(x) -Rastrigin(x[1],x[2]),
         lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
         popSize = 50, maxiter = 100,
         monitor = monitor)
```
## Setting some members of the initial population

The *suggestions* argument to *ga()* function call can be used to provide a matrix of solutions to be included in the initial population.

For example, consider the optimisation of the Rastrigin function introduced above:

```
  suggestedSol <- matrix(c(0.2,1.5,-1.5,0.5), nrow = 2, ncol = 2, byrow = TRUE)
  
  GA1 <- ga(type = "real-valued",
            fitness = function(x) -Rastrigin(x[1],x[2]),
            lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
            suggestions = suggestedSol,
            popSize = 50, maxiter = 1)
 head(GA1@population)
  
```
As it can be seen, the first two solutions considered are those provided, whereas the rest is filled randomly as usual.
A full search can be obtained as follows:

```
  GA <- ga(type = "real-valued",
          fitness = function(x) -Rastrigin(x[1],x[2]),
          lower = c(-5.12,-5.12), upper = c(5.12,5.12),
          suggestions = suggestedSol,
          popSize = 50, maxiter = 100)
 summary(GA)
```

## Constrained optimisation

This example shows how to minimize an objective function subject to nonlinear inequality constraints and bounds using GAs.
Source : http://www.mathworks.it/it/help/gads/examples/constrained-minimization-using-the-genetic-algorithm.html

We want to minimize a simple function of two variables x1 and x2
    
          min f(x) = 100(x1^2-x2)^2 + (1-x1)^2;
           x 
subject to the following nonlinear inequality constraints and bounds:

* x1x2 + x1 - x2 + 1.5 <= 0 (inequality constraint),
* 10 - x1x2 <= 0 (inequality constraint),
* 0 <= x1 <= 1 (bounds), and
* 0 <= x2 <= 13 (bounds).

The above fitness function is known as "cam" as described in L.C.W. Dixon and G.P. Szego(eds.), __Towards Global Optimisation 2__, North-Holland, Amsterdam, 1978.

```
    f <- function(x) 
    { 100 * (x[1]^2 - x[2])^2 + (1 - x[1])^2}
    
    c1 <- function(x) 
    {x[1]*x[2] + x[1] - x[2] + 1.5}
    c2 <- function(x) 
    { 10 - x[1]*x[2]}
```

Plot the function and the feasible regions(coloured areas):

```
ngrid <- 250
x1 <- seq(0,1, length = ngrid)
x2 <- seq(0, 13, length = ngrid)
x12 <- expand.grid(x1,x2)

col <- adjustcolor(bl12gr.colors(4)[2:3], alpha = 0.2)
plot(x1, x2, type ="n", xaxs = "i", yaxs = "i")
image(x1, x2, matrix(ifelse(apply(x12,1,c1) <= 0,0,NA), ngrid, ngrid), col = col[1], add= TRUE)
image(x1,x2, matrix(ifelse(apply(x12,1,c2) <= 0,0,NA), ngrid, ngrid), col = col[2], add = TRUE)
contour(x1, x2, matrix(apply(x12,1,f), ngrid,ngrid), nlevels = 21, add = TRUE)
```

A GA solution can be obtained by defining a penalised fitness function:

``` 
  fitness <- function(x) 
  { 
    f <- -f(x) # we need to maximise -f(x)
    pen <- sqrt(.Machine$double.xmax) # penalty term
    penalty1 <- max(c1(x), 0)*pen # penalisation for 1st inequality constraint
    penalty2 <- max(c2(x), 0)*pen # penalisation for 2nd inequality constraint
    f - penalty1 - penalty2   # fitness function value
    
  }
```

Then 

```
GA <- ga("real-valued", fitness = fitness, lower = c(0,0), upper= c(1,13),
          # selection = GA:::gareal_lsSelection_R,
          maxiter = 1000, run = 200, seed = 123)
          
summary(GA)          

fitness(GA@solution)
f(GA@solution)
c1(GA@solution)
c2(GA@solution)

```

A graph showing the solution found is obtained as:

```
 plot(x1, x2, type = "n", xaxs = "i", yaxs = "i")
image(x1, x2, matrix(ifelse(apply(x12, 1, c1) <= 0, 0, NA), ngrid, ngrid), 
      col = col[1], add = TRUE)
image(x1, x2, matrix(ifelse(apply(x12, 1, c2) <= 0, 0, NA), ngrid, ngrid), 
      col = col[2], add = TRUE)
contour(x1, x2, matrix(apply(x12, 1, f), ngrid, ngrid), 
        nlevels = 21, add = TRUE)
points(GA@solution[1], GA@solution[2], col = "dodgerblue3", pch = 3)  # GA solution
```
  
  
