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
    penalty1 <- max(c1(x), 0)*pen # penalisation for 1st inequality constraint    penalty2 <- max(c2(x), 0)*pen # penalisation for 2nd inequality constraint
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
`      col = col[2], add = TRUE)
contour(x1, x2, matrix(apply(x12, 1, f), ngrid, ngrid), 
        nlevels = 21, add = TRUE)
points(GA@solution[1], GA@solution[2], col = "dodgerblue3", pch = 3)  # GA solution
```
  
  
## Hybrid GAs
Hybrid Genetic Algorithms (HGAs) incorporate efficient local search algorithms into GAs. In case of real-valued optimisation problems, the GA package provides a simple way to start local searches from GA solutions after a certain number of iterations, so that, once a promising region is identified, the convergence to the global optimum can be speed up.

`The use of HGAs is controlled by the optional argument optim = TRUE (by default is set to FALSE). Local searches are executed using the base R function optim(), which makes available general-purpose optimisation methods, such as Nelder–Mead, quasi-Newton with and without box constraints, and conjugate-gradient algorithms. The local search method to be used and other parameters are controlled with the optional argument optimArgs, which must be a list with the following structure and defaults:

```
optimArgs = list(method = "L-BFGS-B",
                 poptim = 0.05,
                 pressel = 0.5,
                 control = list(fnscale = -1, maxit = 100))

GA <- ga(type = "real-valued",
         fitness = function(x) -Rastrigin(x[1], x[2]),
         lower = c(-5.12,-5.12), upper = c(5.12, 5.12),
         popSize = 50, maxiter = 1000, run = 100,
         optim = TRUE)
summary(GA)

plot(GA)
```
For more details see help(ga).

Consider again the two-dimensional Rastrigin function defined previously. A HGA search is obtained as follows:

```
GA <- ga(type = "real-valued",

         fitness = function(x) -Rastrigin(x[1], x[2]),
         lower = c(-5.12,-5.12), upper = c(5.12, 5.12),
         popSize = 50, maxiter = 1000, run = 100,
         optim = TRUE)

summary(GA)

plot(GA)
```
Note the improved solution obtained.

## Parallel computing 

By default searches performed using the GA package occour sequentially. In some cases, particularly when the evaluation of the fitness function is time consuming, parallelisation of the search algorithm may be able to speedup computing time. Starting with version 2.0, the GA package provides facilities for implementing parallelisation of genetic algorithms.
Parallel computing with GA requires the following packages to be installed: parallel (available in base R), doParallel, foreach, and iterators.
To use parallel computing with the GA package on a single machine with multiple cores is simple as manipulating the optional argument parallel in the ga() function call.
The argument parallel can be a logical argument specifying if parallel computing should be used (TRUE) or not (FALSE, default) for evaluating the fitness function. This argument could also be used to specify the number of cores to employ; by default, this is taken from detectCores() function in parallel package.
Two types of parallel functionality are implemented depending on system OS: on Windows only snow type functionality is available, while on POSIX operating systems, such as Unix, GNU/Linux, and Mac OSX, both snow and multicore (default) functionalities are available. In the latter case a string can be used to specify which parallelisation method should be used.
In all cases described above, at the end of GA iterations the cluster is automatically stopped by shutting down the workers.
Consider the following simple example where a pause statement is introduced to simulate an expensive fitness function.

```
fitness <- function(x)
{
  Sys.sleep(0.01)
  x*runif(1)
} 
install.packages("parallel")
library(parallel)
install.packages("doParallel")
library(doParallel)
install.packages("rbenchmark")
library(rbenchmark)

out <- benchmark(GA1 = ga(type = "real-valued",
                          fitness = fitness, lower = 0, upper = 1,
                          popSize = 50, maxiter = 100, monitor = FALSE,
                          seed = 12345),
                 
                 GA2 = ga(type = "real-valued", 
                          fitness = fitness, lower = 0, upper = 1,
                          popSize = 50, maxiter = 100, monitor = FALSE,
                          seed = 12345, parallel = TRUE),
                 GA3 = ga(type = "real-valued", 
                          fitness = fitness, lower = 0, upper = 1,
                          popSize = 50, maxiter = 100, monitor = FALSE,
                          seed = 12345, parallel = 2),
                 GA4 = ga(type = "real-valued", 
                          fitness = fitness, lower = 0, upper = 1,
                          popSize = 50, maxiter = 100, monitor = FALSE,
                          seed = 12345, parallel = "snow"),          


                 columns = c("test", "replications", "elapsed", "relative"),
                 order = "test",
                 replications = 10)
                          
out$average <- with(out, average <- elapsed/replications)
out[,c(1:3,5,4)]

```
If a cluster of multiple machines is available, ga() can be executed in parallel using all, or a subset of, the cores available to the machines belonging to the cluster. However, this option requires more work from the user, who needs to set up and register a parallel back end.
For instance, suppose that we want to create a cluster of two computers having IP addresses 141.250.100.1 and 141.250.105.3, respectively. For each computer we require 8 cores, so we aim at having a cluster of 16 cores evenly distributed on the two machines. Note that comunication between the master worker and the cluster nodes is done via SSH, so you should configure ssh to use password-less login. For more details see McCallum and Weston (2011, Chapter 2).

```
library(doParallel)
workers <- rep(c("141.250.100.1", "141.250.105.3"), each = 8)
cl <- makeCluster(workers, type = "PSOCK")
registerDoParallel(cl)

```

## Island evolution

GAs can be designed to evolve using an Island evolution approach. Here the population is partitioned in a set of subpopulations (islands) in which isolated GAs are executed on separated processor runs. Occasionally, some individuals from an island migrate to another island, thus allowing subpopulations to share genetic material
This approach is implemented in the gaisl() function, which has the same input arguments as the ga() function, with the addition of the following argument:

* numIslands : an integer value specifying the number of islands to use (by default is set to 4)
* migrationRate : a value in the range (0,1) which gives the proportion of individuals that undergo migration between islands in every exchange (by default equal to 0.10)
* migrationInterval : an integer value specifying the number of iterations at which exchange of individuals takes place (by default set at 10).

Parallel computing is used by default in the Island evolution approach. Hybridisation by local search is also available as discussed previously.
As an example, consider again the two-dimensional Rastrigin function. An Island GA search is obtained as follows:

```
GA <- gaisl(type = "real-valued",
            fitness = function(x) -Rastrigin(x[1], x[2]),
            lower = c(-5.12,-5.12), upper = c(5.12,5.12),
            popSize = 100,
            maxiter = 1000, run = 100,
            numIslands = 4,
            migrationRate = 0.2,
            migrationInterval = 50)

summary(GA)
plot(GA, log= "x") 
```
## Memoization 

In certain circumstances, particularly with binary GAs, memoization can be used to speed up calculations by using cached results. This is easily obtained using the memoise package

```
install.packages("memoise")
library(memoise)
install.packages("UsingR")
library(UsingR)

data(fat)
summary(fat)
str(fat)

mod <- lm(body.fat.siri ~ age + weight + height + neck + chest + abdomen + hip + thigh + knee+ ankle + bicep + forearm + wrist, data = fat)

head(mod)
summary(mod)

x <- model.matrix(mod)[,-1]
y <- model.response(mod$model)

fitness <- function(string){
  
mod <- lm(y ~ x[,string ==1])
-BIC(mod)
}

library(memoise)
mfitness <- memoise(fitness)

is.memoised(fitness)

is.memoised(mfitness)

library(rbenchmark)
tab <- benchmark(GA1 = ga("binary", fitness = fitness, nBits = ncol(x),
                          popSize = 100, maxiter = 100, seed = 1, monitor = FALSE),
                 GA2 = ga("binary", fitness = mfitness, nBits = ncol(x), 
                          popSize = 100, maxiter = 100, seed = 1, monitor = FALSE),
                 columns = c("test", "replications", "elapsed", "relative"),
                 replications = 10)
tab$average <- with(tab, elapsed/replications)
tab

# To clear cache use

forget(mfitness)
```
