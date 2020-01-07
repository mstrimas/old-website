---
layout: post
title: "Field Guide to ILP Solvers in R for Conservation Prioritization"
published: true
excerpt: >
  A field guide to all the open-source integer linear programming solvers that
  have R packages. The focus is on finding an open-source alternative to Gurobi
  for conservation prioritization and proteced area design.
category: prioritization
tags: r gurobi optimization marxan
---





In this post I'll compare alternative integer linear programming (ILP) solvers for conservation planning. The goal is to develop a tool to solve the Marxan reserve design problem using ILP rather than simulated annealing. Unlike simulated annealing, ILP can find the true global optimum of an optimization problem or, if time constraints are an issue, it can return a solution that is within a specified distance form the optimum. This ability to evaluate the quality of a solution, makes IP an excellent candidate for conservation planning.

In a [previous post](http://strimas.com/r/gurobi/), I demonstrated how to use the commercial optimization software [Gurobi](http://www.gurobi.com/). Gurobi is great, but expensive, so in this post I'll explore solving the reserve design problem using free and open-source solvers.

# Problem formulation

In general, the goal of an optimization problem is to minimize an **objective function** over a set of **decision variables**, subject to a series of **constraints**. The decision variables are what we control, while the constraints can be thought of as rules that need to be followed. In the particular case of Marxan, the reserve design problem is formulated, for \\( n \\) planning units and \\( m \\) conservation features, as:

$$
\text{Minimize} \sum_{i=1}^{n}x_i c_i + 
b \sum_{i=1}^{n} \sum_{j=1}^{n}x_i (1-x_j)\nu_{ij} +
b \sum_{i=1}^{n} x_i\nu_{ii}
\text{ subject to } \sum_{j=1}^{n}x_j r_{ij}
\geq T_i \space \forall \space i 
$$

where \\( x_i \\) is a binary decision variable specifying whether planning unit \\( i \\) has been selected (1) or not (0), \\( c_i \\) is the cost of planning unit \\( i \\), \\( r_{ij} \\) is the representation level of feature \\( j \\) in planning unit \\( i \\), and \\( T_i \\) is the target for feature \\( i \\). \\( \nu_{ij} \\) is a matrix where off diagonal components are the lengths of the shared boundaries between planning units and diagonal components are external boundaries of planning units (i.e. those that are not shared with another planning unit). Finally, \\( b \\) is known as the **Boundary Length Modifier (BLM)** and determines how much emphasis should be placed on producing compact solutions relative to meeting targets and minimizing cost.

Since the decision variable is binary, this problem falls into the class of optimization problems known as integer programs. Furthermore, since the objective function is quadratic in the decision variables, this is an Integer Quadratic Program (IQP). When the BLM is zero, the objective function is linear and this is an Integer Linear Program (ILP). A variety of specialized  tools, both commercial and open source, exist to solve this class of problems.

To keep things simple, I'll start by focusing on the ILP form of the reserve design problem, which is a also known as the minimum set cover problem.

$$
\text{Minimize} \sum_{i=1}^{n}x_i c_i
\text{ subject to } \sum_{j=1}^{n}x_j r_{ij}
\geq T_i \space \forall \space i 
$$

# Solvers

The [CRAN Task View for Optimization](https://cran.r-project.org/web/views/Optimization.html#MathematicalProgrammingSolvers) list several open source linear programming solvers, and their R package interfaces:

- [**lp_solve**](http://lpsolve.sourceforge.net/5.5/) with R packages [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html) and [lpSolveAPI](https://cran.r-project.org/web/packages/lpSolveAPI/index.html).
- [**COIN-OR SYMPHONY**](https://projects.coin-or.org/SYMPHONY) with R packages [Rsymphony](https://cran.r-project.org/web/packages/Rsymphony/index.html), from CRAN, and [lpsymphony](https://www.bioconductor.org/packages/3.3/bioc/html/lpsymphony.html), from Bioconductor.
- [**COIN-OR Clp**](http://projects.coin-or.org/Clp) with R package [clpAPI](https://cran.r-project.org/web/packages/clpAPI/index.html).
- The [**GNU Linear Programming Kit**](http://www.gnu.org/software/glpk/) with R packages [glpkAPI](https://cran.r-project.org/web/packages/glpkAPI/index.html) and [Rglpk](https://cran.r-project.org/web/packages/Rglpk/index.html).

These solvers and packages can be challenging to install and configure correctly, and the specific steps required to do so are platform dependent. To get around this, I've created a [Docker](https://www.docker.com/) image that includes R Studio with all these R packages and their dependencies installed. To run a Docker container based on this image, visit the [repository on Docker Hub](https://hub.docker.com/r/mstrimas/optimizr/). The code in this post can be run inside this Docker container if you have trouble installing the R packages yourself.

[Gurobi](http://www.gurobi.com/) is a powerful commercial linear and quadratic programming solver. It is significantly more efficient than any of the above open source solvers, but it does require a license. I'll use it as a standard against which to measure the open source solvers. Gurobi is included in the Docker image, but requires configuration to add a license before use. Further details for how to do this are in the [Docker Hub](https://hub.docker.com/r/mstrimas/optimizr/) repository.

## Bounds and gap

As I understand it, most of these optimizers use some variation of an algorithm called [branch and bound](http://www.gurobi.com/resources/getting-started/mip-basics). One of the key features of this algorithm is that upper and lower bounds on the objective function are calculated at each step. The difference between the upper and lower bounds is known as the **gap** and gives an estimate of how close the current best solution is to the true global optimum. As the algorithm progresses, and the solution is refined, this gap becomes smaller until eventually is becomes zero when the global optimum is found. However, the algorithm can also stop at any point and return the current best solution along with an estimate of the quality of the solution (i.e. the gap). Finding the true global optimum is often too time consuming, so this ability to assess the quality of a non-optimal solution is one of the key benefits of integer programming.

## Packages


```r
library(dplyr)
library(sp)
library(raster)
library(rasterVis)
library(viridis)
library(slam)
library(protectr) # devtools::install_github("mstrimas/protectr")
# solvers
library(gurobi)
library(lpSolve)
library(lpSolveAPI)
library(Rsymphony)
library(lpsymphony)
library(clpAPI)
library(glpkAPI)
library(Rglpk)
set.seed(1)
```

# Preparation

First, I'll set up the problem and prepare some data.

## Data generation

To test the various methods, I'll generate 9 species distributions and a cost layer over a 10x10 grid of planning units (100 total). I've intentionally chosen an extremely simplified problem to start with to ensure that all the solvers will be able to handle it; later we'll look at different problem sizes. Although fabricated, these layers have some spatial auto-correlation built it to make them semi-realistic.


```r
# raster template
r <- extent(0, 100, 0, 100) %>% 
  raster(nrows = 10, ncols = 10, vals = 1)

# generate 9 feature distributions with different spatial scales and rarities
species <- mapply(function(x, y, r) gaussian_field(r = r, range = x, prop = y),
                  rep(c(5, 15, 25), each = 3),
                  rep(c(0.1, 0.25, 0.5), times = 3),
                  MoreArgs = list(r = r)) %>% 
  stack %>% 
  setNames(., letters[1:nlayers(.)])
levelplot(species, main = 'Species Distributions', layout = c(3, 3),
          scales = list(draw = FALSE),
          col.regions = c("grey20", "#fd9900"), colorkey = FALSE)
```

<img src="/figures//ilp-field-guide_species-1.png" title="plot of chunk species" alt="plot of chunk species" style="display: block; margin: auto;" />

```r
# genrate cost layer
cost_raster <- gaussian_field(r, 20, mean = 1000, variance = 500) %>% 
  setNames("cost")
levelplot(cost_raster, main = "Cost", margin = FALSE, col.regions = viridis)
```

<img src="/figures//ilp-field-guide_species-2.png" title="plot of chunk species" alt="plot of chunk species" style="display: block; margin: auto;" />

## Pre-processing

The various components of the optimization problem need to be prepared. Where possible, I use sparse matrices from the `slam` package to save memory. First, the representation matrix \\( r_{ij} \\) stores the representation level of feature \\( j \\) in planning unit \\( i \\).


```r
rij <- as.simple_triplet_matrix(t(unname(species[])))
```

I arbitrarily set targets to 30% of the total level of representation across the whole study area.


```r
targets <- 0.3 * cellStats(species, "sum")
```

Finally, I convert the cost raster to a numeric vector.


```r
cost <- cost_raster[[1]][]
```

# Solvers

For each solver I'll define a wrapper function that will take the same set of inputs and return the same outputs.

## General format

Wherever possible I'll try to conform to a standard format for the wrapper functions for each solver:

```
function(cost, rij, targets, gap, time_limit, first_feasible, bound)
```

The first three arguments will be the cost vector, representation matrix, and targets respectively, which specify the optimization model. The remaining arguments will specify the stopping condition:

- `gap`: the relative gap to optimality at which to stop. For example, setting `gap = 0.01` will ensure that the returned solution is at worst within 1% of optimality.
- `time_limit`: the amount of time in seconds to allow the algorithm to run. Using this argument will cause the solver to return the best solution after the given amount of time.  
- `first_feasible`: if this argument is `TRUE` the solver will return the first feasible solution found, i.e. the first solution found that meets all the constraints.

These are important because, in general, it will be too time consuming to find the true global optimum. Note that not all solvers will allow all these stopping conditions to be set and, in these cases, the wrapper function will only accept a subset of these arguments.

In some cases, the solver doesn't return the objective function bounds and therefore the gap to optimality can't be calculated. In these cases, the wrapper function accepts a final argument, `bound`, which gives the lower bound of the objective function as calculated by some other source (e.g. Gurobi).

## Gurobi

[Gurobi](http://gurobi.com) is a cutting edge commercial optimization solver. It's extremely fast, has an easy to use R interface, and can solve both linear and quadratic programs. The downside is that for many conservation applications it's prohibitively expensive at $12,000 for a single license, hence there is a need a viable free open source alternative. I treat Gurobi first since it will be the gold standard against which the remaining free alternatives will be measured.

The function `gurobi(model, params)` takes two arguments: `model` contains the various elements that define the optimization model and `params` is a named list of components specifying Gurobi parameters. The Gurobi documentation describes the [components of the model object](https://www.gurobi.com/documentation/6.5/refman/solving_models_with_the_gu.html) and contains the full list of [possible parameters](https://www.gurobi.com/documentation/6.5/refman/parameters.html#sec:Parameters).


```r
msc_gurobi <- function(cost, rij, targets, gap = 1e-4,
                       time_limit = Inf,
                       first_feasible = FALSE) {
  # construct model
  model <- list()
  # goal is to minimize objective function
  model$modelsense <- "min"
  # binary decision variables
  model$vtype <- "B"
  # objective function
  model$obj <- cost
  # structural constraints
  model$A <- rij
  model$rhs <- targets
  model$sense <- rep(">=", length(targets))

  # stopping conditions
  # gap to optimality
  params <- list(Presolve = -1, MIPGap = gap)
  # stop after specified number of seconds
  if (is.finite(time_limit)) {
    params <- c(params, TimeLimit = time_limit)
  }
  # first feasible solution
  if (first_feasible) {
    params <- c(params, SolutionLimit = 1)
  }
  
  # solve
  t <- system.time(
    results <- gurobi::gurobi(model, params)
  )
  # get rid of log file
  if (file.exists("gurobi.log")) {
    unlink("gurobi.log")
  }
  
  # prepare return object
  list(time = summary(t)[["user"]],
       x = results$x,
       objval = results$objval,
       objbound = results$objbound,
       gap = (results$objval / results$objbound - 1))
}
results_gurobi <- msc_gurobi(cost, rij, targets, gap = 0)
```

Gurobi easily finds the true global optimum almost instantly; no surprise since this was an intentionally simple problem. Let's look at the results.


```r
# objective function value for returned solution
results_gurobi$objval
#> [1] 19959.27
# gap to optimality
results_gurobi$gap
#> [1] 2.220446e-16
# time to solve
results_gurobi$time
#> [1] 0.021
# plot
plot_selection(cost_raster, results_gurobi$x, title = "Gurobi")
```

<img src="/figures//ilp-field-guide_gurobi-results-1.png" title="plot of chunk gurobi-results" alt="plot of chunk gurobi-results" style="display: block; margin: auto;" />

Now that we have this optimal solution, we can test out the other solvers to see if they produce similar results.

## lp_solve

`lp_solve` is a open source ILP solver with two different R packages that can access it.

### lpSolve package

The `lpSolve` package is a high-level interface to `lp_solve`. Since it doesn't interact with the low-level API functions of `lp_solve`, it has limited functionality and tends to be slow. The function `lp()` is the main interface to `lp_solve` and the arguments to this function define the optimization problem in a similar fashion to the components of the `model` object for `gurobi()`. As far as I can tell there is no means of specifying any of the stopping conditions with `lpSolve`.


```r
msc_lpsolve <- function(cost, rij, targets, bound = NA) {
  # convert rij to full matrix form, lpSolve can't handle sparse matrices
  rij <- as.matrix(rij)

  # solve
  t <- system.time(
    results <- lpSolve::lp(
      # goal is to minimize objective function
      direction = "min",
      # binary decision variables
      all.bin = TRUE,
      # objective function
      objective.in = cost,
      # structural constraints
      const.mat = rij,
      const.rhs = targets,
      const.dir = rep(">=", length(targets))
    )
  )
  
  # prepare return object
  list(time = summary(t)[["user"]],
       x = results$solution,
       objval = results$objval,
       objbound = bound,
       gap = (results$objval / bound - 1))
}
results_lpsolve <- msc_lpsolve(cost, rij, targets, bound = results_gurobi$objbound)
```

`lpSolve` finds the same optimal solution as Gurobi, and does so fairly quickly.


```r
# check that correct optimal solution was found
all.equal(results_lpsolve$x, results_gurobi$x)
#> [1] TRUE
# gap to optimality
results_lpsolve$gap
#> [1] 7.194245e-14
# time to solve
results_lpsolve$time
#> [1] 14.025
```

One strength of `lpSolve` is that it's extremely easy to install directly from CRAN; there are no dependencies or external libraries required. However after some testing, I've concluded it isn't a viable option for conservation planning. It doesn't provide a bound on the objective function so it's impossible to assess solution quality. More importantly, for problems only slightly bigger than this extremely simple example, `lpSolve` takes prohibitively long to produce a solution. Finally, there is no ability to specify a stopping condition.

### lpSolveAPI package

An alternative for using `lp_solve` is the `lpSolveAPI` package, which provides a low-level API interface for building and solving linear programs with `lp_solve`. With this package the optimization model is built up with a series of functions calls. Note also that the model object isn't an R object, R just stores a pointer to an external C object. Thus the model has to be built, and the solution accessed, via the set and get methods provided by the package since R can't directly access the model object.

Unlike `lpSolve`, `lpSolveAPI` provides access to all three stopping conditions discussed above. However, `lp_solve` does not explicitly return a lower bound on the objective function or the gap to optimality. It does print the objective function of the **relaxed solution** to screen though. The relaxed solution is simply the solution ignoring the constraints that the decision variables must be binary, i.e. it lets them be real numbers between 0 and 1 instead. This optimization problem is much easier to solve and acts as a starting point for finding the true solution. However, the relaxed solution is also a valid lower bound for the optimization problem because adding constraints can only function to make the minimum value of the objective function larger. Therefore, I extract this value from the screen output and use it as the lower bound.


```r
msc_lpsolveapi <- function(cost, rij, targets,
                           gap = 1e-4,
                           time_limit = Inf,
                           first_feasible = FALSE,
                           bound = NA) {
  # construct model with given number of constraints (i.e. features)
  # and decision variables (i.e. planning units)
  model <- lpSolveAPI::make.lp(nrow = nrow(rij), ncol = ncol(rij))
  # goal is to minimize objective function
  lpSolveAPI::lp.control(model, sense = "min")
  # binary decision variables
  lpSolveAPI::set.type(model, columns = seq_along(cost), type = "binary")
  # objective function
  lpSolveAPI::set.objfn(model, obj = cost)
  # structural constraints
  # set non-zero elements of constraint matrix
  for (k in seq_along(rij$v)) {
    set.mat(model, i = rij$i[k], j = rij$j[k], value = rij$v[k])
  }
  lpSolveAPI::set.rhs(model, b = targets)
  lpSolveAPI::set.constr.type(model, types = rep(">=", length(targets)))

  # the % gap to optimality at which to terminate
  lpSolveAPI::lp.control(model, 
                         verbose = "normal",
                         #presolve = c("rows", "cols"),
                         # first feasible solution
                         break.at.first = first_feasible,
                         # gap to optimality
                         mip.gap = c(1e-11, gap),
                         # stop after specified number of seconds
                         timeout = ifelse(is.finite(time_limit), time_limit, 0))

  # solve
  t <- system.time({
    screen_out <- capture.output(lpSolveAPI::solve.lpExtPtr(model))
  })
  
  # extract lower bound from screen output
  if (is.na(bound)) {
    bound <- stringr::str_subset(screen_out, "Relaxed solution")
    bound <- stringr::str_match(bound, "Relaxed solution\\s+([-+.e0-9]+)")
    bound <- ifelse(nrow(bound) > 0, as.numeric(bound[1, 2]), NA)
  }
  # prepare return object
  results <- list(time = summary(t)[["user"]],
                  x = lpSolveAPI::get.variables(model),
                  objval = lpSolveAPI::get.objective(model),
                  objbound = bound,
                  gap = (lpSolveAPI::get.objective(model) / bound - 1))
  #lpSolveAPI::delete.lp(model)
  return(results)
}
results_lpsolveapi <- msc_lpsolveapi(cost, rij, targets)
```

Again, `lpSolveAPI` finds the correct optimal solution, and does so fairly quickly.


```r
# check that correct optimal solution was found
all.equal(results_lpsolveapi$x, results_gurobi$x)
#> [1] TRUE
# gap to optimality
results_lpsolveapi$gap
#> [1] 0.03085664
# time to solve
results_lpsolveapi$time
#> [1] 2.967
# compare lower bound with Gurobi
c(gurobi = results_gurobi$objbound, lpSolveAPI = results_lpsolveapi$objbound)
#>     gurobi lpSolveAPI 
#>   19959.27   19361.83
```

This interface to `lp_solve` is also super easy to install, but there are significant improvements over `lpSolve`. First, it appears to be faster and able to handle larger problems. `lpSolveAPI` also allows for much more control over the solving algorithm, including the ability to set the stopping conditions. Setting the time limit and requesting the first feasible solution both work as expected. However, I haven't been able to get the gap to optimality condition to work. In my testing the solver will typically keep running even if the first feasible solution is already within the specified gap to optimality. I'm unclear why this is happening.

The lower bound taken from the screen output is just the relaxed solution objective function and isn't very precise. Gurobi does a much better job here.

`lpSolveAPI` also struggles with larger problems. For example, even with a 100x100 grid (10,000 planning units), it won't find the true global optima in a reasonable time frame. Furthermore, I've found that `lpSolveAPI` typically returns a first feasible solution quickly, but typically won't improve upon that solution even when left to run for a long time. Again, not sure what's going on here, but my takeaway has been that `lpSolveAPI` is mostly useful for finding that first feasible solution, which fortunately is often quite good.

## SYMPHONY

[SYMPHONY](https://projects.coin-or.org/SYMPHONY) is another open-source integer programming solver. It's part of the [Computational Infrastructure for Operations Research](http://www.coin-or.org/) (COIN-OR) project, an initiative to promote development of open-source tools for operations research (a field that includes linear programming).

Two different R packages have been developed to interact with SYMPHONY. On Mac OS I found it challenging to install both SYMPHONY and the R package interfaces. This is where the Docker image comes in handy.

### RSymphony package

`RSymphony` is the original R SYMPHONY interface and it appears on CRAN. It took me several hours to figure out how to get it working on Mac OS and there are many posts online with similar installation issues. Over on StackOverflow I outlined [exactly how I eventually got things working](http://stackoverflow.com/questions/32129191/osx-installing-rsymphony-linking-headers-and-libs/37599406#37599406). 

The key function for this package is `Rsymphony_solve_LP()` and the first several arguments specify the optimization model in much the same way as the components of the `model` object in `gurobi`. The remaining arguments can be used to set parameters for the solver, including the stopping conditions.

Note that `Rsymphony` requires an absolute gap to optimality, while the other solvers use a gap relative to the optimum, which is more intuitive. To address this by first solving the relaxed problem (i.e. with no constraint on decision variables being binary). The objective function for the relaxed solution is a lower bound on the objective function for the fully constrained problem, so multiplying it by the relative gap gives an estimate of the absolute gap. Ideally, I'd like to avoid this workaround since it does add to the execution time, but I see no other solution. 


```r
# find the relaxed solution
relaxed_rsymphony <- function(cost, rij, targets) {
  # bounded between 0 and 1
  n_pu <- length(cost)
  bounds <- list(lower = list(ind = seq.int(n_pu), val = rep(0, n_pu)),
                 upper = list(ind = seq.int(n_pu), val = rep(1, n_pu)))
  results <- Rsymphony::Rsymphony_solve_LP(
      # objective function
      obj = cost,
      # structural constraints
      mat = rij,
      rhs = targets,
      dir = rep(">=", length(targets)),
      # decision variables between 0 and 1
      types = "C",
      bounds = bounds,
      # goal is to minimize objective function
      max = FALSE
    )
  list(x = results$solution, objval = results$objval)
}
# solve the actual problem
msc_rsymphony <- function(cost, rij, targets,
                          gap = 1e-4,
                          time_limit = Inf,
                          first_feasible = FALSE,
                          bound = NA) {
  # convert relative to absolute gap
  relaxed <- relaxed_rsymphony(cost, rij, targets)
  gap <- gap * relaxed$objval
  rm(relaxed)
  t <- system.time(
    results <- Rsymphony::Rsymphony_solve_LP(
      # objective function
      obj = cost,
      # structural constraints
      mat = rij,
      rhs = targets,
      dir = rep(">=", length(targets)),
      # binary decision variables
      types = "B",
      # goal is to minimize objective function
      max = FALSE,
      # gap to optimality
      gap_limit = gap,
      # stop after specified number of seconds
      time_limit = ifelse(is.finite(time_limit), time_limit, -1),
      # first feasible solution
      first_feasible = first_feasible
    )
  )
  # prepare return object
  list(time = summary(t)[["user"]],
       x = results$solution,
       objval = results$objval,
       objbound = bound,
       gap = (results$objval / bound - 1))
}
results_rsymphony <- msc_rsymphony(cost, rij, targets, 
                                   bound = results_gurobi$objbound)
```

`Rsymphony` finds the correct optimal solution fairly quickly.


```r
# check that correct optimal solution was found
all.equal(results_rsymphony$x, results_gurobi$x)
#> [1] "Mean relative difference: 2"
# gap to optimality
results_rsymphony$gap
#> [1] 0.001121979
# time to solve
results_rsymphony$time
#> [1] 0.169
```

Overall, I like the simple interface that `Rsymphony` uses, however, the installation problems are a major deterrent. Furthermore, there appears to be no means of determining the optimality gap or lower bound.

### lpsymphony package

`lpsymphony` is almost identical to `Rsymphony`, however, the package ostensibly includes SYMPHONY itself, which is meant to ease installation. On Linux, I didn't find it any easier to install, but on Mac and Windows it appears to work much better. However, installing SYMPHONY directly ensures you get the most recent version. Also, `lpsymphony` is on [Bioconductor](https://www.bioconductor.org/packages/3.3/bioc/html/lpsymphony.html), not CRAN, so can't be installed with `install.packages()`. For these reasons I'm going to stick with `rsymphony` for now, but I include it here for completeness. On the up side, it does have a nice [vignette](https://www.bioconductor.org/packages/3.3/bioc/vignettes/lpsymphony/inst/doc/lpsymphony.pdf), which `Rsymphony` doesn't have.


```r
# find the relaxed solution
relaxed_lpsymphony <- function(cost, rij, targets) {
  # bounded between 0 and 1
  n_pu <- length(cost)
  bounds <- list(lower = list(ind = seq.int(n_pu), val = rep(0, n_pu)),
                 upper = list(ind = seq.int(n_pu), val = rep(1, n_pu)))
  results <- lpsymphony::lpsymphony_solve_LP(
      # objective function
      obj = cost,
      # structural constraints
      mat = rij,
      rhs = targets,
      dir = rep(">=", length(targets)),
      # decision variables between 0 and 1
      types = "C",
      bounds = bounds,
      # goal is to minimize objective function
      max = FALSE
    )
  list(x = results$solution, objval = results$objval)
}
# solve the actual problem
msc_lpsymphony <- function(cost, rij, targets,
                          gap = 1e-4,
                          time_limit = Inf,
                          first_feasible = FALSE,
                          bound = NA) {
  # convert relative to absolute gap
  relaxed <- relaxed_lpsymphony(cost, rij, targets)
  gap <- gap * relaxed$objval
  rm(relaxed)
  t <- system.time(
    results <- lpsymphony::lpsymphony_solve_LP(
      # objective function
      obj = cost,
      # structural constraints
      mat = rij,
      rhs = targets,
      dir = rep(">=", length(targets)),
      # binary decision variables
      types = "B",
      # goal is to minimize objective function
      max = FALSE,
      # gap to optimality
      gap_limit = gap,
      # stop after specified number of seconds
      time_limit = ifelse(is.finite(time_limit), time_limit, -1),
      # first feasible solution
      first_feasible = first_feasible,
      write_lp = TRUE,
      write_mps = TRUE
    )
  )
  
  # clean up files
  if(file.exists("_prep.MPS")) {
    unlink("_prep.MPS")
  }
  if(file.exists("_prep.LPT")) {
    unlink("_prep.LPT")
  }
  
  # prepare return object
  list(time = summary(t)[["user"]],
       x = results$solution,
       objval = results$objval,
       objbound = bound,
       gap = (results$objval / bound - 1))
}
results_lpsymphony <- msc_lpsymphony(cost, rij, targets,
                                     bound = results_gurobi$objbound)
```

`lpsymphony` finds the correct optimal solution.


```r
# check that correct optimal solution was found
all.equal(results_lpsymphony$x, results_gurobi$x)
#> [1] TRUE
# gap to optimality
results_lpsymphony$gap
#> [1] 2.220446e-16
# time to solve
results_lpsymphony$time
#> [1] 0.108
```

My comments on `lpsymphony` are the same as for `Rsymphony` since they're essentially the same package.

## Clp

[Clp](https://projects.coin-or.org/Clp) is another open-source linear programming solver that's part of the COIN-OR project. Unlike all the other solvers covered in this post, Clp does not allow for decision variables to be set as integer or binary. Therefore it isn't appropriate for reserve design, but I'll cover it anyway for completeness. It does allow bounds to be set on the decision variables though, so I constrain all decision variables to be between 0 and 1, which makes this as close to the actual reserve design problem as possible. This is known as the relaxation of the integer linear program because the constraints are relaxed.

The R package for this solver is `clpAPI`, which interacts with the low-level Clp functions. This package has a [short vignette](https://cran.r-project.org/web/packages/clpAPI/vignettes/clpAPI.pdf) outlining the basic usage. The interface for this package is similar to `lpSolveAPI`: the optimization model object is stored in memory and built up with a series of functions calls.

As with SYMPHONY, I found both the solver and package challenging to install on Mac OS. I was able to compile Clp from source, but struggled to install the R package because it could locate the Clp libraries and shared objects.


```r
msc_clpapi <- function(cost, rij, targets, bound = NA) {
  # construct model with given number of constraints (i.e. features)
  # and decision variables (i.e. planning units)
  model <- clpAPI::initProbCLP()
  # goal is to minimize objective function, max = -1, min = 1
  clpAPI::setObjDirCLP(model, lpdir = 1)
  # in vector of non-zero constraint matrix elements need indices for where
  # new columns start
  new_row <- diff(rij$j)
  new_row_i <- which(new_row != 0)
  col_starts <- c(0, rep(new_row_i, times = new_row[new_row_i]), length(rij$i))
  # load optimization problem
  clpAPI::loadProblemCLP(model,
                         # constraint matrix dimensions
                         ncols = ncol(rij), nrows = nrow(rij),
                         # row indices for non-zero elements, 0-indexed
                         ia = (rij$i - 1), 
                         # start indices of new columns, 0-indexed
                         ja = col_starts,
                         ra = rij$v, 
                         # bounds on decision variables
                         lb = rep(0, length(cost)), ub = rep(1, length(cost)),
                         # objective function coefficients
                         obj_coef = cost,
                         # bounds on contraints, i.e. rhs of contraint
                         rlb = targets, rub = NULL)
  # solve
  t <- system.time(
    clpAPI::solveInitialCLP(model)
  )
  # prepare return object
  results <- list(time = summary(t)[["user"]],
                  x = clpAPI::getColPrimCLP(model),
                  objval = clpAPI::getObjValCLP(model),
                  objbound = bound,
                  gap = (clpAPI::getObjValCLP(model) / bound - 1))
  clpAPI::delProbCLP(model)
  return(results)
}
results_clpapi <- msc_clpapi(cost, rij, targets,
                             bound = results_gurobi$objbound)
```

The results:


```r
# objective function value for returned solution
results_clpapi$objval
#> [1] 19361.83
# gap to optimality
results_clpapi$gap
#> [1] -0.02993301
# time to solve
results_clpapi$time
#> [1] 0.003
# plot
clp_sol <- cost_raster
clp_sol[] <- results_clpapi$x
levelplot(clp_sol, main = "clpAPI", margin = FALSE, col.regions = viridis)
```

<img src="/figures//ilp-field-guide_clpapi-results-1.png" title="plot of chunk clpapi-results" alt="plot of chunk clpapi-results" style="display: block; margin: auto;" />

Since Clp solves the relaxation, with continuous decision variables, we see that fractional protection is allowed here. The solution does appear structurally similar to the proper ILP solution though. Also, note that the gap to optimality is negative, indicating the solution is better than optimal. This occurs because the relaxed solution has fewer constraints and therefore can find a better solution than the proper ILP. This isn't what we're looking for, so I won't consider Clp any further in this post. 

## GLPK

The [GNU Linear Programming Kit](https://www.gnu.org/software/glpk/) is an open-source package for solving linear and integer linear programs. As with `lp_solve`, there are two distinct R packages for this solver.

### Rglpk package

`Rglpk` provides a simple, high-level interface to GLPK. The main function, `Rglpk_solve_LP()`, is almost identical to `Rsymphony_solve_LP()` from the `Rsymphony` package: individual components of the optimization model are passed as arguments to this function. There are few additional parameters to control the optimization algorithm, and no means of applying any of the stopping conditions.

CRAN has working binaries for this package, making it easy to install on Mac OS.


```r
msc_rglpk <- function(cost, rij, targets, bound = NA) {
  # prepare model and solve
  t <- system.time(
    screen_out <- capture.output({
      results <- Rglpk::Rglpk_solve_LP(
        # objective function
        obj = cost,
        # structural constraints
        mat = rij,
        rhs = targets,
        dir = rep(">=", length(targets)),
        # binary decision variables
        types = "B",
        # goal is to minimize objective function
        max = FALSE,
        # output to screen
        verbose = TRUE
      )
    })
  )
  
  # extract lower bound from screen output
  if (is.na(bound)) {
    bound <- stringr::str_subset(screen_out, "^\\+")
    bound <- stringr::str_subset(bound, ">=\\s+(-?[+.e0-9]+)")
    bound <- stringr::str_match(bound, ">=\\s+(-?[+.e0-9]+)")
    bound <- ifelse(nrow(bound) > 0, as.numeric(bound[nrow(bound), 2]), NA)
  }
  
  # prepare return object
  list(time = summary(t)[["user"]],
       x = results$solution,
       objval = results$optimum,
       objbound = bound,
       gap = (results$optimum / bound - 1))
}
results_rglpk <- msc_rglpk(cost, rij, targets)
```

`Rglpk` finds the correct optimal solution, and does so quickly.


```r
# check that correct optimal solution was found
all.equal(results_rglpk$x, results_gurobi$x)
#> [1] TRUE
# gap to optimality
results_rglpk$gap
#> [1] 0.02218146
# time to solve
results_rglpk$time
#> [1] 0.247
# compare lower bound with Gurobi
c(gurobi = results_gurobi$objbound, glpkAPI = results_rglpk$objbound)
#>   gurobi  glpkAPI 
#> 19959.27 19526.16
```

`Rglpk` benefits from a simple interface and an easy install process, however, it's lacking some important features. In particular, there's no way to set stopping conditions or to assess solution quality. However, lower bounds at each step are printed to screen as the algorithm progresses, so I've captured the output and parsed out the last bound listed. Unfortunately, the lower bound *after* the algorithm finishes isn't listed, hence the bound I extract is only an approximation and will generally be more lower (i.e. more conservative) than the true lower bound.

### glpkAPI package

`glpkAPI` provides a more feature rich interface to the low-level GLPK API. The interface for this package is similar to `lpSolveAPI`: the optimization model object is stored in memory and built up with a series of functions calls. This model is a C object and most of the functions in `glpkAPI` are directly equivalent to corresponding C functions in the GLPK API. There are various functions to create and modify an optimization model, to retrieve model information, and solve an optimization model.

This package is also easy to install from the conveniently provided binary on CRAN.  There are two nice vignettes for `glpkAPI`: a [quick start guide](https://cran.r-project.org/web/packages/glpkAPI/vignettes/glpkAPI.pdf) and a [full-featured vignette](https://cran.r-project.org/web/packages/glpkAPI/vignettes/glpk-gmpl-intro.pdf). However, neither of them explain how to solve integer linear programs. It turns out this isn't trivial and I had to consult the [GLPK](http://www.gnu.org/software/glpk/glpk.html#documentation) documentation to understand how it's done with the C API, then convert that over to the corresponding `glpkAPI` function calls.


```r
msc_glpkapi <- function(cost, rij, targets,
                        gap = 1e-4,
                        time_limit = Inf,
                        bound = NA) {
  # initialize an empty model
  model <- glpkAPI::initProbGLPK()
  glpkAPI::setProbNameGLPK(model, "reserve-design")
  # goal is to minimize objective function
  glpkAPI::setObjDirGLPK(model, glpkAPI::GLP_MIN)
  # initialize decision variables
  glpkAPI::addColsGLPK(model, ncols = length(cost))
  # objective function
  # also specify no bounds on decision variables
  glpkAPI::setColsBndsObjCoefsGLPK(model, j = seq_along(cost),
                                 lb = NULL, ub = NULL,
                                 obj_coef = cost,
                                 type = rep(glpkAPI::GLP_FR, length(cost)))
  # binary decision variables
  glpkAPI::setColsKindGLPK(model, j = seq_along(cost),
                           kind = rep(glpkAPI::GLP_BV, length(cost)))
  # structural constraints
  # initialize
  glpkAPI::addRowsGLPK(model, nrows = length(targets))
  # set non-zero elements of constraint matrix
  glpkAPI::loadMatrixGLPK(model, ne = length(rij$v),
                          ia = rij$i, ja = rij$j, ra = rij$v)
  # rhs
  glpkAPI::setRowsBndsGLPK(model, i = seq_along(targets),
                           lb = targets, ub = NULL,
                           type = rep(glpkAPI::GLP_LO, length(targets)))
  # presolve and automatically calculate relaxed solution
  # otherwise glpkAPI::solveSimplexGLPK(model) must be called first
  glpkAPI::setMIPParmGLPK(PRESOLVE, GLP_ON)
  glpkAPI::setMIPParmGLPK(MSG_LEV, GLP_MSG_ALL)
  # gap to optimality
  glpkAPI::setMIPParmGLPK(MIP_GAP , gap)
  # stop after specified number of seconds, convert to milliseconds
  if (is.finite(time_limit)) {
    glpkAPI::setMIPParmGLPK(TM_LIM, 1000 * time_limit)
  }
  
  # solve
  t <- system.time({
    screen_out <- capture.output(glpkAPI::solveMIPGLPK(model))
  })
  print(screen_out)
  # extract lower bound from screen output
  if (is.na(bound)) {
    bound <- stringr::str_subset(screen_out, "^\\+")
    bound <- stringr::str_subset(bound, ">=\\s+(-?[+.e0-9]+)")
    bound <- stringr::str_match(bound, ">=\\s+(-?[+.e0-9]+)")
    bound <- ifelse(nrow(bound) > 0, as.numeric(bound[nrow(bound), 2]), NA)
  }
  # prepare return object
  results <- list(time = summary(t)[["user"]],
                  x = glpkAPI::mipColsValGLPK(model),
                  objval = glpkAPI::mipObjValGLPK(model),
                  objbound = bound,
                  gap = (glpkAPI::mipObjValGLPK(model) / bound - 1))
  glpkAPI::delProbGLPK(model)
  return(results)
}
results_glpkapi <- msc_glpkapi(cost, rij, targets)
#> [1] "[1] 14"
```

`glpkAPI` finds the correct optimal solution, and does so quickly. As with `Rglpk`, the lower bound from the screen output is output to screen, however, unlike `Rglpk` it goes to the C standard out not the R standard out so can't easily be captured.


```r
# check that correct optimal solution was found
all.equal(results_glpkapi$x, results_gurobi$x)
#> [1] TRUE
# gap to optimality
results_glpkapi$gap
#> [1] NA
# time to solve
results_glpkapi$time
#> [1] 0.26
```

`glpkAPI` has a lot going for it: it's easy to install, feature rich, and quick. The interface can be confusing to use and it took me some sleuthing to figure out how to get everything working as desired. However, this is the benefit of the wrapper function: all that confusion can be abstracted away from the end user. The main thing it's missing it a good method for directly accessing the objective function lower bound.

## Conclusions

In my next post, I'll apply some of these open-source solvers to more realistic problems and do some bench-marking to get a sense for their performance. However, to conclude this post, I'll provide a summary of my observations thus far. Keep in mind that this is all based on an extremely simplified problem. Furthermore, some of these solvers have a wide range of settings that can be tweaked to improve performance. I'm still quite new to this field, so I haven't explored these more advanced options yet.

I'll break down my comments into the R packages I'm including and excluding from the next round of tests. Within these groups the packages are listed in order of my subjective ranking of how much I like them.

### General comments

There seems to a general trend of many of these packages being hard to install. Since they all rely on the libraries of external solvers, they need to link to these libraries when being compiled. Unfortunately, this is an extremely non-trivial process. Presumably for hardcore software developers it wouldn't be an issue, but for the average R user it's a major deterrent. I consider myself to be fairly computer savvy, so if this stuff is challenging for me, it's likely impossible for many conservation researchers. I created the [Docker image](https://hub.docker.com/r/mstrimas/optimizr/) mentioned previously to address this issue.

After installation challenges, my biggest complaint is that none of these packages provide a direct means of estimating the gap to optimality for the returned solution. I'm frankly baffled that this is the case. All the solvers need to calculate the gap internally for their algorithms to work, so why not return it with the solution? Assessing solution quality seems like the most important feature. I know very little about this field, so I must be missing something. Perhaps most users solve their problems all the way to optimality and therefore don't need to estimate quality. Or, maybe there is a way to access this information that I'm missing.

Finally, it's imporant to note that, unlike Gurobi, none of these packages can handle quadratic objective functions. This precludes all these open source packages from solving the full Marxan objective function, instead they can only be applied to the linear version with no boundary length modifier. It is possible to linearize the full Marxan objective function, however, this adds considerable theoretical and computational complexity.

### Exclude

Even from the simple testing I've done thus far, it's clear these packages aren't going to cut it.

- `Rglpk`: easy to install, an intuitive interface, fast, and able to handle larger problems. This is an excellent package for a excellent solver. If you're looking to solve small- to medium-sized problems all the way to optimality, and you don't want to worry about any of the more advanced settings, this is your best option. Unfortunately, conservation prioritization problems can be huge, so I need a package that provides access to these more advanced settings, especially the stopping conditions.
- `lpSolve`: this package is extremely easy to install and use, and it could be a good option for solving simple optimization problems. Unfortunately, it can't handle larger problems and doesn't provide access to any stopping conditions.
- `lpsymphony`: excluded because it's redundant since it's essentially identical to `Rsymphony`. If it was on CRAN not Bioconductor, I go with this package instead.
- `clpAPI`: doesn't allow for integer constraints, so of limited use for conservation prioritization. This package could be a good option if you're trying to solve linear programs with continuous variables. I did find it hard to install though, and the interface is somewhat confusing because it formulates the optimization model slightly differently than all the other solvers.

### Include

In the next post, I'll look into the performance of these packages and consider them as open-source alternatives to Gurobi for conservation prioritization.

- `glpkAPI`: easy to install, fast, able to handle larger problems, and very feature rich! The interface isn't very intuitive and to really harness the full power of the package you'll likely have to consult the GLPK documentation. I think this is the package with the most promise as an open-source alternative to Gurobi.
- `Rsymphony`: a nice user-friendly interface giving access to the most important settings (e.g. stopping conditions), but not cluttering things up with less frequently used settings. The biggest weakness is that, at least on Mac OS X, this package is extremely challenging to install. I wasted several hours and still couldn't get it to work. Regardless, this package has potential as an open-source alternative to Gurobi.
- `lpSolveAPI`: easy to install, fairly easy to use, and access to many more features than `lpSolve`, including stopping conditions. As with `lpSolve`, in preliminary testing it struggles with larger problems, but I'll keep it in the running for the time being.
