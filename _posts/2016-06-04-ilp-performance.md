---
layout: post
title: "Open-source ILP Solver Performance Comparison"
published: true
excerpt: >
  Comparing the performance of open-source integer linear programming solvers.
  The focus is on finding an open-source alternative to Gurobi for conservation
  prioritization and proteced area design.
category: prioritization
tags: r gurobi optimization marxan
---





Conservation prioritization is fundamentally an optimization problem: what actions should we take to get the best possible conservation outcome for the lowest possible cost. In the context of protected area design, the goal is to identify the suite of sites that meets a set of biodiversity targets for the lowest possible cost. This is known as an integer linear progamming (ILP) problem in the language of operations research and mathematical optimization. A wide variety of tools (ranging from free to quite expensive) exist to solve such problems. In this series of posts I am invesitgating which tools are best for solving conservation prioritization problems in R.

In my [first post](http://strimas.com/prioritization/gurobi/), I demonstrated how to use the commercial optimization software [Gurobi](http://www.gurobi.com/) to solve the Marxan reserve design problem in R. Gurobi is great, but expensive, so in the [next post](http://strimas.com/prioritization/ilp-field-guide/) I explored seven R packages that interface with 4 open-source ILP solvers. That post left me with three R packages that are potential candidates for a free alternative to Gurobi. In this post, I'll focus on comparing the performance of these solvers.

## Packages


```r
library(dplyr)
library(tidyr)
library(sp)
library(raster)
library(rasterVis)
library(viridis)
library(slam)
library(protectr) # devtools::install_github("mstrimas/protectr")
library(ggplot2)
# solvers
library(gurobi)
library(lpSolveAPI)
library(Rsymphony)
library(glpkAPI)
set.seed(1)
```

# Preparation

First, I'll set up the problem and prepare some data.

## Data generation

To test the various methods, I'll generate 9 species distributions and a cost layer over a 100x100 grid of planning units (10,000 total). Although fabricated, these layers have some spatial auto-correlation built it to make them semi-realistic.


```r
# raster template
r <- extent(0, 100, 0, 100) %>% 
  raster(nrows = 100, ncols = 100, vals = 1)

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

<img src="/figures//ilp-performance_species-1.png" title="plot of chunk species" alt="plot of chunk species" style="display: block; margin: auto;" />

```r
# genrate cost layer
cost_raster <- gaussian_field(r, 20, mean = 1000, variance = 500) %>% 
  setNames("cost") %>% 
  {.[[1]]}
levelplot(cost_raster, main = "Cost", margin = FALSE, col.regions = viridis)
```

<img src="/figures//ilp-performance_species-2.png" title="plot of chunk species" alt="plot of chunk species" style="display: block; margin: auto;" />

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
cost <- cost_raster[]
```

# Solvers

In [my previous post](http://strimas.com/prioritization/ilp-field-guide/), I looked at four linear programming solvers, accessed via seven R packages. For a [variety of reasons](http://strimas.com/prioritization/ilp-field-guide/#conclusions) outlined in that post, I excluded several packages, leaving me with three potential candidate R packages for an open-source alternative to Gurobi:

- **Solver:** `lp_solve`; **R Package:** `lpSolveAPI`
- **Solver:** SYMPHONY; **R Package:** `Rsymphony`
- **Solver:** GLPK; **R Package:** `glpkAPI`

In [that post](http://strimas.com/prioritization/ilp-field-guide/), I went through the process of building up standardized wrapper functions for solving the reserve design problem using these ILP solvers. Consult that post, or the [source code](https://github.com/mstrimas/mstrimas.github.io/blob/master/_source/2016-05-27-ilp-performance.rmd) for this post, for the function definitions. Here I just give the signatures of the wrapper functions.



## Gurobi

[Gurobi](http://gurobi.com) is a cutting edge commercial optimization solver. It's extremely fast, has an easy to use R interface, and can solve both linear and quadratic programs. The downside is that for many conservation applications it's prohibitively expensive at $12,000 for a single license, hence there is a need a viable free open source alternative. I treat Gurobi first since it will be the gold standard against which the remaining free alternatives will be measured.

```r
msc_gurobi(cost, rij, targets,
           gap = 1e-4, time_limit = Inf, first_feasible = FALSE,
           bound = NA)
```

## lp_solve

`lp_solve` is a open source ILP solver that is best accessed with the `lpSolveAPI` package, which provides a low-level API interface for building and solving linear programs .

```r
msc_lpsolve(cost, rij, targets,
            gap = 1e-4, time_limit = Inf, first_feasible = FALSE,
            bound = NA)
```

## SYMPHONY

[SYMPHONY](https://projects.coin-or.org/SYMPHONY) is another open-source integer programming solver. It's part of the [Computational Infrastructure for Operations Research](http://www.coin-or.org/) (COIN-OR) project, an initiative to promote development of open-source tools for operations research (a field that includes linear programming). The `RSymphony` package provides a convenient R interface to this solver.

```r
msc_symphony(cost, rij, targets,
             gap = 1e-4, time_limit = Inf, first_feasible = FALSE,
             bound = NA)
```

## GLPK

The [GNU Linear Programming Kit](https://www.gnu.org/software/glpk/) is an open-source package for solving linear and integer linear programs. The R package `glpkAPI` provides a more feature rich interface to the low-level GLPK API.

```r
msc_glpk(cost, rij, targets,
         gap = 1e-4, time_limit = Inf,
         bound = NA)
```

# Stopping conditions

For simple problems, these solvers will all find the true optimal solution. However, in general, this will not be possible due to computational constraints and time limitations. Therefore, it's usually necessary to tell the solvers to stop early and I've built in three stopping conditions to do this. Before looking at performance, I'll see how these solvers deal with the stopping conditions.

I'll start by finding the exact solution with Gurobi. This will facilitate evaluating how well each solver does.


```r
exact <- msc_gurobi(cost, rij, targets, gap = 0)
exact_val <- exact$objval
```

## Time limit

The simplest approach is to give the solver a time limit and ask it to return the best solution obtained after it's run for that amount of time. I'll give each solver 5 seconds to work on the problem defined above.


```r
tl_lpsolve <- msc_lpsolve(cost, rij, targets, time_limit = 5, bound = exact_val)
tl_symphony <- msc_symphony(cost, rij, targets, time_limit = 5, bound = exact_val)
tl_glpk <- msc_glpk(cost, rij, targets, time_limit = 5, bound = exact_val)
```

The run times are all about 5 seconds as desired.


```r
tl_lpsolve$time
#> [1] 4.608
tl_symphony$time
#> [1] 1.871
tl_glpk$time
#> [1] 5.443
```

All the solvers got fairly close to optimality.


```r
# % gap to optimality
100 * tl_lpsolve$gap
#> [1] 0.02992902
100 * tl_symphony$gap
#> [1] 0.01071333
100 * tl_glpk$gap
#> [1] 0.008289119
```

Usually time is the limiting factor in an analysis, and specifying the run time allows you to control this. Therefore it's usually the best place to start. Run the solver for a reasonable amount of time, then inspect the solution to see how close to optimality it is. Depending on the results you may want to run the solver for a longer time or with a gap specified.

## First feasible

These solvers all start by solving the relaxed problem, which is the optimization problem without the constraint that all the decision variables are integer, then work to find solutions that do meet the integrality constraints. The first solution that meets all the constraints is called the first feasible solution. All the packages, except `glpkAPI`, give the option of returning this solution. Unless you're lucky this won't be the true optimum, but it's often pretty good because it's not just any feasible solution, it's one derived from the optimal relaxed solution.


```r
ff_lpsolve <- msc_lpsolve(cost, rij, targets, first_feasible = TRUE,
                          bound = exact_val)
ff_symphony <- msc_symphony(cost, rij, targets, first_feasible = TRUE,
                            bound = exact_val)
```

This approach has the benefit of being a the quickest way to get a feasible solution.


```r
ff_lpsolve$time
#> [1] 1.658
ff_symphony$time
#> [1] 1.822
```

However, note that `lpSolveAPI` does no better than this first feasible solution when allowed to run longer. In fact, even allowing it to run *much* longer has no impact, the returned solution is still the first feasible solution.


```r
# ratio of objective function values, first feasible / 5 second run
ff_lpsolve$objval / tl_lpsolve$objval
#> [1] 1
ff_symphony$objval / tl_symphony$objval
#> [1] 1
```

Requesting the first feasible solution is a good approach if you're more concerned with time that getting a solution close to optimality. This will get you a solution that satisfies all the constraints quickest.

## Gap to optimality

These solvers all use some variation of an algorithm called [branch and bound](http://www.gurobi.com/resources/getting-started/mip-basics). One of the key features of this algorithm is that it calculates and constantly refines upper and lower bounds on the objective function. The difference between the upper and lower bounds is known as the **gap** and gives an estimate of how close the current best solution is to the true global optimum. As the algorithm progresses, and the solution is refined, this gap becomes smaller until eventually is becomes zero when the global optimum is found. However, the algorithm can also stop at any point and return the current best solution along with an estimate of the quality of the solution (i.e. the gap).

Here I'll request that the solvers return when they get within 0.1% from the optimal solution.


```r
gap_lpsolve <- msc_lpsolve(cost, rij, targets, gap = 0.001, time_limit = 60,
                           bound = exact_val)
gap_symphony <- msc_symphony(cost, rij, targets, gap = 0.001, bound = exact_val)
gap_glpk <- msc_glpk(cost, rij, targets, gap = 0.001, bound = exact_val)
```

Note that for all three solvers, the 5 second run above was sufficient to get a solution within 1% of optimality. Despite this, when I ran `lpSolveAPI` for several minutes it never returned a solution, so I had to enforce a 60 second time limit. The gap that these solvers calculate is an upper bound, which is refined as the algorithm progresses. In contrast, I know exactly how far these solutions are from optimality since I found the optimum with Gurobi. Given that `lp_solve` isn't stopping execution despite having a solution within the specified gap, it suggests that something is wrong. Either `lp_solve` isn't doing a great job of estimating the objective function bounds, or there's an issue with the interface.

In addition, the gap that SYMPHONY takes is an absolute gap, rather than a relative gap like the other solvers. I've gotten around this by first solving the relaxed problem (i.e. with no constraint on decision variables being binary) with [Clp](https://projects.coin-or.org/Clp). The objective function for the relaxed solution is a lower bound on the objective function for the fully constrained problem, so multiplying it by the relative gap gives an estimate of the absolute gap. Ideally, I'd like to avoid this workaround since it does add to the execution time, but I see no other solution. 


```r
gap_lpsolve$time # set to stop at 60 seconds
#> [1] 59.017
gap_symphony$time
#> [1] 1.887
gap_glpk$time
#> [1] 0.568
```

All the solvers do much better than the desired 0.1% gap.


```r
100 * gap_lpsolve$gap
#> [1] 0.02992902
100 * gap_symphony$gap
#> [1] 0.006386605
100 * gap_glpk$gap
#> [1] 0.01525037
```

## Bound

While all these solvers must internally calculate a lower bound on the objective function, none of the open-source solvers return these bounds in an easily accessible form. Knowing this lower bound is critical for calculating the quality of the solution (i.e. the gap to optimality). I've had to deal with getting a bound in a different manner for each solver:

- `gurobi`: conveniently returned as a component of the Gurobi results object (`result$objbound`).
- `lpSolveAPI`: printed to screen while lp_solve runs, however, it only allows for precision up to a tenth of a percent, anything below 0.1% appears as 0.0%. I extract this from the screen output with regular expressions.
- `Rsymphony`: not provided at all.
- `glpkAPI`: printed to screen when GLPK finishes, however, not accessible because it's sent to C standard output not R standard output.

So, really none of the open-source solvers do a great job in this department, either the bound is missing or returned in a simplified form.

## Summary

Here are some general takeaways from the above exploration with respect to stopping condition:

- `gurobi` can be used with a gap, a time limit, or to return the first feasible solution.
- `lpSolveAPI` can only be used to return the first feasible solution. Setting a gap doesn't appear to work and, even with a long time limit, `lp_solve` doesn't improve upon this first feasible solution.
- `Rsymphony` can be used with a gap, a time limit, or to return the first feasible solution, however, requires a workaround to convert a relative to absolute gap.
- `glpkAPI` works with either a gap or a time limit. There is no option to return the first feasible solution.

# Performance comparison

Now let's look at how these solvers stack up in terms of performance. I'll dis-aggregate the study grid into successively smaller grid cells (i.e. more decision variables) and solve the reserve design problem to within 0.1% of optimality using all the solvers except `lpSolveAPI` since it's quite slow.


```r
# data frame template to store results
comp_perf <- data_frame(solver = character(),
                        dimensions = integer(),
                        n = integer(),
                        time = numeric(),
                        objval = numeric(),
                        objbound = numeric(),
                        gap = numeric())
# iteratively disaggregate and solve
for (i in seq(1, 10, by = 3)) {
  message(sprintf("i = %i", i))
  dims <- 100 * i
  # disaggregate
  cost_tmp <- disaggregate(cost_raster, i)
  species_tmp <- disaggregate(species, i)
  
  # prepare inputs
  rij <- as.simple_triplet_matrix(t(unname(species_tmp[])))
  targets <- 0.3 * cellStats(species_tmp, "sum")
  cost <- cost_tmp[]
  rm(cost_tmp, species_tmp)
  
  # solve
  # gurobi
  message("solve gurobi")
  sol_gurobi <- msc_gurobi(cost, rij, targets, gap = 0.001)
  gurobi_bound <- sol_gurobi$objbound
  sol_gurobi$x <- NULL
  comp_perf <- list(solver = "gurobi",
                    dimensions = dims,
                    n = dims^2) %>%
    c(sol_gurobi) %>% 
    bind_rows(comp_perf, .)
  rm(sol_gurobi); gc()
  # symphony
  message("solve symphony")
  sol_symphony <- msc_symphony(cost, rij, targets, gap = 0.001,
                               bound = gurobi_bound)
  sol_symphony$x <- NULL
  comp_perf <- list(solver = "symphony",
                    dimensions = dims,
                    n = dims^2) %>% 
    c(sol_symphony) %>% 
    bind_rows(comp_perf, .)
  rm(sol_symphony); gc()
  # glpk
  message("solve glpk")
  fname <- sprintf("data/ilp-performance/model-dim-%i.mps", dims)
  sol_glpk <- msc_glpk(cost, rij, targets, gap = 0.001,
                       bound = gurobi_bound,
                       filename = fname)
  sol_glpk$x <- NULL
  comp_perf <- list(solver = "glpk",
                    dimensions = dims,
                    n = dims^2) %>% 
    c(sol_glpk) %>% 
    bind_rows(comp_perf, .)
  rm(sol_glpk); gc()
}
saveRDS(comp_perf, "data/ilp-performance/comp-perf.rds")
```



Let's take a look at the results:


```r
comp_perf %>% 
  mutate(solver = reorder(solver, time)) %>% 
  ggplot(aes(x = n/1000, y = time/60, colour = solver)) +
    geom_point() +
    geom_line() +
    scale_x_log10(breaks = c(10, 100, 1000), labels = scales::comma) +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = scales::comma) +
    scale_color_brewer("ILP Solver", palette = "Set1") +
    labs(x = "# of planning units (1000's)",
         y = "Time (minutes)",
         title = "ILP Performance Comparison (Gap = 0.1%)") +
    theme(text = element_text(family = "Helvetica Neue Light"),
          legend.position = c(0.15, 0.85),
          legend.background = element_rect(colour = "black"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16, face = "bold"),
          legend.key.size = unit(2, units = "lines"))
```

<img src="/figures//ilp-performance_performance-results-1.png" title="plot of chunk performance-results" alt="plot of chunk performance-results" style="display: block; margin: auto;" />

As expected, Gurobi is the clear winner here and, of the open source solvers, SYMPHONY does the best. Note that this is on a log-log scale, so the differences between solvers are actually huge. Looking at the solve times in **minutes** in tabular form:


```r
comp_perf %>% 
  mutate(solver = reorder(solver, time), time = time / 60) %>% 
  dplyr::select(solver, time, n) %>% 
  spread(solver, time) %>% 
  knitr::kable(digits = 1, format.args = list(big.mark = ","))
```



|         n| gurobi| symphony| glpk|
|---------:|------:|--------:|----:|
|    10,000|    0.0|      0.0|  0.0|
|   160,000|    0.0|      0.3|  2.2|
|   490,000|    0.1|      1.7| 23.3|
| 1,000,000|    0.1|      6.8| 94.1|

For a problem with a million planning units, GLPK takes 1.5 hours, SYMPHONY takes 7 minutes, and Gurobi just a few seconds. One interesting observation not reflected in these data is that GLPK spends the majority of its time solving the relaxed problem. In contrast, SYMPHONY solves the relaxed problem quite quickly.

# Conclusion

This posts confirms that Gurobi is clearly much better than the open-source alternatives. It solves large problems in a matter of seconds that the other solves either can't solve or take several orders of magnitude longer to solve. In addition, we've been looking at a simplified linear objective function, but the real Marxan objective function is quadratic. None of the open source solvers I've looked at thus far can handle quadratic problems, however, Gurobi can.

The impressive performance of Gurobi aside, conservation planners may need or want to use free, open-source software. In this post, and the [previous](http://strimas.com/prioritization/ilp-field-guide/), I explored seven R packages that interface with four open-source solvers. I've reduced that field of seven to two possible candidates for a Gurobi alternative for conservation prioritization:

- The R package `Rsymphony`, which interfaces with the SYMPHONY solver, can handle large prioritization problems reasonably quickly. The interface is intuitive and easy to use, however, the ability to customize algorithm parameters is limited. The downside is that it can be challenging to install. In some cases, `lpsymphony` on Bioconductor may be an easier-to-install alternative.
- In contrast, the `glpkAPI` package is easier to install and provides fine scale control of algorithm parameters. It can handle large prioritization problems, but the run time is slow above a couple hundred thousands planning units.

Neither solver returns gap to optimality or a lower bound on the solution.
## Next steps

Moving forward, there are three things I'd like to investigate further:

1. Can tweaking the parameters of the solving algorithms lead to significant performance benefits?
2. I haven't covered the [SCIP](http://scip.zib.de) solver because it doesn't have an R interface. Several sources have suggested to me that this is the best free solver. I'd like to test SCIP on the same set of prioritization problems as in this post. If it does prove to be much faster, I'll look into how hard it would be to create an R interface.
3. Finally, I hope to build an R package that provides access to a variety of solvers for conservation prioritization.
