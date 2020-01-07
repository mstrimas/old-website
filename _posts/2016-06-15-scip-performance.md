---
layout: post
title: "ILP Solver Comparison Part II: SCIP"
published: true
excerpt: >
  Comparing the performance of open-source integer linear programming solvers.
  In this post, I look at the SCIP solver, which doesn't have an R interface.
category: prioritization
tags: r gurobi optimization marxan
---





In the last few posts I've been investigating the use of integer linear programming (ILP) to solve conservation prioritization problems. First, I [demonstrated how to use the commercial optimization software Gurobi](http://strimas.com/prioritization/gurobi/) to solve the Marxan reserve design problem in R. Gurobi is great, but expensive, so in the [next post](http://strimas.com/prioritization/ilp-field-guide/) I explored seven R packages that interface with 4 open-source ILP solvers. Then I looked in detail at 3 of those R packages to [compare their performance](http://strimas.com/prioritization/ilp-performance/) in solving conservation prioritization problems of different sizes. 

In this post, I'm going to look at one more ILP solver. [SCIP](http://scip.zib.de/) (Solving Constraint Integer Programs) is a free, though not totally open source solver, that has performed quite well in benchmarks against other free solvers. In addition, it is the only free solver I'm aware of that can handle problems with quadratic objective functions, such as the full Marxan objective function with the boundary length term.

The catch is that there is no R package for SCIP. Instead I've written a shell scripts to run all four solvers (Gurobi, SYMPHONY, GLPK, and SCIP) from the command line for the same set problem as outlined in the previous post, which I've exported in MPS format, a standard format for specifying optimization problems. In particular, I have 9 species defined on grids of increasing sizes: 100x100, 400x400, 700x700, 1,000x1,000.

# Shell script

Here's the shell script to run all four solvers for each of the four problems. As in the previous post, I've chosen to solve each problem to within 0.1% of optimality (`gap = 0.001`).

```bash
#!/bin/bash

# set up output file
echo "solver,n,time" > times.csv

# output time in seconds
TIMEFORMAT=%R

for i in *.mps; do
  printf "Processing file: %s\n" $i
  
  # grid dimension and number of planning units
  dim=$(echo $i | grep -o "[0-9]*\.mps$" | cut -d . -f 1)
  ((n = dim * dim))
  
  ## gurobi
  # location: /Library/gurobi651/mac64/bin/gurobi_cl
  printf ">Gurobi: "
  fname="output/gurobi-sol-$dim.sol"
  # solve and capture time
  t=$( { time gurobi_cl LogToConsole=0 ResultFile=$fname \
    MIPGap=0.001 $i; } 2>&1 )
  printf "gurobi,%i,%.4f\n" $n $t >> times.csv
  echo "$t seconds"
  # extract objective function value to calculate absolute gap for symphony
  # need to convert from scientific notation
  obj=$(grep "^# Objective" $fname | cut -d = -f 2)
  obj=$(printf "%.f" $obj)
  ((gap = obj / 1000)) 
  
  ## scip
  # location: /usr/local/scipoptsuite-3.2.1/scip-3.2.1/bin/
  printf ">SCIP: "
  fname="output/scip-log-$dim.txt"
  # solve and capture time
  t=$( { time scip -s scip-gap-0.1.set -l $fname -q -f $i; } 2>&1 )
  printf "scip,%i,%.4f\n" $n $t >> times.csv
  echo "$t seconds"
  
  ## symphony
  # location: /usr/local/Cellar/symphony/5.6.6/bin/
  printf ">SYMPHONY: "
  fname="output/symphony-log-$dim.txt"
  # solve and capture time
  t=$( { time symphony -g $gap -F $i > $fname; } 2>&1 )
  printf "scip,%i,%.4f\n" $n $t >> times.csv
  echo "$t seconds"
  
  ## glpk
  # location: /usr/local/Cellar/glpk/4.60/bin/
  printf ">GLPK: "
  fname="output/glpk-out-$dim.txt"
  lfile="output/glpk-log-$dim.txt"
  # solve and capture time
  t=$( { time glpsol --freemps $i --min --presol --primal --mipgap 0.001 \
    -o $fname > $lfile; } 2>&1 )
  printf "glpk,%i,%.4f\n" $n $t >> times.csv
  echo "$t seconds"
  
  printf "\n"
done
```

# Packages


```r
library(readr)
library(dplyr)
library(ggplot2)
```

# Results

I import the times for each run and plot them.


```r
times <- read_csv("data/ilp-performance/times.csv") %>% 
  mutate(via = "shell")
times_r <- readRDS("data/ilp-performance/comp-perf.rds") %>% 
  select(solver, n, time) %>% 
  mutate(via = "R")
times <- rbind(times, times_r)
# plot
mutate(times, solver = reorder(solver, time)) %>% 
  ggplot(aes(x = n/1000, y = time/60, colour = solver, linetype = via)) +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = c(10, 100, 1000), labels = scales::comma) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = scales::comma) +
  scale_color_brewer("ILP Solver", palette = "Set1") +
  scale_linetype("Via") +
  labs(x = "# of planning units (1000's)",
       y = "Time (minutes)",
       title = "ILP Performance Comparison (Gap = 0.1%)") +
  theme(text = element_text(family = "Helvetica Neue Light", size = 14),
        title = element_text(family = "Helvetica Neue Bold"),
        legend.position = c(0.25, 0.75),
        legend.background = element_rect(colour = "black"),
        legend.box = "horizontal",
        legend.key.size = unit(2, units = "lines"))
```

<img src="/figures//scip-performance_times-1.png" title="plot of chunk times" alt="plot of chunk times" style="display: block; margin: auto;" />

# Conclusions

First, there are a couple other interesting features visible on this plot. SYMPHONY does much better when run from the executable compared to via R (33 seconds vs. 400 seconds for a million planning units). Makes me feel like I've done something wrong. In addition, I'm very suprised by how poorly GLPK does compard to the other solvers. It's almost 100 times slower than SCIP and 300 times slower than SYMPHONY. Again, this makes me suspicious that I've done something wrong.

The main takeway is that SCIP performs well, however, it isn't quite as quick as SYMPHONY. For a problem with a million planning units, SCIP takes about 100 seconds, while SYMPHONY only about 30. However, SCIP does have the ability to solve quadratic problems, which makes it a promising candidate for a Gurobi alternative. If only there were an R package...
