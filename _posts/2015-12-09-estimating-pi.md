---
layout: post
title: Estimating pi with GIS and Monte Carlo methods in R
published: true
excerpt: >
  Using Monte Carlo methods, randomly placed spatial objects, and topological 
  operations in R to estimate pi.
category: gis
tags: r spatial gis monte-carlo
---

While reading the Wikipedia page on [Monte Carlo methods](https://en.wikipedia.org/wiki/Monte_Carlo_method) recently I came across two interesting spatial approaches for estimating \\( \\pi \\). Both involve randomly placing a large number of spatial objects on a plane (points and lines, respectively), performing some topological operations (testing for containment and intersection, respectively), and summarizing over all spatial objects to estimate \\( \\pi \\). While both approaches are probably more easily and efficiently done outside of a GIS, I thought it would be fun to try them using some of the spatial objects and spatial analysis tools in R.  

# Monte Carlo Methods  

The fact that I was reading about Monte Carlo on Wikipedia is a pretty good indication of the amount I know about these techniques (not much!). So, I'll resort to the concise definition from [Wolfram](http://mathworld.wolfram.com/MonteCarloMethod.html):  

> Any method which solves a problem by generating suitable random numbers and observing that fraction of the numbers obeying some property or properties.  

This definition is a bit vague because the term encompasses a very broad class of algorithms spanning many disciplines. The key seems to be the use of repeated draws of random numbers from some probability distribution to solve a complex problem or simulate a complex system. I think the following examples will make this more concrete.  

## Required packages  


```r
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
set.seed(1)
```

# Incribed Circle Approach  

[Wikipedia](https://en.wikipedia.org/wiki/Monte_Carlo_method#Introduction) gives the following method for estimating \\( \\pi \\) as a simple example of a Monte Carlo method:  

1.  Draw a square with a circle inscribed within it. Note that the ratio of the areas of the circle and square is \\( \\pi / 4 \\).  
2.  Randomly distribution a large number of points over the square, such that there is a uniform probability of a given point being anywhere within the square.  
3.  Count the number of points inside the circle and the total number of points.  
4.  The ratio of the numbers of points will equal the ratio of the areas (i.e. \\( \\pi / 4 \\)), so multiplying this ratio by 4 will give an estimate of  \\( \\pi \\).  

To see where the area ratio of \\( \\pi / 4 \\) comes from, note that a square with side length \\( s \\) has area \\( s^2 \\), and a circle inscribed within this square will have radius \\( s / 2 \\) and hence area \\(s ^ 2 / 4 \\). The Wikipedia page has a nice animated visual representation of this technique:  

<p align='center'>
  <img src='https://upload.wikimedia.org/wikipedia/commons/8/84/Pi_30K.gif'
  alt='Inscribed circle approach' />  
</p>

## R Implementation  

First I draw a square with a circle inscribed.  


```r
square <- readWKT('POLYGON((-1 -1,-1 1,1 1,1 -1, -1 -1))')
circle <- readWKT('POINT(0 0)') %>% 
  gBuffer(width = 1, quadsegs = 25)
plot(square, axes = F, border = '#FA6900', lwd = 2)
plot(circle, add = T, col = '#69D2E7', border = 'transparent')
```

<img src="/figures//estimating-pi_square-circle-1.svg" title="plot of chunk square-circle" alt="plot of chunk square-circle" style="display: block; margin: auto;" />

[Well-known Text](https://en.wikipedia.org/wiki/Well-known_text) is a simple markup language for representing vector geometries in plain text, and `rgeos::readWKT()` creates spatial objects based on a WKT strings. It's the easiest way I know of to quickly create simple geometries in R. To generate the circle I put a 1 unit buffer around a point. The `quadsegs` parameter in `rgeos::gBuffer()` sets the number of line segments to use to approximate a quarter circle, which is a measure of the smoothness of the buffering polygon. The default value of `quadsegs = 5` results in a circle with noticeable edges, but `quadsegs = 25` seems to create a nice smooth circle.  

Next I randomly sample a large number points within this square using `sp::spsample()`.  


```r
n_pts <- 1000
pts <- spsample(square, n_pts, type = 'random')
plot(square, axes = F, border = '#FA6900', lwd = 2)
plot(circle, add = T, col = '#69D2E7', border = 'transparent')
plot(pts, add = T, pch = 21, cex = 0.25, col = '#333333')
```

<img src="/figures//estimating-pi_spsample-1.svg" title="plot of chunk spsample" alt="plot of chunk spsample" style="display: block; margin: auto;" />

There are a few ways to determine which of these points are within the circle. First, the function `over(x, y)` in the `sp` package gives the indexes of spatial object `y` at the spatial locations of object `x`. If a given feature in `x` is not contained within `y` `NA` is returned.  


```r
over(pts, circle) %>% 
  {!is.na(.)} %>% 
  sum
#> [1] 768
```

Alternatively, the `sp` package offers an idiom of sorts to extract the features in one geometry that are within another geometry.  


```r
(pts_within <- pts[circle, ])
#> class       : SpatialPoints 
#> features    : 768 
#> extent      : -0.9738448, 0.9855762, -0.9905714, 0.981201  (xmin, xmax, ymin, ymax)
#> coord. ref. : NA
length(pts_within)
#> [1] 768
```

In addition to these methods in the `sp` package, the `rgeos` package provides a more complete set of topological operations. `gIntersects()` tests if two geometries overlap, either overall (with `byid = F`) or at the level of individual features (with `byid = T`).  


```r
n_within <- gIntersects(pts, circle, byid = T) %>% 
  sum
n_within
#> [1] 768
```

Finally, taking the ratio of points within the circle to total points and multiplying by 4, gives an estimate of \\( \\pi \\).  


```r
(pi_est <- 4 * n_within / n_pts)
#> [1] 3.072
round(abs(100 * (pi_est / pi - 1)), 2)
#> [1] 2.22
```

So, with 1000 points I get an error of about 2.22%; not too bad!  

## Convergence  

To see how the estimates from this approach converge to the true value of \\( \\pi \\), I repeat the process for different numbers of points. First, I'll write a simple function to estimate pi based on a given number of points. Out of curiosity, I also time it.  


```r
estimate_pi <- function(n_pts) {
  t <- system.time({
    square <- readWKT('POLYGON((-1 -1,-1 1,1 1,1 -1, -1 -1))')
    circle <- readWKT('POINT(0 0)') %>% 
      gBuffer(width = 1, quadsegs = 25)
    n_within <- spsample(square, n_pts, type = 'random') %>% 
      gIntersects(circle, byid = T) %>% 
      sum
    pi_estimate <- (4 * n_within / n_pts)
  })
  return(data.frame(pi_estimate, t = t[['user.self']]))
}
estimate_pi(1e4)
#>   pi_estimate     t
#> 1      3.1408 0.083
```

And I run this for a range of values of `n_pts`, estimating \\( \\pi \\) multiple times at each to get a sense of the variability.  


```r
n_values <- 10^(3:6)
n_reps <- 10
estimates <- expand.grid(n_values, 1:n_reps) %>% 
  setNames(c('n_pts', 'run')) %>% 
  tbl_df %>% 
  group_by(n_pts, run) %>% 
  do(estimate_pi(.$n_pts))
```

Plotting the results:  


```r
se <- function(x) {
  sd(x) / sqrt(length(x))
}
estimate_summary <- estimates %>% 
  mutate(pct_error = abs((pi_estimate / pi - 1))) %>% 
  group_by(n_pts) %>% 
  summarize(pi_mean = mean(pi_estimate),
            pi_se = se(pi_estimate),
            error_mean = mean(pct_error),
            error_se = se(pct_error),
            t_mean = mean(t))
# estimate
ggplot(estimate_summary, aes(x = n_pts, y = pi_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = pi_mean - pi_se, ymax = pi_mean + pi_se), width = 0.1) +
  geom_line() +
  geom_hline(aes(yintercept=pi), color = 'red', linetype = 'dashed') +
  scale_x_log10() +
  labs(x = 'Number of Points', y = expression(pi ~ estimate))
```

<img src="/figures//estimating-pi_pi-convergence-1.svg" title="plot of chunk pi-convergence" alt="plot of chunk pi-convergence" style="display: block; margin: auto;" />

```r
# error
ggplot(estimate_summary, aes(x = n_pts, y = error_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = error_mean - error_se, ymax = error_mean + error_se), width = 0.1) +
  geom_line() +
  scale_x_log10() +
  scale_y_continuous(labels = percent, breaks = 0:6 * 0.01) +
  labs(x = 'Number of Points', y = 'Error')
```

<img src="/figures//estimating-pi_pi-convergence-2.svg" title="plot of chunk pi-convergence" alt="plot of chunk pi-convergence" style="display: block; margin: auto;" />

```r
# times
ggplot(estimate_summary, aes(x = n_pts, y = t_mean)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(labels = comma) +
  labs(x = 'Number of Points', y = 'Execution Time (s)')
```

<img src="/figures//estimating-pi_pi-convergence-3.svg" title="plot of chunk pi-convergence" alt="plot of chunk pi-convergence" style="display: block; margin: auto;" />

So, with a million points, I get an error of 0.0813%; however, there is clearly an effect of diminishing returns: execution time increases linearly with number of points, but gains in precision are decreasing exponentially.  

## Alternate method using distance  

In the above approach, the call to `gIntersects()` is the bottle neck for large numbers of points. An alternate approach that avoids this is to calculate the distance of each point from the origin and check is that distance is less that the radius of the circle (in this case 1). This is facilitated by the `gDistance()` function from the `rgeos` package, which calculates a pairwise distance matrix between features in two geometries when `byid = T`.  


```r
estimate_pi_distance <- function(n_pts) {
  t <- system.time({
    square <- readWKT('POLYGON((-1 -1,-1 1,1 1,1 -1, -1 -1))')
    centre <- readWKT('POINT(0 0)')
    n_within <- spsample(square, n_pts, type = 'random') %>% 
      gDistance(centre, byid = T) %>% 
      {. <= 1} %>% 
      sum
    pi_estimate <- (4 * n_within / n_pts)
  })
  return(data.frame(pi_estimate, t = t[['user.self']]))
}
estimate_pi_distance(1e5)
#>   pi_estimate     t
#> 1      3.1438 0.643
```

For comparison, I estimate \\( \\pi \\) at a range of numbers of points as above.  


```r
estimates_d <- expand.grid(n_values, 1:n_reps) %>% 
  setNames(c('n_pts', 'run')) %>% 
  tbl_df %>% 
  group_by(n_pts, run) %>% 
  do(estimate_pi_distance(.$n_pts))
```

Plotting the results:  


```r
estimate_d_summary <- estimates_d %>% 
  mutate(pct_error = abs((pi_estimate / pi - 1))) %>% 
  group_by(n_pts) %>% 
  summarize(pi_mean = mean(pi_estimate),
            pi_se = se(pi_estimate),
            error_mean = mean(pct_error),
            error_se = se(pct_error),
            t_mean = mean(t))
# estimate
ggplot(estimate_d_summary, aes(x = n_pts, y = pi_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = pi_mean - pi_se, ymax = pi_mean + pi_se), width = 0.1) +
  geom_line() +
  geom_hline(aes(yintercept=pi), color = 'red', linetype = 'dashed') +
  scale_x_log10() +
  labs(x = 'Number of Points', y = expression(pi ~ estimate))
```

<img src="/figures//estimating-pi_pi-convergence-distance-1.svg" title="plot of chunk pi-convergence-distance" alt="plot of chunk pi-convergence-distance" style="display: block; margin: auto;" />

This method does seem to converge a bit better, with a million points, we get an error of 0.026%.  

# Buffon's Needle  

Another approach for estimating \\( \\pi \\) that Wikipedia mentions is known as [Buffon's Needle](https://en.wikipedia.org/wiki/Buffon%27s_needle). Incidentally, this method is named after [Georges-Louis Leclerc, Comte de Buffon](https://en.wikipedia.org/wiki/Georges-Louis_Leclerc,_Comte_de_Buffon) who, among many other contributions, was one of the first to think about evolution scientifically prior to Darwin.  

<p align='center'>
  <img src='https://upload.wikimedia.org/wikipedia/commons/f/f6/Buffon_needle.gif'
  alt="Buffon's Needle" />  
</p>

This approach is more mathematical than the above circle method. It is based on the following scenario: given a needle (i.e. a line) of length \\( l \\) randomly placed onto a plane ruled with parallel lines a distance \\( t \\) apart (see diagram), the probability that the needle will cross a line, given that \\( l < t \\), is:  

$$
P(l, t) = \frac{2l}{t\pi}
$$

The derivation involves some calculus and is given [on Wikipedia](https://en.wikipedia.org/wiki/Buffon%27s_needle). Given this result, if \\( N \\) needles are dropped onto a ruled plane, and \\( n \\) intersect lines, then \\( \\pi \\) can be estimated as:  

$$
\pi \approx {\frac{2lN}{tn}}
$$

To simulate this experiment in R, I set \\( l = 0.5 \\) and \\( t = 1 \\), which simplifies the estimate of \\( \\pi \\) to the ratio \\( N / n \\). First, I generate a minimal ruled plane and a needle dropping function.  


```r
plane <- readWKT('POLYGON((-1 -1,-1 1,1 1,1 -1, -1 -1))')
ruled_lines <- readWKT('MULTILINESTRING((-0.5 -1, -0.5 1),(0.5 -1, 0.5 1))')
drop_needle <- function(n_needles, plane, l = 0.5) {
  midpoints <- spsample(plane, n_needles, type = 'random') %>% 
    coordinates
  angles <- 2 * pi * runif(n_needles)
  shift <- cbind(l * cos(angles) / 2, l * sin(angles) / 2)
  cbind(midpoints - shift, midpoints + shift) %>% 
    {dimnames(.) <- NULL; .} %>% 
    alply(1, function(x) {Line(rbind(x[1:2], x[3:4]))}) %>% 
    {SpatialLines(list(Lines(., ID = 'a')))} %>% 
    disaggregate
}
plot(plane, axes = T, border = '#F38630', col = 'transparent', lwd = 1, lty = 2)
plot(ruled_lines, axes = T, col = '#FA6900', lwd = 2, add = T)
drop_needle(50, plane, l = 0.5) %>% 
  plot(add = T, col = '#69D2E7', lwd = 1.5)
```

<img src="/figures//estimating-pi_needle-setup-1.svg" title="plot of chunk needle-setup" alt="plot of chunk needle-setup" style="display: block; margin: auto;" />

Creating these needles (i.e. random line segments) is messier than I'd ideally like, but it gets the job done for now. Based on these spatial objects, I estimate \\( \\pi \\) for 1000 needles.  


```r
n_needles <- 1000
needles <- drop_needle(n_needles, plane, l = 0.5)
(n_cross <- sum(gIntersects(needles, ruled_lines, byid = T)))
#> [1] 289
(pi_est <- n_needles / n_cross)
#> [1] 3.460208
round(abs(100 * (pi_est / pi - 1)), 2)
#> [1] 10.14
```

With 1000 needles I get an error of about 10.14%. So, this method is less precise than the circle method and it's also much slower due to the convoluted way I create the needles.  

## Convergence  

As before, I wrap all this into a function for estimating \\( \\pi \\).  


```r
estimate_pi_needle <- function(n_needles) {
  t <- system.time({
    plane <- readWKT('POLYGON((-1 -1,-1 1,1 1,1 -1, -1 -1))')
    ruled_lines <- readWKT('MULTILINESTRING((-0.5 -1, -0.5 1),(0.5 -1, 0.5 1))')
    needles <- drop_needle(n_needles, plane, l = 0.5)
    n_cross <- sum(gIntersects(needles, ruled_lines, byid = T))
    pi_estimate <- n_needles / n_cross
  })
  return(data.frame(pi_estimate, t = t[['user.self']]))
}
estimate_pi_needle(1000)
#>   pi_estimate     t
#> 1    3.154574 1.409
```

And estimate \\( \\pi \\) for a range of parameters.  


```r
n_values <- 10^(3:5)
n_reps <- 10
estimates_n <- expand.grid(n_values, 1:n_reps) %>% 
  setNames(c('n_pts', 'run')) %>% 
  tbl_df %>% 
  group_by(n_pts, run) %>% 
  do(estimate_pi_needle(.$n_pts))
```

Then I plot the resulting estimates.  


```r
estimate_n_summary <- estimates_n %>% 
  mutate(pct_error = abs((pi_estimate / pi - 1))) %>% 
  group_by(n_pts) %>% 
  summarize(pi_mean = mean(pi_estimate),
            pi_se = se(pi_estimate),
            error_mean = mean(pct_error),
            error_se = se(pct_error),
            t_mean = mean(t))
# estimate
ggplot(estimate_n_summary, aes(x = n_pts, y = pi_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = pi_mean - pi_se, ymax = pi_mean + pi_se), width = 0.1) +
  geom_line() +
  geom_hline(aes(yintercept=pi), color = 'red', linetype = 'dashed') +
  scale_x_log10() +
  labs(x = 'Number of Needles', y = expression(pi ~ estimate))
```

<img src="/figures//estimating-pi_pi-convergence-needle-1.svg" title="plot of chunk pi-convergence-needle" alt="plot of chunk pi-convergence-needle" style="display: block; margin: auto;" />

# Comparison of Methods  

Lastly, I combine the results from all three approaches for comparison.  


```r
comparison <- bind_rows(
  mutate(estimate_summary, method = 'Inscribed Circle'),
  mutate(estimate_d_summary, method = 'Distance'),
  mutate(estimate_n_summary, method = 'Needle'))

# estimate
ggplot(comparison, aes(x = n_pts, y = pi_mean, color = method)) +
  geom_point() +
  geom_errorbar(aes(ymin = pi_mean - pi_se, ymax = pi_mean + pi_se), width = 0.1) +
  geom_line() +
  geom_hline(aes(yintercept=pi), color = 'red', linetype = 'dashed') +
  scale_x_log10() +
  scale_color_brewer(name = 'Method', palette = 'Set1') +
  labs(x = 'Number of Needles', y = expression(pi ~ estimate)) +
  theme(legend.position=c(0.75, 0.75))
```

<img src="/figures//estimating-pi_comparison-1.svg" title="plot of chunk comparison" alt="plot of chunk comparison" style="display: block; margin: auto;" />

```r
# error
ggplot(comparison, aes(x = n_pts, y = error_mean, color = method)) +
  geom_point() +
  geom_errorbar(aes(ymin = error_mean - error_se, ymax = error_mean + error_se), width = 0.1) +
  geom_line() +
  scale_x_log10() +
  scale_color_brewer(name = 'Method', palette = 'Set1') +
  labs(x = 'Number of Points', y = 'Error') +
  theme(legend.position=c(0.75, 0.75))
```

<img src="/figures//estimating-pi_comparison-2.svg" title="plot of chunk comparison" alt="plot of chunk comparison" style="display: block; margin: auto;" />

```r
# times
ggplot(comparison, aes(x = n_pts, y = t_mean, color = method)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(labels = comma) +
  scale_color_brewer(name = 'Method', palette = 'Set1') +
  labs(x = 'Number of Points', y = 'Execution Time (s)') +
  theme(legend.position=c(0.75, 0.75))
```

<img src="/figures//estimating-pi_comparison-3.svg" title="plot of chunk comparison" alt="plot of chunk comparison" style="display: block; margin: auto;" />

Clearly the needle method is way less efficient!  
