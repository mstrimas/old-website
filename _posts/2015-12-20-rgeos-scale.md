---
layout: post
title: Scale and numerical precision in RGEOS
published: true
excerpt: >
  Exploring the unexpected results that can arise in GEOS topology operations 
  from problems with numerical precision or registration in coordinates.
category: spatial
tags: r spatial gis
---



In this post, I explore the unexpected results that can arise in GEOS topology operations from problems with numerical precision or registration in coordinates. It is inspired by a discussion on [this R-sig-Geo thread](http://r-sig-geo.2731867.n2.nabble.com/gUnaryUnion-Not-Dissolving-Correctly-td7589145.html), especially the comments provided by Roger Bivand.  

Many of these issues arise from challenges associated with comparing floating point numbers, a topic discussed in [R FAQ 7.31](https://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-doesn_0027t-R-think-these-numbers-are-equal_003f). It is noted there that "the only numbers that can be represented exactly in R’s numeric type are integers...as a result, two floating point numbers will not reliably be equal unless they have been computed by the same algorithm". They provide the following example to illustrate this:  


```r
a <- sqrt(2)
a * a == 2
#> [1] FALSE
a * a - 2
#> [1] 4.440892e-16
```

Similar issues can arise when comparing spatial coordinates stored as floating point numbers. To facilitate comparisons, GEOS, and consequently `rgeos`, shifts all coordinates to an integer grid after multiplying by a scale factor. This scale factor determines the precision at which differences in the coordinates will be detected.  

## Example Data  

Load packages required for this document. The `maptools` packages provides the `elide()` function, which can be used to translate the coordinates of spatial objects.  


```r
library(sp)
library(rgeos)
library(maptools)
# the effects of setting these rgeos parameters will be explored later
set_RGEOS_polyThreshold(0)
set_RGEOS_warnSlivers(TRUE)
set_RGEOS_dropSlivers(FALSE)
```

Generate two adjacent unit squares and demonstrate the effects of `elide()`.  


```r
p1 <- readWKT("POLYGON((0 0,0 1,1 1,1 0,0 0))")
row.names(p1) <- '1'
p2 <- readWKT("POLYGON((1 0,1 1,2 1,2 0,1 0))")
row.names(p2) <- '2'

plot(p1, col='red', axes=T, xlim=c(0,2.5), ylim=c(0, 1))
plot(p2, col='blue', add=T)
plot(elide(p2, shift=c(0.25, 0)), border='orange', add=T, lty=2, lwd=2)
```

<img src="/figures//rgeos-scale_example-data-1.svg" title="plot of chunk example-data" alt="plot of chunk example-data" style="display: block; margin: auto;" />

These `SpatialPolygons` objects share a boundary whose coordinates are identical. In what follows, the polygons will be translated by small amounts to examine how slight coordinate mismatches affect the topology operations of `rgeos`.  

## `rgeos` Scale Factor  

By default the `rgeos` scale factor is \\( 10^8 \\), corresponding to a precision of \\( 10^{-8} \\), and points will be treated as equal if differences in their coordinates are smaller than this precision. The `rgeos` functions `setScale()` and `getScale()` can be used to change the scale or return the current scale, respectively.  


```r
setScale(1e8)
getScale()
#> [1] 1e+08
```

### Test of intersection  

`gIntersects(x,y)` tests if geometries `x` and `y` have at least one point in common. Since the two square polygons under consideration share a boundary they do intersect. If the right square is shifted further to the right, the squares no longer overlap; however, `gIntersects()` only considers the geometries non-intersecting if this shift is detectable within the precision set by the scale factor.  


```r
gIntersects(p1, p2)
#> [1] TRUE
gIntersects(p1, elide(p2, shift=c(0.1, 0))) # shift > precision => FALSE
#> [1] FALSE
gIntersects(p1, elide(p2, shift=c(1e-8, 0))) # shift = precision => FALSE
#> [1] FALSE
gIntersects(p1, elide(p2, shift=c(1e-9, 0))) # shift < precision => TRUE
#> [1] TRUE
```

Note that in the final case the coordinates only differ by \\( 10^{-9} \\), less than the precision (\\( 10^{8} \\)), hence GEOS treats these coordinates as equal.  

### Union  

`gUnion(x,y)` merges geometries `x` and `y` if they are intersecting. As above, moving the polygons apart by more than the precision results in two non-intersecting polygons that are not merged.  


```r
plot(gUnion(p1, elide(p2, shift=c(0.1, 0))), col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_gunion-nomerge-1.svg" title="plot of chunk gunion-nomerge" alt="plot of chunk gunion-nomerge" style="display: block; margin: auto;" />

```r
plot(gUnion(p1, elide(p2, shift=c(1e-8, 0))), col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_gunion-nomerge-2.svg" title="plot of chunk gunion-nomerge" alt="plot of chunk gunion-nomerge" style="display: block; margin: auto;" />

However, if the difference in coordinates is too small to detect at the given precision, then the polygons are considered to be intersecting and are merged.  


```r
plot(gUnion(p1, elide(p2, shift=c(1e-9, 0))), col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_gunion-merge-1.svg" title="plot of chunk gunion-merge" alt="plot of chunk gunion-merge" style="display: block; margin: auto;" />

### Intersection  

`gIntersection(x,y)` returns the intersection between the geometries `x` and `y`, or `NULL` if the geometries do not intersect. Since the original polygons share an edge, their intersection is a line.  


```r
class(gIntersection(p1, p2))
#> [1] "SpatialLines"
#> attr(,"package")
#> [1] "sp"
plot(p1, col='lightblue', border='transparent', axes=T, xlim=c(0,2), 
     ylim=c(0, 1))
plot(p2, col='transparent', border='black', add=T)
plot(gIntersection(p1, p2), col='orange', add=T, lty=2, lwd=3)
```

<img src="/figures//rgeos-scale_gi-line-1.svg" title="plot of chunk gi-line" alt="plot of chunk gi-line" style="display: block; margin: auto;" />

If the two polygons overlap, then a polygon geometry results from their intersection.   


```r
p <- elide(p2, shift=c(-0.1, 0))
class(gIntersection(p1, p))
#> [1] "SpatialPolygons"
#> attr(,"package")
#> [1] "sp"
plot(p1, col='lightblue', border='transparent', axes=T, xlim=c(0,2), 
     ylim=c(0, 1))
plot(p, col='transparent', border='black', add=T)
plot(gIntersection(p1, p), border='orange', add=T, lty=2, lwd=3)
```

<img src="/figures//rgeos-scale_gi-poly-1.svg" title="plot of chunk gi-poly" alt="plot of chunk gi-poly" style="display: block; margin: auto;" />

However, if the amount of overlap is small enough that the difference in coordinates is below the precision (\\( 10^{8} \\)), then the polygons are treated as having a shared boundary and the intersection returns a line.  


```r
p <- elide(p2, shift=c(-1e-8, 0)) # overlap > precision => polygon
class(gIntersection(p1, p)) 
#> [1] "SpatialPolygons"
#> attr(,"package")
#> [1] "sp"
p <- elide(p2, shift=c(-1e-9, 0)) # overlap < precision => line
class(gIntersection(p1, p)) 
#> [1] "SpatialLines"
#> attr(,"package")
#> [1] "sp"
```

In contrast, if the two polygons are separated, such that they no longer intersect, the intersection operation returns `NULL`.  


```r
p <- elide(p2, shift=c(0.1, 0))
class(gIntersection(p1, p))
#> [1] "NULL"
plot(p1, col='lightblue', border='transparent', axes=T, xlim=c(0,2.1), 
     ylim=c(0, 1))
plot(p, col='transparent', border='black', add=T)
```

<img src="/figures//rgeos-scale_gi-null-1.svg" title="plot of chunk gi-null" alt="plot of chunk gi-null" style="display: block; margin: auto;" />

Again, the difference between coordinates relative to the precision determines the outcome of the intersection.  


```r
p <- elide(p2, shift=c(1e-8, 0)) # separation > precision => no overlap => NULL
class(gIntersection(p1, p)) 
#> [1] "NULL"
p <- elide(p2, shift=c(1e-9, 0)) # separation < precision => shared edge => line
class(gIntersection(p1, p)) 
#> [1] "SpatialLines"
#> attr(,"package")
#> [1] "sp"
```

### Change of scale factor  

Changing the scale to \\( 10^4 \\), so that precision is lower, shows how the behaviour demonstrated above is influenced by the scale factor.  


```r
setScale(1e4)
gIntersects(p1, elide(p2, shift=c(1e-4, 0)))
#> [1] FALSE
gIntersects(p1, elide(p2, shift=c(1e-5, 0)))
#> [1] TRUE

plot(gUnion(p1, elide(p2, shift=c(1e-4, 0))), col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_low-precision-1.svg" title="plot of chunk low-precision" alt="plot of chunk low-precision" style="display: block; margin: auto;" />

```r
plot(gUnion(p1, elide(p2, shift=c(1e-5, 0))), col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_low-precision-2.svg" title="plot of chunk low-precision" alt="plot of chunk low-precision" style="display: block; margin: auto;" />

```r

class(gIntersection(p1, elide(p2, shift=c(-1e-4, 0)))) # overlap > precision => polygon
#> [1] "SpatialPolygons"
#> attr(,"package")
#> [1] "sp"
class(gIntersection(p1, elide(p2, shift=c(-1e-5, 0)))) # overlap < precision => line
#> [1] "SpatialLines"
#> attr(,"package")
#> [1] "sp"
class(gIntersection(p1, elide(p2, shift=c(1e-5, 0)))) # separation < precision => line
#> [1] "SpatialLines"
#> attr(,"package")
#> [1] "sp"
class(gIntersection(p1, elide(p2, shift=c(1e-4, 0)))) # separation > precision => NULL
#> [1] "NULL"
```

## Slivers  

These issues of numerical precision or mis-registration can lead to slivers being generated by the topology operations. These small area polygons are artifacts of the computational geometry. The `rgeos` function `set_RGEOS_polyThreshold()` can be used to set an area threshold for valid polygons, such that any polygon whose area is below this threshold will be considered a sliver. The threshold is 0 by default, which will ignore slivers, but can be set to a small number to detect slivers.  


```r
setScale(1e4)
set_RGEOS_polyThreshold(1e-2)
```

Furthermore, to report on detected slivers `warnSlivers` must be set to `TRUE` using `set_RGEOS_warnSlivers(TRUE)`. By default, this value is `FALSE`.  


```r
set_RGEOS_warnSlivers(TRUE)
```

Now, polygons with small, but increasing amounts of overlap are intersected with sliver detection and reporting turned on. With the current parameters, precision is set to \\( 10^{-4} \\) and `polyThreshold` to \\( 10^{-2} \\).  

For overlap below the level of detection with the current precision, a line is returned and no sliver warnings are raised.  


```r
class(gIntersection(p1, elide(p2, shift=c(-1e-5, 0))))
```

```
[1] "SpatialLines"
attr(,"package")
[1] "sp"
```

For overlap at or above the current precision, a polygon is returned; however, if the area of this polygon is below the `polyThreshold` a sliver warning is raised.  

```r
gi <- gIntersection(p1, elide(p2, shift=c(-1e-4, 0)))
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : 1: Polygon object 1 area 0.0001
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : Exterior ring 0 of object 1 area 0.0001
class(gi)
#> [1] "SpatialPolygons"
#> attr(,"package")
#> [1] "sp"
gArea(gi)
#> [1] 1e-04
gi <- gIntersection(p1, elide(p2, shift=c(-1e-3, 0)))
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : 1: Polygon object 1 area 0.001
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : Exterior ring 0 of object 1 area 0.001
class(gi)
#> [1] "SpatialPolygons"
#> attr(,"package")
#> [1] "sp"
gArea(gi)
#> [1] 0.001
```

And, with sufficient overlap, the resulting polygon will have an area greater than the `polyThreshold`. In this case, the polygon will no longer be considered a sliver and no warning will be given.  


```r
gi <- gIntersection(p1, elide(p2, shift=c(-1e-2, 0)))
class(gi)
#> [1] "SpatialPolygons"
#> attr(,"package")
#> [1] "sp"
gArea(gi)
#> [1] 0.01
```

### Change of threshold


```r
set_RGEOS_polyThreshold(1e-3)
```

With the threshold lowered from \\( 10^{-2} \\) to \\( 10^{-3} \\), `rgeos` is less sensitive to slivers. A horizontal shift of \\( 10^{-4} \\) still raises a sliver warning.  


```r
gi <- gIntersection(p1, elide(p2, shift=c(-1e-4, 0)))
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : 1: Polygon object 1 area 0.0001
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : Exterior ring 0 of object 1 area 0.0001
class(gi)
#> [1] "SpatialPolygons"
#> attr(,"package")
#> [1] "sp"
gArea(gi)
#> [1] 1e-04
```

However, with a shift of \\( 10^{-3} \\), `rgoes` no longer treats the resulting polygon as a sliver since the area is above the new threshold. No warning is raised.  


```r
gi <- gIntersection(p1, elide(p2, shift=c(-1e-3, 0)))
class(gi)
#> [1] "SpatialPolygons"
#> attr(,"package")
#> [1] "sp"
gArea(gi)
#> [1] 0.001
```

### Threshold is area based  

Note that it isn't the linear overlap that triggers the warning, it is that the area of the resulting polygons are below the threshold. In the above examples the size of the linear shift was equal to the area of the resulting polygon because the original geometries were unit squares. 


```r
gi1 <- gIntersection(p1, elide(p2, shift=c(-1e-3, 0)))
gArea(gi1)
#> [1] 0.001
gArea(gi1) / get_RGEOS_polyThreshold()
#> [1] 1
```

However, this need not be the case. Now a warning is raised because a slight shift in the vertical direction has caused the polygon resulting from the intersection to have area just less than the \\( 10^{-3} \\) threshold.  


```r
gi2 <- gIntersection(p1, elide(p2, shift=c(-1e-3, -1e-3)))
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : 1: Polygon object 1 area 0.000999
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : Exterior ring 0 of object 1 area 0.000999
gArea(gi2)
#> [1] 0.000999
gArea(gi2) / get_RGEOS_polyThreshold()
#> [1] 0.999
```

### Dropping slivers  

`rgeos` can also be set to automatically drop slivers resulting from topology operations. This is accomplished by setting `dropSlivers` to `TRUE` using `set_RGEOS_dropSlivers()`. By default `dropSlivers` is `FALSE`.  

In the next example, the intersection yields two polygons: one valid and one sliver resulting from a slight misalignment.  


```r
p_a <- rbind(p1, elide(p2, shift=c(0.5, 0)))
p_b <- elide(p2, shift=c(-1e-3, 0))
plot(p_a, col='lightgrey', axes=T, xlim=c(0,2.5), ylim=c(0,1))
plot(p_b, border='orange', add=T, lty=2, lwd=3)
```

<img src="/figures//rgeos-scale_drop-slivers-ex-1.svg" title="plot of chunk drop-slivers-ex" alt="plot of chunk drop-slivers-ex" style="display: block; margin: auto;" />

With `dropSlivers` set to `FALSE`, both are returned.  


```r
set_RGEOS_polyThreshold(1e-2)
set_RGEOS_dropSlivers(FALSE)
gi <- gIntersection(p_a, p_b, byid=T)
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : 1: Polygon object 1 2 area 0.001
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : Exterior ring 0 of object 1 2 area 0.001
gArea(gi, byid=T)
#>   1 2   2 2 
#> 0.001 0.499
```

However, `dropSlivers` set to `TRUE`, the small area sliver is removed from the resulting geometry.  


```r
set_RGEOS_dropSlivers(TRUE)
gi <- gIntersection(p_a, p_b, byid=T)
#> Warning in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td,
#> unaryUnion_if_byid_false, : 1: Polygon object 1 2 area 0.001
gArea(gi, byid=T)
#>   2 2 
#> 0.499
set_RGEOS_dropSlivers(FALSE)
```

### Slivers from union operations  

Slivers can also arise as a result of union operations.  


```r
setScale(1e4)
set_RGEOS_polyThreshold(1e-2)
set_RGEOS_warnSlivers(TRUE)
set_RGEOS_dropSlivers(FALSE)

p3 <- readWKT("POLYGON((0 1,0 2,2 2,2 1,0 1))")
row.names(p3) <- '3'
p4 <- readWKT("POLYGON((0 -1,0 0,2 0,2 -1,0 -1))")
row.names(p4) <- '4'
plot(rbind(p1, p3, p4), axes=T)
plot(p2, add=T, col='red')
```

<img src="/figures//rgeos-scale_reset-params-1.svg" title="plot of chunk reset-params" alt="plot of chunk reset-params" style="display: block; margin: auto;" />

Now the the middle right (i.e. red) square is shifted to the right by increasing amounts. If the shift is below the precision, the misalignment of the middle edge is not picked up.  


```r
pshift <- elide(p2, shift=c(1e-5, 0))
pp <- rbind(p1, p3, p4, pshift)
guu <- gUnaryUnion(pp)
plot(guu, col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_sliver-union-no-1.svg" title="plot of chunk sliver-union-no" alt="plot of chunk sliver-union-no" style="display: block; margin: auto;" />

For a shift within the limits of precision, the misalignment of the middle edge is picked up and a very narrow hole appears in the resulting geometry. A warning is raised since this interior ring has area below the `polyThreshold`.    


```r
pshift <- elide(p2, shift=c(1e-4, 0))
pp <- rbind(p1, p3, p4, pshift)
guu <- gUnaryUnion(pp)
#> Warning: Interior ring 0 of Polygon 0 of object 1 area 0.0001
plot(guu, col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_sliver-union-warn-1.svg" title="plot of chunk sliver-union-warn" alt="plot of chunk sliver-union-warn" style="display: block; margin: auto;" />

However, for a larger shift, the hole persists, but no warning is raised since the area is now above the `polyThreshold`.  


```r
pshift <- elide(p2, shift=c(1e-2, 0))
pp <- rbind(p1, p3, p4, pshift)
plot(gUnaryUnion(pp), col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_sliver-union-no-warn-1.svg" title="plot of chunk sliver-union-no-warn" alt="plot of chunk sliver-union-no-warn" style="display: block; margin: auto;" />

The fact that this is a hole and not a vertical line becomes apparent when the shift is larger.  


```r
pshift <- elide(p2, shift=c(0.1, 0))
pp <- rbind(p1, p3, p4, pshift)
plot(gUnaryUnion(pp), col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_sliver-union-big-1.svg" title="plot of chunk sliver-union-big" alt="plot of chunk sliver-union-big" style="display: block; margin: auto;" />

Finally, `set_RGEOS_dropSlivers()` can be used to repair the geometry by removing these interior slivers.  


```r
set_RGEOS_dropSlivers(TRUE)
pshift <- elide(p2, shift=c(1e-4, 0))
pp <- rbind(p1, p3, p4, pshift)
guu <- gUnaryUnion(pp)
#> Warning: Interior ring 0 of Polygon 0 of object 1 area 0.0001
plot(guu, col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_sliver-union-drop-1.svg" title="plot of chunk sliver-union-drop" alt="plot of chunk sliver-union-drop" style="display: block; margin: auto;" />

### Inward dangles  

It is possible that an inward dangle (a zero area line in from a edge) will escape detection even when `warnSlivers` or `dropSlivers` are `TRUE`.  


```r
set_RGEOS_warnSlivers(TRUE)
set_RGEOS_dropSlivers(TRUE)
pshift <- elide(p2, shift=c(1e-4, 0))
pp <- rbind(p1, p4, pshift)
guu <- gUnaryUnion(pp)
plot(guu, col='lightgrey', axes=T)
```

<img src="/figures//rgeos-scale_dangle-1.svg" title="plot of chunk dangle" alt="plot of chunk dangle" style="display: block; margin: auto;" />

Note that this dangle has no impact on the area of the resulting geometry, suggesting that it has zero area itself, which explains how it escapes detection.  


```r
cat(gArea(guu, byid=T))
#> 4
```
