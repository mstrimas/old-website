---
layout: post
title: "Using the new Equal Earth projection in R"
published: true
excerpt: >
  A short post demonstrating how to use the new Equal Earth projection in R.
category: gis
tags: r gis 
editor_options: 
  chunk_output_type: console
---

The [Equal Earth map projection](http://equal-earth.com/equal-earth-projection.html) is a new global equal-area projection developed by Bojan Šavrič, Tom Patterson, and Bernhard Jenny. It's been designed to be aesthetically pleasing, with a shape similar to the Robinson projection, while also showing all countries, north or south, correctly sized. It's been blowing up on Twitter, and has been implemented in [PROJ4](https://proj4.org/operations/projections/eqearth.html), but I haven't been able to get on the bandwagon. Until now!

To get this working in R, start by upgrading PROJ4 to version 5.2.0 using homebrew with `brew upgrade`, then update the `sf` package making sure to compile it from source with `install.package("sf", type = "source")`. Now, when you load `sf` you should see that it's linking to `proj 5.2.0`. Unfortunately, although PROJ4 has Equal Earth support, the current version of GDAL on homebrew doesn't, so using `st_transform(crs = "+proj=eqearth")` won't work. I found the solution in [this GitHub issue](https://github.com/OSGeo/gdal/issues/870): use `"+proj=eqearth +wktext"` to pass the coordinates unmodified to PROJ. Hopefully GDAL will be updated soon, making this work around obsolute, but if you want to get on this Equal Earth train now, here's a simple example.

First, to plot a world map with `sf`, including country borders and graticules from [Natural Earth](https://www.naturalearthdata.com/):


```r
library(sf)
library(rnaturalearth)

# natural earth data
countries <- ne_countries(returnclass = "sf") %>% 
  st_transform(crs = "+proj=eqearth +wktext") %>% 
  st_geometry()
graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>% 
  st_transform(crs = "+proj=eqearth +wktext") %>% 
  st_geometry()
bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>% 
  st_transform(crs = "+proj=eqearth +wktext") %>% 
  st_geometry()

# sf plot
par(mar = c(0, 0, 0, 0))
plot(graticules, col = "grey20", lwd = 0.3)
plot(bb, col = NA, border = "grey20", lwd = 1.2, add = TRUE)
plot(countries, col = "grey80", border = "grey40", lwd = 0.8, add = TRUE)
```

<img src="/figures//equal-earth_sf-1.png" title="plot of chunk sf" alt="plot of chunk sf" style="display: block; margin: auto;" />

To plot a similar map with `ggplot2`:


```r
library(rnaturalearth)
library(ggplot2)

countries <- ne_countries(returnclass = "sf")
graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf")
bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf")
ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data = graticules, col = "grey20", lwd = 0.1) +
  geom_sf(data = countries, fill = "grey80", col = "grey40", lwd = 0.3) +
  coord_sf(crs = "+proj=eqearth +wktext") +
  theme_minimal() +
  theme(axis.text = element_blank())
```

<img src="/figures//equal-earth_ggplot-1.png" title="plot of chunk ggplot" alt="plot of chunk ggplot" style="display: block; margin: auto;" />

