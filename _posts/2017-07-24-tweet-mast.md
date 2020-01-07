---
layout: post
title: "Twitter citizen science: mapping spruce masting events"
published: true
excerpt: >
  In this post, I the use the Twitter API to access geocoded tweets on spruce 
  masting events in North America. Based on this crowd sourced data, I produce a
  map of spruce masting.
category: r
tags: r gis twitter
leaflet: true
---

One of the fondest memories of my life was spending the summer of 2012 and winter of 2013 working in the Yukon on the [Kluane Red Squirrel Project (KRSP)](http://redsquirrel.biology.ualberta.ca/). KRSP is a long-term, large-scale field experiment studying the effects of food availability on the evolution and ecology of red squirrels. While working on the project, I spent each day outside, walking around the boreal forest monitoring red squirrels on one of the project's study grids. The location was beautiful, empty, and wild, I came to love the fiesty little red squirrels, and I worked with a great crew of fellow squirrelers!

<div style="text-align: center">
  <img src="/img/twitter-mast/kluane.jpg"/>
</div>

The squirrels mostly feed on the cones produced by white spruce trees, and this system is particularly interesting because spruce trees exhibit masting behaviour: they produce almost no cones in most years, but every few years they produce huge numbers of cones. This masting is often synchronized over large areas. The result is that squirrels experience massive inter-annual fluctuations in the availability of food.

One of the KRSP PIs, [Andrew McAdam](http://www.mcadamlab.ca/), recently had a cool idea: use Twitter to crowd source data on where spruce trees are masting this year. He suggested using the hashtags [#SpruceMast2017Yes ](https://twitter.com/hashtag/SpruceMast2017Yes?src=hash) and [#SpruceMast2017No ](https://twitter.com/hashtag/SpruceMast2017No?src=hash), combined with geo tagging tweets, to identify whether spruce trees are masting or not in your area. I think this is a cool idea for quickly engaging citizen scientists.

In this post, I'll use R to extract data from the Twitter API on spruce masting, then make a couple of quick maps of the results.

## Required packages


```r
library(dplyr)
library(ggplot2)
library(stringr)
library(rtweet)
library(sf)
library(rnaturalearth)
library(ggmap)
library(leaflet)
```

## Extract Twitter data

The R package `rtweet` provides access to the Twitter API, making it easy to grab tweets based on a search query, such as a hashtag. Fortunately, the data frame returned by `search_tweets()` contains the locataion assigned to the tweet in the form of the city (e.g. Guelph, ON). I then geocode these locations using `ggmap::geocode()` to turn these city names into coordinates.

Note that there is some initial set up here to get a Twitter API key and store it as an R environment variable. This is covered nicely in the `rtweet` [documentation](http://rtweet.info/articles/auth.html).


```r
# masting
tweet_mast_yes <- search_tweets("#SpruceMast2017Yes") %>% 
  # exclude retweets
  filter(is.na(retweet_status_id)) %>% 
  # exlude tweets without location
  filter(!is.na(place_full_name)) %>% 
  # exlcude tweets with no and yes
  filter(!str_detect(hashtags, "SpruceMast2017No")) %>% 
  select(screen_name, status_id, text, place = place_full_name) %>% 
  mutate(mast = "Yes") %>% 
  bind_cols(geocode(.$place))
# not masting
tweet_mast_no <- search_tweets("#SpruceMast2017No") %>% 
  filter(is.na(retweet_status_id)) %>% 
  filter(!is.na(place_full_name)) %>% 
  filter(!str_detect(hashtags, "SpruceMast2017Yes")) %>% 
  select(screen_name, status_id, text, place = place_full_name) %>% 
  mutate(mast = "No") %>% 
  bind_cols(geocode(.$place))
# comine
tweet_mast <- bind_rows(tweet_mast_yes, tweet_mast_no)
```

Next I use the latitude and longitude to turn these data into a spatial object with the `sf` package.


```r
tweet_mast <- st_as_sf(tweet_mast, coords = c("lon", "lat"), crs = 4326)
```

## Map 

Now that I have the data prepared, I put them on a couple maps, one static and one interactive. First I grab boundaries of North American States and provinces.


```r
na_boundaries <- ne_states(country = c("united states of america", "canada"),
                           returnclass = "sf") %>% 
  filter(name != "Hawaii")
```

### Static map (ggplot2)

Here I'll use the new `ggplot2` functionality for plotting `sf` objects. Note that there is currently [an issue](https://github.com/tidyverse/ggplot2/issues/2037) when plotting points, which prevents including a legend in the plot. I've used green for locations where masting has been observed, and red for locations where it hasn't, which I think is intuitive without the legend.


```r
proj <- paste0("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 ",
               "+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
# dummy data frame for legend
ggplot() +
  geom_point() +
  geom_sf(data = na_boundaries, size = 0.1) +
  geom_sf(data = tweet_mast, aes(color = mast), size = 0.1,
          show.legend = FALSE) +
  scale_color_manual(values = c("No" = "#e41a1c",  "Yes" = "#4daf4a")) +
  coord_sf(crs = proj, xlim = c(-2700000, 3000000), ylim = c(0, 4000000)) +
  labs(title = "Spruce masting events reported on Twitter (2017)",
       subtitle = paste("Geotagged tweets with hashtags #SpruceMast2017Yes",
                        "and #SpruceMast2017No")) +
  theme(panel.background = element_rect(fill = "lightblue"), 
        panel.grid.major = element_line(color = "white", size = 0.5, 
                                        linetype = 2),
        axis.ticks = element_blank())
```

<a href="/figures//tweet-mast_gg-1.png"><img src="/figures//tweet-mast_gg-1.png" title="plot of chunk gg" alt="plot of chunk gg" style="display: block; margin: auto;" /></a>

### Interactive map (Leaflet)

Finally, I use the `leaflet` package to make an interactive map of the same data. First, I write a function that generates HTML for the popups.


```r
popup_html <- function(text, place) {
  text <- stringr::str_replace_all(text, "\\s+", " ")
  paste("<strong>Tweet:</strong>", text, "<br />", 
        "<strong>Location:</strong>", place)
}
tweet_mast <- mutate(tweet_mast, popup = popup_html(text, place))
```

Now I produce the map.


```r
pal <- colorFactor(c("#e41a1c", "#4daf4a"), domain = c("No", "Yes"))
leaflet(tweet_mast) %>% 
  addProviderTiles("OpenMapSurfer.Roads") %>% 
  addCircleMarkers(color = ~ pal(mast), radius = 5, popup = ~ popup) %>% 
  addLegend(pal = pal, values = ~ mast, opacity = 1,
            title = "Mast observed?")
```



<iframe src="/assets/leaflet/tweet-mast.html" style="border: none; width: 800px; height: 500px"></iframe>

The fullscreen map is viewable <a href="/assets/leaflet/tweet-mast.html" target="_blank"><strong>here</strong></a>.
