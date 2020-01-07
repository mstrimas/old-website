---
layout: post
title: "Under eBirded counties: filling the gaps in eBird data"
published: true
excerpt: >
  Highlighting the least eBirded counties in the US by state.
category: r
tags: r spatial gis ebird
leaflet: true
---

After seeing my [previous post](/r/ebird-county/), one of my colleagues suggested it would be interesting to produce a map of the least eBirded counties in the US. Highlighting these counties might encourage more eBirders to visit them, thereby filling in the geographic gaps in eBird data. So, here it is, the counties within each state having the fewest checklists. I'll also throw on the county with the most checklists for comparison. At the suggestion of another colleague, I've normalized the number of checklists by county to account for variation in county size.


```r
library(dplyr)
library(stringr)
library(sf)
library(leaflet)
library(leaflet.extras)
library(here)

# load the data generated from the previous post
ebird_counties <- here("_source", "data", "ebird-county", 
                       "ebird-county_sf.rds") %>% 
  readRDS() %>% 
  # filter to counties in lower 48
  filter(str_detect(region_code, "^US"), !state %in% c("AK", "DC")) %>% 
  mutate(checklist_density = 1e6 * n_checklists / area)

# add popups
popup_maker <- function(code, name, n_species, n_checklists, 
                        checklist_density) {
  link_style <- "target='_blank' style='text-decoration: none;'"
  region_url <- paste0(
    "<strong><a href='http://ebird.org/ebird/subnational2/",
    code, "?yr=all' ", link_style, ">", name, "</a></strong>")
  counts <- paste0("<strong>", format(n_checklists, big.mark = ","), 
                   " checklists</strong> (", round(checklist_density, 3), 
                   " per km<sup>2</sup>) | ", 
                   n_species, " species")
  paste(region_url, counts, sep = "<br/>")
}
ebird_counties <- ebird_counties %>% 
  mutate(region_name = paste(region_name, state, sep = ", "),
         popup = popup_maker(region_code, region_name, 
                             n_species, n_checklists, checklist_density)) 

# bottom county by state
ebird_counties_least <- ebird_counties %>% 
  group_by(state) %>% 
  top_n(-1, checklist_density) %>% 
  ungroup() %>% 
  mutate(label = "Bottom counties", grp = "Bottom county") %>% 
  select(popup, label, grp)
# bottom 5% by state
ebird_counties_b5 <- ebird_counties %>% 
  group_by(state) %>% 
  filter(checklist_density <= quantile(checklist_density, 0.05)) %>% 
  ungroup() %>% 
  mutate(label = "Bottom counties", grp = "Bottom 5%") %>% 
  select(popup, label, grp)
# bottom 10% by state
ebird_counties_b10 <- ebird_counties %>% 
  group_by(state) %>% 
  filter(checklist_density <= quantile(checklist_density, 0.1)) %>% 
  ungroup() %>% 
  mutate(label = "Bottom counties", grp = "Bottom 10%") %>% 
  select(popup, label, grp)
# most checklists by state
ebird_counties_most <- ebird_counties %>% 
  group_by(state) %>% 
  top_n(1, checklist_density) %>% 
  ungroup() %>% 
  mutate(label = "Top county", grp = "Top county") %>% 
  select(popup, label, grp)

# leaflet map
p <- colorFactor(c("#e41a1c", "#4daf4a"), 
                 domain = c("Bottom counties", "Top county"))
leaflet() %>%
  addProviderTiles("OpenMapSurfer.Roads") %>% 
  addFullscreenControl() %>% 
  addPolygons(data = ebird_counties_least, group = "Bottom county",
              color = ~ p(label), opacity = 1.0,
              weight = 2, smoothFactor = 0.5,
              fillColor = "#555555", fillOpacity = 0.2,
              popup = ~ popup) %>% 
  addPolygons(data = ebird_counties_b5, group = "Bottom 5%",
              color = ~ p(label), opacity = 1.0,
              weight = 2, smoothFactor = 0.5,
              fillColor = "#555555", fillOpacity = 0.2,
              popup = ~ popup) %>% 
  addPolygons(data = ebird_counties_b10, group = "Bottom 10%",
              color = ~ p(label), opacity = 1.0,
              weight = 2, smoothFactor = 0.5,
              fillColor = "#555555", fillOpacity = 0.2,
              popup = ~ popup) %>% 
  addPolygons(data = ebird_counties_most, group = "Top county",
              color = ~ p(label), opacity = 1.0,
              weight = 2, smoothFactor = 0.5,
              fillColor = "#555555", fillOpacity = 0.2,
              popup = ~ popup) %>% 
  addLegend("bottomright", pal = p,
            values = c("Bottom counties", "Top county"),
            title = "Counties with highest/lowest eBird checklist density",
            opacity = 1) %>% 
  addLayersControl(
    baseGroups = c("Bottom county", "Bottom 5%", "Bottom 10%", "Top county"),
    overlayGroups = "Top county",
    options = layersControlOptions(collapsed = FALSE))
```



<iframe src="/assets/leaflet/under-ebirded_density.html" style="border: none; width: 800px; height: 600px"></iframe>
<a href="/assets/leaflet/under-ebirded_density.html" target="_blank"><strong>Fullscreen</strong></a> | Data source: [eBird](http://ebird.org/)

If you're interested in the unnormalized version, ranking counties by the raw number of checklists, take a [look at this map](/assets/leaflet/under-ebirded.html). More importantly, if you're a birder and live near one of these counties, get out there and submit some checklists, you might find something interesting!
