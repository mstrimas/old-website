---
layout: post
title: "County eBirding: web scraping and web mapping in R"
published: true
excerpt: >
  Scraping the ebird website to find the top hotspot in each county. Covers 
  scraping data from websites with rvest, manipulating spatial data with sf, and 
  making interactive maps with leaflet.
category: r
tags: r spatial gis ebird
leaflet: true
---

[eBird](http://ebird.org) is my go to tool for keeping track of my sightings when I'm out birding. The eBird website is also a great tool for deciding where to go birding. You can find out where particular target birds have been seen recently, or see what other birders have seen in different geographical regions or birding hotspots. Recently, when I've been on a road trip and I pass into a new county, I'll try to stop a least once to go birding. Usually I'll choose a location by looking at nearby [eBird hotspots](https://ebird.org/ebird/hotspots) and picking the one with the most species observed. Part of the fun of this is filling in as many counties as possible on my [eBird profile](https://ebird.org/ebird/profile/MTU0MjQz/US-NY) for my home state, New York.

<div style="text-align: center">
  <img src="/img/ebird-county/county-map.png" width = "500"/>
</div>

With this county birding in mind, I decided to make a map of the top eBird hotspot in each county in the US. As always, this is really just an excuse to mess around in R, and this post will cover scraping data from websites with `rvest` and making interactive web maps with `leaflet`. I'll also be using the new R package, `sf`, for working with spatial data.

## Required packages


```r
library(tidyverse)
library(rvest)
library(stringr)
library(sf)
library(leaflet)
library(leaflet.extras)
library(USAboundaries)
library(rnaturalearth)
library(viridis)
library(here)
```

## North American regional stats

For each country in North America, I drill down to the state and county level and scrape eBird to get the number of species seen and the total number of checklists submitted.

### Country data

Country data needs to be treated separately, since there is no North America page summarizing the countries, so I manually extract it here. `rvest` makes scraping this sort of data from websites almost trivial. `read_html()` is used to read the HTML source of a page and `html_nodes()` is used extract specific components of the page using CSS selectors. The only slightly tricky part of scraping pages is figuring out what [CSS selector](https://www.w3schools.com/cssref/css_selectors.asp) to use to target a given page component. I typically use Chrome's *Inspect* tool; however, you can also use the [SelectorGadget](http://selectorgadget.com/) Chrome plugin.


```r
ebird_regions <- data_frame(
  region_code = c("CA", "US", "MX"),
  region_name = c("Canada", "United States", "Mexico"),
  region_level = "country")
# function to get number of counts and checklists per country
country_data <- function(country) {
  # species and checklist  counts
  base_url <- "http://ebird.org/ebird/country/%s?yr=all"
  page <- sprintf(base_url, country) %>% 
    read_html()
  counts <- html_nodes(page, ".hs-section-count") %>% 
    html_text() %>% 
    str_extract("[0-9,]+") %>% 
    parse_number()
  data_frame(n_species = counts[1], n_checklists = counts[2])
}
# apply this to our country list
ebird_regions <- ebird_regions %>% 
  mutate(country_df = map(region_code, country_data)) %>% 
  unnest()
```

### Extraction function

Now I define a function that extracts all the sub-region data within a region. So I don't hit the eBird website with a bunch of requests all at once, I've used `Sys.sleep()` to pause between each page load.


```r
extract_subregion_data <- function(region_code, region_level, sleep = 5) {
  Sys.sleep(sleep)
  # species and checklist counts
  base_url <- "http://ebird.org/ebird/%s/%s/regions?yr=all"
  page <- sprintf(base_url, region_level, region_code) %>% 
    read_html()
  region_stats <- html_table(page)
  # skip this region if no data at sub-region
  if (length(region_stats) == 0) {
    return(data_frame())
  }
  region_stats <- region_stats[[1]] %>% 
    set_names(c("rank", "region_name", "n_species", "n_checklists"))
  # region codes
  region_urls <- page %>% 
    html_nodes("td a") %>% 
    html_attr("href")
  region_codes <- region_urls %>% 
    str_extract("[-A-Z0-9]+\\?") %>% 
    str_replace("\\?", "")
  region_levels <- region_urls %>% 
    str_extract("ebird/[a-z0-9]+") %>% 
    str_replace("ebird/", "")
  mutate(region_stats, 
         region_code = region_codes, 
         region_level = region_levels) %>% 
    select(region_code, region_name, region_level,
           n_species, n_checklists) %>% 
    as_tibble()
}
```

### State data

I use the above function to grab all the state/province data.


```r
state_df <- map_df(ebird_regions$region_code, extract_subregion_data,
                   region_level = "country")
```

### County data

And the same for county data. Note that Mexico has no data below state.


```r
county_df <- map_df(state_df$region_code, extract_subregion_data,
                    region_level = "subnational1")
```

### Combine

Combine the country, state, and county data together into one data frame.


```r
ebird_regions <- bind_rows(ebird_regions, state_df, county_df)
```

## Top hotspots

Next I want to know the top hotspot in each of these regions. Again, I'll write a function to scrape the eBird website.


```r
extract_top_hotspot <- function(region_code, region_level, sleep = 5) {
  Sys.sleep(sleep)
  # species and checklist  counts
  base_url <- "http://ebird.org/ebird/%s/%s/hotspots?yr=all"
  th <- sprintf(base_url, region_level, region_code) %>% 
    read_html() %>% 
    html_node("td a") %>% 
    html_attr("href") %>% 
    str_extract("L[0-9]+")
  if (!is.character(th) || length(th) != 1) {
    return(NA_character_)
  } else {
    return(th)
  }
}
```

I call this function for every country, state, and county. This takes awhile...


```r
ebird_regions <- ebird_regions %>% 
  mutate(top_hotspot = map2_chr(region_code, region_level, extract_top_hotspot))
```

### Hotspot details

Each hotspot has its own page from which I extract the number of species seen, number of checklists, and location. The following function scrapes a hotspot page. Note here (and earlier in the post), that once you've selected a specific component of the page using `read_nodes()`, you extract data using either `html_text()` or `html_attr()`. `html_text()` grabs the text between a set of tags, while `html_attr()` grabs the value for a particular attribute (e.g. the `href` value of an `<a>` tag).


```r
extract_hotspot_data <- function(hotspot_id, sleep = 5) {
  Sys.sleep(sleep)
  if (is.na(hotspot_id)) {
    return(data_frame())
  }
  # species and checklist counts
  base_url <- "http://ebird.org/ebird/hotspot/%s"
  page <- sprintf(base_url, hotspot_id) %>% 
    read_html()
  # hotspot name
  hotspot_name <- page %>% 
    html_node(".hotspot--name") %>% 
    html_text() %>% 
    trimws()
  # number of species and checklists
  hotspot_stats <- page %>% 
    html_nodes("span.hs-section-count") %>% 
    html_text() %>% 
    parse_number()
  # location
  hotspot_loc <- page %>% 
    html_nodes("div.sub-nat a") %>% 
    html_attr("href") %>% 
    str_subset("google") %>% 
    str_match("ll=([-.0-9]+),([-.0-9]+)$") %>% 
    `[`(1, 2:3) %>% 
    parse_number()
  data_frame(hotspot_id = hotspot_id,
             hotspot_name = hotspot_name,
             lat = hotspot_loc[1],
             lng = hotspot_loc[2],
             hotspot_species = hotspot_stats[1],
             hotspot_checklists = hotspot_stats[2])
}
top_hotspots <- ebird_regions$top_hotspot %>% 
  unique() %>% 
  map_df(extract_hotspot_data) %>% 
  mutate(hotspot_name = str_replace(hotspot_name, "\\([ .A-z]+ Co\\.\\)", ""),
         hotspot_name = trimws(hotspot_name))
```

Now I bring the hotspot name into the region data and convert the hotspot data frame into a spatial object using the latitude and longitude.


```r
class(top_hotspots) <- c("tbl_df", "tbl", "data.frame")
# regions
ebird_regions <- top_hotspots %>% 
  select(top_hotspot = hotspot_id, hotspot_name) %>% 
  left_join(ebird_regions, ., by = "top_hotspot")
# hotspots
top_hotspots <- top_hotspots %>% 
  st_as_sf(coords = c("lng", "lat")) %>% 
  st_set_crs(4326)
```





## Maps

Now that I have the data from eBird, I'll make a few maps. First, I'll produce static `ggplot2` choropleth maps of the number of species seen in each state and county. Then I'll produce an interactive `leaflet` map of the top eBird hotspot in each county.

### ggplot2 state map

Before diving into the county data, I'll map the number of species within each state or province in North America. First, I use `rnaturalearth` to access the geographical boundaries of the states.


```r
# lookup table to fix region codes
mx_codes <- ebird_regions %>% 
  filter(region_level == "subnational1", str_detect(region_code, "^MX")) %>% 
  select(region_name, region_code) %>% 
  deframe()
# state boundaries
na_states <- ne_states(iso_a2 = c("CA", "US", "MX"), returnclass = "sf") %>% 
  mutate(region_code = if_else(iso_a2 == "MX", mx_codes[name], iso_3166_2)) %>% 
  select(region_code, country = iso_a2, state = postal)
# join attribute and spatial data
ebird_states <- ebird_regions %>% 
  filter(region_level == "subnational1", region_name != "Hawaii") %>% 
  select(region_code, region_name, n_species) %>% 
  inner_join(na_states, ., by = "region_code")
```

Now, I produce a choropleth map of species seen in each state using the new `ggplot2::geom_sf()` ggplot2 geom for working with `sf` objects.


```r
# lambers conformal conic projection
proj <- paste0("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 ",
               "+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
ebird_states_aea <- st_transform(ebird_states, crs = proj) %>% 
  st_simplify(preserveTopology = TRUE, dTolerance = 1000)
# internal country boundary lines
bl <- ne_download(scale = 110, type = "admin_0_boundary_lines_land", 
                  returnclass = "sf") %>% 
  st_transform(crs = proj)
#> OGR data source with driver: ESRI Shapefile 
#> Source: "/var/folders/mg/qh40qmqd7376xn8qxd6hm5lwjyy0h2/T//RtmpHIrCxz", layer: "ne_110m_admin_0_boundary_lines_land"
#> with 185 features
#> It has 4 fields
#> Integer64 fields read as strings:  scalerank
bl <- bl[ebird_states_aea, ]
# plot
ggplot(ebird_states_aea) + 
  geom_sf(aes(fill = n_species), color = "white", size = 0.1) +
  geom_sf(data = bl, color = "white", size = 0.5) +
  scale_fill_viridis("# species", 
                     limits = c(0, 750), breaks = seq(0, 750, by = 150)) +
  ggtitle("Number of species reported on eBird by state") +
  labs(caption = "Data source: ebird.org") +
  theme_void() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.2, 0.05), legend.justification = c(0, 0), 
        legend.direction = "vertical",
        legend.key.width = unit(0.02, "npc"),
        legend.key.height = unit(0.09, "npc"))
```

<a href="/figures//ebird-county_choropleth-species-st-1.png"><img src="/figures//ebird-county_choropleth-species-st-1.png" title="plot of chunk choropleth-species-st" alt="plot of chunk choropleth-species-st" style="display: block; margin: auto;" /></a>

In the US, California (683 species) and Texas (660 species) are, unsurprisingly, the big winners here. Mexico's Oaxaca is overall the most diverse state in North America, with an impressive 742 species reported; more than all of Canada. At the other end, Canada's Nunavut has a paltry 230 species. 

### County boundaries

First, I'll need to get the boundaries for the counties. The pain here is that county definitions seem to be fluid and eBird doesn't exactly match with any spatial data I've found online. So, I need to combine various sources of data, each of which needs to be manually modified to match eBird. This next bit it a total mess and is probably best skipped over.

For the US county boundaries, I use the `USAboundaries` package.


```r
us_counties <- us_counties() %>% 
  st_as_sf() %>% 
  filter(!state_name %in% c("Hawaii", "Puerto Rico")) %>% 
  mutate(state = state.abb[match(state_name, state.name)],
         state = if_else(state_name == "District of Columbia", "DC", state),
         region_code = paste("US", state, countyfp, sep = "-"),
         area = st_area(.) %>% as.numeric()) %>% 
  select(region_code, state, area) %>% 
  # adjust counties to match ebird
  mutate(region_code = recode(region_code,
                              "US-MD-510" = "US-MD-005", 
                              "US-SD-102" = "US-SD-113",
                              "US-AK-230" = "US-AK-232", 
                              "US-AK-105" = "US-AK-232",
                              "US-AK-275" = "US-AK-280", 
                              "US-AK-195" = "US-AK-280",
                              "US-AK-198" = "US-AK-201",
                              "US-AK-158" = "US-AK-102")) %>% 
  # combine some counties together
  group_by(region_code, state) %>% 
  # note that summarise calls st_union() on the geometries
  summarise(area = sum(area)) %>% 
  ungroup()
```

For Canadian boundaries, I use the [GADM](http://gadm.org/) level 2 administrative regions, which seems to most closely match eBird. Again, I need to apply some manual modifications, especially for Québec where GADM is outdated.


```r
ca_counties <- here("_source", "data", "ebird-county", "CAN_adm.gpkg") %>% 
  read_sf("CAN_adm2") %>% 
  mutate(region_code = str_replace_all(HASC_2, "\\.", "-"),
         state = str_match(region_code, "-([A-Z]{2})-")[, 2],
         area = st_area(.) %>% as.numeric()) %>% 
  # remove great lakes
  filter(region_code != "CA-ON-WB") %>% 
  # remove counties needing splitting
  filter(!region_code %in% c("CA-QC-SR", "CA-QC-FS", "CA-QC-MB", 
                             "CA-QC-NQ", "CA-QC-FR", "CA-ON-HN")) %>% 
  # fix county names
  mutate(
    region_code = str_replace(region_code, "-NF-", "-NL-"),
    region_code = if_else(NAME_2 == "Division No. 11" & 
                            region_code == "CA-NL-TE",
                          "CA-NL-EL", region_code),
    region_code = recode(region_code, 
                         "CA-QC-CC" = "CA-QC-DE",
                         "CA-QC-AM" = "CA-QC-APP")) %>% 
  select(region_code, state, area) %>% 
  # combine split counties
  group_by(region_code, state) %>% 
  summarize(area = sum(area)) %>% 
  ungroup()
# split counties using more granular data
can_adm3 <- here("_source", "data", "ebird-county", "CAN_adm.gpkg") %>%
  read_sf("CAN_adm3")
sr_split <- can_adm3 %>%
  filter(NAME_2 %in% "Sept-Rivières--Caniapiscau") %>%
  mutate(region_code = if_else(ID_3 %in% c(4436, 4437, 4439, 4440, 4441, 4442,
                                           4445, 4448, 4451),
                               "CA-QC-CAN", "CA-QC-SR"),
         state = "QC",
         area = st_area(.) %>% as.numeric()) %>%
  group_by(region_code, state) %>%
  summarize(area = sum(area)) %>%
  ungroup()
fs_split <- can_adm3 %>%
  filter(NAME_2 %in% "Le Fjord-du-Saguenay") %>%
  mutate(region_code = if_else(ID_3 %in% c(3792, 3794, 3796, 3797, 3801, 3812,
                                           3813),
                               "CA-QC-SAG", "CA-QC-FS"),
         state = "QC",
         area = st_area(.) %>% as.numeric()) %>%
  group_by(region_code, state) %>%
  summarize(area = sum(area)) %>%
  ungroup()
mb_split <- can_adm3 %>%
  filter(NAME_2 %in% "Minganie--Basse-Côte-Nord") %>%
  mutate(region_code = if_else(ID_3 %in% c(4183, 4184, 4185, 4186, 4188,
                                           4194, 4195, 4198),
                               "CA-QC-GSL", "CA-QC-MB"),
         state = "QC",
         area = st_area(.) %>% as.numeric()) %>%
  group_by(region_code, state) %>%
  summarize(area = sum(area)) %>%
  ungroup()
nq_split <- can_adm3 %>%
  filter(NAME_2 %in% "Nord-du-Québec") %>%
  mutate(region_code = if_else(ID_3 %in% c(4249:4255, 4267, 4269:4274,
                                           4284:4291),
                               "CA-QC-JAM", "CA-QC-NQ"),
         state = "QC",
         area = st_area(.) %>% as.numeric()) %>%
  group_by(region_code, state) %>%
  summarize(area = sum(area)) %>%
  ungroup()
fr_split <- can_adm3 %>%
  filter(NAME_2 %in% "Francheville") %>%
  mutate(region_code = if_else(ID_3 %in% c(3446, 3448:3450, 3459, 3460),
                               "CA-QC-FR", "CA-QC-LCH"),
         state = "QC",
         area = st_area(.) %>% as.numeric()) %>%
  group_by(region_code, state) %>%
  summarize(area = sum(area)) %>%
  ungroup()
hn_split <- can_adm3 %>%
  filter(NAME_2 %in% "Haldimand-Norfolk") %>%
  mutate(region_code = if_else(ID_3 == 2495, "CA-ON-NF", "CA-ON-HD"),
         state = "ON",
         area = st_area(.) %>% as.numeric()) %>%
  group_by(region_code, state) %>%
  summarize(area = sum(area)) %>%
  ungroup()
# bring everything together
ca_counties <- rbind(ca_counties, sr_split, fs_split, mb_split, nq_split,
                     fr_split, hn_split) %>%
  st_set_geometry("geometry") %>%
  select(-geom)
```





Mexico is fortunately easy since eBird only goes to state level, which I already have from above.


```r
mx_states <- na_states %>% 
  filter(country == "MX") %>% 
  mutate(area = st_area(.) %>% as.numeric()) %>% 
  select(region_code, state, area)
```

And I join everything together.


```r
na_counties <- rbind(us_counties %>% st_transform(crs = 4326), 
                     ca_counties %>% st_transform(crs = 4326),
                     mx_states %>% st_transform(crs = 4326))
```

Well, that was a total pain! I now have some sort of Frankenstein map of North American counties, and I can bring in the ebird data.


```r
ebird_counties <- ebird_regions %>% 
  filter(
    (region_level == "subnational2" & str_detect(region_code, "^US|CA")) | 
      (region_level == "subnational1" & str_detect(region_code, "^MX"))) %>% 
  inner_join(na_counties, ., by = "region_code") %>% 
  select(region_code, state, region_name, area, n_species, n_checklists,
         top_hotspot, hotspot_name)
```



### ggplot2 county map

Finally, I produce a choropleth map of species seen in each county using `ggplot2`.


```r
ebird_counties_aea <- st_transform(ebird_counties, crs = proj) %>% 
  # simplify to get smaller file size
  st_simplify(preserveTopology = TRUE, dTolerance = 1000)
# plot
ggplot(ebird_counties_aea) + 
  geom_sf(aes(fill = n_species), color = "white", size = 0) +
  geom_sf(data = bl, color = "white", size = 0.5) +
  scale_fill_viridis("# species", 
                     limits = c(0, 750), breaks = seq(0, 750, by = 150)) +
  ggtitle("Number of species reported on eBird",
          "County eBird totals for US/Canada, state totals for Mexico") +
  labs(caption = "Data source: ebird.org") +
  theme_void() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = c(0.2, 0.05), legend.justification = c(0, 0), 
        legend.direction = "vertical",
        legend.key.width = unit(0.02, "npc"),
        legend.key.height = unit(0.09, "npc"))
```

<a href="/figures//ebird-county_choropleth-species-1.png"><img src="/figures//ebird-county_choropleth-species-1.png" title="plot of chunk choropleth-species" alt="plot of chunk choropleth-species" style="display: block; margin: auto;" /></a>

As you'd expect, Mexico is the winner here, both because there is higher bird diversity, but also because the Mexican data are at the state level while the US and Canada data are at the (smaller) county level. In fact, this map is misleading in general because the mapped regions are so variable in size. Within the US, the southern states, and coastal states, have the most species reported.

### Leaflet choropleth map

Now I reproduce the same choropleth map, but using Leaflet to map it interactively. To simplify things, I'll subset the data to just include the contiguous US. Also, I'll throw on a point for the top hotspot in each county.

I want to have popups appear on the map for each county. To do this, I'll create an additional column in the data frame with the HTML for the popup. Here I define a function to generate the popup HTML from data.


```r
# polygon popups
popup_maker <- function(code, name, hs_code, hs_name, n) {
  link_style <- "target='_blank' style='text-decoration: none;'"
  region_url <- paste0("<strong><a href='http://ebird.org/ebird/subnational2/",
                       code, "?yr=all' ", link_style, ">", 
                       name, "</a></strong> ", n, " species")
  hs_url <- paste0("<strong>Top Hotspot: </strong><a href='", 
                   "http://ebird.org/ebird/hotspot/", hs_code, "?yr=all' ",
                   link_style, ">", hs_name, "</a>")
  paste(region_url, hs_url, sep = "<br/>")
}
# just contiguous us
ebird_counties_us <- ebird_counties %>% 
  filter(str_detect(region_code, "^US"), state != "AK") %>% 
  # define popup
  mutate(region_name = paste(region_name, state, sep = ", "),
         popup = popup_maker(region_code, region_name, 
                             top_hotspot, hotspot_name, 
                             n_species))
# point popups
popup_maker_hs <- function(code, name, n) {
  link_style <- "target='_blank' style='text-decoration: none;'"
  hs_url <- paste0("<strong><a href='http://ebird.org/ebird/hotspot/",
                   code, "?yr=all' ", link_style, ">", name, "</a></strong>")
  n_species <- paste(n, " species")
  paste(hs_url, n_species, sep = "<br/>")
}
# just us hotspots
hotspots_us <- top_hotspots %>% 
  filter(hotspot_id %in% ebird_counties_us$top_hotspot) %>% 
  mutate(hotspot_name = str_replace(hotspot_name, "\\([ .A-z]+ Co\\.\\)", ""),
         hotspot_name = trimws(hotspot_name)) %>% 
  mutate(popup = popup_maker_hs(hotspot_id, hotspot_name, hotspot_species))
```

Then I generate the Leaflet map.


```r
# color palette
p <- colorNumeric("viridis", domain = ebird_counties_us$n_species)
leaflet(ebird_counties_us) %>%
  addTiles() %>% 
  addFullscreenControl() %>% 
  addPolygons(color = "#444444", weight = 0.25, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.7,
              fillColor = ~ p(n_species),
              popup = ~ popup) %>% 
  addLegend("bottomright", pal = p, values = ~ n_species,
            title = "# species", opacity = 1) %>% 
  addCircles(data = hotspots_us, color = "#444444",
             radius = 2000, stroke = FALSE, fillOpacity = 1,
             group = "Hotspots", popup = ~ popup) %>% 
  addLayersControl(overlayGroups = "Hotspots",
                   options = layersControlOptions(collapsed = FALSE))
```


<iframe src="/assets/leaflet/ebird-county_counties.html" style="border: none; width: 800px; height: 500px"></iframe>
<a href="/assets/leaflet/ebird-county_counties.html" target="_blank"><strong>Fullscreen</strong></a> | Data source: [eBird](http://ebird.org/)

### Leaflet hotspot map

Finally, I produce a similar map, this time highlighting the hotspots rather than the counties. For this one, I'll zoom in on the NY, otherwise the points overlap and look messy.


```r
# color palette
p <- colorNumeric("inferno", domain = hotspots_us$hotspot_species)
leaflet(hotspots_us) %>%
  addProviderTiles("Stamen.Terrain") %>% 
  # focus on ne us
  fitBounds(-79, 41, -72, 45) %>%
  addFullscreenControl() %>% 
  # county boundaries
  addPolygons(data = ebird_counties_us,
              color = "#444444", weight = 0.5, smoothFactor = 0.5,
              opacity = 1, fillColor = "#AAAAAA", fillOpacity = 0.1) %>%
  # hotspots
  addCircleMarkers(data = hotspots_us, color = ~ p(hotspot_species),
                   radius = 5, stroke = FALSE, fillOpacity = 1,
                   popup = ~ popup) %>% 
  addLegend("bottomright", pal = p, values = ~ hotspot_species,
            title = "# species", opacity = 1)
```



The fullscreen map is viewable <a href="/assets/leaflet/ebird-county_hotspots.html" target="_blank"><strong>here</strong></a>.
