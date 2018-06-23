
# Define a the study area boundary along the North American Pacific rim,
# to roughly mirror the domain of Chile in South America.
# Based on watershed boundaries, following the crest of the cordillera 
# from Cabo through Glacier Bay.

# Matthew Kling, June 2018


library(tidyverse)
library(raster)
library(rgeos)
library(rgdal)

setwd("e:/chilefornia")


####### NORTH AMERICA #######

# watersheds -- data from http://www.cec.org/tools-and-resources/map-files/watersheds
w <- readOGR("E:/chilefornia/Watersheds_Shapefile/NA_Watersheds/data/NA_Watersheds", 
             "watershed_p_v2")

# boundary polygon manually drawn in google earth
b <- readOGR("e:/chilefornia/chilefornia_v2.kml") %>%
      spTransform(crs(w))

# watersheds intersecting boundary
o <- over(w, b)
wb <- w[which(!is.na(o$Name)),]
writeOGR(wb, dsn="watersheds_chilefornia_intersect.kml", layer="watersheds", 
         driver="KML", overwrite=T)

# elevation data
altm <- getData("alt", country="MEX")
altu <- getData("alt", country="USA") %>% Reduce("merge", .)
altc <- getData("alt", country="CAN")
alt <- altu %>% merge(altm) %>% merge(altc)
alt <- crop(alt, spTransform(wb, crs(alt)))

# projections
b <- spTransform(b, crs(alt))
w <- spTransform(w, crs(alt))
wb <- spTransform(wb, crs(alt))

# final study area: watersheds fully contained within boundary
o <- over(wb, as(b, "SpatialLines"))
ww <- wb[which(is.na(o)),]
writeOGR(ww, dsn="watersheds_chilefornia_contained.kml", layer="watersheds", 
         driver="KML", overwrite=T)
sa <- gUnaryUnion(ww) %>%
      SpatialPolygonsDataFrame(data.frame(x=1))
writeOGR(sa, dsn="chilefornia_study_area_north.kml", layer="chilefornia", 
         driver="KML", overwrite=T)
saveRDS(sa, "study_area_north.rds")

# plots
plot(alt, col=colorRampPalette(c("darkgreen", "yellow", "white"))(20))
plot(wb, add=T)
plot(ww, add=T, border="red")
plot(b, add=T, border="blue")

plot(sa)
plot(w, add=T, border="gray80")
plot(sa, add=T, border="darkred")



####### SOUTH AMERICA #######

chile <- getData("GADM", country="CHL", level=0)
#ext <- drawExtent()
ext <- extent(c(-77.87009, -65.24649, -56.5389, -16.84332))
chile <- crop(chile, ext)
saveRDS(chile, "data/study_area_south.rds")





############################

# side-by-side plots

s <- readRDS("data/study_area_south.rds")
n <- readRDS("data/study_area_north.rds") %>%
      spTransform(crs(s))

d <- broom::tidy(n) %>% mutate(region="north") %>%
      rbind(broom::tidy(s) %>% mutate(region="south"))

p <- ggplot(d %>% mutate(lat=abs(lat),
                         long=ifelse(region=="south", long-35, long)), 
            aes(long, abs(lat), group=paste(region, group))) +
      geom_polygon(fill="gray50", color=NA) +
      theme_minimal() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank()) +
      coord_map("stereographic")
ggsave("e:/chilefornia/chilefornia_map2.png", width=10, height=10, units="in")

p <- ggplot(d %>% mutate(lat=abs(lat),
                         long=ifelse(region=="south", long-40, long),
                         lat=ifelse(region=="south", abs(lat)+5, lat)), 
            aes(long, lat, group=paste(region, group))) +
      geom_polygon(fill="gray50", color=NA) +
      theme_minimal() +
      theme(axis.title=element_blank(),
            axis.text=element_blank()) +
      coord_map("stereographic")
ggsave("e:/chilefornia/chilefornia_map3.png", width=10, height=10, units="in")
