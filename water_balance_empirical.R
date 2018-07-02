

library(raster)
library(tidyverse)

source("e:/chilefornia/chilefornia/water_balance.R")

# organize input rasters
f <-  list.files("F:/chelsa/monthly48", full.names=T) %>%
      ecoclim::parseMetadata(is.dir=F, skips="2_land",
                             keys=list(var=c("prec", "temp10", "tmin10", "tmax10"),
                                       mo=1:12)) %>%
      mutate(var=sub("10", "", var)) %>%
      arrange(var, mo)

# generate latitude raster
#lat <- latitude(f$path[1], already_latlong=T)
#writeRaster(lat, "f:/chelsa/derived/latitude.tif")

# compute water balance variables
r <- stack(c(f$path, "f:/chelsa/derived/latitude.tif"))
names(r) <- c(paste0(f$var, f$mo), "latitude")
wb_chelsa <- water_balance(r, temp_scalar=0.1, ncores=7)
writeRaster(wb_chelsa, "f:/chelsa/derived/water_balance.tif")


########## plots #########


# plot input layers to diagnose artic issue
for(i in 1:nlayers(r)){
      png(paste0("f:/chelsa/charts/chelsa", i, ".png"), width=500, height=500)
      plot(r[[i]], main=names(r)[i], colNA="black")
      dev.off()
}

#pt <- drawExtent()
pt <- data.frame(x=-129, y=68)
coordinates(pt) <- c("x", "y")
x <- raster::extract(r, pt)


td <- expand.grid(month=1:12, latitude=-90:90) %>%
      rowwise() %>%
      mutate(S0=monthly_S0(month, latitude))
ggplot(td, aes(latitude, S0, color=factor(month))) +
      geom_line() +
      scale_x_continuous(limits=c(-90,90), breaks=seq(-90,90,90)) +
      theme_minimal()



txt <- td %>%
      filter(latitude %in% seq(-90, 90, 10)) %>%
      group_by(latitude) %>%
      filter(month==max(month[is.finite(S0)], na.rm=T))
p <- ggplot(filter(td, latitude %in% seq(-90, 90, 10)),
       aes(month, S0, color=factor(latitude))) +
      geom_line() + geom_point() +
      geom_text(data=txt, aes(month, S0, label=latitude), hjust=-.25) +
      scale_x_continuous(breaks=1:12) +
      scale_color_manual(values=colorRampPalette(c("magenta", "red", "black", "blue", "cyan"))(length(seq(-90, 90, 10)))) +
      theme_minimal() +
      theme(legend.position="none") +
      labs(y="S0 (extraterrestrial solar radiation)")
ggsave("e:/chilefornia/S0_month_latitude.png", p, width=8, height=6, units="in")

# the method fails above the artic circle (66.3 degrees)
ggplot(filter(td, latitude %in% c(66,67)),
       aes(month, S0, color=factor(latitude))) +
      geom_line() + geom_point() +
      #geom_text(data=txt, aes(month, S0, label=latitude), hjust=-.25) +
      scale_x_continuous(breaks=1:12) +
      theme_minimal() +
      theme(legend.position="none")

png("e:/chilefornia/water_balance_chelsa.png", width=500, height=1500)
plot(wb_chelsa, nc=1)
dev.off()



############### calculate WB for CHELSA #####################

library(raster)
library(tidyverse)

f <-  x %>%
      ecoclim::parseMetadata(is.dir=F, skips="2_land",
                             keys=list(var=c("prec", "temp10", "tmin10", "tmax10"),
                                       mo=1:12)) %>%
      mutate(var=sub("10", "", var)) %>%
      arrange(var, mo)

ext <- extent(-124, -119, 36, 41)
r <- f %>%
      select(path) %>%
      unlist() %>%
      lapply(raster) %>%
      lapply(crop, y=ext) %>%
      stack()
names(r) <- paste0(f$var, f$mo)

# correct units by removing multiplier from temperature
for(i in 1:nlayers(r)) if(!grepl("prec", names(r)[i])) r[[i]] <- r[[i]] / 10


# calculate water balance
start <- Sys.time()
wb_chelsa <- water_balance(r, enforce_latlong=F)
elapsed <- Sys.time() - start

total <- elapsed * ncell(raster(f$path[1])) / ncell(wb[[1]])
6000 / 60 / 24



########### calculate WB for CWNA  ###########

# project extent to CNA CRS
ext_cna <- as(ext, "SpatialPolygons")
crs(ext_cna) <- crs(raster("F:/chelsa/monthly48/CHELSA_prec_1_V1.2_land.tif"))
crs_cna <- crs(raster("F:/ClimateNA/1km_Hamman_raw/27_biovars/ClimateNA_Reference/ClimateNA_DEM.asc"))
ext_cna <- spTransform(ext_cna, crs_cna)

r <- list.files("F:/ClimateNA/NA_NORM_6190_Monthly_ASCII", full.names=T) %>%
      lapply(raster) %>%
      lapply(function(x){crs(x) <- crs_cna; return(x)}) %>%
      lapply(crop, y=ext_cna) %>%
      stack()

# correct units by removing multiplier from temperature
#for(i in 1:nlayers(r)) if(!grepl("prec", names(r)[i])) r[[i]] <- r[[i]] / 10


# calculate water balance
wb_cna <- water_balance(r)

# compare to native CNA CWD variable
cna <- raster("F:/ClimateNA/NA_NORM_6190_Bioclim_ASCII/CMD.asc") %>%
      crop(ext_cna)

wb_cna$native <- cna
#pairs(wb_cna)
plot(wb_cna[[c("native", "CWD")]])

dcna <- as.data.frame(rasterToPoints(stack(wb_cna, lat)))

cor(dcna$native, dcna$CWD)
apply(dcna, 2, range)

dcna %>%
      sample_n(1000) %>%
      ggplot(aes(native, CWD, color=PPT01)) +
      geom_point() +
      viridis::scale_color_viridis()




############ compare to other datasets ############

d <- as.data.frame(rasterToPoints(stack(wb, lat)))

pts <- d
coordinates(pts) <- c("x", "y")
crs(pts) <- crs(wb)

# bcm cwd, for comparison
bcm <- raster("F:/BCM/cwd1981_2010_ave_HST_1486167838/cwd1981_2010_ave_HST_1486167838.tif")
#bcm <- raster("F:/BCM/ppt1981_2010_ave_HST_1486148746/ppt1981_2010_ave_HST_1486148746.tif")
ext_bcm <- as(ext, "SpatialPolygons")
crs(ext_bcm) <- crs(r)
ext_bcm <- spTransform(ext_bcm, crs(bcm))
bcm <- crop(bcm, ext_bcm)
bcm <- aggregate(bcm, 3)
d$bcm <- raster::extract(bcm, spTransform(pts, crs(bcm)))

# climateNA cwd, for comparison
cna <- raster("F:/ClimateNA/NA_NORM_6190_Bioclim_ASCII/CMD.asc")
#cna <- raster("F:/ClimateNA/NA_NORM_6190_Bioclim_ASCII/MAP.asc")
prj <- raster("F:/ClimateNA/1km_Hamman_raw/27_biovars/ClimateNA_Reference/ClimateNA_DEM.asc") %>%
      crs()
crs(cna) <- prj
ext_cna <- spTransform(ext_bcm, crs(cna))
cna <- crop(cna, ext_cna)
d$cna <- raster::extract(cna, spTransform(pts, crs(cna)))



cor(select(d, cna, bcm, CWD), use="pairwise.complete.obs")

ggplot(sample_n(d, 5000), 
       aes(cna, bcm, color=prec1)) + 
      geom_point() +
      geom_abline(slope=1, intercept=0, color="red") +
      viridis::scale_color_viridis()


md <- d %>%
      select(x, y, bcm, cna, CWD) %>%
      gather(dataset, cwd, bcm, cna, CWD)

ggplot(md, aes(x, y, fill=cwd)) +
      geom_raster() +
      viridis::scale_fill_viridis() +
      facet_grid(.~dataset) +
      theme_void()






par(mfrow=c(1,3))
plot(bcm)
plot(cna)
plot(wb$CWD)

library(microbenchmark)
microbenchmark(w(x))
