
library(tidyverse)
library(raster)
library(rgeos)
library(fastcluster)
library(FNN)
library(colormap)
library(ggplot2)


setwd("e:/chilefornia/chilefornia")

# load climate rasters
r <- list.files("f:/chelsa/bio19", full.names=T) %>% stack() %>% 
      subset(paste0("CHELSA_bio10_", c(2,7,10,11,13,14,18,19)))

# load study area polygons
s <- readRDS("data/study_area_south.rds") %>% spTransform(crs(r))
n <- readRDS("data/study_area_north.rds") %>% spTransform(crs(r))
b <- gUnion(n, s)
saveRDS(b, "../chilefornia_shapefiles.rds")
b <- readRDS("../chilefornia_shapefiles.rds")

# climate within study areas
r <- r %>% crop(b) %>% mask(b)
saveRDS(r, "../chilefornia_climate.rds")
r <- readRDS("../chilefornia_climate.rds")

# convert raster to matrix
v <- cbind(coordinates(r), values(r)) %>% na.omit()

# log-transform ppt vars; scale all vars
alog <- function(x){
      x[x==0] <- min(x[x!=0])
      log10(x)
}
for(i in 7:10) v[,i] <- alog(v[,i])
for(i in 3:ncol(v)) v[,i] <- scale(v[,i])


# subsample pixels for speed
px <- sample(nrow(v), 20000) # change this to 100k for production run

# find sampled pixel most similar to each non-sampled pixel
nn <- get.knnx(v[px,3:ncol(v)], v[,3:ncol(v)], k=1)

# pca for colorspace
pc <- prcomp(v[,3:ncol(v)])$x[,1:3]
col3d <- colors3d(pc) %>%
      col2rgb() %>%
      t()

# plotting data frame
res <- resolution(v[,2][v[,2]>0])
pd <- as.data.frame(v) %>% 
      mutate(x=plyr::round_any(ifelse(y<0, x-35, x), res, ceiling),
             y=ifelse(y<0, 
                      plyr::round_any(y, res, floor),
                      plyr::round_any(y, res, ceiling)),
             y=abs(y))

# continuous color plot
p <- ggplot(pd, aes(x, y)) + 
      geom_raster(fill=rgb(col3d, maxColorValue=255)) +
      theme_minimal() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank()) +
      coord_fixed() +
      labs(y="degrees poleward")
png(paste0("climate_figures/climate_map_continuous.png"), width=8, height=8, units="in", res=1000)
plot(p)
dev.off()


# fit clusters
tree <- hclust.vector(v[px,3:ncol(v)], method="ward")

# cut tree into specified number of clusters
for(k in c(5, 10, 15)){
      clust <- cutree(tree, k)
      
      # transfer cluster identities to non-sampled pixels
      cluster <- clust[nn$nn.index]
      
      # back-convert to raster format and export
      kr <- r[[1]]
      kr[!is.na(values(kr))] <- cluster
      
      # visualize w distant colors
      clrs <- distant_colors(k)[cluster]
      p <- ggplot(pd, aes(x, y)) + 
            geom_raster(fill=clrs) +
            theme_minimal() +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank()) +
            coord_fixed() +
            labs(y="degrees poleward")
      png(paste0("climate_figures/climate_map_clusters_distcol_", k, ".png"), width=8, height=8, units="in", res=1000)
      plot(p)
      dev.off()
      
      # visualize w hierarchical colors
      hclrs <- as.data.frame(cbind(cluster, col3d)) %>%
            group_by(cluster) %>%
            mutate_each(funs(mean)) %>%
            mutate(hex=rgb(red, green, blue, maxColorValue=255))
      p <- ggplot(pd, aes(x, y)) + 
            geom_raster(fill=hclrs$hex) +
            theme_minimal() +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank()) +
            coord_fixed() +
            labs(y="degrees poleward")
      png(paste0("climate_figures/climate_map_clusters_hiercol_", k, ".png"), width=8, height=8, units="in", res=1000)
      plot(p)
      dev.off()
}



##################################################################

stop("old code, for worldclim scatterplots")


library(raster)
library(dplyr)


setwd("e:/chilefornia")

# get data

tmin <- getData("worldclim", var="tmin", res=5)
tmax <- getData("worldclim", var="tmax", res=5)
ppt <- getData("worldclim", var="prec", res=5)
bio <- getData("worldclim", var="bio", res=5)


cali <- getData("GADM", country="USA", level=1)
cali <- cali[cali$NAME_1 == "California",]
chile <- getData("GADM", country="CHL", level=0)

bio <- bio[[c(10,11,18,19)]]

ca <- bio %>% crop(cali) %>% mask(cali)
ch <- bio %>% crop(chile) %>% mask(chile)

cad <- ca %>% rasterToPoints %>% as.data.frame() %>% mutate(region="california")
chd <- ch %>% rasterToPoints %>% as.data.frame() %>% mutate(region="chile")

d <- rbind(cad, chd) %>%
      mutate_at(vars(bio18, bio19), log10)

library(ecoclim)
pd <- d %>% 
      pairsData(c("bio10", "bio11", "bio18", "bio19"),
                c("x", "y", "region")) %>%
      mutate(x_variable=sub("of", "of\n", 
                            translate(as.character(x_var), "words")),
             y_variable=sub("of", "of\n", 
                            translate(as.character(y_var), "words")))

ord <- sub("of", "of\n",  translate(c("bio10", "bio11", "bio18", "bio19"), "words"))
pd <- mutate(pd,
             x_variable = factor(x_variable, levels=ord),
             y_variable = factor(y_variable, levels=ord))

p <- ggplot(sample_n(pd, nrow(pd)), 
            aes(x_value, y_value, color=region)) +
      geom_point(size=.5) +
      facet_grid(y_variable ~ x_variable, scales="free") +
      theme_minimal() +
      theme(legend.position="top") +
      labs(x="degC*10 or log10mm",
           y="degC*10 or log10mm")
ggsave("climate_scatter.png", p, width=8, height=8, units="in")





library(caret)

fit <- d %>%
      filter(is.finite(bio18+bio19)) %>%
      group_by(region) %>%
      sample_n(5000) %>%
      ungroup() %>%
      dplyr::select(bio10:bio19) %>%
      preProcess(method="pca")

pc <- filter(d, is.finite(bio18+bio19)) %>%
      predict(fit, .)
pc$color <- colormap::colors3d(dplyr::select(pc, PC1:PC3))
pc <- filter(pc, region=="california" | x > -80)

ggplot(pc, aes(x, y)) +
      geom_raster(fill=pc$color) +
      facet_wrap(~region, scales="free")



library(rgl)
plot3d(dplyr::select(pc, PC1:PC3), col=c("black", "red")[factor(pc$region)])





##########################################################################








