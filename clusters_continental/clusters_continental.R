
library(tidyverse)
library(raster)
library(rgeos)
library(fastcluster)
library(FNN)
library(colormap)
library(ggplot2)
library(caret)
library(grid)
library(gridExtra)


# load climate rasters
r <- list.files("f:/chelsa/bio19", full.names=T) %>%
      stack()
wb <- stack("f:/chelsa/derived/water_balance.tif")
names(wb) <- c("PPT", "PET", "AET", "CWD", "RAR")
r <- stack(r[[5:6]], wb[[c(1,4)]])



############# prep continental climate data #######

# crop to continental bounding boxes, convert to equal area
bb_s <- extent(-82.77317, -29.61833, -56, 12.5)
bb_n <- extent(-169.0535, -18.83329, 12.5, 66)
cea <- CRS("+proj=cea +lon_0=0 +lat_ts=37.5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
clim_s <- crop(r, bb_s) %>% projectRaster(crs=cea)
clim_n <- crop(r, bb_n) %>% projectRaster(crs=cea)

# transform variables
mat_s <- coordinates(clim_s) %>% cbind(values(clim_s)) %>% na.omit()
mat_n <- coordinates(clim_n) %>% cbind(values(clim_n)) %>% na.omit()
xy_s <- mat_s[,1:2]
xy_n <- mat_n[,1:2]
mat <- rbind(mat_s[,3:6], mat_n[,3:6])
trans <- preProcess(mat[sample(nrow(mat), 100000),], 
                    c("YeoJohnson", "center", "scale"))
mat_s <- predict(trans, mat_s[,3:6])
mat_n <- predict(trans, mat_n[,3:6])

# subsample pixels for speed -- change n to 100k for production run
# and find sampled pixel most similar to each non-sampled pixel
npx <- 100000
px_s <- sample(nrow(mat_s), npx)
px_n <- sample(nrow(mat_n), npx)
nn_s <- get.knnx(mat_s[px_s,], mat_s, k=1)
nn_n <- get.knnx(mat_n[px_n,], mat_n, k=1)

# fit trees
tree_s <- hclust.vector(mat_s[px_s,], method="ward")
tree_n <- hclust.vector(mat_n[px_n,], method="ward")



########## prep climate within study areas ##########

# load study area polygons
sa_s <- readRDS("data/study_area_south.rds") %>% spTransform(cea) %>% gBuffer(width=0)
sa_n <- readRDS("data/NA_April.rds") %>% spTransform(cea)

# climate within study areas
clim_sa_s <- clim_s %>% crop(sa_s) %>% mask(sa_s)
clim_sa_n <- clim_n %>% crop(sa_n) %>% mask(sa_n)

# transform variables
mat_sa_s <- coordinates(clim_sa_s) %>% cbind(values(clim_sa_s)) %>% na.omit()
mat_sa_n <- coordinates(clim_sa_n) %>% cbind(values(clim_sa_n)) %>% na.omit()
xy_sa_s <- mat_sa_s[,1:2]
xy_sa_n <- mat_sa_n[,1:2]
mat_sa_s <- predict(trans, mat_sa_s[,3:6])
mat_sa_n <- predict(trans, mat_sa_n[,3:6])

# nearest neighbors in the continental training dataset
nn_sa_s <- get.knnx(mat_s[px_s,], mat_sa_s, k=1)
nn_sa_n <- get.knnx(mat_n[px_n,], mat_sa_n, k=1)



######### slice and dice #########

for(k in c(5, 10, 15, 20)){
      
      # cut tree and assign cluster ids to all cells in both datasets
      clust_s <- cutree(tree_s, k)
      clust_n <- cutree(tree_n, k)
      cluster_s <- clust_s[nn_s$nn.index]
      cluster_n <- clust_n[nn_n$nn.index]
      cluster_sa_s <- clust_s[nn_sa_s$nn.index]
      cluster_sa_n <- clust_n[nn_sa_n$nn.index]
      
      # continental extent of each cluster
      freq_s <- data.frame(cluster = cluster_s) %>% count(cluster)
      freq_n <- data.frame(cluster = cluster_n) %>% count(cluster)
      
      d_s <- cbind(xy_s, cluster_s) %>% as.data.frame()
      d_n <- cbind(xy_n, cluster_n) %>% as.data.frame()
      d_sa_s <- cbind(xy_sa_s, cluster_sa_s) %>% as.data.frame()
      d_sa_n <- cbind(xy_sa_n, cluster_sa_n) %>% as.data.frame()
      
      col <- distant_colors(k)
      
      map_s <- ggplot(d_s, aes(x, y, fill=factor(cluster_s, levels=1:k))) +
            geom_raster() +
            scale_fill_manual(values=col) +
            coord_fixed() +
            theme_void() +
            theme(legend.position="none")
      map_sa_s <- ggplot(d_sa_s, aes(x, y, fill=factor(cluster_sa_s, levels=1:k))) +
            geom_raster() +
            scale_fill_manual(values=col) +
            coord_fixed() +
            theme_void() +
            theme(legend.position="none")
      bar_s <- ggplot(freq_s, aes(cluster, n, fill=factor(cluster, levels=1:k))) +
            geom_bar(stat="identity") +
            scale_fill_manual(values=col) +
            ylim(0, max(c(freq_s$n, freq_n$n))) +
            theme_minimal() +
            theme(legend.position="none",
                  axis.text.x=element_blank()) +
            labs(x="cluster", y="continental land area")
      p <- arrangeGrob(bar_s, map_s, ncol=1)
      p <- arrangeGrob(map_sa_s, p, nrow=1, widths=c(1, 2))
      ggsave(paste0("clusters_continental/south_k", k, ".png"), p, 
             width=8, height=8, units="in")
      
      
      map_n <- ggplot(d_n, aes(x, y, fill=factor(cluster_n, levels=1:k))) +
            geom_raster() +
            scale_fill_manual(values=col) +
            coord_fixed() +
            theme_void() +
            theme(legend.position="none")
      map_sa_n <- ggplot(d_sa_n, aes(x, y, fill=factor(cluster_sa_n, levels=1:k))) +
            geom_raster() +
            scale_fill_manual(values=col) +
            coord_fixed() +
            theme_void() +
            theme(legend.position="none")
      bar_n <- ggplot(freq_n, aes(cluster, n, fill=factor(cluster, levels=1:k))) +
            geom_bar(stat="identity") +
            scale_fill_manual(values=col) +
            ylim(0, max(c(freq_s$n, freq_n$n))) +
            theme_minimal() +
            theme(legend.position="none",
                  axis.text.x=element_blank()) +
            labs(x="cluster", y="continental land area")
      p <- arrangeGrob(bar_n, map_n, ncol=1)
      p <- arrangeGrob(map_sa_n, p, nrow=1, widths=c(1, 2))
      ggsave(paste0("clusters_continental/north_k", k, ".png"), p, 
             width=8, height=8, units="in")
      
}