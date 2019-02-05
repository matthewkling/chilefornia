
library(tidyverse)
library(raster)
library(rgeos)
library(fastcluster)
library(FNN)
library(colormap)
library(ggplot2)


setwd("e:/chilefornia/chilefornia")

# load climate rasters
r <- list.files("f:/chelsa/bio19", full.names=T) %>%
      stack()
wb <- stack("f:/chelsa/derived/water_balance.tif")
names(wb) <- c("PPT", "PET", "AET", "CWD", "RAR")

# load study area polygons
s <- readRDS("data/study_area_south.rds") %>% spTransform(crs(r)) %>% gBuffer(width=0)
n <- readRDS("data/study_area_north.rds") %>% spTransform(crs(r))
b <- gUnion(n, s)
saveRDS(b, "../chilefornia_shapefiles.rds")
b <- readRDS("../chilefornia_shapefiles.rds")

# climate within study areas
r <- r %>% crop(b) %>% mask(b)
writeRaster(r, "../chilefornia_climate.tif", overwrite=T)

wb <- wb %>% crop(b) %>% mask(b)
writeRaster(wb, "../chilefornia_water_balance.tif", overwrite=T)

stack("../chilefornia_climate.tif") %>%
      stack(wb) %>%
      writeRaster("../chilefornia_climate_all.tif", overwrite=T)


r <- stack("../chilefornia_climate_all.tif")
template <- r[[1]]

# convert raster to matrix
a <- !is.na(r[[1]][])
v <- coordinates(r)[a,]
for(i in 1:nlayers(r)){
      message(paste("layer", i))
      v <- cbind(v, r[[i]][][a])
}
v0 <- v
a0 <- apply(v0, 1, function(x) !is.na(sum(x)))
v <- na.omit(v)

ar <- a
ar[ar] <- a0

vars <- list.files("f:/chelsa/bio19") %>%
      gsub("CHELSA_bio10_|\\.tif", "", .) %>%
      paste0("bio", .) %>%
      c(c("PPT", "PET", "AET", "CWD", "RAR"))
colnames(v) <- c("x", "y", vars)

saveRDS(v, "chilefornia_climate_matrix_allvars.rds")
v <- readRDS("chilefornia_climate_matrix_allvars.rds")
v <- v[,colnames(v) != "PPT"] # eliminate duplicate variable

# log-transform ppt vars; scale all vars
alog <- function(x){
      x[x==0] <- min(x[x!=0])
      log10(x)
}
for(i in which(colnames(v) %in% c(paste0("bio", 12:18)))) v[,i] <- alog(v[,i])
for(i in 3:ncol(v)) v[,i] <- scale(v[,i])


# correlation plot
r <- cor(v[,3:ncol(v)], method="spearman")
row.names(r) <- ecoclim::translate(row.names(r), "words")
row.names(r)[row.names(r)=="NULL"] <- colnames(r)[row.names(r)=="NULL"]
png("climate_figures/climate_r2.png", width=800, height=500)
corrplot::corrplot(r^2, order="hclust", is.corr=F, 
                   col=colorRampPalette(c("white", "white",  "white", "white", "white",
                                          "white", 
                                          "yellow", "orange", "red", "darkmagenta", "black"))(72))
dev.off()





# experiment with algorithmically finding the lowest-correlation variable set
r <- cor(v[,3:ncol(v)], method="spearman")
cmb <- combn(1:ncol(r), 5) %>% t()

cors <- function(x){
      #x <- cmb[1234,]
      rx <- r[x,]
      rxi <- rx[,x] ^ 2
      #rxe <- rx[,-x] ^ 2
      return(c(sum(rxi), max(rxi[rxi!=1])))
}

rd <- t(apply(cmb, 1, cors)) %>%
      as.data.frame() %>%
      cbind(t(apply(cmb, 1, function(x) colnames(r)[x])), .)
colnames(rd) <- c(paste0("v", 1:5), "ri", "mi")

rd <- rd %>% mutate(id=1:nrow(.),
                    ntemp = apply(rd[,1:5], 1, function(x) sum(paste0("bio", 1:11) %in% x)),
                    nprecip = apply(rd[,1:5], 1, function(x) sum(paste0("bio", 12:19) %in% x)),
                    ncwd = apply(rd[,1:5], 1, function(x) sum("CWD" %in% x)))


rd %>%
      filter(ncwd==1, ntemp==2, nprecip==2) %>%
      filter(ri==min(ri))

rd %>%
      filter(ncwd==1, ntemp==2, nprecip==2) %>%
      filter(mi==min(mi)) %>%
      filter(ri==min(ri)) %>%
      gather(v, var, v1:v5) %>%
      mutate(words = ecoclim::translate(var, "words"))

rd %>%
      filter(mi==min(mi)) %>%
      filter(ri==min(ri)) %>%
      gather(v, var, v1:v5) %>%
      mutate(words = ecoclim::translate(var, "words"))




### final selection of 5 variables ###

variables <- c("bio5", "bio6", "bio15", "AET", "CWD")
variables <- c("bio5", "bio6", "AET", "CWD")
variables <- c("bio1", "bio11", "bio17", "CWD")
v <- v[,match(c("x", "y", variables), colnames(v))]





# subsample pixels for speed
px <- sample(nrow(v), 100000) # change this to 100k for production run

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
c3d <- apply(col3d, 2, function(x) ecdf(x)(x))
p <- ggplot(pd, aes(x, y)) + 
      #geom_raster(fill=rgb(col3d, maxColorValue=255)) +
      geom_raster(fill=rgb(1-c3d[,c(3,2,1)], maxColorValue=1)) +
      theme_minimal() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank()) +
      coord_fixed() +
      labs(y="degrees poleward")
png(paste0("climate_figures/climate_continuous.png"), width=8, height=8, units="in", res=1000)
plot(p)
dev.off()
png(paste0("climate_figures/climate_continuous_v2.png"), width=8, height=8, units="in", res=1000)
plot(p %+% mutate(pd, y=-y))
dev.off()


# fit clusters
tree <- hclust.vector(v[px,3:ncol(v)], method="ward")

# cut tree into specified number of clusters
for(k in c(5:10)){
      clust <- cutree(tree, k)
      
      # transfer cluster identities to non-sampled pixels
      cluster <- clust[nn$nn.index]
      
      
      # raster of clusters
      cr <- template
      cr[] <- NA
      cr[ar] <- cluster
      writeRaster(cr, paste0("cluster_rasters/clusters_k", k, "_",
                             paste(variables, collapse="+"), ".tif"),
                  overwrite=T)
      #next()
      
      # visualize w hierarchical colors
      hclrs <- as.data.frame(cbind(cluster, col3d)) %>%
            group_by(cluster) %>%
            mutate_all(funs(mean)) %>%
            mutate(hex=rgb(red, green, blue, maxColorValue=255))
      p <- ggplot(pd, aes(x, y)) + 
            geom_raster(fill=hclrs$hex) +
            theme_minimal() +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank()) +
            coord_fixed() +
            labs(y="degrees poleward")
      png(paste0("climate_figures/climate_clusters_", k, "_hiercol.png"), width=8, height=8, units="in", res=1000)
      plot(p)
      dev.off()
      png(paste0("climate_figures/climate_clusters_", k, "_hiercol_v2.png"), width=8, height=8, units="in", res=1000)
      plot(p %+% mutate(pd, y=-y))
      dev.off()
      
      
      # visualize w distant colors
      clrz <- distant_colors(k)
      clrs <- clrz[cluster]
      p <- ggplot(pd, aes(x, y)) + 
            geom_raster(fill=clrs) +
            theme_minimal() +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank()) +
            coord_fixed() +
            labs(y="degrees poleward")
      png(paste0("climate_figures/climate_clusters_", k, "_distcol.png"), width=8, height=8, units="in", res=1000)
      plot(p)
      dev.off()
      png(paste0("climate_figures/climate_clusters_", k, "_distcol_v2.png"), width=8, height=8, units="in", res=1000)
      plot(p %+% mutate(pd, y=-y))
      dev.off()
      
      x <- pd %>%
            mutate(cluster=cluster) %>%
            gather(var, value, -x, -y, -cluster) %>%
            mutate(var=gsub("CHELSA_bio10_", "", var),
                   vr=var)
      xv <- x$var
      xv[!xv %in% c("AET", "CWD")] <- unlist(ecoclim::translate(xv[!xv %in% c("AET", "CWD")], "words"))
      x$var <- xv
      
      p <- ggplot(x %>% mutate(cluster=factor(cluster)), 
                  aes(value, color=cluster, fill=cluster)) +
            geom_density(alpha=.3) +
            scale_color_manual(values=clrz) +
            scale_fill_manual(values=clrz) +
            facet_wrap(~var, scales="free_y", ncol=1) +
            theme_minimal() +
            theme(legend.position="none",
                  axis.text.y=element_blank())
      png(paste0("climate_figures/climate_clusters_", k, "_histograms.png"), width=8, height=8, units="in", res=1000)
      plot(p)
      dev.off()
}


############################################

# assess which 5 variables best predict known vegetation formations
# using random forests + recursive feature elimination

library(rgdal)
library(randomForest)

#### data prep: climate and veg values for random points in chile

r <- stack("../chilefornia_climate_all.tif")

vars <- list.files("f:/chelsa/bio19") %>%
      gsub("CHELSA_bio10_|\\.tif", "", .) %>%
      paste0("bio", .) %>%
      c(c("PPT", "PET", "AET", "CWD", "RAR"))
names(r) <- vars



#vegN <- readOGR("e:/chilefornia/veg_NA", "Ecoregions_NA")
#veg <- readOGR("e:/chilefornia/Veg_Form_Chile", "Veg_Form_Chile")
#veg <- spTransform(veg, crs(r))


vegS <- raster("e:/chilefornia/vegformations/chile_form.tif")
vegN <- raster("e:/chilefornia/vegformations/na_form.tif")

imp <- data.frame()
for(veg in list(vegS, vegN)){
      
      pts <- r[[1]] %>% 
            #crop(extent(-137.5668, -66.41681, -55.98347, 0)) %>%
            crop(veg) %>%
            trim() %>%
            rasterToPoints() %>%
            as.data.frame() %>%
            sample_n(100000)
      coordinates(pts) <- c("x", "y")
      crs(pts) <- crs(r)
      
      veg_pts <- extract(veg, pts)
      clim_pts <- extract(r, pts)
      
      d <- cbind(veg_pts, clim_pts) %>% na.omit()
      veg_pts <- d[,1]
      clim_pts <- d[,2:ncol(d)]
      
      vars <- colnames(clim_pts)
      vars <- vars[!vars %in% c("PPT", "bio3", "bio7", "RAR")] # duplicated & unwanted variables
      
      nreps <- 1000
      rfe <- expand.grid(rep=1:nreps,
                         var=vars,
                         elim=NA,
                         imp=NA,
                         stringsAsFactors=F)
      
      
      library(doParallel)
      cl <- makeCluster(detectCores()-2)
      registerDoParallel(cl)
      rfe <- foreach(rp = 1:nreps, .combine="rbind") %dopar% {
            
            require(randomForest)
            require(tidyverse)
            
            rfd <- filter(rfe, rep==rp)
            ss <- sample(nrow(clim_pts), 1000)
            
            for(i in 1:(length(vars)-1)){
                  
                  # fit a random forest model
                  vrs <- colnames(clim_pts) %in% setdiff(vars, rfd$var[!is.na(rfd$elim)])
                  fit <- randomForest(x=clim_pts[ss,vrs], 
                                      y=factor(veg_pts[ss]), 
                                      importance=TRUE)
                  
                  # identify the least important variable 
                  # and add it to list of eliminated variables
                  imp <- fit$importance
                  imp <- imp[,ncol(imp)-1]
                  v <- imp[which.min(imp)]
                  rfd$elim[rfd$var==names(v)] <- i
                  rfd$imp[rfd$var==names(v)] <- v
                  #message(paste("rep", rep, "-- eliminating", names(v)))
            }
            
            # document the last variable that didn't get eliminated
            imp <- fit$importance
            imp <- imp[,ncol(imp)-1]
            v <- imp[which.max(imp)]
            rfd$elim[rfd$var==names(v)] <- i+1
            rfd$imp[rfd$var==names(v)] <- v
            
            return(rfd)
      }
      stopCluster(cl)
      
      rfe$variable <- NA
      
      rfe$variable <- ifelse(rfe$var %in% c("PPT", "RAR", "CWD", "AET", "PET"),
                             as.character(rfe$var),
                             ecoclim::translate(as.character(rfe$var)))
      rfe$variable <- unlist(rfe$variable)
      
      rfe$region <- names(veg)[1]
      
      imp <- rbind(imp, rfe)
      
}





labs <- imp %>%
      group_by(variable) %>%
      summarize(elim=mean(elim), imp=mean(imp))

rfed <- imp %>%
      group_by(variable) %>%
      mutate(rank=mean(elim)) %>%
      arrange(rank) %>%
      ungroup() %>%
      mutate(variable=factor(.$variable, levels=unique(variable))) %>%
      group_by(variable, elim) %>%
      mutate(rank=mean(rank),
             n=length(rank))

p <- ggplot(rfed, aes(elim, variable, size=n)) +
      geom_point() +
      theme_minimal() +
      theme_minimal() +
      theme(plot.title=element_text(hjust=.5)) +
      labs(x="elimination order (high values are most important variables)",
           title="Climate variable importance in predicting Chilean vegetation formations\n(estimated using random forest models and recursive feature elimination)")
ggsave("climate_figures/variable_importance.png", p, width=8, height=6, units="in")





labs <- imp %>%
      group_by(variable, region) %>%
      summarize(elim=mean(elim), imp=mean(imp))

rfed <- imp %>%
      group_by(variable) %>%
      mutate(rank=mean(elim)) %>%
      arrange(rank) %>%
      ungroup() %>%
      mutate(variable=factor(.$variable, levels=unique(variable))) %>%
      group_by(variable, elim, region) %>%
      mutate(rank=mean(rank),
             n=length(rank))

p <- ggplot(rfed, aes(elim, variable, size=n)) +
      geom_point() +
      facet_wrap(~region) +
      theme_minimal() +
      theme(plot.title=element_text(hjust=.5)) +
      labs(x="elimination order (high values are most important variables)",
           title="Climate variable importance in predicting vegetation formations\n(estimated using random forest models and recursive feature elimination)")
ggsave("climate_figures/variable_importance.png", p, width=14, height=6, units="in")
p <- ggplot(rfed %>% arrange(desc(n)), 
            aes(elim, variable, size=n, color=region)) +
      geom_point() +
      theme_minimal() +
      theme(plot.title=element_text(hjust=.5)) +
      labs(x="elimination order (high values are most important variables)",
           title="Climate variable importance in predicting Chilean vegetation formations\n(estimated using random forest models and recursive feature elimination)")
#ggsave("climate_figures/variable_importance.png", p, width=8, height=6, units="in")


rd <- imp %>%
      group_by(variable, region) %>%
      summarize(rank=mean(elim),
                min=quantile(elim, .05),
                max=quantile(elim, .95)) %>%
      arrange(rank) %>%
      ungroup() %>%
      mutate(variable=factor(.$variable, levels=unique(variable))) %>%
      gather(stat, value, rank, min, max) %>%
      unite(region_stat, region, stat) %>%
      spread(region_stat, value)

p <- ggplot(rd, aes(chile_form_rank, na_form_rank)) +
      geom_abline(slope=1, intercept=0, linetype=2, color="gray80") +
      geom_segment(aes(x=chile_form_min, xend=chile_form_max,
                       y=na_form_rank, yend=na_form_rank), color="gray") +
      geom_segment(aes(x=chile_form_rank, xend=chile_form_rank,
                       y=na_form_min, yend=na_form_max), color="gray") +
      geom_point() +
      geom_text(aes(label=variable, x=chile_form_rank+.15,
                    y=na_form_rank-.2), hjust=0) +
      theme_minimal() +
      theme(panel.grid=element_blank()) +
      labs(x="Chile variable importance (high values are more important)",
           y="North America variable importance")
ggsave("climate_figures/variable_importance_scatter.png", p, width=8, height=8, units="in")


############################################



#dtc <- list.files("f:/chelsa/bio19", full.names=T)[1] %>% raster() %>%
#      reclassify(c(NA,NA,1,-Inf,Inf,NA)) %>%
#      distance()
#saveRDS(dtc, "../dtc_chelsa.rds")

elev <- raster("f:/chelsa/elevation/mn30_grd")



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








