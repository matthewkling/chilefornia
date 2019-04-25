

library(raster)
library(tidyverse)
library(ecoclim)
library(geometry)
library(ecospat)
library(grid)
library(gridExtra)
library(corrplot)
extract <- raster::extract
select <- dplyr::select


#### data setup

hydro <- stack("F:/chelsa/chilefornia/historic/derived/historic.gri")
names(hydro) <- c("ppt", "pet", "aet", "cwd", "rar")

bio <- stack(list.files("F:/chelsa/bio19", full.names=T)[c(1, 4, 14:16, 19)])
names(bio) <- paste0("bio", c(1, 12, 4, 5, 6, 9))
bio <- crop(bio, hydro)

pwp <- raster("f:/chelsa/chilefornia/historic/derived/historic_pwp.grd") %>%
      crop(bio)
names(pwp) <- "pwp"

climate <- stack(hydro, bio, pwp)
coords <- coordinates(climate)
climate$x <- climate$y <- climate[[1]]
climate$x[] <- coords[,1]
climate$y[] <- coords[,2]

bs <- readRDS("E:/chilefornia/chilefornia/data/SA_new_january19.rds") %>% spTransform(crs(bio))
bn <- readRDS("E:/chilefornia/chilefornia/data/study_area_north.rds") %>% spTransform(crs(bio))

ds <- extract(climate, bs) %>% as.data.frame() %>% na.omit() %>% mutate(region="south")
dn <- extract(climate, bn) %>% as.data.frame() %>% na.omit() %>% mutate(region="north")
d <- rbind(dn, ds)
saveRDS(d, "E:/chilefornia/ns_data.rds")



#### analysis

d <- readRDS("E:/chilefornia/ns_data.rds") %>% 
      mutate(bio7 = bio6 - bio5) %>% 
      select(-bio12)

cm <- cor(select(d, -x, -y, -region) %>% sample_n(10000))
png("E:/chilefornia/chilefornia/overlap_figures/corr.png")
corrplot(cm^2, order="hclust", is.corr=F)
dev.off()


sets <- list("set1" = c("cwd", "aet", "bio5", "bio6"),
             "set2" = c("cwd", "ppt", "bio5", "bio6"),
             "set3" = c("cwd", "aet", "bio1", "bio7"),
             "set4" = c("cwd", "ppt", "bio1", "bio7"),
             "set5" = c("cwd", "aet", "bio5", "bio6", "pwp"),
             "set6" = c("cwd", "ppt", "bio5", "bio6", "pwp"),
             "set7" = c("cwd", "aet", "bio1", "bio7", "pwp"),
             "set8" = c("cwd", "ppt", "bio1", "bio7", "pwp"),
             "set9" = c("cwd", "bio1", "bio7", "pwp"),
             "set10" = c("cwd", "bio5", "bio6", "pwp"))


for(sn in names(sets)[c(1:4, 7:10)]){
      message(sn)
      set <- sets[[sn]]
      
      s <- d %>% 
            #mutate(ppt = log10(ppt + 1)) %>%
            na.omit() %>%
            distinct() %>%
            arrange(region) %>%
            split(.$region)
      
      ## convex hull method ##
      
      hd <- s %>%
            map(select, set) %>%
            map(as.matrix)
      h <- do.call("rbind", hd)
      
      hull <- map(hd, convhulln) %>%
            map(inhulln, p=h) %>%
            do.call("cbind", .) %>%
            as.data.frame()
      names(hull) <- paste0("hull_", names(hull))
      
      
      ## MESS method ##
      
      hd <- s %>%
            map(select, c("x", "y", set))
      h <- do.call("rbind", hd)
      
      mess <- map(hd, function(x) ecospat.mess(proj=h, cal=x)) %>%
            map(as.data.frame) %>%
            map(select, MESS) %>%
            do.call("cbind", .)
      names(mess) <- c("mess_north", "mess_south")
      
      
      ## ExDet method ##
      
      hd <- s %>%
            map(select, set)
      h <- do.call("rbind", hd)
      
      exdet <- map(hd, ecospat.climan, p=h) %>%
            do.call("cbind", .) %>%
            as.data.frame()
      names(exdet) <- c("exdet_north", "exdet_south")
      
      
      
      
      ## plot ##
      
      p <- s %>%
            bind_rows() %>%
            bind_cols(hull, mess, exdet) %>%
            mutate(hull = hull_north & hull_south,
                   mess = exdet_north > 0 & exdet_south > 0,
                   exdet = exdet_north > 0 & exdet_north < 1 &
                         exdet_south > 0 & exdet_south < 1) %>%
            select(x, y, region, hull:exdet) %>%
            gather(metric, overlap, hull:exdet) %>%
            mutate(metric = factor(metric, levels=c("mess", "exdet", "hull"))) %>%
            
            ggplot(aes(x, y, fill=overlap)) +
            geom_raster() +
            scale_fill_manual(values=c("darkred", "gray75")) +
            facet_wrap(~region+metric, scales="free", nrow=1) +
            theme_void() +
            labs(title=paste0(sn, ": ", paste(set, collapse=", "), "\n")) +
            theme(plot.title=element_text(hjust=.5),
                  legend.position="top")
      ggsave(paste0("E:/chilefornia/chilefornia/overlap_figures/map_", sn, ".png"), 
             width=12, height=6, units="in")
      
}



d %>%
      sample_n(2000) %>%
      ggplot(aes(aet, bio6, color=pwp)) +
      geom_point() +
      scale_color_viridis_c() +
      facet_wrap(~region)

