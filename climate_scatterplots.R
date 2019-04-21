

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
d <- readRDS("E:/chilefornia/ns_data.rds")


#### analysis


d <- mutate(d, bio7 = bio6 - bio5) %>% select(-bio12)

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
             "set8" = c("cwd", "ppt", "bio1", "bio7", "pwp"))

for(sn in names(sets)[2]){
      message(sn)
      set <- sets[[sn]]
      
      s <- d %>%
            mutate(ppt = log10(ppt + 1)) %>%
            select(c(set, "region")) %>%
            group_by(region) %>%
            sample_n(10000) %>%
            ungroup() %>%
            sample_n(nrow(.))
      pd <- pairsData(s, set, "region", mirror=T) %>%
            mutate(x_var = factor(x_var, levels=set),
                   y_var = factor(y_var, levels=set))
      
      # LDA
      #formula <- paste("region ~", paste(set, collapse=" + "))
      #fit <- MASS::lda(as.formula(formula), data=d)
      #s$fit <- predict(fit, s)$class
      #message(mean(s$region==s$fit, na.rm=T))
      
      ### 2D overlap ###
      if(F){
            v <- c("bio1", "bio4")
            if(sn=="set2") v <- c("bio5", "bio6")
            
            hd <- s %>%
                  mutate_at(vars(set), signif, digits=5) %>%
                  split(.$region) %>%
                  map(select, v) %>%
                  map(na.omit) %>%
                  map(distinct) %>%
                  map(as.matrix)
            h <- do.call("rbind", hd)
            
            ih <- map(hd, convhulln) %>%
                  map(inhulln, p=h) %>%
                  do.call("cbind", .) %>%
                  as.data.frame()
            
            h <- cbind(h, ih)
            
            
            p2 <- ggplot(h %>% mutate(clr=paste(north, south)), 
                         aes_string(v[1], v[2], color="clr")) +
                  geom_point() +
                  scale_color_manual(values=c("red", "blue", "purple")) +
                  xlim(min(d[,v[1]]), max(d[,v[1]])) +
                  ylim(min(d[,v[2]]), max(d[,v[2]])) +
                  theme_minimal() +
                  theme(legend.position="none") +
                  labs(title=paste(signif(mean(ih$north & ih$south), 3)*100, "% overlap"))
            ggsave(paste0("E:/chilefornia/chilefornia/overlap_figures/pairs_", sn, "_2D.png"), 
                   width=8, height=8, units="in")
      }
      
      ### full-D overlap ###
      
      hd <- s %>%
            mutate_at(vars(set), signif, digits=5) %>%
            split(.$region) %>%
            map(select, set) %>%
            map(na.omit) %>%
            map(distinct) %>%
            map(as.matrix)
      h <- do.call("rbind", hd)
      
      ih <- map(hd, convhulln) %>%
            map(inhulln, p=h) %>%
            do.call("cbind", .) %>%
            as.data.frame()
      
      p6 <- ggplot(pd, aes(x_value, y_value, color=region)) +
            geom_point(size=.2) +
            scale_color_manual(values=c("dodgerblue", "darkred")) +
            facet_grid(y_var ~ x_var, scales="free") +
            theme_minimal() +
            labs(title=paste0(sn, ": ", paste(set, collapse=", "),"\n", 
                             signif(mean(ih$north & ih$south), 3)*100, "% overlap"))
      ggsave(paste0("E:/chilefornia/chilefornia/overlap_figures/pairs_", sn, ".png"), 
             width=8, height=8, units="in")
      
      
      if(sn=="set2"){
            s <- d %>% 
                  mutate(ppt = log10(ppt + 1)) %>%
                  na.omit()
            
            hd <- s %>%
                  mutate_at(vars(set), signif, digits=5) %>%
                  arrange(region) %>%
                  split(.$region) %>%
                  map(select, set) %>%
                  map(as.matrix)
            h <- do.call("rbind", hd)
            
            ih <- map(hd, convhulln) %>%
                  map(inhulln, p=h) %>%
                  do.call("cbind", .) %>%
                  as.data.frame() %>%
                  bind_cols(s)
            
            p <- ggplot(ih %>% mutate(overlap = north & south) %>%
                              group_by(region) %>%
                              mutate(po = signif(mean(overlap)*100, 3),
                                     label = paste0(region, ": ", po, "% overlap")), 
                   aes(x, y, fill=overlap)) +
                  geom_raster() +
                  facet_wrap(~label, scales="free") +
                  theme_void() +
                  labs(title=paste0(sn, ": ", paste(set, collapse=", "), "\n"))
            ggsave(paste0("E:/chilefornia/chilefornia/overlap_figures/map_", sn, ".png"), 
                   width=6, height=8, units="in")
      }
      
      
      next()
      ### ecospat.mess ###
      
      md <- d %>%
            split(.$region) %>%
            map(select, c("x", "y", v)) %>%
            map(na.omit)
      
      mn <- ecospat.mess(md[[1]], md[[2]]) %>% 
            as.data.frame() %>% 
            cbind(select(md[[1]], v)) 
      ms <- ecospat.mess(md[[2]], md[[1]]) %>% 
            as.data.frame() %>% 
            cbind(select(md[[2]], v))
      
      pn <- ggplot() +
            geom_point(data=ms %>% sample_n(5000), 
                       aes_string(v[1], v[2]), color="gray80") +
            geom_point(data=mn %>% sample_n(5000), 
                       aes_string(v[1], v[2], color="MESSw")) +
            scale_color_gradientn(colours=c("darkred", "orange",  "dodgerblue", "darkblue"),
                                  values=c(0, .49, .51, 1),
                                  limits=max(abs(range(mn[,"MESSw"])))*c(-1, 1)) +
            xlim(min(d[,v[1]]), max(d[,v[1]])) +
            ylim(min(d[,v[2]]), max(d[,v[2]])) +
            theme_minimal() +
            theme() +
            labs(title="MESS north")
      ps <- ggplot() +
            geom_point(data=mn %>% sample_n(5000), 
                       aes_string(v[1], v[2]), color="gray80") +
            geom_point(data=ms %>% sample_n(5000), 
                       aes_string(v[1], v[2], color="MESSw")) +
            scale_color_gradientn(colours=c("darkred", "orange", "dodgerblue", "darkblue"),
                                  values=c(0, .49, .51, 1),
                                  limits=max(abs(range(ms[,"MESSw"])))*c(-1, 1)) +
            xlim(min(d[,v[1]]), max(d[,v[1]])) +
            ylim(min(d[,v[2]]), max(d[,v[2]])) +
            theme_minimal() +
            theme() +
            labs(title="MESS south")
      #ggsave(paste0("E:/chilefornia/chilefornia/overlap_figures/pairs_", sn, "_2D_mess.png"), 
      #       width=8, height=8, units="in")
      
      p <- arrangeGrob(p2, pn, ps, nrow=1)
      png(paste0("E:/chilefornia/chilefornia/overlap_figures/", sn, ".png"), 
          width=1200, height=400)
      grid.draw(p)
      dev.off()
      
}

