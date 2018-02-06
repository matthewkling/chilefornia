


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
