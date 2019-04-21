
library(dismo)
library(ecospat)
library(tidyverse)
library(grid)
library(gridExtra)

# a simulated niche in climate space
d <- expand.grid(a=1:100, b=1:100) %>%
      mutate(suit = case_when(abs(a - b) < 5 & a > 20 & a < 80 & b > 20 & b < 80 ~ 1,
                              TRUE ~ 0))
p1 <- ggplot(d, aes(a, b, fill=suit)) + geom_raster() + 
      scale_fill_viridis_c() + coord_fixed() + theme_minimal() +
      labs(title="simulated niche (training dataset)") +
      theme(legend.position="top")

# fit maxent model on entire dataset
fit <- maxent(select(d, a, b), d$suit)
d$pred <- predict(fit, d)
p2 <- ggplot(d, aes(a, b, fill=pred)) + geom_raster() + 
      scale_fill_viridis_c() + coord_fixed() + theme_minimal() +
      labs(title="maxent fit based on full training dataset") +
      theme(legend.position="top")

# MESS does not correctly identify extrapolation
mess <- ecospat.mess(d %>% select(a, b) %>% mutate(x=1, y=1) %>% select(x, y, a, b), 
                     d %>% select(a, b) %>% mutate(x=1, y=1) %>% select(x, y, a, b)) %>% 
      as.data.frame() %>%
      cbind(d)
p3 <- ggplot(mess, aes(a, b, fill = MESS>=0)) + geom_raster() + 
      coord_fixed() + theme_minimal() +
      labs(title="MESS says no extrapolation") +
      theme(legend.position="top")

# use a subset of climate space for model fitting
md <- mutate(d, training = a>=b & b > 10)
p4 <- ggplot(md, aes(a, b, fill=training)) + geom_raster() + 
      coord_fixed() + theme_minimal() +
      xlim(0, 100) +
      labs(title="filtered training dataset to test for extrapolation") +
      theme(legend.position="top")

# see how maxent fits and extrapolates
fit <- maxent(md %>% filter(training) %>% select(a, b), 
              md %>% filter(training) %>% pull(suit))
d$pred <- predict(fit, d)
p5 <- ggplot(d, aes(a, b, fill=pred)) + geom_raster() + 
      scale_fill_viridis_c() + coord_fixed() + theme_minimal() +
      labs(title="maxent fit based on filtered training dataset") +
      theme(legend.position="top")

# MESS does not correctly identify extrapolation
mess <- ecospat.mess(d %>% select(a, b) %>% mutate(x=1, y=1) %>% select(x, y, a, b), 
                     md %>% filter(training) %>% select(a, b) %>% mutate(x=1, y=1) %>% select(x, y, a, b)) %>% 
      as.data.frame() %>%
      cbind(d)
p6 <- ggplot(mess, aes(a, b, fill = MESS>=0)) + geom_raster() + 
      coord_fixed() + theme_minimal() +
      labs(title="MESS misses extrapolation") +
      theme(legend.position="top")

# save plots
p <- arrangeGrob(p1, p2, p3, p4, p5, p6, nrow=2)
ggsave("E:/chilefornia/chilefornia/overlap_figures/maxent_mess_simulation.png", 
       p, width=12, height=8, units="in")
