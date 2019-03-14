
library(tidyverse)

setwd("~/desktop/chilefornia/GCM_selection")

files <- list.files(full.names=T, recursive=T, pattern="DifPre_")
files <- files[grepl("bio1_bio14|bio11_bio17", files)]
      
d <- files %>%
      lapply(read_csv)
for(i in 1:length(d)) d[[i]]$scenario <- files[i] %>% substr(nchar(.)-19, nchar(.)-15)

fx <- function(x) na.omit(x)[1]
d <- Reduce("full_join", d) %>%
      group_by(GCM, scenario) %>%
      summarize_all(fx) %>%
      filter(! GCM %in% c("ENSEMBLE", "BASELINE"))
d <- d[, c(1, 2, 5, 6, 9, 10)]
names(d) <- c("GCM", "scenario", "bio1", "bio14", "bio11", "bio17")

# remove models that don't have all 3 scenarios
d <- d %>% 
      group_by(GCM) %>% 
      add_tally() %>%
      filter(n==3)

# highlight models that are consistent across scenarios
# in being above- or below-average for all 4 climate variabes
dd <- d %>% 
      group_by(scenario) %>%
      mutate_at(vars(bio1:bio17), scale, scale=F) %>%
      mutate(quadrant = paste(bio1 > 0, bio14 > 0,
                              bio11 > 0, bio17 > 0)) %>%
      group_by(GCM) %>%
      mutate(consistent = length(unique(quadrant))==1) %>%
      arrange(scenario) 

p <- ggplot(dd, aes(bio11, bio17, 
                    color=scenario, group=GCM, alpha=consistent)) +
      geom_vline(xintercept=0)  +
      geom_hline(yintercept=0)  +
      geom_line(color="black") +
      geom_point() +
      geom_text(data= . %>% filter(scenario=="rcp45"),
                aes(label=GCM),
                color="black") +
      labs(title="anomalies relative to ensemble, by scenario")
p <- ggsave("gcm_selection_bio11_bio17.png", p, width=8, height=6, units="in")

p <- ggplot(dd, aes(bio1, bio14, 
                    color=scenario, group=GCM, alpha=consistent)) +
      geom_vline(xintercept=0)  +
      geom_hline(yintercept=0)  +
      geom_line(color="black") +
      geom_point() +
      geom_text(data= . %>% filter(scenario=="rcp45"),
                aes(label=GCM),
                color="black") +
      labs(title="anomalies relative to ensemble, by scenario")
p <- ggsave("gcm_selection_bio1_bio14.png", p, width=8, height=6, units="in")

