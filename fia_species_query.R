

library(data.table)
library(dplyr)
library(tidyr)
library(raster)
library(ggplot2)
select <- dplyr::select


#### load, subset, and join FIA data

# this is the "ENTIRE" Phase 2 dataset, downloaded August 2017

tree <- fread("E:/fia/raw_data/ENTIRE/TREE.csv", stringsAsFactors=F)
tree <- select(tree, PLT_CN, PLOT, SUBP, TREE, SPCD, DIA, HT, WDLDSTEM)
gc()

plot <- fread("E:/fia/raw_data/ENTIRE/PLOT.csv", stringsAsFactors=F)
subplot <- fread("E:/fia/raw_data/ENTIRE/SUBPLOT.csv", stringsAsFactors=F)
species <- fread("E:/fia/raw_data/REF_SPECIES.csv", stringsAsFactors=F)

plot <- select(plot, CN, PLOT, LAT, LON, ELEV, STATECD) %>%
      rename(PLT_CN=CN)
subplot <- select(subplot, PLT_CN, SUBP, SLOPE, ASPECT)
species <- select(species, SPCD, COMMON_NAME, GENUS, SPECIES)
gc()

d <- left_join(tree, subplot) %>%
      left_join(plot) %>%
      left_join(species)
colnames(d) <- tolower(colnames(d))

rm(tree)
rm(plot)
rm(subplot)
rm(species)
gc()



ca <- d %>%
      filter(statecd %in% c(6, 53, 41)) %>%
      mutate(woodland = !is.na(wdldstem)) %>%
      group_by(genus, species, common_name) %>%
      summarize(dia_max_cm = max(dia*2.54, na.rm=T),
                ht_max_m = max(ht/3.28, na.rm=T),
                woodland = any(woodland),
                n_records = n()) %>%
      mutate(dia_10cm = dia_max_cm > 10,
             dia_20cm = dia_max_cm > 20,
             ht_3m = ht_max_m > 3,
             dia_20cm_ht_3m = dia_max_cm > 20 & ht_max_m > 3)

write.csv(ca, "e:/chilefornia/california_fia_species.csv", row.names=F)



source("e:/chilefornia/jepson_parse.r")

jepson <- filter(ca, species=="spp.")$genus %>%
      genus_urls() %>%
      lapply(species_urls) %>%
      unlist() %>%
      lapply(species_data) %>%
      Reduce("full_join", .)
      
write.csv(jepson, "e:/chilefornia/jepson_traits.csv", row.names=F)
