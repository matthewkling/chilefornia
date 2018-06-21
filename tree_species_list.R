

library(data.table)
library(dplyr)
library(tidyr)
library(raster)
library(ggplot2)
select <- dplyr::select


#### FIA tree species and traits ####

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


# restrict to CA, OR, WA, AK species
ca <- d %>%
      filter(statecd %in% c(2, 6, 41, 53)) %>%
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
write.csv(ca, "data/fia_species.csv", row.names=F)




#### JEPSON: get trait data for the species in FIA list ####

source("e:/chilefornia/chilefornia/jepson_parse.r")

jepson <- filter(ca, species=="spp.")$genus %>%
      genus_urls() %>%
      lapply(species_urls) %>%
      unlist() %>%
      lapply(species_data) %>%
      Reduce("full_join", .)

write.csv(jepson, "data/jepson_traits.csv", row.names=F)



#### BIEN: get list of species in study area, and get data to classify trees ####

library(BIEN)


## spatial species query

boundary <- readRDS("study_area_north.rds")
boundary <- readOGR("e:/chilefornia/chilefornia_v2.kml")
spp <- BIEN_list_spatialpolygons(boundary)
saveRDS(spp, "study_area_north_species.rds")


## trait query

traits <- BIEN_trait_list()
tr <- BIEN_trait_trait(c("diameter at breast height (1.3 m)", 
                         "maximum whole plant height",
                         "whole plant growth form", 
                         "whole plant height", 
                         "whole plant woodiness"))
saveRDS(tr, "bien_traits_raw.rds")
tr <- readRDS("bien_traits_raw.rds")


trees <- tr %>%
      filter(trait_name == "whole plant growth form",
             trait_value %in% c("tree", "Tree")) %>%
      select(scrubbed_species_binomial) %>%
      distinct()

trait <- tr %>%
      filter(!trait_name %in% c("whole plant growth form", "whole plant woodiness")) %>%
      mutate(trait_value=as.numeric(trait_value)) %>%
      group_by(scrubbed_species_binomial, trait_name) %>%
      summarise(mean=mean(trait_value, na.rm=T),
                max=max(trait_value, na.rm=T)) %>%
      ungroup() 


## intersect study area species and tree species








###########################

# alternative to polygon query: state/province species lists
pol <- BIEN_metadata_list_political_names()
us_spp <- BIEN_list_state(country="United States", 
                          state=c("California", "Oregon", "Washington", "Alaska"))
canada_spp <- BIEN_list_state(country="Canada", 
                          state=c("British Columbia"))
mexico_spp <- BIEN_list_state(country="Mexico", 
                              state=c("Estado de Baja California", "Estado de Baja California Sur"))
