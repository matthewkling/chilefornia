

#### scrape jepson eflora ####

# get characteristics for all species in target genera


library(rvest)
library(stringr)
library(dplyr)
library(tidyr)


genus_urls <- function(genera){
      #genera <- c("Pinus", "Prunus")
      
      html <- read_html("http://ucjeps.berkeley.edu/eflora/toc.html")
      
      gen <- html %>%
            html_nodes(".bodyText a") %>%
            html_text()
      
      tids <- html %>%
            html_nodes(".bodyText a") %>%
            html_attrs() %>%
            unlist() %>%
            str_split("tid=") %>%
            sapply(function(x) x[2])
      
      base_url <- "http://ucjeps.berkeley.edu/eflora/eflora_display.php?tid="
      
      gen <- data.frame(genus=gen, 
                           tid=tids, 
                           url=paste0(base_url, tids),
                           stringsAsFactors=F) %>%
            filter(genus %in% genera)
      
      return(gen$url)
      
}


species_urls <- function(genus_url){
      #genus_url <- "http://ucjeps.berkeley.edu/eflora/eflora_display.php?tid=10029"
      
      html <- read_html(genus_url)
      
      tids <- html %>%
            html_node("select") %>%
            html_children() %>%
            html_attrs() %>%
            unlist() %>%
            str_split("tid=") %>%
            sapply(function(x) x[2]) %>%
            na.omit()
      
      base_url <- "http://ucjeps.berkeley.edu/eflora/eflora_display.php?tid="
      
      urls <- paste0(base_url, tids)
      
      return(urls)
      
}


species_data <- function(species_url){
      #url <- "http://ucjeps.berkeley.edu/eflora/eflora_display.php?tid=11520"
      
      html <- read_html(species_url)
      
      vars <- html %>%
            html_nodes("table+ table div") %>%
            html_children()
      
      taxon <- vars[1] %>%
            html_children() %>%
            html_text() %>%
            paste(collapse=" var. ")
      
      vars <- vars %>% html_text()
      vars <- vars[nchar(vars)>0 & 
                         substr(vars, nchar(vars), nchar(vars))==":"]
      
      text <- html %>%
            html_nodes("table+ table div") %>%
            html_text(trim=T)
      text <- text[[2]]
      
      for(var in vars){
            text <- text %>%
                  str_split(var) %>%
                  unlist()
      }
      
      d <- data.frame(attribute = "Taxon:",
                      value = taxon,
                      stringsAsFactors = F)
      
      for(var in c("Status:", vars)){
            d <- rbind(d, c(var, sub(": ", "", text[1])))
            text <- text[2:length(text)]
      }
      
      d %>%
            mutate(attribute=sub(":", "", attribute)) %>%
            column_to_rownames("attribute") %>%
            t() %>% 
            as.data.frame() %>%
            mutate(url=species_url) %>%
            return()
}



#d <- c("Pinus", "Quercus") %>%
#      genus_urls() %>%
#      lapply(species_urls) %>%
#      unlist() %>%
#      lapply(species_data) %>%
#      Reduce("full_join", .)
