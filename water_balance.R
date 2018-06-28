

# Extraterrestrial solar radiation, per Shuttleworth 1993
ETSR <- function(psi, # latitude
                 J){ # julian day
      
      # Shuttleworth's test: 
      # the following should output c(15.0, 15.1, 11.2)
      # ETSR(psi=c(30,0,-30), J=105)
      
      # convert degrees latitude to radians
      psi <- psi * pi / 180
      
      # solar declination
      delta <- 0.4093 * sin((J*2*pi/365) - 1.405)
      
      # sunset hour angle
      omega <- acos(-tan(psi) * tan(delta))
      
      # relative earth-to-sun distance
      dr <- 1 + 0.033 * cos(J*2*pi/365)
      
      # extraterrestrial solar radiation (mm / day)
      S0 <- 15.392 * dr * (omega * sin(psi) * sin(delta) + cos(psi) * cos(delta) * sin(omega))
      
      return(S0)
}

# mean S0 for a month of the year
monthly_S0 <- function(month, latitude){
      
      # sequence of months corresponding to julian dates for a non-leap year
      jmonths <- as.Date(0:364, format="%j", origin=as.Date("2018-01-01"))
      jmonths <- as.integer(substr(as.character(jmonths), 6, 7))
     
      # julian days for target month
      jdays <- which(jmonths==month)
      
      # average ETSR across all days of the month
      mean(ETSR(latitude, jdays))
}

# Hargreaves equation for evapotranspiration, per Shuttleworth 1993
hargreaves <- function(S0, tmean, tmin, tmax){
      #Erc, in mm/day
      0.0023 * S0 * (tmax - tmin) * (tmean - 17.8)
}


hydro <- function(latitude, # integer
                  ppt, tmean, tmax, tmin){ # vectors of length 12
      #browser()
      if(is.na(tmean[1])) return(rep(NA, 3))
      
      # get S0 values for all 12 months
      S0 <- sapply(1:12, monthly_S0, latitude=latitude)
      
      pet_mo <- hargreaves(S0, tmean, tmin, tmax)
      aet_mo <- pmin(pet_mo, ppt)
      cmd_mo <- pmax(0, pet_mo - ppt)
      rr_mo <- pmax(0, ppt - pet_mo)
      
      # annual sums
      return(c(CMD=sum(cmd_mo), AET=sum(aet_mo), RR=sum(rr_mo)))
}


function(rasters){ #stack of 48 rasters: ppt1-12, tmean1-12, tmax1-12, tmin1-12
      
      # genrate a latitude raster
      lat <- rasters[[1]]
      lat[] <- coordinates(lat)[,2]
      rasters <- stack(rasters, lat)
      
      # compute the 3 annual water balance variables
      w <- function(x, ...) hydro(latitude=x[49], 
                                  ppt=x[1:12], 
                                  tmean=x[13:24],
                                  tmax=x[25:36],
                                  tmin=x[37:48])
      water <- calc(rasters, w)#, forceapply=T)
      
}

####################################

library(raster)
library(tidyverse)

ext <- extent(-123.6866, -120.7557, 36.95263, 39.87957)
f <- list.files("f:/chelsa/monthly48", full.names=T) %>%
      ecoclim::parseMetadata(is.dir=F, skips="2_land",
                             keys=list(var=c("prec", "temp10", "tmin10", "tmax10"),
                                       mo=1:12)) %>%
      mutate(var=sub("10", "", var)) %>%
      arrange(var, mo)

r <- f %>%
      select(path) %>%
      unlist() %>%
      lapply(raster) %>%
      lapply(crop, y=ext) %>%
      stack()

names(r) <- paste0(f$var, f$mo)

# correct units by removing multiplier from temperature
for(i in 1:nlayers(r)) if(!grepl("prec", names(r)[i])) r[[i]] <- r[[i]] / 10

rasters <- r

x <- rasters[100,100]
