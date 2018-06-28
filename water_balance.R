
# Extraterrestrial solar radiation, per Shuttleworth 1993
ETSR <- function(psi, # latitude, in degrees
                 J){ # julian day
      
      # Shuttleworth's test: 
      # the following should output c(15.0, 15.1, 11.2)
      # round(ETSR(psi=c(30,0,-30), J=105), 1)
      
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
           
      # julian days for target month
      jdays <- which(jmonths==month)
      
      # average ETSR across all days of the month
      mean(ETSR(latitude, jdays))
}

# Hargreaves equation for evapotranspiration
# per Hargreaves & Samani 1985
# (Shuttleworth 1993 seems to incorrectly omit the exponent)
hargreaves <- function(S0, tmean, tmin, tmax){
      
      #ETP, in mm/day
      0.0023 * S0 * (tmax - tmin)^.5 * (tmean + 17.8)
}

# calculate monthly and annual water balance variables for one site
hydro <- function(latitude, # integer
                  ppt, tmean, tmax, tmin){ # vectors of length 12
      #browser()
      if(is.na(tmean[1])) return(rep(NA, 5))
      
      # get S0 values for all 12 months
      S0 <- sapply(1:12, monthly_S0, latitude=latitude)
      
      # monthly water balance variables
      pet <- hargreaves(S0, tmean, tmin, tmax) * dayspermonth
      pet[tmean<0] <- 0 # per Wang et al 2012
      aet <- pmin(pet, ppt)
      cmd <- pmax(0, pet - ppt)
      rr <- pmax(0, ppt - pet)
      
      # annual sums
      return(c(PPT=sum(ppt), PET=sum(pet), AET=sum(aet), CWD=sum(cmd), RAR=sum(rr)))
}

# calculate annual water balance variables for a raster stack
water_balance <- function(rasters, #stack of 48 rasters: ppt1-12, tmean1-12, tmax1-12, tmin1-12
                          enforce_latlong=TRUE){ # leave as T unless rasters are already in lat-long projection
      
      # create latitude raster
      coords <- as.data.frame(coordinates(rasters))
      if(enforce_latlong){
            coordinates(coords) <- c("x", "y")
            crs(coords) <- crs(rasters)
            coords <- spTransform(coords, crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
            coords <- coordinates(coords)
      }
      lat <- rasters[[1]]
      lat[] <- coords[,2] # note: Wang et al 2012 page 21 use a latitude correction that we could decide to implement here
      rasters <- stack(rasters, lat)
      
      # sequence of months corresponding to julian dates for a non-leap year
      jmonths <- as.Date(0:364, format="%j", origin=as.Date("2018-01-01"))
      jmonths <- as.integer(substr(as.character(jmonths), 6, 7))
      dayspermonth <- as.vector(table(jmonths))
      
      # compute annual water balance variables
      w <- function(x, ...) hydro(latitude=x[49], 
                                  ppt=x[1:12], 
                                  tmean=x[13:24],
                                  tmax=x[25:36],
                                  tmin=x[37:48])
      wb <- calc(rasters, w, progress='text')
      names(wb) <- c("PPT", "PET", "AET", "CWD", "RAR")
      return(wb)
}

