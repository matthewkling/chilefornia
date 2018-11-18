
# This script includes a set of functions used to calculate annual
# hydrologic variables from monthly temperature and precipitation data.

# References:
#     Wang et al 2012 -- ClimateWNA -- Journal of Applied Meteorology and Climatology
#     Shuttleworth 1993 -- Evaporation (chapter 4 in Handbook of Hydrology)
#     Hargreaves and Samani 1985 -- Reference Crop Evaporation from Ambient Air Temperature

# Matthew Kling, June 2018



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
      #j <- as.Date(0:364, format="%j", origin=as.Date("2018-01-01"))
      #j <- as.integer(substr(as.character(j), 6, 7))
      j <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
             3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
             5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
             6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
             7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
             8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
             9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 
             10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
             11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 
             12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12)
      jdays <- which(j==month)
      
      # average ETSR across all days of the month
      mean(ETSR(latitude, jdays))
}


# Hargreaves equation for evapotranspiration, per Hargreaves & Samani 1985
# (Shuttleworth 1993 seems to incorrectly omit the exponent)
hargreaves <- function(S0, tmean, tmin, tmax){
      
      #ETP, in mm/day
      0.0023 * S0 * (tmax - tmin)^.5 * (tmean + 17.8)
}


# calculate monthly and annual water balance variables for one site
hydro <- function(latitude, # integer
                  ppt, tmean, tmax, tmin){ # vectors of length 12
      
      if(is.na(tmean[1])) return(rep(NA, 5))
      
      # [note: Wang et al 2012 page 21 use a latitude correction;
      # this note is a placeholder for that calculation if we decide to use it]
      
      # get S0 values for all 12 months
      S0 <- sapply(1:12, monthly_S0, latitude=latitude)
      
      # monthly water balance variables
      pet <- hargreaves(S0, tmean, tmin, tmax) * 
            c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      pet[tmean<0] <- 0 # per Wang et al 2012
      aet <- pmin(pet, ppt)
      cmd <- pmax(0, pet - ppt)
      rr <- pmax(0, ppt - pet)
      
      # annual sums
      return(c(PPT=sum(ppt), PET=sum(pet), AET=sum(aet), CWD=sum(cmd), RAR=sum(rr)))
}


# calculate annual water balance variables for a raster stack
water_balance <- function(rasters, # stack of 49 rasters: ppt1-12, tmean1-12, tmax1-12, tmin1-12, latitude
                          temp_scalar=1, # multiplier to convert input temperatures to deg C
                          ppt_scalar=1, # multiplier to convert input precipitation to mm
                          ncores=1){ # number of computing cores to use for parallel processing
      
      # compute annual water balance variables
      w <- function(x, ...) hydro(latitude=x[49], 
                                  ppt=x[1:12], 
                                  tmean=x[13:24],
                                  tmax=x[25:36],
                                  tmin=x[37:48])
      if(temp_scalar != 1 | ppt_scalar != 1) w <- function(x, ...) hydro(latitude=x[49], 
                                                       ppt=x[1:12] * ppt_scalar, 
                                                       tmean=x[13:24] * temp_scalar,
                                                       tmax=x[25:36] * temp_scalar,
                                                       tmin=x[37:48] * temp_scalar)
      
      beginCluster(ncores, type="SOCK")
      wb <- clusterR(rasters, calc, args=list(fun=w), 
                     export=c("ETSR", "hargreaves", "hydro", "monthly_S0"))
      endCluster()
      
      names(wb) <- c("PPT", "PET", "AET", "CWD", "RAR")
      return(wb)
}


# generate a raster layer with call values representing degrees latitude
latitude <- function(template, # file path to a raster layer
                     already_latlong=FALSE){ # leave as F unless rasters are already in lat-long projection
      
      lat <- raster(template)
      lat <- lat * 1 # force into RAM
      coords <- coordinates(lat)
      
      if(!already_latlong){
            coords <- as.data.frame(coords)
            coordinates(coords) <- c("x", "y")
            crs(coords) <- crs(lat)
            coords <- spTransform(coords, crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
            coords <- coordinates(coords)
      }
      
      lat[] <- coords[,2]
      return(lat)
}

