# This script serves to prepare Viirs data for THREDDS uplaoding

viirs2thredds <- function(inPath, inFile, outPath, verbose = FALSE, ...){

  # initial time stamp
  print(Sys.time())	
  
  # declare product lists for extraction
  prodList <- c("chl_oc3", "a_410_qaa", "a_443_qaa", "a_486_qaa", "a_551_qaa", "a_671_qaa", "a_745_qaa", "a_862_qaa", "bb_551_qaa", 
                 "aph_443_qaa", "adg_410_qaa", "c_551_qaa", "Rrs_410", "Rrs_443", "Rrs_486", "Rrs_551", "Rrs_671", "Rrs_745", "Rrs_862",
                 "ndvi", "evi", "pic", "poc", "class_34k_w_owmc", "sst", "qual_sst")
  productUnits <- c("mg m^-3","m^-1", "m^-1", "m^-1", "m^-1", "m^-1", "m^-1", "m^-1", "m^-1", "m^-1", "m^-1", "m^-1", "sr^-1", "sr^-1", 
                    "sr^-1", "sr^-1", "sr^-1", "sr^-1", "sr^-1", "NA", "NA", "mol m^-3", "mg m^-3", "NA", "NA", "NA", "Celsius", "Celsius")
  productName <- c("Chlorophyll Concentration, OC3 Algorithm","Total absorption at 410 nm QAA algorithm", "Total absorption at 443 nm QAA algorithm", 
                   "Total absorption at 486 nm, QAA algorithm", "Total absorption at 551 nm, QAA algorithm", 
                   "Total absorption at 671 nm, QAA algorithm", "Total absorption at 745 nm, QAA algorithm", "Total absorption at 862 nm, QAA algorithm",
                   "Total backscattering at 551 nm, QAA algorithm", "Absorption due to phytoplankton at 443 nm, QAA algorithm", 
                   "Absorption due to gelbstoff and detrital material at 410 nm, QAA algorithm", "Beam attenuation at 551 nm, QAA algorithm",
                   "Remote sensing reflectance at 410 nm", "Remote sensing reflectance at 443 nm", "Remote sensing reflectance at 486 nm", 
                   "Remote sensing reflectance at 551 nm", "Remote sensing reflectance at 671 nm", "Remote sensing reflectance at 745 nm", 
                   "Remote sensing reflectance at 862 nm", "Normalized Difference Vegetation Index", "Enhanced Vegetation Index", 
                   "Calcite Concentration, Balch and Gordon", "Particulate Organic Carbon, D. Stramski, 2007 (443/555 version)", 
                   "OWMC merged Wards and K_means (34K+W) Classification", "Sea Surface Temperature", "Qual Sea Surface Temperature") 

  #--------DATA EXTRACTION & DIMENSIONS-------------
  
  # open the inFile.nc if it exists and return error if it doesn't
  if (file.exists(file.path(inPath, inFile))){
    viirsFile <- ncdf4::nc_open(file.path(inPath,inFile))
  } else {
    stop("inPath/inFile does not exist!")
  }

  # Extract the current file time and write to a dimension in posixct time. Note: this only works for the current storage syntax from basin
  ncyear  <- substr(inFile,2,5)
  jday    <- as.numeric(substr(inFile,6,8))
  ncmon   <- format(as.POSIXct(strptime(jday, "%j"), tz = "GMT", origin = paste0(ncyear, "01-01")), "%m")
  ncday   <- format(as.POSIXct(strptime(jday, "%j"), tz = "GMT", origin = paste0(ncyear, "01-01")), "%d")
  nchr    <- substr(inFile,9,10)
  ncmin   <- substr(inFile,11,12)
  ncsec   <- substr(inFile,13,14)
  ncposix <- as.POSIXct(strptime(paste0(ncyear,"-",ncmon,"-",ncday," ", nchr, ":", ncmin, ":",ncsec), "%Y-%m-%d %H:%M:%S", tz="GMT"))
  
  # Create the now_time dimension
  file_time <- ncdf4::ncdim_def(name='time', units="seconds since 1970-01-01T00:00:00Z", 
                                   vals=as.numeric(ncposix), create_dimvar=TRUE, unlim=TRUE)
  
  # write dimensions
  dim <- viirsFile$dim
  dim$time <- file_time
  
  # Initialize the lists where we are going to store stuff
  dat.list <- list()
  var.list <- list()

  # Extract the data from the input file
  for (j in seq_along(prodList)){
    dat.list[[j]] <- ncdf4::ncvar_get(viirsFile, as.character(prodList[j]))
    var.list[[j]] <- ncdf4::ncvar_def(name     = as.character(prodList[j]), 
                                      units    = as.character(productUnits[j]),
                                      longname = as.character(productName[j]),
                                      dim      = dim,
                                      missval  = -999, 
                                      verbose  = verbose)
  }

  # close the viirsFile
  ncdf4::nc_close(viirsFile)
  
  # Replace 0s with NAs for sst
  for (j in seq_along(prodList)){
    if (prodList[j] == "sst" | prodList[j] == "qual_sst"){
      dat.list[[j]][dat.list[[j]] == 0] <- -999
    }
  }
  
  #--------WRITE NEW .NC FILE-------------
  
  #creating the nc4 output file
  loc.file <- paste0(outPath,substr(inFile,1,14),"_", substr(inFile,19,21),".nc")
  
  #writing all we need to the output file
  loc <- ncdf4::nc_create(filename=loc.file, vars=var.list)
  for(j in seq_along(prodList)){
    ncdf4::ncvar_put(nc=loc, varid=as.character(prodList[j]), vals=dat.list[[j]])
  }
  
  # adding metadata information
  ncatt_put(loc, 0, "Conventions", "CF-1.0")
  ncatt_put(loc, 0, "creator_name", "James Simkins")
  ncatt_put(loc, 0, "creator_email", "simkins@udel.edu")
  ncatt_put(loc, 0, "institution", "University of Delaware Ocean Exploration, Remote Sensing and Biogeography Group (ORB)")
  ncatt_put(loc, 0, "url", "http://orb.ceoe.udel.edu/")
  ncatt_put(loc, 0, "source", "satellite observation NASA MODIS-Aqua instrument")
  ncatt_put(loc, 0, "groundstation", "University of Delaware, Newark, Center for Remote Sensing")
  ncatt_put(loc, 0, "software", 'SeaDas')
  ncatt_put(loc, 0, "createAgency", 'NASA')
  ncatt_put(loc, 0, "product_list", prodList)
  ncatt_put(loc, 0, "summary", "VIIRS ocean color and sst calculation by SeaDas; Regridded to Mercator lon/lat projection. Processed at the University of Delaware")
  ncatt_put(loc, 0, "history", " ")
  ncatt_put(loc, "lon", "long_name", "Longitude")
  ncatt_put(loc, "lat", "long_name", "Latitude")
  ncatt_put(loc, "time", "long_name", "Time")
  
  ncdf4::nc_close(loc)
}

