
# This program has been developed by: -----------------------------------

# Begum Rushi, Regional Associate, HKH Region -----------------------------

# For any query, please feel free to contact: begumrabeya.rushi@nasa.gov -------------


rm(list = ls())

library("optparse")

#!/usr/bin/env Rscript
option_list = list(
  make_option(c("-o", "--corrected_chirps"), type="character", default="C:/biased_correction/corrected_chirps/", 
              help="input folder directory [default= %default]", metavar="character"),
   make_option(c("-i", "--persian"), type="character", default="C:/biased_correction/persian/", 
              help="observed folder directory [default= %default]", metavar="character"),
  make_option(c("-Y", "--Year"), type="integer", default="2015", 
              help=" starting year[default= %default]", metavar="integer"),
  make_option(c("-N", "--End_Year"), type="integer", default="2015", 
              help=" End year[default= %default]", metavar="integer"),
  make_option(c("-c", "--corrected_output"), type="character", default="C:/biased_correction/corrected_persian/", 
              help="output folder directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if (is.null(opt$Year)){
  print_help(opt_parser)
  stop("All arguments must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$End_Year)){
  print_help(opt_parser)
  stop("All arguments must be supplied (input file).n", call.=FALSE)
}



Year<-c(opt$Year:opt$End_Year)

mon<-c(1:12) #1:12 # Number of Months

month_names<-c("Jan","Feb","Mar","Apr","May","Jun","July","Aug","Sep","Nov","Dec")

# Step 3: Make mean monthly matrix for chirps and SPP --------------------------

data_names<-c(opt$corrected_chirps,opt$persian)


for (d in 1:length(data_names)) {
  
  for (m in 1:length(mon)){
    myList <- list()
    for (y in 1:length(Year)) {
      
      
      inputfolder<- paste(data_names[d],Year[y],"/",sep = "")
      
      setwd(inputfolder)      
      
      files <- list.files(path=".",pattern = paste(Year[y],".",sprintf("%02d", mon[m]),".*",sep="")) #path=inputfolder,
      
      for(i in files) { 
        
        
        require(raster)
        
        directory<- paste(inputfolder,i,sep = "")
        
        myList[[length(myList)+1]]<-as.matrix(raster(directory)) }
      
      
      
    }
    MonthlyMean<-Reduce("+", myList) / length(myList)
    mean_monthly<-as.matrix(MonthlyMean)

    rb <- raster(mean_monthly)
    class(rb)

    # replace with correct coordinates
    extent(rb) <- c(82,98,23.75,31.5)
    
    assign(paste(mon[m],data_names[d],"monthly_mean.tif", sep = "_"), rb)
    
    monthly_folder<-data_names[d]
    
    setwd(monthly_folder)
    
    
    monthly_tif<-paste(monthly_folder,mon[m],"_monthly_mean.tif",sep="")
    
    writeRaster(rb, filename=monthly_tif, format="GTiff", overwrite=TRUE)
    
    
  }
  
  
}



# Step 2: -----------------------------------------------------------------




for (y in 1:length(Year)) {
  
  for (m in 1:length(mon)) {
    
    inputfolder_spp<- paste(opt$persian,Year[y],"/",sep = "")
    setwd(inputfolder_spp)      
    
    
    
    files_spp <- list.files(path=".",pattern = paste(Year[y],".",sprintf("%02d", mon[m]),".*",sep="")) #path=inputfolder,
    
    
    
    for(i in files_spp) { 
      
      require(raster)
    
      
      
      directory_spp<-paste(inputfolder_spp,i,sep="")
      spp_daily<- raster(directory_spp)
      
      
  
      monthly_chirps_factor<-raster(paste(opt$corrected_chirps,mon[m],"_monthly_mean.tif",sep=""))
      monthly_spp_factor<-raster(paste(opt$persian,mon[m],"_monthly_mean.tif",sep=""))
      
      
      monthly_bias_factor<-overlay(monthly_chirps_factor, monthly_spp_factor, fun=function(x, y){ x/y} )
      
      corrected_spp<-overlay(spp_daily, monthly_bias_factor,fun=function(x, y){ x*y})
      
      
      correct_folder<-paste(opt$corrected_output,Year[y],"/",sep="")

      if(!file.exists(correct_folder))dir.create(correct_folder)
      setwd(correct_folder)


      corrected_file<-paste(correct_folder,i,sep="")
      
      corrected_file
      
      
      writeRaster(corrected_spp, filename=corrected_file, format="GTiff", overwrite=TRUE)

      
      
      
      
    }
    
  }
  
}







