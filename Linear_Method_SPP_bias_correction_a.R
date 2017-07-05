
# This program has been developed by: -----------------------------------

# Begum Rushi, Regional Associate, HKH Region -----------------------------

# For any query, please feel free to contact: begumrabeya.rushi@nasa.gov -------------



rm(list = ls())

library("optparse","rgdal")

#!/usr/bin/env Rscript
option_list = list(
  make_option(c("-i", "--chirps"), type="character", default="C:/biased_correction/chirps/", 
              help="input folder directory [default= %default]", metavar="character"),
  make_option(c("-l", "--lat_lon"), type="character", default="C:/biased_correction/NCDC_point_available.csv", 
              help="observed folder directory [default= %default]", metavar="character"),
  make_option(c("-o", "--observed"), type="character", default="C:/biased_correction/observed/", 
              help="observed folder directory [default= %default]", metavar="character"),
  make_option(c("-Y", "--Year"), type="integer", default="2015", 
              help=" starting year[default= %default]", metavar="integer"),
  make_option(c("-N", "--End_Year"), type="integer", default="2015", 
              help=" End year[default= %default]", metavar="integer"),
  make_option(c("-p", "--comparison"), type="character", default="C:/biased_correction/chirps_obs_daily_comparison/", 
              help="output comparison directory [default= %default]", metavar="character"),
  make_option(c("-c", "--output"), type="character", default="C:/biased_correction/corrected_chirps/", 
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
# Necessary Library -------------------------------------------------------

library(raster)
library(rgdal)
library(sp)



# Declare: Time Extent ----------------------------------------------------

Year<- c(opt$Year:opt$End_Year)

##"Ganges","Indus","Meghna")



# Step1: Loading NCDC_point for River Basin ---------------------------------------


NCDC_points<- read.csv(opt$lat_lon,header = T) # read the station names

NCDC_co_ordinate<-data.frame(NCDC_points$Lon,NCDC_points$Lat) # Only Lat and lon have been kept. 

colnames(NCDC_co_ordinate)<- c("Lon","Lat")

NCDC_co<-data.matrix(NCDC_co_ordinate,rownames.force = T)

station_point<-nrow(NCDC_co)

NCDC_points<-as.matrix(NCDC_points) 



# Step 2: Extracting cell information of chirps for observed value comparison -------------------------------

for (y in 1:length(Year)) {
  
  date_begin<-paste(Year[y],"-1-1",sep = "")
  date_end<-paste(Year[y],"-12-31",sep = "")
  
  date1<- seq(as.Date(date_begin), as.Date(date_end), "days")
  date_seq<- gsub("-",".",date1)
  head(date_seq)
    
  for(d in 1:length(date_seq)) {
    
    # station_chirps<- matrix(0,no_days,4)
    st_chirps<-data.frame(matrix(NA, nrow = station_point, ncol = 5))
    colnames(st_chirps)<-c("StationName","Lon","Lat","CHIRPSPrcp","ObsPrcp")
    st_chirps$StationName<-NCDC_points[,1]
    
    
    
    inputfile<- paste(opt$chirps,Year[y],"/",date_seq[d],".tif",sep = "") 
    
    chirps_global_crop<- raster(inputfile)
    
 
    
  st_chirps[,2:4]<-data.frame(coordinates(NCDC_co),extract(chirps_global_crop,NCDC_co))
    
    
    for(p in 1:station_point) {   
      

     obsinfolder<-paste(opt$observed,NCDC_points[p,1],"/",Year[y],sep="")
    
    
     inputfile<-paste(obsinfolder,"/",Year[y],"_daily.csv",sep="")
    
     if (file.exists(inputfile)){
     obs_data<- read.csv(paste(inputfile,sep = ""),header = T)
    
     if (length(which(obs_data[,1]==gsub("[.]","",date_seq[d])))>0){
     
     st_chirps[p,5]<-obs_data[which(obs_data[,1]==gsub("[.]","",date_seq[d])),4]
     } else {print(paste("Observed data for ",NCDC_points[p,1]," for day ",date_seq[d]," is not available",sep=""))}
     }
     
     
     st_chirps[is.na(st_chirps)] <- 0
     
     chirps_obs_prcp<- as.data.frame(st_chirps)
     
     
     chirps_obs_prcp_final<-data.frame(chirps_obs_prcp$CHIRPSPrcp,chirps_obs_prcp$ObsPrcp) 
     
  
     chirps_obs_prcp_final$chirps_obs_prcp.CHIRPSPrcp <- ifelse(chirps_obs_prcp_final$chirps_obs_prcp.CHIRPSPrcp<1,0,chirps_obs_prcp_final$chirps_obs_prcp.CHIRPSPrcp)
     
     chirps_obs_prcp_final$chirps_obs_prcp.ObsPrcp <- ifelse(chirps_obs_prcp_final$chirps_obs_prcp.ObsPrcp<1,0,chirps_obs_prcp_final$chirps_obs_prcp.ObsPrcp)
     
     colnames(chirps_obs_prcp_final) <-c("chirps_prcp","obs_prcp")
     
    } 
  
    chirps_obs_outfolder<-paste(opt$comparison,Year[y],"/",sep="")
    
    if(!file.exists(chirps_obs_outfolder))dir.create(chirps_obs_outfolder)

    chirps_obs_outfile<-paste(chirps_obs_outfolder,date_seq[d],"_chirps_obs_prcp.csv",sep="")


    write.table(chirps_obs_prcp_final, chirps_obs_outfile,row.names=FALSE, na="",col.names=T, sep=",")
    

  }
}
    


# Step 3: Monthly Mean Component ----------------------------------------------



mon<-c(1:12) #1:12 # Number of Months

month_names<-c("Jan","Feb","Mar","Apr","May","Jun","July","Aug","Sep","Nov","Dec")

for (m in 1:length(mon)) {
  myList <- list()
  

  for (y in 1:length(Year)) {
    
    
    inputfolder<- paste(opt$comparison,Year[y],sep = "")
    setwd(inputfolder)      
    
    files <- list.files(path=".",pattern = paste(Year[y],".",sprintf("%02d", mon[m]),".*",sep="")) #path=inputfolder,
    
    for(i in files) { myList[[length(myList)+1]]<-as.matrix(read.csv(i)) }
    
    
    
  }
  
  MonthlyMean<-Reduce("+", myList)/length(myList)
  
  mean_monthly<-as.data.frame(MonthlyMean)
  
  mean_monthly$Ratio<- mean_monthly$obs_prcp/mean_monthly$chirps_prcp
  
  
  mean_monthly_garbage_removal <- mean_monthly[!is.infinite(mean_monthly$Ratio),]
  
  mean_monthly_final <- mean_monthly_garbage_removal[!is.na(mean_monthly_garbage_removal$Ratio),]
  
  monthly_bias_factor<- as.numeric(mean(mean_monthly_final$Ratio))
  
  
  
  assign(paste(mon[m],"monthly_obs_bias_factor", sep = "_"), monthly_bias_factor)
  
}



# Step 4: Correction of CHIRPS with observed Precipitation ----------------


for (y in 1:length(Year)) {
  
  for (m in 1:length(mon)) {
  
    inputfolder_chirps<- paste(opt$chirps,Year[y],"/",sep = "")
    setwd(inputfolder_chirps)      
    
  
    
    files_chirps <- list.files(path=".",pattern = paste(Year[y],".",sprintf("%02d", mon[m]),".*",sep="")) #path=inputfolder,
    
   
    
     for(i in files_chirps) { 
       
       require(raster)
      
       directory<- paste(inputfolder_chirps,i,sep = "")
       
       chirps_daily<- as.matrix(raster(directory))
    
    
      
      monthly_factor<- as.numeric(lapply(paste(mon[m],"monthly_obs_bias_factor", sep = "_"),get))
      
      corrected_chirps<-data.matrix(chirps_daily*monthly_factor,rownames.force = T)
      

    rb <- raster(corrected_chirps)
    class(rb)
    
    # replace with correct coordinates
    extent(rb) <- c(82,98,23.75,31.5)
    
    
    correct_folder<-paste(opt$output,Year[y],"/",sep="")

    if(!file.exists(correct_folder))dir.create(correct_folder)
    setwd(correct_folder)


    correct_chirps_file<-paste(correct_folder,i,sep="")

    writeRaster(rb, filename=correct_chirps_file, format="GTiff", overwrite=TRUE)
    
    
    
              }
    
    
    
  }}




