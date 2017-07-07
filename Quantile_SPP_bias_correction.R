# This program has been developed by: -----------------------------------

# Begum Rushi, Regional Associate, HKH Region -----------------------------

# For any query, please feel free to contact: begumrabeya.rushi@nasa.gov -------------


rm(list = ls())

library("optparse","rgdal")

#!/usr/bin/env Rscript
option_list = list(
  make_option(c("-o", "--corrected_chirps"), type="character", default="C:/bias_correction_SPP/Test_Data_Set/Nyando/chirps/", 
              help="input folder directory [default= %default]", metavar="character"),
  make_option(c("-i", "--persian"), type="character", default="C:/bias_correction_SPP/Test_Data_Set/Nyando/persian/", 
              help="observed folder directory [default= %default]", metavar="character"),
  make_option(c("-Y", "--Year"), type="integer", default="2015", 
              help=" starting year[default= %default]", metavar="integer"),
  make_option(c("-N", "--End_Year"), type="integer", default="2015", 
              help=" End year[default= %default]", metavar="integer"),
  make_option(c("-c", "--corrected_output"), type="character", default="C:/bias_correction_SPP/Test_Data_Set/Nyando/corrected_persian_quantile/", 
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


library(raster)
require(abind)
library(MASS)
library(pscl)
library(EDISON)
library(MCMCpack)
library(invgamma)


years<-c(opt$Year:opt$End_Year)
mon<-c(1)
month_names<-c("Jan","Feb","Mar","Apr","May","Jun","July","Aug","Sep","Nov","Dec")


   
      
for (p in 1:length(mon)){
        myList <- list()
        counter<-0
        counter2<-0
        for (y in 1:length(years)) {

 # Loading CHIRPS and SPP --------------------------------------------------     
          
ObsInputfolder<- paste(opt$corrected_chirps,years[y],"/",sep = "")


SatInputfolder<- paste(opt$persian,years[y],"/",sep = "")
    


      
files_chirps <- list.files(path=ObsInputfolder,pattern = paste(years[y],".",sprintf("%02d", mon[p]),".*",sep="")) #path=inputfolder,

files_sat <- list.files(path=SatInputfolder,pattern = paste(years[y],".",sprintf("%02d", mon[p]),".*",sep="")) #path=inputfolder,




for(i in files_chirps) { 
  counter<-counter+1
  print(counter)
  ObsDirectory<- paste(ObsInputfolder,i,sep = "")
  if (counter==1){
    
    
    chirps<-as.matrix(raster(ObsDirectory))
  } else {
    
    chirps<-cbind(chirps,as.matrix(raster(ObsDirectory)))
  }
  #myList[[length(myList)+1]]<-as.matrix(raster(i)) 

}


for(j in files_sat) { 
  counter2<-counter2+1
  print(counter2)
  SatDirectory<- paste(SatInputfolder,j,sep = "")
  if (counter2==1){
    sat_prcp<-as.matrix(raster(SatDirectory))
    Sat_names<- as.list(j)
  } else {
    
    sat_prcp<-cbind(sat_prcp,as.matrix(raster(SatDirectory)))
    
    Sat_names<- cbind(Sat_names,as.list(j))
    
    
  }
  
}

} 


# Making 3-D Matrix for Each Month Considering All Years ------------------


dim(chirps) <- c(dim(raster(SatDirectory))[1], dim(raster(SatDirectory))[2], counter) # 3-D arrray formation

#dim(chirps) <- c(dim(raster(SatDirectory))[1], dim(raster(SatDirectory))[2], length(files_chirps)) # converting it to 3d matrix 
Drizzle<-1 # less than 1 mm rain is considered drizzle
chirps[which(chirps<Drizzle)]<-0

chirps

dim(sat_prcp) <- c(dim(raster(paste(SatInputfolder,j,sep = "")))[1], dim(raster(paste(SatInputfolder,j,sep = "")))[2], counter2) # converting it to 3d matrix 
Drizzle<-1 # less than 1 mm rain is considered drizzle
sat_prcp[which(sat_prcp<Drizzle)]<-0

#CHIRPS<-array(3:63, dim=c(3,4,5))


#z <- array(1:60, dim=c(3,4,5))
#z[1,1,1]<-0 ## for testing

#GammaCDF <- array(0, dim=c(3,4,5))
#BCz<- array(0, dim=c(3,4,5)) # for storing bias corrected data

GammaCDF_chirps <- array(0, dim=c(dim(raster(ObsDirectory))[1], dim(raster(ObsDirectory))[2],counter))
GammaCDF_sat<- array(0, dim=c(dim(raster(SatDirectory))[1], dim(raster(SatDirectory))[2], counter2)) # for storing bias corrected data


#GammaParms<-matrix(0,dim(raster(i))[1],dim(raster(i))[2])
#library(fitdistrplus)

###CHIRPS param ###
CHIRPSParmsLambda<-matrix(0,dim(raster(ObsDirectory))[1],dim(raster(ObsDirectory))[2])
CHIRPSParmsTheta<-matrix(0,dim(raster(ObsDirectory))[1],dim(raster(ObsDirectory))[2])
### 
GammaParmsLambda<-matrix(0,dim(raster(SatDirectory))[1],dim(raster(SatDirectory))[2])
GammaParmsTheta<-matrix(0,dim(raster(SatDirectory))[1],dim(raster(SatDirectory))[2])



 for (m in seq_len(dim(sat_prcp)[1])) {
      for (n in seq_len(dim(sat_prcp)[2])) {
        #print(i)
        ### CHIRPS
        IndexNonZeroCHIRPS<-which(chirps[m,n,]>0)
        NonZeroCHIRPS<-chirps[m,n,][which(chirps[m,n,]>0)]
        #### Satellite
        IndexNonZeroSat<-which(sat_prcp[m,n,]>0)
        NonZeroSat<-sat_prcp[m,n,][which(sat_prcp[m,n,]>0)]
       
        
        if (length(IndexNonZeroCHIRPS)>5 & length(IndexNonZeroSat)>5 & length(unique(NonZeroCHIRPS)) >5 & length(unique(NonZeroSat))>5 )
        {
        CHIRPSParmsLambda[m,n]<-fitdistr(NonZeroCHIRPS, "gamma")$estimate[1] #lambda OR SHAPE
        CHIRPSParmsTheta[m,n]<-fitdistr(NonZeroCHIRPS, "gamma")$estimate[2] #theta or rate
      

        
        #GammaParmsLambda[i,j]<-fitdistr(z[i,j,], "gamma")$estimate[1] #lambda
        #GammaParmsTheta[i,j]<-fitdistr(z[i,j,], "gamma")$estimate[2] #theta
        GammaParmsLambda[m,n]<-fitdistr(NonZeroSat, "gamma")$estimate[1] #lambda
        GammaParmsTheta[m,n]<-fitdistr(NonZeroSat, "gamma")$estimate[2] #theta

        
        #GammaCDF[i,j,]<-pgamma(z[i,j,], GammaParmsLambda[i,j], rate = GammaParmsTheta[i,j], log = FALSE)
        NonZeroGammaCDF<-pgamma(NonZeroSat, GammaParmsLambda[m,n], rate = GammaParmsTheta[m,n], log = FALSE)
        #print(NonZeroGammaCDF)
        
        #BCz[i,j,] <-qgamma(GammaCDF[i,j,],CHIRPSParmsLambda[i,j], CHIRPSParmsTheta[i,j])
        GammaCDF_sat[m,n,IndexNonZeroSat]<-qgamma(NonZeroGammaCDF,CHIRPSParmsLambda[m,n], CHIRPSParmsTheta[m,n]) #inverse
        }else {
          print(NonZeroSat)
          GammaCDF_sat[m,n,IndexNonZeroSat]<- NonZeroSat ## no bias correction is done if only 2 points are available 
          
       
        }
        

        
      }
 }

for (kk in 1:counter2){
  
corrected_spp_daily <- as.matrix(GammaCDF_sat[,,kk])

 rb <- raster(corrected_spp_daily)
 class(rb)

 # replace with correct coordinates
 extent(rb) <- c(34.75,36,-0.5,0.1)
 
 
 correct_folder<-paste(opt$corrected_output,sep="")

if(!file.exists(correct_folder))dir.create(correct_folder)

 correct_spp_file<-paste(correct_folder,Sat_names[kk],sep="")
 
 
print(correct_spp_file)
 
 writeRaster(rb, filename=correct_spp_file, format="GTiff", overwrite=TRUE)
 
  }
 
}
 


