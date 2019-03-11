##############################################
#training code file
#this uses data downloaded from the Atlas of Living Australia related to the amphibians there - all data downloaded from Australia and then cleaned so that only the last ten years are left - specific species used is Litoria aureas
#for environmental layers the 19 bicolimatic layers from bioclim have been downloaded as well as a landuse, soil, population density (humans), elevation and NDVI layers. 
#the environmental layers are first cut into shape and cell size and extent is set
#dropbox -> Maxent_training (on my own computer - users will have to change directories and file names to the relevant ones)


#tells you where the dismo package is located and where you need to place the maxent file
system.file("java", package="dismo")

 #MAXENT run for single species including set up of training region
#packages needed
library(rtiff)
library(rgdal)
library(raster)
library(rgeos)
library(sm)
library(maptools)
library(sp)
library(dismo)
library(rJava)
library(vegan)
library(usdm)
library(MASS)
library(magrittr)


##########################################################################################
###########################Preparation of the raster layers###############################
##########################################################################################
#filenames <- list.files("Full_bioclim/", pattern="*.tif", full.names=TRUE)
#ldf <- lapply(filenames, raster)

#before running loop read in Australian boundary
boundary <- readShapePoly("Australia_boundary.shp")
evi <- raster('raster_files/evi.asc')

k <-seq(from = 1, to = 19, by = 1)
k<-data.frame("bioclim"=k)
p<-c(39)
p<-rep(p,times=19,each=1)
p<-data.frame("tile1" = p)
q<-c(49)
q<-rep(q,times=19,each=1)
q<-data.frame("tile2" = q)
t<-c(310)
t<-rep(t,times=19,each=1)
t<-data.frame("tile3" = t)
s<-c(410)
s<-rep(s,times=19,each=1)
s<-data.frame("tile4" = s)
o<-c(311)
o<-rep(o,times=19,each=1)
o<-data.frame("tile4" = o)
v<-c(411)
v<-rep(v,times=19,each=1)
v<-data.frame("tile4" = v)
i = 1

#loop to automatically read in, merge, clip and write out the ascii files

for (i in 1:19){
  r <- raster(paste("Full_bioclim/bio",k[i,],"_",p[i,],".tif",sep=""))#read tile 1
  r2 <-raster(paste("Full_bioclim/bio",k[i,],"_",q[i,],".tif",sep=""))#read tile 2
  r3 <-raster(paste("Full_bioclim/bio",k[i,],"_",t[i,],".tif",sep=""))#read tile 3
  r4 <-raster(paste("Full_bioclim/bio",k[i,],"_",s[i,],".tif",sep=""))#read tile 4
  r5 <-raster(paste("Full_bioclim/bio",k[i,],"_",o[i,],".tif",sep=""))#read tile 5
  r6 <-raster(paste("Full_bioclim/bio",k[i,],"_",v[i,],".tif",sep=""))#read tile 6
  rc1 <- merge(r,r2,r3,r4,r5,r6)#merge the tiles
  rc1 <- crop(rc1, extent(evi))#crop according to boundary file
  clip <-mask(rc1,boundary)#clip along the boundary (Austrlian coastline)
  projection(clip) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # coodinate reference system for WGS84
  rf <-writeRaster (clip, filename = paste("Ascii/bio",k[i,],sep=""), format = "ascii", overwrite = TRUE)
}


########################################################################
######FOLLOW THESE STEPS IF THE CELL SIZE NEEDS TO BE CHANGED######
#Set extent and size
#specify the number of pixels for an empty raster file (use the data from whatever file is already in the correct format)
r <- raster(ncol = 5204, nrow = 4100)
#set extent again
bb <- extent(111.3640000000000043,154.7307000000000130,-44.2250000000000014,-10.0583070714834619)
extent(r) <- extent(bb)

#loop to automate changing the bioclimatic files
for (i in 1:19){
  #read in the ascii file that needs to be changed
  rp2 <- raster(paste("Ascii_2/bio",k[i,],".asc",sep=""))
  #Specify extent of raster
  #resample the file that needs the pixel size changed to the empty raster
  res3 <- resample(rp2, r)
  #set extent of file
  extent(res3) <- extent(bb)
  rf <-writeRaster (res3, filename = paste("Ascii/bio",k[i,],sep=""), format = "ascii", overwrite = TRUE)
}



############################################################################################
#########################Collinearity analysis of raster layers#############################
############################################################################################
#read in raster layers into a stack

PCAstack <- stack(
  raster('Ascii/bio1.asc'),
  raster('Ascii/bio2.asc'),
  raster('Ascii/bio3.asc'),
  raster('Ascii/bio4.asc'),
  raster('Ascii/bio5.asc'),
  raster('Ascii/bio6.asc'),
  raster('Ascii/bio7.asc'),
  raster('Ascii/bio8.asc'),
  raster('Ascii/bio9.asc'),
  raster('Ascii/bio10.asc'),
  raster('Ascii/bio11.asc'),
  raster('Ascii/bio12.asc'),
  raster('Ascii/bio13.asc'),
  raster('Ascii/bio14.asc'),
  raster('Ascii/bio15.asc'),
  raster('Ascii/bio16.asc'),
  raster('Ascii/bio17.asc'),
  raster('Ascii/bio18.asc'),
  raster('Ascii/bio19.asc'),
  raster('raster_files/gridsoil.asc'),
  raster('raster_files/evi.asc'),
  raster('raster_files/landuse.asc'),
  raster('raster_files/popd10.asc'),
  raster('raster_files/altitude_F.asc')
)

#Calculate collinearity


#shows all of the vif values
vif(PCAstack)
#vifcor, first find a pair of variables which has the maximum linear correlation (greater than th), and exclude one of them which has greater VIF. The procedure is repeated untill no variable with a high corrrelation coefficient (grater than threshold) with other variables remains
v1 <- vifcor (PCAstack, th = 0.9)
#vifstep calculate VIF for all variables, exclude one with highest VIF (greater than threshold), repeat the procedure untill no variables with VIF greater than th remains.
v2 <- vifstep (PCAstack, th = 10)

#update the raster stack to only include the layers which are not collinear using the vif method you have decided upon. Base the usage on ecological intuition rather than solely relying on the statistics (i.e. if you think that altitude is highly important and one of the methods excluded it it might be better to pick the one that included it)
environstack <- exclude (PCAstack, v1)

############################################################################################
#########################Creation of kernel density for training#############################
############################################################################################

#create a master file for kernel density, set the working directory to the folder of the files first 
#we use all of the CSV files downloaded from all amphibian species in Australia, these files have been cleaned and only the records from the last ten years are used
files = list.files(pattern="*.csv")
myfile = do.call(rbind, lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE)))

write.csv(myfile,"master_occurence.csv")

#as it is likely that people went looking for frogs would record more than just their own species we use all of the records we have to create a kernel density layer for sampling effort, it will weigh the pseudoabsence points as being more likely to be true absences if they are in an area with a high sampling effort value

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#put the occurences on a grid, it is improtant that they have the same resolution as the climate data that we are using!
#read in one of the raster files
evi <- raster("Current_Large/evi.asc")

#prepare the master file for use
master <-myfile[,c(3,2)]

#rasterize the occurence records according to the cimate layer (giving > or = 1 the value of 1 and anything below 0)
occur.ras<-rasterize(master,evi,1)
plot(occur.ras)

#needed if you have downloaded the data fresh online
#read in the polygon of the shapefile outlining the area we are either projecting to
#states<-readShapePoly("H:/Shapefiles/US_States/states.shp",proj4string=occur.ras@crs)
#mask the area (this will remove any points that are outside the shape file thatw e are not interested in)
#the code above is not needed in this case as we already read in our file earlier and it is called "boundary" 

occur.boundary<-mask(occur.ras,boundary) %>% crop(.,boundary)



#perform the kernel density 
presences<-which(values(occur.boundary)==1)
pres.locs<-coordinates(occur.boundary)[presences,]
dens<-kde2d(pres.locs[,1],pres.locs[,2],n=c(nrow(occur.boundary),ncol(occur.boundary)))
dens.ras<-raster(dens)
plot(dens.ras)



############################################################################################
#########################Creation of buffer zone for training###############################
############################################################################################

#Creation of prediction region, in this case the whole of Borneo
#First read in one of the created layers that you know all values are above 0 for the region you want to use as mask
rasmask <- raster ("raster_files/bio7.asc")
#reclassify all values above 0 as 1 and all NA's as NA
rc2 <- reclassify(rasmask, c(NA,NA,NA, 0.0,Inf,1))
#check that the file plotted correctly
rc2
#write out the file in its own folder to separate prediction region and training region
dir.create(paste("Prediction",sep=""), recursive=TRUE)
setwd(paste("Prediction",sep=""))
#please note that the reason for the sperate folders is that the training and prediction region has to have the SAME name and come in the same order in the raster stack, as they have the same name ensure you are able to separate between the two adequately in another manner
rf <-writeRaster (rc2, filename = "mask", format = "ascii", overwrite = TRUE)
#set the directory back to main folder
setwd("../..")

#Now we will create the training region, this is a smaller area and can be based on areas where the species is known to occur and have been looked for extensively (so that background points that the species is not found in are more likely to actually be absences) or it can be an area that we just "make up" as we have so few points for the species and it hasn't been looked for enough so we will create a buffer region around the known points of occurence
#Read in the occurrence points
species <- read.csv("Litoria_area_areas.csv", head = TRUE) #this needs to be the CSV file of the species we are interested in predicting for, in our case we use Litora aurea from Australia within the last decade
head(species)

#Create dataframe
n1 <- length(species$LONG)

df1 <- data.frame(ID = 1:n1, Value = 1)

#bind together the original species data into the new dataframe 
speciesdf <- cbind(species,df1)

#specify the coordinates to be used
coordinates(speciesdf) <- c("LONG","LAT")

#Specify CRS otherwise we cannot perform the transformation to UTM
projection(speciesdf) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # coodinate reference system for WGS84

#Transform to UTM CRS
speciesdf <- spTransform(speciesdf, CRS("+init=epsg:29850"))

################################################################################
#Create buffers 
################################################################################
#
#specify the width of the buffer, currently set to 25 km
buff <- gBuffer(speciesdf, width = 25000) 
#you can check the class of the file you are working with by
class(buff)
#to change into spatialpolygons the following code would be used (its alreay there though)
buff <- spTransform(buff, CRS("+proj=longlat +datum=WGS84"))

################################################################################
#Transform to spatial polygon to SP data frame 
################################################################################


#Extract polygon ID's
buff_ID <- sapply(slot(buff, "polygons"), function(x) slot(x, "ID"))
#Create dataframe with correct rownames
buff_df <- data.frame( ID=1:length(buff), row.names = buff_ID)
#Coerce into sSPDF
buff <- SpatialPolygonsDataFrame(buff, data=buff_df)
#Check that it is the correct class
class(buff)
#Add column giving all the buffer areas a value of 1
buff@data$Value <- 1 
#Write out as shapefile
writeOGR(buff,paste("buff25",sep=''),paste("250km",sep=""),driver="ESRI Shapefile",overwrite=TRUE)

################################################################################
#Transform file to raster
################################################################################
#Specify the number of pixels in the raster file ####these were already set earlier for some other files so r and bb are already correct, this is just highlighting that they are getting used again
r <- raster(ncol = 5204, nrow = 4100)
#set extent again
bb <- extent(111.3640000000000043,154.7307000000000130,-44.2250000000000014,-10.0583070714834619)
extent(r) <- extent(bb)


#read in shapefile
buff_sp <- readOGR(paste("buff250",sep=''),paste("250km",sep=""))
#Convert shape file to raster
buff_raster <- rasterize(buff_sp, r,"Value")

################################################################################
#Clip file along the Australian coast
################################################################################

#we will use the boundary file that was read in earlier again

#Identifying pixels that lie within the borders
clip <-mask(buff_raster,boundary)

#reclassify file to 1 within border and NA without
rc2 <- reclassify(clip, c(NA,NA,NA, 0.00,Inf,1))
#set extent of the file
extent(rc2) <- bb
rc3 <- setExtent(rc2, bb, keepres=TRUE)
#write out raster ensuring it goes in the correct folder
#dir.create(paste("Training",sep=""), recursive=TRUE)
#setwd(paste("Training",sep=""))
rf <-writeRaster (rc3, filename = "mask", format = "ascii", overwrite = TRUE)

###############################################################################
###############################################################################
#################### M A X E N T ##############################################
###############################################################################
###############################################################################


# get species records read in
Current <- read.csv("Litoria_aurea_areas.csv", head = T)
#Ensure the values needed go in where they need to be
Current <- Current[,c(3,2)]



##############################################################################
# Training stacks
#############################################################################

# using all layers
# this stack has a small masking for training (It has to be first in the list!!!!!). Use the layers that the collinearity analysis told you to keep

trainingStack <- stack(
  raster(paste("Training/mask.asc",sep="")),
  environstack
)

# set the coordinate reference system for the raster stacks... this is not absolutely necessary if the rasters are unprojected (e.g., WGS84), but we'll do it to avoid warning messages below

projection(trainingStack) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # coodinate reference system for WGS84


##############################################################################
# Prediction stack
#############################################################################

# using all layers
# this stack has a large masking for prediction (It has to be first in the list!!!!!). Use the layers that the collinearity analysis told you to keep

predStack <- stack(
  raster(paste("Prediction/mask.asc",sep="")),
  environstack
)

# set the coordinate reference system for the raster stacks... this is not absolutely necessary if the rasters are unprojected (e.g., WGS84), but we'll do it to avoid warning messages below

projection(predStack) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # coodinate reference system for WGS84


##############################################################################
#Random point creation
#############################################################################

#read in the training mask
mask <- raster(paste("Training/mask.asc",sep=""))
projection(mask)<- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # coodinate reference system for WGS84
#Random background points to base the bankground predictions on, the value can be set to what you want but the default and minimum value should be 10000. This is a different way to using the kernel density mask
pointCurrent <- randomPoints(mask, 10000)

#if you call pointCurrent it will tell you the coordinates of all of the background points that will be used
pointCurrent

##############################################################################
#Maxent training
#############################################################################

#might need to run the code .jinit() if the Java Virtual Machine does not start

# train 
trainMod <- maxent(
  x=trainingStack,
  p=Current,
  a=pointCurrent,
  path=paste("25kmtrain",sep=""),
  biasfile=dens.ras,
  args=c(
    'randomtestpoints=30',
    'betamultiplier=1',
    'linear=true',
    'quadratic=true',
    'product=true',
    'threshold=true',
    'hinge=true',
    'threads=2',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false'
  )
)



##############################################################################
#Prediction maps
#############################################################################



# predict small extent 25 model to entire extent 
mapExtent <- predict(
  object=trainMod,
  x=predStack,
  filename=paste("25kmtrain",sep=""),
  na.rm=TRUE,
  format='GTiff',
  overwrite=TRUE
)

predmap <-raster("25kmtrain.tif")


