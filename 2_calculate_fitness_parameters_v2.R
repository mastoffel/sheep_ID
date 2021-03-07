#
# Calculate Fitness Parameters
# Script by Susan Johnston 
# with slight modifications 
# Feb 2019
#

# This scripts processes files from the sheep database into annual fitness data

# Note from MAS: This script results in the file 
# 1_Annual_Fitness_Measures_April_20190501.txt
# which can be found in the Zenodo data repository 
# accompanying our manuscript.

library(ggplot2)
library(plyr)
library(beepr)
library(reshape)
library(lubridate)
library(tidyr)
library(dplyr)
library(lubridate)

options(error = function(){    # Beep on error
  beepr::beep()
  Sys.sleep(1)
}
)

#
library(RJDBC)
#### database connection ######
dbname <- "../sheep/data/db/StKilda_Data.accdb"
driver <- "net.ucanaccess.jdbc.UcanloadDriver"
driverpath <- "../sheep/data/db/UCanAccess/loader/ucanload.jar"
options <- paste0("jdbc:ucanaccess://", dbname, ";memory=false")

con <- DBI::dbConnect(JDBC(driver, driverpath), options)
# src <- src_dbi(con)

# names of tables
tbls <- dbGetTables(con)
# names of variables in a table
flds <- dbGetFields(con, "Sheep")
# get a table 
Sheep <- dbGetQuery(con, "Select * from Sheep")
Census <- dbGetQuery(con, "Select * from CensusData")
Capture <- dbGetQuery(con, "Select * from CaptureData")
tblPreg <- dbGetQuery(con, "Select * from tblPregnancies")

tblSheepCore <- dbGetQuery(con, "Select * from tblSheepCoreNotes")
# censusdata formatting to look like old file
censusdata_new <- dbGetQuery(con, "Select * from CensusData") 
censusdata <- censusdata_new %>% 
                  as_tibble() %>% 
                  dplyr::select(ID:Shelter) %>% 
                  mutate(Date = ymd_hms(Date)) %>% 
                  mutate(Date = date(Date)) %>% 
                  mutate(Date = format(Date, format = "%d/%m/%Y")) %>% 
                  mutate_at(c("Veg", "Act"), replace_na, "") %>% 
                  mutate_if(is.double, as.integer)

consortdata_new <- dbGetQuery(con, "Select * from Consorts")
consortdata <- consortdata_new %>% 
                as_tibble() %>% 
                dplyr::select(-RowRef) %>% 
                mutate(Date = ymd_hms(Date)) %>% 
                mutate(Date = format(Date, format = "%d/%m/%Y")) %>% 
                mutate(Time = ymd_hms(Time)) %>% 
                mutate(Time = format(Time, format = "%H:%M")) %>% 
                mutate_at(c("UnknownTup", "UnknownEwe"), replace_na, "") %>% 
                mutate_if(is.double, as.integer)

capdata_new <- dbGetQuery(con, "Select * from CaptureData")
capdata <- capdata_new %>% 
                as_tibble() %>% 
                dplyr::select(one_of(names(capdata_old)))
write_delim(capdata, "~/Desktop/SoayCaptureData_DB,txt", delim = "\t")

basedata_new <- dbGetQuery(con, "Select * from Sheep") %>% 
                as_tibble()

dbDisconnect(con)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   1. Read in data for calculating fitness.      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

censusdata_old  <- read.table("../sheep/data/db_tables/20190208/20190208_CensusData.txt",       sep = "\t", header = T, stringsAsFactors = F) %>% as_tibble()
consortdata_old <- read.table("../sheep/data/db_tables/20190208/20190208_ConsortData.txt",      sep = "\t", header = T, stringsAsFactors = F)%>% as_tibble()
capdata_old     <- read.table("../sheep/data/db_tables/20190208/20190208_SoayCaptureData.txt",  sep = "\t", quote = "\"", header = T, stringsAsFactors = F)%>% as_tibble()

basedata_old    <- read.table("../sheep/data/db_tables/20190208/20190208_SoayBaseData.txt",     sep = "\t", header = T, stringsAsFactors = F)%>% as_tibble()

oldpink_old     <- read.table("../sheep/data/db_tables/20190208/OPAges_ReproEdit.txt",          sep = "\t", header = T, stringsAsFactors = F)%>% as_tibble()

foetusIDs_old   <- read.table("../sheep/data/db_tables/20190208/20190208_FoetusTable.txt",      sep = "\t", header = T, stringsAsFactors = F)[,1] %>% as_tibble()
pedigree_old    <- read.table("../sheep/data/db_tables/20190208/20190208_Full_Pedigree.txt",    sep = "\t", header = T, stringsAsFactors = F)%>% as_tibble()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   2. Format data                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Census information for each sheep

names(censusdata) <- toupper(names(censusdata))

censusdata <- subset(censusdata, !is.na(ID))
censusdata$DATE <- as.Date(censusdata$DATE, "%d/%m/%Y")

locationdata <- subset(censusdata, select = c(ID, DATE, EASTING, NORTHING))
locationdata$EASTING <- as.numeric(locationdata$EASTING)
locationdata$NORTHING <- as.numeric(locationdata$NORTHING)

censusdata <- subset(censusdata, select = c(ID, PLACE, DATE))


#~~ Consort information

names(consortdata) <- toupper(names(consortdata))

consortdata <- subset(consortdata, select = c(TUPID, PLACE, DATE))
names(consortdata)[1] <- "ID"
consortdata <- subset(consortdata, !is.na(ID))
consortdata$DATE <- as.Date(consortdata$DATE, "%d/%m/%Y")
consortdata$PLACE <- paste0(consortdata$PLACE, ".Consort")

#~~ Add Consort information where males were seen to the census data

censusdata <- rbind(censusdata, consortdata)
rm(consortdata)

#~~ Base data for each sheep

names(basedata) <- toupper(names(basedata))
basedata$Censused <- basedata$ID %in% censusdata$ID

basedata$BIRTHDATE <- apply(basedata[,c("BIRTHYEAR", "BIRTHMONTH", "BIRTHDAY")], 1, function(x) paste(x, collapse = "-"))
basedata$BIRTHDATE[grep("NA", basedata$BIRTHDATE)] <- NA
basedata$BIRTHDATE <- as.Date(basedata$BIRTHDATE, "%Y-%m-%d")

basedata$DEATHMONTH[which(basedata$DEATHMONTH == 0)] <- NA
basedata$DEATHDAY[which(basedata$DEATHDAY == 0)] <- NA

basedata$DEATHDATE <- apply(basedata[,c("DEATHYEAR", "DEATHMONTH", "DEATHDAY")], 1, function(x) paste(x, collapse = "-"))
basedata$DEATHDATE[grep("NA", basedata$DEATHDATE)] <- NA
basedata$DEATHDATE <- as.Date(basedata$DEATHDATE, "%Y-%m-%d")

basedata$SEX[which(basedata$SEX == 1)] <- "F"
basedata$SEX[which(basedata$SEX == 2)] <- "M"
basedata$SEX[which(basedata$SEX == 3)] <- "Cas"


table(basedata$SIBCOUNT)
basedata$TWIN <- ifelse(basedata$SIBCOUNT > 0, 1, 0)

basedata <- unique(subset(basedata, select = -c(SIBCOUNT, SIBNOOBSERVED, TWINID, TRIPLETID)))

#~~ Add pedigree information to basedata

basedata$ID <- as.character(basedata$ID)

basedata <- join(basedata, pedigree)

#~~ Old Pink Ages

for(i in which(basedata$ID %in% oldpink$ID)){
  if(is.na(basedata$BIRTHYEAR[i])){
    
    x <- oldpink[oldpink$ID == basedata$ID[i],]$ProbCohort
    x <- unique(x)
    
    if(any(is.na(x))) x <- x[-which(is.na(x))]
    if(length(x) > 1) basedata$BIRTHYEAR[i]  <- min(x)
    if(length(x) == 1) basedata$BIRTHYEAR[i] <- x
    if(length(x) == 0) basedata$BIRTHYEAR[i] <- NA
  }
}

rm(x, i)

basedata$OldPink <- basedata$ID %in% oldpink$ID

#~~ Capture Phenotype Information

head(capdata)
capdata$CAPDATE <- apply(capdata[,c("CapYear", "CapMonth", "CapDay")], 1, function(x) paste(x, collapse = "-"))
capdata$CAPDATE[grep("NA", capdata$CAPDATE)] <- NA
capdata$CAPDATE <- as.Date(capdata$CAPDATE, "%Y-%m-%d")
capdata$Sex[which(capdata$Sex == 1)] <- "F"
capdata$Sex[which(capdata$Sex == 2)] <- "M"
capdata$Sex[which(capdata$Sex == 3)] <- "Cas"
#capdata$Comments <- NULL


#~~ Foetus information

basedata$Foetus <- basedata$ID %in% foetusIDs

rm(oldpink, foetusIDs)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Add Census Information to the Base Data File #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Add live capture data to the census information

sub.capdata <- subset(capdata, select = c(ID, LiveMeasure, CAPDATE))
names(sub.capdata) <- c("ID", "PLACE", "DATE")
sub.capdata <- subset(sub.capdata, PLACE == "Live" & !is.na(DATE))

censusdata <- subset(censusdata, !is.na(DATE))
head(censusdata)

censusdata <- rbind(censusdata, sub.capdata)
censusdata <- unique(censusdata)
rm(sub.capdata)

#~~ If the individual sired an offspring (even a foetus), add a date to the censusdata at 01/11/(Birthyear - 1)

transped <- melt(pedigree, id.vars = "ID")
names(transped) <- c("ID", "PARENT", "PARENT.ID")
transped <- join(transped, basedata)

transped <- subset(transped, !is.na(PARENT.ID) & !is.na(BIRTHYEAR))
transped <- subset(transped, select = c(PARENT.ID, BIRTHYEAR))

names(transped) <- c("ID", "YEAR")
transped$DATE <- as.Date(paste(transped$YEAR - 1, "11", "01", sep = "-"), "%Y-%m-%d")
transped$PLACE <- "RUT"
transped <- subset(transped, select = -YEAR)

censusdata <- rbind(censusdata, transped)
censusdata <- unique(censusdata)
rm(transped)

#~~ Add the date of birth of IDs as an observation of it's mother.

transped <- melt(pedigree, id.vars = "ID")
names(transped) <- c("ID", "PARENT", "PARENT.ID")
transped <- join(transped, basedata)

transped <- subset(transped, !is.na(PARENT.ID) & !is.na(BIRTHYEAR) & PARENT == "MOTHER")

transped$BIRTHDATE2 <- paste(transped$BIRTHYEAR, transped$BIRTHMONTH, transped$BIRTHDAY, sep = "-")

tempvec <- which(is.na(transped$BIRTHDATE) & !is.na(transped$BIRTHYEAR) & !is.na(transped$BIRTHMONTH))
transped$BIRTHDATE2[tempvec] <- paste(transped$BIRTHYEAR[tempvec], transped$BIRTHMONTH[tempvec], "01", sep = "-")
rm(tempvec)

tempvec <- which(is.na(transped$BIRTHDATE) & !is.na(transped$BIRTHYEAR) & is.na(transped$BIRTHMONTH))
transped$BIRTHDATE2[tempvec] <- paste(transped$BIRTHYEAR[tempvec], "04", "01", sep = "-")
rm(tempvec)
transped$BIRTHDATE2[grep("NA", transped$BIRTHDATE2)]

transped <- subset(transped, select = c(PARENT.ID, BIRTHDATE2))

names(transped) <- c("ID", "DATE")
transped$DATE <- as.Date(transped$DATE, "%Y-%m-%d")
transped$PLACE <- "Mum.BIRTH"

censusdata <- rbind(censusdata, transped)
censusdata <- unique(censusdata)
rm(transped)

#~~ Finally add the birth date to the census data

sub.basedata <- subset(basedata, select = c(ID, BIRTHDATE))
names(sub.basedata) <- c("ID", "DATE")
sub.basedata$PLACE <- "BIRTH"
sub.basedata <- subset(sub.basedata, !is.na(DATE))

censusdata <- rbind(censusdata, sub.basedata)
rm(sub.basedata)

#~~ Get rid of information from wrongly observed IDs

censusdata <- join(censusdata, subset(basedata, select = c(ID, DEATHDATE)))
censusdata <- censusdata[-which(censusdata$DEATHDATE < censusdata$DATE & !censusdata$PLACE %in% c("Live", "RUT", "BIRTH")),]

censusdata <- join(censusdata, subset(basedata, select = c(ID, DEATHMONTH, DEATHYEAR)))
censusdata$DEATHMONTH <- ifelse(!is.na(censusdata$DEATHYEAR) & is.na(censusdata$DEATHMONTH), 12, censusdata$DEATHMONTH)
censusdata$DEATHDAY <- ifelse(!is.na(censusdata$DEATHYEAR) & censusdata$DEATHMONTH %in% c(1, 3, 5, 7, 8, 10, 12), 31, NA)
censusdata$DEATHDAY <- ifelse(!is.na(censusdata$DEATHYEAR) & censusdata$DEATHMONTH %in% c(9, 4, 6, 11), 30, censusdata$DEATHDAY)
censusdata$DEATHDAY <- ifelse(!is.na(censusdata$DEATHYEAR) & censusdata$DEATHMONTH %in% c(2), 28, censusdata$DEATHDAY)
censusdata$TempDeathDate <- as.Date(paste(censusdata$DEATHYEAR, censusdata$DEATHMONTH, censusdata$DEATHDAY, sep = "-"), "%Y-%m-%d")

censusdata <- censusdata[-which(censusdata$TempDeathDate < censusdata$DATE & !censusdata$PLACE %in% c("Live", "RUT", "BIRTH")),]

censusdata <- censusdata[,1:3]

head(censusdata)

#~~ Determine the dates first seen and last seen from capture and census info

maxdate <- censusdata %>% group_by(ID) %>% summarise(max(DATE))
names(maxdate) <- c("ID", "LastSeen")

mindate <- censusdata %>% group_by(ID) %>% summarise(min(DATE))
names(mindate) <- c("ID", "FirstSeen")

basedata <- join(basedata, maxdate)
basedata <- join(basedata, mindate)

rm(maxdate, mindate)

#~~ Find the maximum death date

maxdeaddate <- subset(capdata, LiveMeasure == "Dead") %>% group_by(ID) %>% summarise(max(CAPDATE, na.rm = T))
names(maxdeaddate) <- c("ID", "MaxDeadDate")
maxdeaddate$ID <- as.character(maxdeaddate$ID)

basedata <- join(basedata, maxdeaddate)
rm(maxdeaddate)

#~~ Edit maximum death date to be 01/Month/Year

tempvec <- which(!is.na(basedata$DEATHMONTH) & !is.na(basedata$DEATHYEAR))

basedata$MaxDeadDate <- paste0("01-", month(basedata$MaxDeadDate), "-", year(basedata$MaxDeadDate))
basedata$MaxDeadDate[tempvec] <- paste0("01-", basedata$DEATHMONTH[tempvec], "-", basedata$DEATHYEAR[tempvec])

basedata$MaxDeadDate[grep("NA", basedata$MaxDeadDate)] <- NA
basedata$MaxDeadDate <- as.Date(basedata$MaxDeadDate, "%d-%m-%Y")

rm(tempvec)

#~~ Find the last date the sheep was alive

basedata$LastAliveDate <- as.Date(NA)
basedata$LastAliveDate[which(!is.na(basedata$DEATHDATE))] <-  basedata$DEATHDATE[which(!is.na(basedata$DEATHDATE))]
basedata$LastAliveDate[which(is.na(basedata$DEATHDATE) & !is.na(basedata$LastSeen))] <-  basedata$LastSeen[which(is.na(basedata$DEATHDATE) & !is.na(basedata$LastSeen))]

basedata$DeathAgeDays <- basedata$DEATHDATE - basedata$BIRTHDATE
basedata$DeathAgeDays <- ifelse(is.na(basedata$DEATHDATE) & !is.na(basedata$MaxDeadDate), basedata$MaxDeadDate - basedata$BIRTHDATE, basedata$DeathAgeDays)
basedata$DeathAgeDays[which(basedata$DeathAgeDays < 0)] <- 0

basedata$LastSeenAgeDays <- basedata$LastSeen - basedata$BIRTHDATE
basedata$LastSeenAgeDays[which(!is.na(basedata$DeathAgeDays))] <- basedata$DeathAgeDays[which(!is.na(basedata$DeathAgeDays))] 

basedata$LastSeenFirstSeenDiff <- basedata$LastSeen - basedata$FirstSeen

#~~ Identify sheep that were last seen at birth

basedata$LastSeenAtBirth <- ifelse(basedata$LastSeen - basedata$BIRTHDATE <= 3, "yes", NA)

#~~ Get location centroids

locationdata <- na.omit(locationdata)
locationdata$YEAR <- ifelse(month(locationdata$DATE) < 5, year(locationdata$DATE) - 1, year(locationdata$DATE))
locationdata <- subset(locationdata, NORTHING > 10)

head(locationdata)
studycentroids <- unique(subset(locationdata, select = c(EASTING, NORTHING)))

ggplot(studycentroids, aes(EASTING, NORTHING)) + 
  geom_point() + 
  geom_point(x = mean(locationdata$EASTING), y = mean(locationdata$NORTHING), colour = "red")

basecentroids = data.frame(MeanLifeEASTING = tapply(locationdata$EASTING, locationdata$ID, mean),
                           MeanLifeNORTHING = tapply(locationdata$NORTHING, locationdata$ID, mean),
                           CountLifeEASTING = tapply(locationdata$EASTING, locationdata$ID, length),
                           CountLifeNORTHING = tapply(locationdata$NORTHING, locationdata$ID, length))

basecentroids$ID <- row.names(basecentroids)

basecentroids <- tbl_df(basecentroids)

basecentroids$LifeCartesianDistance <- sqrt(
  (basecentroids$MeanLifeEASTING - mean(locationdata$EASTING))^2 +
    (basecentroids$MeanLifeNORTHING - mean(locationdata$NORTHING))^2
  )

ggplot(basecentroids, aes(LifeCartesianDistance)) + geom_histogram(binwidth = 0.2, col = "grey")

basecentroidsann = data.frame(MeanEASTING = tapply(locationdata$EASTING,   list(locationdata$ID, locationdata$YEAR), mean),
                              MeanNORTHING = tapply(locationdata$NORTHING, list(locationdata$ID, locationdata$YEAR), mean),
                              CountEASTING = tapply(locationdata$EASTING,   list(locationdata$ID, locationdata$YEAR), length),
                              CountNORTHING = tapply(locationdata$NORTHING, list(locationdata$ID, locationdata$YEAR), length))
head(basecentroidsann)
basecentroidsann$ID <- row.names(basecentroidsann)

basecentroidsann <- melt(basecentroidsann)
basecentroidsann <- na.omit(basecentroidsann)
basecentroidsann <- separate(basecentroidsann, col = "variable", into = c("variable", "SheepYear"), sep = "\\.", remove = T)
basecentroidsann <- cast(basecentroidsann, formula = ID + SheepYear ~ variable)

basecentroidsann$CartesianDistance <- sqrt(
  (basecentroidsann$MeanEASTING - mean(locationdata$EASTING))^2 +
    (basecentroidsann$MeanNORTHING - mean(locationdata$NORTHING))^2
)

head(basedata)

head(basecentroidsann)
temp <- basecentroidsann
temp$SheepYear <- as.numeric(temp$SheepYear) - 1
names(temp) <- c("ID", "SheepYear", "CountEASTING.p1", "CountNORTHING.p1", "MeanEASTING.p1", 
                 "MeanNORTHING.p1", "CartesianDistance.p1")
temp$SheepYear <- as.character(temp$SheepYear)

basecentroidsann <- join(basecentroidsann, temp)

basedata <- join(basedata, basecentroids)

locationdata2 <- locationdata
locationdata2 <- subset(locationdata2, month(locationdata2$DATE) %in% 3:5)

basecentroidsspring = data.frame(MeanSpringEASTING = tapply(locationdata2$EASTING,   list(locationdata2$ID, locationdata2$YEAR), mean),
                                 MeanSpringNORTHING = tapply(locationdata2$NORTHING, list(locationdata2$ID, locationdata2$YEAR), mean),
                                 CountSpringEASTING = tapply(locationdata2$EASTING,   list(locationdata2$ID, locationdata2$YEAR), length),
                                 CountSpringNORTHING = tapply(locationdata2$NORTHING, list(locationdata2$ID, locationdata2$YEAR), length))
head(basecentroidsspring)
basecentroidsspring$ID <- row.names(basecentroidsspring)

basecentroidsspring <- melt(basecentroidsspring)
basecentroidsspring <- na.omit(basecentroidsspring)
basecentroidsspring <- separate(basecentroidsspring, col = "variable", into = c("variable", "SheepYear"), sep = "\\.", remove = T)
basecentroidsspring <- cast(basecentroidsspring, formula = ID + SheepYear ~ variable)

basecentroidsspring$SpringCartesianDistance <- sqrt(
  (basecentroidsspring$MeanSpringEASTING - mean(locationdata2$EASTING))^2 +
    (basecentroidsspring$MeanSpringNORTHING - mean(locationdata2$NORTHING))^2
)

head(basecentroidsspring)

basecentroidsspring <- tbl_df(basecentroidsspring)


locationdata2 <- locationdata
locationdata2 <- subset(locationdata2, month(locationdata2$DATE) %in% 10:12)

basecentroidsrut = data.frame(MeanRutEASTING = tapply(locationdata2$EASTING,   list(locationdata2$ID, locationdata2$YEAR), mean),
                                 MeanRutNORTHING = tapply(locationdata2$NORTHING, list(locationdata2$ID, locationdata2$YEAR), mean),
                                 CountRutEASTING = tapply(locationdata2$EASTING,   list(locationdata2$ID, locationdata2$YEAR), length),
                                 CountRutNORTHING = tapply(locationdata2$NORTHING, list(locationdata2$ID, locationdata2$YEAR), length))
head(basecentroidsrut)
basecentroidsrut$ID <- row.names(basecentroidsrut)

basecentroidsrut <- melt(basecentroidsrut)
basecentroidsrut <- na.omit(basecentroidsrut)
basecentroidsrut <- separate(basecentroidsrut, col = "variable", into = c("variable", "SheepYear"), sep = "\\.", remove = T)
basecentroidsrut <- cast(basecentroidsrut, formula = ID + SheepYear ~ variable)

basecentroidsspring$SpringCartesianDistance <- sqrt(
  (basecentroidsspring$MeanSpringEASTING - mean(locationdata2$EASTING))^2 +
    (basecentroidsspring$MeanSpringNORTHING - mean(locationdata2$NORTHING))^2
)

head(basecentroidsrut)

basecentroidsrut <- tbl_df(basecentroidsrut)

#~~ need to change the year as it will match the previous year...
basecentroidsrut$SheepYear <- as.numeric(basecentroidsrut$SheepYear) + 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Calculate Early survival at different resolutions    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(basedata)
max(censusdata$DATE)

# ~~ Did they survive past April the Year after they were born? NB. If death
# year the year after the birth year, but with no month/day and was last seen in
# birth year, then say it did not survive the winter.

basedata$OverWinterSurvival <- NA

# TRUE

# seen or died age 2 upwards
basedata$OverWinterSurvival[which((year(basedata$LastSeen) - basedata$BIRTHYEAR) >= 2)] <- T 
basedata$OverWinterSurvival[which((basedata$DEATHYEAR - basedata$BIRTHYEAR) >= 2)]      <- T 
basedata$OverWinterSurvival[which((year(basedata$MaxDeadDate) - year(basedata$FirstSeen)) >= 2)] <- T
basedata$OverWinterSurvival[which((year(basedata$LastSeen) - year(basedata$FirstSeen)) >= 2)] <- T


# seen or died age 1 after March
basedata$OverWinterSurvival[which((year(basedata$LastSeen) - basedata$BIRTHYEAR) == 1 & # which sheep were seen on or after April the year after birth
                                    month(basedata$LastSeen) > 4)] <- T
basedata$OverWinterSurvival[which((basedata$DEATHYEAR - basedata$BIRTHYEAR) == 1 & basedata$DEATHMONTH > 4) ] <- T # which sheep died on or after April the year after birth
basedata$OverWinterSurvival[which((year(basedata$MaxDeadDate) - basedata$BIRTHYEAR) == 1 & month(basedata$MaxDeadDate) > 4)] <- T
basedata$OverWinterSurvival[which((year(basedata$LastSeen) - year(basedata$FirstSeen)) == 1 & month(basedata$LastSeen) > 4)] <- T


# FALSE (will overwrite some MaxDeadDate IDs from January)

# which sheep died the year they were born?
basedata$OverWinterSurvival[which(basedata$DEATHYEAR == basedata$BIRTHYEAR)] <- F 
basedata$OverWinterSurvival[which(year(basedata$MaxDeadDate) == basedata$BIRTHYEAR)] <- F # which sheep died the year they were born?

basedata$OverWinterSurvival[which((basedata$DEATHYEAR - basedata$BIRTHYEAR) == 1 & basedata$DEATHMONTH <= 4) ] <- F # which sheep died within a year of birth (before April?)
basedata$OverWinterSurvival[which((basedata$MaxDeadDate - basedata$BIRTHYEAR) == 1 & month(basedata$MaxDeadDate) <= 4) ] <- F # which sheep died within a year of birth (before April?)
basedata$OverWinterSurvival[which((basedata$DEATHYEAR - basedata$BIRTHYEAR) == 1 & month(basedata$MaxDeadDate) <= 4) ] <- F # which sheep died within a year of birth (before April?)


table(basedata$OverWinterSurvival, useNA = "always")


#~~ Did they survive to the October in the year they were born?

basedata$OctSurvival <- NA

# TRUE
basedata$OctSurvival[which(basedata$OverWinterSurvival == T)] <- T
basedata$OctSurvival[which((year(basedata$LastSeen) - basedata$BIRTHYEAR) > 0)] <- T
basedata$OctSurvival[which((basedata$DEATHYEAR - basedata$BIRTHYEAR) > 0)]      <- T
basedata$OctSurvival[which(month(basedata$LastSeen) > 9)] <- T
basedata$OctSurvival[which((basedata$DEATHYEAR - basedata$BIRTHYEAR) == 0 & month(basedata$MaxDeadDate) > 9) ] <- T
basedata$OctSurvival[which(is.na(basedata$BIRTHYEAR) & month(basedata$MaxDeadDate) %in% c(1:3, 10:12))] <- T
basedata$OctSurvival[which(is.na(basedata$BIRTHYEAR) & month(basedata$LastSeen) %in% c(1:3, 10:12))] <- T

# FALSE
basedata$OctSurvival[which(basedata$DEATHYEAR == basedata$BIRTHYEAR & month(basedata$MaxDeadDate) < 10)] <- F

table(basedata$OctSurvival, useNA = "always")



#~~ Did they survive to the August in the year they were born?

basedata$AugSurvival <- NA

# TRUE
basedata$AugSurvival[which(basedata$OverWinterSurvival == T)] <- T
basedata$AugSurvival[which((year(basedata$LastSeen) - basedata$BIRTHYEAR) > 0)] <- T
basedata$AugSurvival[which((basedata$DEATHYEAR - basedata$BIRTHYEAR) > 0)]      <- T
basedata$AugSurvival[which(month(basedata$LastSeen) > 7)] <- T
basedata$AugSurvival[which((basedata$DEATHYEAR - basedata$BIRTHYEAR) == 0 & month(basedata$MaxDeadDate) > 7) ] <- T
basedata$AugSurvival[which(is.na(basedata$BIRTHYEAR) & month(basedata$MaxDeadDate) %in% c(1:3, 8:12))] <- T
basedata$AugSurvival[which(is.na(basedata$BIRTHYEAR) & month(basedata$LastSeen) %in% c(1:3, 8:12))] <- T

# FALSE
basedata$AugSurvival[which(basedata$DEATHYEAR == basedata$BIRTHYEAR & month(basedata$MaxDeadDate) < 8)] <- F

table(basedata$AugSurvival, useNA = "always")

#~~ Did they survive the first three days?

# TRUE
basedata$BirthSurvival <- NA
basedata$BirthSurvival[which(basedata$AugSurvival == T)] <- T
basedata$BirthSurvival[which(basedata$LastSeenAgeDays > 3)] <- T
basedata$BirthSurvival[which(month(basedata$LastSeen) %in% c(1:2, 6:12))] <- T
basedata$BirthSurvival[which(is.na(basedata$BirthSurvival) & basedata$DEATHMONTH > basedata$BIRTHMONTH)] <- T

# FALSE
basedata$BirthSurvival[which(basedata$DeathAgeDays < 4)] <- F
basedata$BirthSurvival[which(basedata$MaxDeadDate - basedata$BIRTHDATE < 4)] <- F


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Calculate Reproductive Success                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

transped <- gather(pedigree, PARENT, PARENT.ID, MOTHER, FATHER)
transped <- join(transped, basedata)

#~~ remove NA parents and Dead Foetuses

transped <- subset(transped, !is.na(PARENT.ID))
transped <- subset(transped, Foetus == FALSE)

transped <- subset(transped, select = c(ID, PARENT.ID, SEX, BIRTHYEAR, OverWinterSurvival, OctSurvival, AugSurvival, BirthSurvival))
head(transped)

#~~ Annual reproductive success

annRS.0 <- data.frame(table(transped$PARENT.ID, transped$BIRTHYEAR))
head(annRS.0)
names(annRS.0) <- c("ID", "BIRTHYEAR", "OffspringBorn")

transped.1 <- subset(transped, BirthSurvival == T)
annRS.1 <- data.frame(table(transped.1$PARENT.ID, transped.1$BIRTHYEAR))
head(annRS.1)
names(annRS.1) <- c("ID", "BIRTHYEAR", "OffspringSurvived")

transped.2 <- subset(transped, AugSurvival == T)
annRS.2 <- data.frame(table(transped.2$PARENT.ID, transped.2$BIRTHYEAR))
head(annRS.2)
names(annRS.2) <- c("ID", "BIRTHYEAR", "AugSurvOffspring")

transped.3 <- subset(transped, OctSurvival == T)
annRS.3 <- data.frame(table(transped.3$PARENT.ID, transped.3$BIRTHYEAR))
head(annRS.3)
names(annRS.3) <- c("ID", "BIRTHYEAR", "OctSurvOffspring")

transped.4 <- subset(transped, OverWinterSurvival == T)
annRS.4 <- data.frame(table(transped.4$PARENT.ID, transped.4$BIRTHYEAR))
head(annRS.4)
names(annRS.4) <- c("ID", "BIRTHYEAR", "OverWinterOffspring")

#~~ Join the tables together

annRS <- join(annRS.0, annRS.1)
annRS <- join(annRS,   annRS.2)
annRS <- join(annRS,   annRS.3)
annRS <- join(annRS,   annRS.4)

rm(annRS.0, annRS.1, annRS.2, annRS.3, annRS.4, transped, transped.1, transped.2, transped.3, transped.4, locationdata, locationdata2, temp)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Calculate individual survival to next April                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cutoffmonth <- 5

basedata$maxyear <- year(basedata$LastSeen)
basedata$maxmonth<- month(basedata$LastSeen)
basedata$minyear <- year(basedata$FirstSeen)

#~~ Create x. IDs will be removed from x when it's annual survival has been determined.

x <- basedata
head(x)
x <- subset(x, select = -c(SEX, COAT, HORN, BIRTHWT, OldPink))
x <- subset(x, Foetus == FALSE)

#~~ Final annual survival values will be stored in newdata.

newdata <- NULL

#~~ 1. Which sheep died in their first year (where DEATHYEAR == BIRTHYEAR)

sampvec <- which(!is.na(x$BIRTHYEAR) & 
                   !is.na(x$DEATHYEAR) & 
                   x$DEATHYEAR == x$BIRTHYEAR)

y1 <- x[sampvec, c("ID", "BIRTHYEAR")]
y1$Age      = 0
y1$Survival = 0
names(y1)[which(names(y1) == "BIRTHYEAR")] <- "SheepYear"
if(length(sampvec) > 0) x <- x[-sampvec,]
newdata <- rbind(newdata, y1)
rm(sampvec, y1)

#~~ 2. Which sheep died before May of any year?     

sampvec <- which(x$DEATHYEAR > x$BIRTHYEAR & month(x$MaxDeadDate) < cutoffmonth)

for(i in sampvec){
  
  if(which(sampvec == i) %in% seq(1, length(sampvec), 100)) print(paste("Analysing row", which(sampvec == i), "of", length(sampvec)))
  
  z <- x[i,]
  maxage <- (z$DEATHYEAR - z$BIRTHYEAR)-1
  y1 <- data.frame(ID = z$ID,
                   SheepYear = seq(z$BIRTHYEAR, (z$DEATHYEAR - 1)),
                   Age = seq(0, maxage, 1),
                   Survival = c(rep(1, times = (maxage )), 0))
  
  newdata <- rbind(newdata, y1)
  rm(y1, maxage, z)
}

x <- x[-sampvec,]
rm(sampvec)

#~~ 3. Which sheep died after April       

sampvec <- which(x$DEATHYEAR > x$BIRTHYEAR & x$DEATHMONTH >= cutoffmonth)

for(i in sampvec){
  
  z <- x[i,]
  maxage <- (z$DEATHYEAR - z$BIRTHYEAR)-1
  y1 <- data.frame(ID = z$ID,
                   SheepYear = seq(z$BIRTHYEAR, z$DEATHYEAR),
                   Age = seq(0, maxage + 1, 1),
                   Survival = c(rep(1, times = (maxage + 1)), 0))
  newdata <- rbind(newdata, y1)
  rm(y1, maxage, z)
}

if(length(sampvec) > 0) x <- x[-sampvec,]
rm(sampvec)

#~~ 4. Which sheep have birthyear, deathyear but no death month, but were seen after April

sampvec <- which(!is.na(x$BIRTHYEAR) &
                   !is.na(x$DEATHYEAR) &
                   is.na(x$DEATHMONTH) &
                   year(x$LastSeen) == x$DEATHYEAR &
                   month(x$LastSeen) >= cutoffmonth)

for(i in sampvec){
  z <- x[i,]
  maxage <- (z$DEATHYEAR - z$BIRTHYEAR)-1
  y1 <- data.frame(ID = z$ID,
                   SheepYear = seq(z$BIRTHYEAR, z$DEATHYEAR),
                   Age = seq(0, maxage + 1, 1),
                   Survival = c(rep(1, times = (maxage + 1)), 0))
  newdata <- rbind(newdata, y1)
  rm(y1, maxage, z)
}

if(length(sampvec) > 0) x <- x[-sampvec,]
rm(sampvec)

#~~ 5. Which sheep died later but the death month is unknown. If last seen the year before, say 0.

newdata$Comment <- NA

sampvec <- which(x$DEATHYEAR > x$BIRTHYEAR & is.na(x$DEATHMONTH) & year(x$LastSeen) < x$DEATHYEAR)

for(i in sampvec){
  
  if(which(sampvec == i) %in% seq(1, length(sampvec), 100)) print(paste("Analysing row", which(sampvec == i), "of", length(sampvec)))
  z <- x[i,]
  
  maxage <- (z$DEATHYEAR - z$BIRTHYEAR)-1
  y1 <- data.frame(ID = z$ID,
                   SheepYear = seq(z$BIRTHYEAR, z$DEATHYEAR - 1),
                   Age = seq(0, maxage, 1),
                   Survival = c(rep(1, times = (maxage )), 0),
                   Comment = c(rep(NA, times = maxage), "Not seen since year before, Likely to be 0"))
  
  newdata <- rbind(newdata, y1)
  rm(y1, maxage, z)
}

if(length(sampvec) > 0) x <- x[-sampvec,]
rm(sampvec)


sampvec <- which(x$DEATHYEAR > x$BIRTHYEAR & is.na(x$DEATHMONTH))

for(i in sampvec){
  
  if(which(sampvec == i) %in% seq(1, length(sampvec), 100)) print(paste("Analysing row", which(sampvec == i), "of", length(sampvec)))
  z <- x[i,]
  
  maxage <- (z$DEATHYEAR - z$BIRTHYEAR)-1
  y1 <- data.frame(ID = z$ID,
                   SheepYear = seq(z$BIRTHYEAR, z$DEATHYEAR - 1),
                   Age = seq(0, maxage, 1),
                   Survival = c(rep(1, times = (maxage )), NA),
                   Comment = c(rep(NA, times = maxage), "Seen Spring, Maybe 0"))
  
  newdata <- rbind(newdata, y1)
  rm(y1, maxage, z)
}

if(length(sampvec) > 0) x <- x[-sampvec,]
rm(sampvec)



#~~ 6. Which sheep have a known Birth year but no know death year, but have been censused after April of a year?

sampvec <- which(!is.na(x$BIRTHYEAR) & is.na(x$DEATHYEAR) & !is.na(x$LastSeen) & month(x$LastSeen) >= cutoffmonth)

for(i in sampvec){
  
  if(which(sampvec == i) %in% seq(1, length(sampvec), 100)) print(paste("Analysing row", which(sampvec == i), "of", length(sampvec)))
  z <- x[i,]
  maxage <- year(z$LastSeen) - z$BIRTHYEAR
  y1 <- data.frame(ID = z$ID,
                   SheepYear = seq(z$BIRTHYEAR, year(z$LastSeen)),
                   Age = seq(0, maxage, 1),
                   Survival = c(rep(1, times = (maxage)), NA),
                   Comment = c(rep(NA, times = maxage), paste("Last seen", month(z$LastSeen, label = T), year(z$LastSeen))))
  newdata <- rbind(newdata, y1)
  rm(y1, maxage, z)
}

if(length(sampvec) > 0) x <- x[-sampvec,]
rm(sampvec)


#~~ 7. Which sheep have a known Birth year, no known death year and were last seen in the birth year or never seen again

sampvec <- sort(unique(c(which(!is.na(x$BIRTHYEAR) & is.na(x$DEATHYEAR) & !is.na(x$LastSeen) & year(x$LastSeen) == x$BIRTHYEAR),
                         which(!is.na(x$BIRTHYEAR) & is.na(x$DEATHYEAR) & is.na(x$LastSeen)))))
for(i in sampvec){
  
  if(which(sampvec == i) %in% seq(1, length(sampvec), 100)) print(paste("Analysing row", which(sampvec == i), "of", length(sampvec)))
  z <- x[i,]
  y1 <- data.frame(ID = z$ID,
                   SheepYear = z$BIRTHYEAR,
                   Age = 0,
                   Survival = NA,
                   Comment = "Last seen at birth")
  newdata <- rbind(newdata, y1)
  rm(y1, z)
}

if(length(sampvec) > 0) x <- x[-sampvec,]
rm(sampvec)


#~~ 8. Which sheep have a known Birth year but no know death year, but have been censused before April of a year?

sampvec <- which(!is.na(x$BIRTHYEAR) & is.na(x$DEATHYEAR) & !is.na(x$LastSeen) & month(x$LastSeen) < cutoffmonth)

for(i in sampvec){
  
  if(which(sampvec == i) %in% seq(1, length(sampvec), 100)) print(paste("Analysing row", which(sampvec == i), "of", length(sampvec)))
  z <- x[i,]
  maxage <- (year(z$LastSeen) - z$BIRTHYEAR) -1
  y1 <- data.frame(ID = z$ID,
                   SheepYear = seq(z$BIRTHYEAR, year(z$LastSeen) -1),
                   Age = seq(0, maxage, 1),
                   Survival = c(rep(1, times = (maxage)), NA),
                   Comment = c(rep(NA, times = maxage), paste("Last seen", month(z$LastSeen, label = T), year(z$LastSeen))))
  newdata <- rbind(newdata, y1)
  rm(y1, maxage, z)
}

if(length(sampvec) > 0) x <- x[-sampvec,]
rm(sampvec)

#~~ 9. Remainder are not added to the table.  Save objects in current form:

# save(newdata, x, annRS, file = "results/1.2_Fitness_b4_processing_v2.Rdata")

#~~ Create a lifetime survival table

annSurv <- newdata
rm(newdata)

lifeSurv<- data.frame(MinSheepYear = tapply(annSurv$SheepYear, annSurv$ID, min),
                      MaxSheepYear = tapply(annSurv$SheepYear, annSurv$ID, max),
                      Longevity    = tapply(annSurv$Survival,  annSurv$ID, sum))
lifeSurv$ID <- row.names(lifeSurv)

hist(lifeSurv$Longevity)

# Survival check

deathage <- subset(basedata, select = c(ID, BIRTHYEAR, DEATHYEAR, DEATHMONTH))
deathage <- left_join(deathage, lifeSurv)
deathage$DEATHAGE <- ifelse(deathage$DEATHMONTH <= 10,
                            deathage$DEATHYEAR - deathage$BIRTHYEAR,
                            (deathage$DEATHYEAR - deathage$BIRTHYEAR) + 1)

head(deathage)
ggplot(deathage, aes(DEATHAGE, Longevity)) + geom_jitter(alpha = 0.1)

deathage[which(deathage$DEATHAGE < deathage$Longevity),]

gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. Consolidate survival and reproductive success information       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(annSurv)
head(annRS)

annRS$SheepYear <- as.numeric(as.character(annRS$BIRTHYEAR))-1   # THIS WAS CHANGED BY ADDING -1
annRS <- subset(annRS, select = -BIRTHYEAR)

annSurv$ID.SheepYear <- paste(annSurv$ID, annSurv$SheepYear)
annRS$ID.SheepYear   <- paste(annRS$ID  , annRS$SheepYear)

annRS$SurvivalMeasureInYear <- annRS$ID.SheepYear %in% annSurv$ID.SheepYear

#~~ Create annfit

annfit <- join(annRS, annSurv, type = "full")
annfit <- join(annfit, basedata[,c("ID", "BIRTHYEAR", "SEX")])

#~~ Remove IDs with unknown Birth Years as fitness measures cannot be accurately determined.

annfit <- subset(annfit, !is.na(BIRTHYEAR))

#~~ Remove lines where Age and Survival are NA

annfit <- subset(annfit, !is.na(Age))

#annfit <- subset(annfit, !is.na(Age) & !is.na(Survival))
#annfit <- subset(annfit, SurvivalMeasureInYear == TRUE)
head(annfit)

#~~ Add a column: SeenInRut the year before

novcensus <- subset(censusdata, month(DATE) %in% c(10, 11, 12))
novcensus$SheepYear <- year(novcensus$DATE)          # THIS WAS CHANGED BY REMOVING +1
novcensus <- unique(novcensus[,c("ID", "SheepYear")])
novcensus$SeenInRut <- "yes"

annfit <- join(annfit, novcensus)
annfit$SeenInRut[which(is.na(annfit$SeenInRut))] <- "no"

table(annfit$SeenInRut, annfit$SheepYear)


#~~ If an ID was seen in the rut but has NA for number of offspring, then change the values to 0

annfit[which(is.na(annfit$OffspringBorn) & annfit$SeenInRut == "yes"),
       c("OffspringBorn",
         "OffspringSurvived",
         "AugSurvOffspring",
         "OctSurvOffspring",
         "OverWinterOffspring")] <- 0

#~~ If an ID was not seen in the rut but has 0 for number of offspring, then change the values to NA

annfit[which(annfit$OffspringBorn == 0 & annfit$SeenInRut == "no" & annfit$SEX == "M"),
       c("OffspringBorn",
         "OffspringSurvived",
         "AugSurvOffspring",
         "OctSurvOffspring",
         "OverWinterOffspring")] <- NA

head(annfit)

#~~ Add the census centroids

annfit$SheepYear <- as.character(annfit$SheepYear)

annfit <- join(annfit, basecentroidsann)
annfit <- join(annfit, basecentroids)
annfit <- join(annfit, basecentroidsspring)
annfit <- join(annfit, basecentroidsrut)


#~~ New redo the lifetime fitness

lifefit <- annfit[which(annfit$Survival == 0),c("ID", "Age")]

#temptab <- subset(annfit, Age != 0)

temptab2 <- data.frame(TotalOffspringBorn = tapply(annfit$OffspringBorn, annfit$ID, sum, na.rm = T),
                       TotalOffspringSurvived = tapply(annfit$OffspringSurvived, annfit$ID, sum, na.rm = T),
                       TotalOffspringAugust = tapply(annfit$AugSurvOffspring, annfit$ID, sum, na.rm = T),
                       TotalOffspringOctober = tapply(annfit$OctSurvOffspring, annfit$ID, sum, na.rm = T),
                       TotalOffspringOverWinter = tapply(annfit$OverWinterOffspring, annfit$ID, sum, na.rm = T))

temptab2$ID <- row.names(temptab2)

lifefit <- join(lifefit, temptab2, type = "full")

rm(temptab2, temptab, novcensus, transped, transped.hold, x)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 8. Sanity Check                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#oldannfit <- read_delim("../Recombination Analysis G/data/NewAnnualFitness20130413.txt", delim = "\t")
oldannfit <- read.table("../20180921_Soay_Recombination_Fitness_Analysis/results/1_Annual_Fitness_Measures_April_20180921.txt", sep = "\t", stringsAsFactors = F, header = T)

head(oldannfit)

oldannfit <- subset(oldannfit, select = c(ID, SheepYear, Survival, OffspringBorn))
names(oldannfit)[3:4] <- c("OldSurvival", "OldOffspringBorn")
oldannfit$ID <- as.character(oldannfit$ID)
oldannfit$SheepYear <- as.character(oldannfit$SheepYear)

test <- join(annfit, oldannfit)
head(test)

ggplot(test, aes(OffspringBorn, OldOffspringBorn)) + geom_point(alpha = 0.1) + stat_smooth(method = "lm")
ggplot(test, aes(Survival, OldSurvival)) + geom_point(alpha = 0.1) + stat_smooth(method = "lm")

rm(test, oldannfit)

oldlifefit <- read.table("../Recombination Analysis G/data/NewLifetimeFitness20130413.txt", sep = "\t", stringsAsFactors = F, header = T)
head(oldlifefit)

oldlifefit <- subset(oldlifefit, select = c(ID, Age, RecSeen, RecCount))
names(oldlifefit)[2:4] <- c("OldAge", "OldRecSeen", "OldRecCount")
oldlifefit$ID <- as.character(oldlifefit$ID)

test <- left_join(lifefit, oldlifefit)
head(test)

ggplot(test, aes(TotalOffspringOctober, OldRecSeen)) + geom_point(alpha = 0.1) + stat_smooth(method = "lm")
ggplot(test, aes(TotalOffspringOctober, OldRecCount)) + geom_point(alpha = 0.1) + stat_smooth(method = "lm")

ggplot(test, aes(Age, OldAge)) + geom_point(alpha = 0.1) + stat_smooth(method = "lm")

rm(test, oldlifefit)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 9. Merge with phenotypic data                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(annfit)
annfit <- subset(annfit, select = -c(ID.SheepYear, SurvivalMeasureInYear))
table(annfit$Age)

#~~ add pedigree info

head(pedigree)
annfit <- join(annfit, pedigree)

#~~ Add birth weight and SheepYear weight

head(basedata)
basedata2 <- subset(basedata, select = c(ID, SEX, TWIN, BIRTHWT))
annfit <- join(annfit, basedata2)

head(capdata)
capdata2 <- subset(capdata, select = c(ID, CapYear, CapMonth, LiveMeasure, Weight, Foreleg, Hindleg, HornLen, HornCirc, BolCirc, BolLen, Horn)) 

head(capdata2)
capdata2 <- subset(capdata2, LiveMeasure == "Live")
capdata2 <- subset(capdata2, select = -LiveMeasure)
capdata2$ID.CapYear <- paste(capdata2$ID, capdata2$CapYear, sep = "_")

#~~ add capdata for males that were not captured in August but were in rut and add an identifier for "Caught during Rut"

capdataNov <- subset(capdata2, CapMonth %in% c(11, 12))
capdata2 <- subset(capdata2, CapMonth == 8)

capdataNov <- subset(capdataNov, !ID.CapYear %in% capdata2$ID.CapYear)
head(capdataNov)

capdata2$RutMeasure <- NA
capdataNov$RutMeasure <- "Caught in Rut Only"

capdata2 <- rbind(capdata2, capdataNov)


#~~ Deal with repeated measures

capdata3 <- data.frame(table(capdata2$ID.CapYear))
head(capdata3)
names(capdata3)[1] <- "ID.CapYear"

capdata2 <- join(capdata2, capdata3)
head(capdata2)
capdata3 <- subset(capdata2, Freq > 1)
capdata2 <- subset(capdata2, Freq == 1)

for(i in unique(capdata3$ID.CapYear)){
  x <- capdata3[which(capdata3$ID.CapYear == i),]
  x1 <- sapply(x[,c("Weight", "Foreleg", "Hindleg", "HornLen", "HornCirc", "BolCirc", "BolLen")], mean, na.rm = T)
  x <- x[1,]
  x[,c("Weight", "Foreleg", "Hindleg", "HornLen", "HornCirc", "BolCirc", "BolLen")] <- x1
  capdata2 <- rbind(capdata2, x)
  rm(x, x1)
}
rm(capdata3, capdataNov)

head(annfit)
head(capdata2)
names(capdata2)[1:2] <- c("ID", "SheepYear")
capdata2$ID <- as.character(capdata2$ID)
capdata2$SheepYear <- as.character(capdata2$SheepYear)

annfit <- join(annfit, capdata2)

rm(basedata2, capdata2)

head(annfit)

#~~ save to file

write.table(annfit, "../sheep/data/1_Annual_Fitness_Measures_April_20190501.txt", sep = "\t", quote = F, row.names = F)


# #~~ LIfetime fitness
# 
# basedata2 <- subset(basedata, select = c(ID, MOTHER, FATHER, SEX, BIRTHYEAR, TWIN, BIRTHWT, DEATHYEAR))
# lifefit <- left_join(lifefit, basedata2)
# capdata4 <- unique(capdata[,c("ID", "Horn")])
# lifefit <- left_join(lifefit, capdata4)
# 
# table(lifefit$BIRTHYEAR)
# 
# 
# annfit0 <- subset(annfit, is.na(OffspringBorn))
# 
# lifefit$YearMissing <- ifelse(lifefit$ID %in% annfit0$ID, "YearsMissing", "AllYearsCounted")
# 
# rm(capdata4, basedata2, annfit0, deathage)
# 
# write.table(lifefit, "results/1_Lifetime_Fitness_Measures_April_20180831.txt", sep = "\t", quote = F, row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 10. Second Sanity Check By Eye                            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


set.seed(2019)

#~~ Sample 20 IDs

basedata <- subset(basedata, SEX != "Cas")

id.samp <- sample(basedata$ID, size = 20, replace = F)

annfit.samp  <- subset(annfit, ID %in% id.samp)
#lifefit.samp <- subset(lifefit, ID %in% id.samp)

write.table(annfit.samp, "results/20190501_annfit_sample.txt", row.names = F, sep = "\t", quote = F)
#write.table(lifefit.samp, "results/20180921_lifefit_sample.txt", row.names = F, sep = "\t", quote = F)
write.table(basedata, "results/20190501_basedata_updated.txt", row.names = F, sep = "\t", quote = F)





