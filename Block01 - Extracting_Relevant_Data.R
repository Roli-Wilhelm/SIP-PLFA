#Ensure all data is in the form of a .csv (note: some .csv formats differ; use the microsoft office formate) and in the working directory. All ".csv" files will be swept into analysis, so becareful not to include any unrelated .csv files.

#Import Mulitple .csv formatted Raw Data (removing header)
temp <- list.files(path = paste("Input\\RAW_Batch\\"),pattern="*.csv")

for (i in temp) {
  x<-read.csv(paste("Input\\RAW_Batch\\",i,sep=""), header=FALSE, stringsAsFactors=FALSE)
  x<-x[-(1:grep(c("CalcMethod"), x[,2], value=FALSE)),]
  assign(i,x) 
  }

x<-do.call(rbind, lapply(temp, as.symbol))
names(x)<-read.csv("Accessory_Files\\names.csv", header=FALSE, stringsAsFactors=FALSE)

#Subset all important data
x<-data.frame(PeakType=x$PeakType, SampleName=x$DataFileName, InjVolume=x$InjVolume, rt=x$RetTime, Peak_Area=x$MajorArea, Delta_value=x$CorrectedDelta1, AcquisitionDate=x$AcquisitionDate)

#remove ghost entries and hexane washes
x<-droplevels(subset(x, PeakType == "S"))
x<-x[-(grep(c("hexane"), x$SampleName, value=FALSE)),]

#attribute standards
x[grep(c("BAME|FAME"), x$SampleName, value=FALSE),ncol(x)+1]<-c("1")
colnames(x)[ncol(x)]<-"Peak_Standard"
x$Peak_Standard<-as.numeric(x$Peak_Standard)
x[grep(c("Standard|Std"), x$SampleName, value=FALSE),ncol(x)+1]<-c("1")
colnames(x)[ncol(x)]<-"Quantitative_Standard"
x$Quantitative_Standard<-as.numeric(x$Quantitative_Standard)
x[grep(c("C20#2|C20#X"), x$SampleName, value=FALSE),ncol(x)+1]<-c("1")
colnames(x)[ncol(x)]<-"Isotope_Standard"
x$Isotope_Standard<-as.numeric(x$Isotope_Standard)
x[is.na(x)] <- 0 

#Get Aquisition Date as a time class
pattern <-"([0-9]+)/([0-9]+)/([0-9]+) ([0-9]+):([0-9]+)"
x$Date<-x$AcquisitionDate
x$Time<-x$AcquisitionDate
newtimeformat <-"\\4:\\5:00"
newdateformat <-"\\1/\\2/\\3"
x$Date<-gsub(pattern, newdateformat, x$Date, perl=TRUE)
x$Time<-gsub(pattern, newtimeformat, x$Time, perl=TRUE)
dts <- dates(x$Date, format = c(dates = "d/m/y"))
tms <- times(x$Time)
x$AcquisitionDate <- chron(dates = dts, times = tms)
x$Time<-NULL

#Add peakID variable, re-order by rt and spit out to user to validate peaks with Y
x$PeakType<-NULL
x$FA<-NA
x$Project<-NA
x<-x[c(2,4,6,7,8,9,5,1,10,3,11,12)]

print("Does your data come from a variety of projects? Do you wish to separate them at this point?")
switch(menu(c("Yes","No")), operator<-c("Yes"), operator<-("No"))

if(operator == "Yes") {
  
  x<-x[order(x$AcquisitionDate),]
  write.csv(x, file="Output\\project_details.csv")
  print("Take the project_details.csv and mark which samples correspond to which project. Leave any unwanted samples as NA. This would also be the appropriate time to adjust your delta values if you have run standards. When complete, copy the file to the *Input* directory. Samples have been ordered by date.")

  } else {
  
    #Reorder according to retention time
    peakID<-read.csv("Accessory_Files\\peakID.csv")
    x<-rbind(x, peakID)
    x$rt<-as.numeric(as.character(x$rt))
    x<-x[order(x$rt, x$AcquisitionDate),]
    x$rt<-x$rt/10 #Puts rt into seconds
    write.csv(x, file="Output\\confirm_peaks.csv")
    print("Take the confirm_peaks.csv and mark peaks which you deem as good with the appropriate peakID and UF# for likely peaks (which are unidentified) and leave as NA for peaks you wish to remove. When complete, copy the file to the *Input* directory.")
}