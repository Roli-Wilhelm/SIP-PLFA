writeLines("\n\nDid you previously chose to input project details?")
switch(menu(c("Yes","No")), operator<-c("Yes"), operator<-("No"))

if(operator == "Yes") {
  x<-read.csv("Input\\project_details.csv")
  x$FA[1]<-"Mark Your Assignments HERE"
  x<-x[,colSums(is.na(x)) != nrow(x)]
  
  writeLines("\n\nWould you like to separate your projects and work with them individually? Note: Recommended if data was collected across wide timeframe?")
  switch(menu(c("Yes","No")), operator<-c("Yes"), operator<-("No"))
  
  if(operator == "Yes") {
    Project_Names<-unique(x$Project)
    Workingdirectory<-getwd()
    peakID<-read.csv("Accessory_Files\\peakID.csv")
    
      for (i in Project_Names) {
      File_name<- toString(i)
      foo <-subset(x, Project == File_name)
      foo$X <- NULL
      foo <-rbind(foo, peakID)
      foo$rt<-as.numeric(as.character(foo$rt))
      foo<-foo[order(foo$rt, foo$AcquisitionDate),]
      foo$rt<-foo$rt/10 #Puts rt into seconds
      assign(i,foo)
      
      dir.create(paste(Workingdirectory,"/Output/",File_name,sep=""))
      write.csv(foo, file=paste(Workingdirectory,"/Output/",File_name,"/","confirm_peaks.csv",sep=""))
      }
  writeLines("\n\nYou will now see a directory for each project containing a file named *confirm_peaks.csv*. From here on out, you will have to move each file to \\Input\\ subdirectory after each output had been modified for the script to proceed")
  writeLines("\nTake the confirm_peaks.csv and mark peaks which you deem as good with the appropriate peakID and UF# for likely peaks (which are unidentified) and leave as NA for peaks you wish to remove")
    
    } else {
  x<-x[,colSums(is.na(x)) != nrow(x)]
  peakID<-read.csv("Accessory_Files\\peakID.csv", stringsAsFactors = FALSE)
  x<-rbind(x, peakID)
  x$X <- NULL
  x$rt<-as.numeric(as.character(x$rt))
  x<-x[order(x$rt, x$AcquisitionDate),]
  x$rt<-x$rt/10 #Puts rt into seconds
  write.csv(x, file="Output\\confirm_peaks.csv")
  writeLines("\n\nTake the confirm_peaks.csv and mark peaks which you deem as good with the appropriate peakID and UF# for likely peaks (which are unidentified) and leave as NA for peaks you wish to remove. It is a good idea to save this file in a different directory in case the analysis is re-run")
  }
}