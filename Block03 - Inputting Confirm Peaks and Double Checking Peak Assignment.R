writeLines("\n\nSelect GO if you have finished identifying your peaks")
switch(menu(c("GO")), operator2<-c("GO"))

if(operator2 == "GO") {
  x<-read.csv("Input\\confirm_peaks.csv", stringsAsFactors = FALSE)
} else {
  
  x<-read.csv("Input\\confirm_peaks.csv", stringsAsFactors = FALSE)
  
}

#Remove columns with random NA
x <- x[,colSums(is.na(x))<nrow(x)]

#Look for misassigned FAs (i.e. samples which have been given multiple peaks for a single FA)

Duplicates<-ddply(x, ~ FA + SampleName, summarize, Duplicated=duplicated(SampleName))
x<-x[order(x$FA, x$SampleName),]
Duplicates<-Duplicates[order(Duplicates$FA, x$SampleName),]
Duplicates<-cbind(x, Duplicates$Duplicated)
Duplicates<-Duplicates[order(Duplicates$rt),]
write.csv(Duplicates, file="Output\\duplicated_FAs.csv")
writeLines("\n\nA file entitled duplicated_FAs has now been created. Samples which have been assigned more than one peak per FA will be marked as TRUE. Please correct duplicate peaks by re-assigning them; those marked *NA* will be removed.")
