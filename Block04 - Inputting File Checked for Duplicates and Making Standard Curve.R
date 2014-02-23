writeLines("Select Complete if you have finished double-checking your peak assignments")
switch(menu(c("Complete")), operator3<-c("Complete"))

if(operator3 == "Complete") {
  x<-read.csv("Input\\duplicated_FAs.csv", stringsAsFactors = FALSE)
  ##Note: did not take into account the slight reduction in delta_value for FAs due to the addition of methyl group during PLFA preparation (according to Lerch, T. Z., M. F. Dignac, E. Barriuso, G. Bardoux, and A. Mariotti.2007. Tracing 2,4-D metabolism inCupriavidus necatorJMP134 with 13C-labelling technique and fatty acid profiling. J. Microbiol. Methods 71:162–174). This is b/c many peaks were not identifiable to their chain length. 
  
##Create levels for FA according to their position in the sequence
  
  #Set-up Quantitative standard quantitation
  standards<-subset(x, Quantitative_Standard == 1)
  #na<-na.exclude(standards$FA)
  #na<-as.numeric(attr(na, 'na.action'))
  #std_names<-standards[-na,]
  standard_values<-as.matrix(read.csv("Accessory_Files\\standards.csv", stringsAsFactors = FALSE, header=FALSE))
  standard_FA<-as.vector(na.omit(unique(standards$FA)))
  
  for(k in standard_FA) {
    
    foo<-subset(standards, FA == k)
    
    for(i in standard_values) {
      foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
    }
    
    #Temporarily display peak_area versus standard dilution to verify the curve looks good
    foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
    foo$Standard_dilution<- 1/foo$SampleName
    
    #Plot curve in ggplot
    graphics.off()
    p <- ggplot(foo, aes(Standard_dilution, Peak_Area))
    p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(k))
    k<-gsub(":", "_", k)
    ggsave(filename=paste("Output\\",k,"_standard_curve.pdf",sep=""), plot = last_plot())
}
  writeLines("\n\nStandard curves have been plotted for each of your standards. Please review these figures and ensure your standards are acceptable for use.")
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
  options(warn=-1)
} 