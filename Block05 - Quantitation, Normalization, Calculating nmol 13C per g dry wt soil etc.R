writeLines("\n\nSelect which standard you would like to quantitate with. If you select *Average,* each peak will be quantitated by each standard separately and these three values averaged.")
switch(menu(c(paste(standard_FA[1]),paste(standard_FA[2]),paste(standard_FA[3]), "Average")), operator<-(paste(standard_FA[1])), operator<-(paste(standard_FA[2])), operator<-(paste(standard_FA[3])), operator<-("Average"))

writeLines("\n\nPlease input required fields based on the preparation of your standards in the spreadsheet entitled *std_conc.csv* in the //Accessory_Files// directory. Once you have done this, hit Proceed.")
switch(menu(c("Proceed")), operator4<-("Proceed"))

if (operator4 == "Proceed"){ 
  
  std_conc<-read.csv("Accessory_Files\\std_conc.csv", stringsAsFactors = FALSE)
  standards<-subset(x, Quantitative_Standard == 1)
  standard_values<-as.matrix(read.csv("Accessory_Files\\standards.csv", stringsAsFactors = FALSE, header=FALSE))
  
  writeLines("\n\nPlease chose whether or not to force regression line through origin.[Note: will not appear in graph, but will be calculated accordingly.]")
  switch(menu(c("Force through origin","Do not force the fit")), force<-("Yes"), force<-("No"))
  
  if (force == "Yes") {
    
    if (operator == paste(standard_FA[1])) {
      
      foo<-subset(standards, FA == paste(operator))
      
      for(i in standard_values) {
        foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
        options(warn=-1)
      }
      
      solo_std<-std_conc[c(grep(paste(operator), std_conc$FA)),]
      
      #Fit peak_area to molC in fatty acids
      foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
      foo$Standard_molC<- 1/foo$SampleName
      foo$Standard_molC<- (((foo$Standard_molC*(solo_std$concentration_in_g_per_mL*1/1000)/solo_std$mol_wt) * solo_std$number_C_per_FA)) ##produces units of mol C per uL
      
      #Get slope equation
      regression_line<-lm(Standard_molC ~ Peak_Area -1, foo)  #y ~ x
      regression_stats<-data.frame(c(regression_line$coefficients))
      row.names(regression_stats)<-c("slope")
      operator<-gsub(":", "_", operator)
      assign(operator, regression_stats)
      
      #Plot curve in ggplot
      graphics.off()
      p <- ggplot(foo, aes(Standard_molC, Peak_Area))
      p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(operator))
      operator<-gsub(":", "_", operator)
      ggsave(filename=paste("Output\\",operator,"mol_per_uL_standard_curve.pdf",sep=""), plot = last_plot())
      
    }
    
    if (operator == paste(standard_FA[2])) {
      
      foo<-subset(standards, FA == paste(operator))
      
      for(i in standard_values) {
        foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
        options(warn=-1)
      }
      
      solo_std<-std_conc[c(grep(paste(operator), std_conc$FA)),]
      
      #Fit peak_area to molC in fatty acids
      foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
      foo$Standard_molC<- 1/foo$SampleName
      foo$Standard_molC<- (((foo$Standard_molC*(solo_std$concentration_in_g_per_mL*1/1000)/solo_std$mol_wt) * solo_std$number_C_per_FA)) ##produces units of mol C per uL    
      
      #Get slope equation
      regression_line<-lm(Standard_molC ~ Peak_Area -1, foo)  #y ~ x
      regression_stats<-data.frame(c(regression_line$coefficients))
      row.names(regression_stats)<-c("slope")
      operator<-gsub(":", "_", operator)
      assign(operator, regression_stats)
      
      #Plot curve in ggplot
      graphics.off()
      p <- ggplot(foo, aes(Standard_molC, Peak_Area))
      p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(operator))
      operator<-gsub(":", "_", operator)
      ggsave(filename=paste("Output\\",operator,"mol_per_uL_standard_curve.pdf",sep=""), plot = last_plot())
      
    }
    
    if (operator == paste(standard_FA[3])) {
      
      foo<-subset(standards, FA == paste(operator))
      
      for(i in standard_values) {
        foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
        options(warn=-1)
      }
      
      solo_std<-std_conc[c(grep(paste(operator), std_conc$FA)),]
      
      #Fit peak_area to molC in fatty acids
      foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
      foo$Standard_molC<- 1/foo$SampleName
      foo$Standard_molC<- (((foo$Standard_molC*(solo_std$concentration_in_g_per_mL*1/1000)/solo_std$mol_wt) * solo_std$number_C_per_FA)) ##produces units of mol C per uL
      
      #Get slope equation
      regression_line<-lm(Standard_molC ~ Peak_Area-1, foo)  #y ~ x
      regression_stats<-data.frame(c(regression_line$coefficients))
      row.names(regression_stats)<-c("slope")
      operator<-gsub(":", "_", operator)
      assign(operator, regression_stats)
      
      #Plot curve in ggplot
      graphics.off()
      p <- ggplot(foo, aes(Standard_molC, Peak_Area))
      p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(operator))
      operator<-gsub(":", "_", operator)
      ggsave(filename=paste("Output\\",operator,"mol_per_uL_standard_curve.pdf",sep=""), plot = last_plot())
      
    }
    
    if (operator == "Average") {
      
      for(m in standard_FA) {
        foo<-subset(standards, FA == paste(m))
        
        for(i in standard_values) {
          foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
          options(warn=-1)
        }
        
        solo_std<-std_conc[c(grep(paste(m), std_conc$FA)),]
        
        #Fit peak_area to molC in fatty acids using data input from *standard_concentrations.csv*
        foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
        foo$Standard_molC<- 1/foo$SampleName
        foo$Standard_molC<- (((foo$Standard_molC*(solo_std$concentration_in_g_per_mL*1/1000)/solo_std$mol_wt) * solo_std$number_C_per_FA)) ##produces units of mol C per uL
        
        #Get slope equation
        regression_line<-lm(Standard_molC ~ Peak_Area-1, foo)  #y ~ x
        regression_stats<-data.frame(c(regression_line$coefficients))
        row.names(regression_stats)<-c("slope")
        m<-gsub(":", "_", m)
        assign(m, regression_stats)
        
        #Plot curve in ggplot
        graphics.off()
        p <- ggplot(foo, aes(Standard_molC, Peak_Area))
        p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(m))
        ggsave(filename=paste("Output\\",m,"mol_per_uL_standard_curve.pdf",sep=""), plot = last_plot())
      }
    }
    
  } else {
    
    if (operator == paste(standard_FA[1])) {
  
      foo<-subset(standards, FA == paste(operator))
      
      for(i in standard_values) {
        foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
        options(warn=-1)
      }
      
      solo_std<-std_conc[c(grep(paste(operator), std_conc$FA)),]
      
      #Fit peak_area to molC in fatty acids
      foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
      foo$Standard_molC<- 1/foo$SampleName
      foo$Standard_molC<- (((foo$Standard_molC*(solo_std$concentration_in_g_per_mL*1/1000)/solo_std$mol_wt) * solo_std$number_C_per_FA)) ##produces units of mol C per uL
        
      #Get slope equation
      regression_line<-lm(Standard_molC ~ Peak_Area, foo)  #y ~ x
      regression_stats<-data.frame(c(regression_line$coefficients))
      row.names(regression_stats)<-c("intercept","slope")
      operator<-gsub(":", "_", operator)
      assign(operator, regression_stats)
      
      #Plot curve in ggplot
      graphics.off()
      p <- ggplot(foo, aes(Standard_molC, Peak_Area))
      p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(operator))
      operator<-gsub(":", "_", operator)
      ggsave(filename=paste("Output\\",operator,"mol_per_uL_standard_curve.pdf",sep=""), plot = last_plot())
      
  }

  if (operator == paste(standard_FA[2])) {
    
    foo<-subset(standards, FA == paste(operator))
    
    for(i in standard_values) {
      foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
      options(warn=-1)
    }
    
    solo_std<-std_conc[c(grep(paste(operator), std_conc$FA)),]
    
    #Fit peak_area to molC in fatty acids
    foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
    foo$Standard_molC<- 1/foo$SampleName
    foo$Standard_molC<- (((foo$Standard_molC*(solo_std$concentration_in_g_per_mL*1/1000)/solo_std$mol_wt) * solo_std$number_C_per_FA)) ##produces units of mol C per uL    
    
    #Get slope equation
    regression_line<-lm(Standard_molC ~ Peak_Area, foo)  #y ~ x
    regression_stats<-data.frame(c(regression_line$coefficients))
    row.names(regression_stats)<-c("intercept","slope")
    operator<-gsub(":", "_", operator)
    assign(operator, regression_stats)
    
    #Plot curve in ggplot
    graphics.off()
    p <- ggplot(foo, aes(Standard_molC, Peak_Area))
    p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(operator))
    operator<-gsub(":", "_", operator)
    ggsave(filename=paste("Output\\",operator,"mol_per_uL_standard_curve.pdf",sep=""), plot = last_plot())
    
  }

  if (operator == paste(standard_FA[3])) {
    
    foo<-subset(standards, FA == paste(operator))
    
    for(i in standard_values) {
      foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
      options(warn=-1)
    }
    
    solo_std<-std_conc[c(grep(paste(operator), std_conc$FA)),]
    
    #Fit peak_area to molC in fatty acids
    foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
    foo$Standard_molC<- 1/foo$SampleName
    foo$Standard_molC<- (((foo$Standard_molC*(solo_std$concentration_in_g_per_mL*1/1000)/solo_std$mol_wt) * solo_std$number_C_per_FA)) ##produces units of mol C per uL
    
    #Get slope equation
    regression_line<-lm(Standard_molC ~ Peak_Area, foo)  #y ~ x
    regression_stats<-data.frame(c(regression_line$coefficients))
    row.names(regression_stats)<-c("intercept","slope")
    operator<-gsub(":", "_", operator)
    assign(operator, regression_stats)
    
    #Plot curve in ggplot
    graphics.off()
    p <- ggplot(foo, aes(Standard_molC, Peak_Area))
    p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(operator))
    operator<-gsub(":", "_", operator)
    ggsave(filename=paste("Output\\",operator,"mol_per_uL_standard_curve.pdf",sep=""), plot = last_plot())
    
  }
  
        
  if (operator == "Average") {
    
    for(m in standard_FA) {
      foo<-subset(standards, FA == paste(m))
               
            for(i in standard_values) {
            foo[c(grep(paste(i), foo$SampleName)),  which(colnames(foo) == "SampleName")] <- paste(i)
            options(warn=-1)
            
            }  
      
      solo_std<-std_conc[grep(paste(m), std_conc$FA),]
    
      #Fit peak_area to molC in fatty acids using data input from *standard_concentrations.csv*
      foo$SampleName<-as.numeric(gsub("1_([0-9]+)_", "\\1", foo$SampleName))
      foo$Standard_molC<- 1/foo$SampleName
      foo$Standard_molC<- (((foo$Standard_molC*(solo_std$concentration_in_g_per_mL*1/1000)/solo_std$mol_wt) * solo_std$number_C_per_FA)) ##produces units of mol C per uL
      
      #Get slope equation
      regression_line<-lm(Standard_molC ~ Peak_Area, foo)  #y ~ x
      regression_stats<-data.frame(c(regression_line$coefficients))
      row.names(regression_stats)<-c("intercept","slope")
      m<-gsub(":", "_", m)
      assign(m, regression_stats)
              
      #Plot curve in ggplot
      graphics.off()
      p <- ggplot(foo, aes(Standard_molC, Peak_Area))
      p + layer(foo, geom = "point") + geom_smooth(method = "lm") + labs(title = paste(m))
      ggsave(filename=paste("Output\\",m,"mol_per_uL_standard_curve.pdf",sep=""), plot = last_plot())
    
  }
  }
  }
}

#Find all places where InjVolume is not equal to standards and scale the Peak_Area accordingly
grep_InjVolume<-x[grep("Standard|Std", x$SampleName), which(colnames(x) == "InjVolume")]
grep_InjVolume<-unique(grep_InjVolume)

for(q in which(x$InjVolume != grep_InjVolume)) {
  scalar<- grep_InjVolume / x[q, match("InjVolume",names(x))] 
  x[q, match("Peak_Area",names(x))] <- x[q, match("Peak_Area",names(x))] * scalar
}

#Clean up standards used for peakID
x<-subset(x, Peak_Standard == 0 & Quantitative_Standard == 0 & Isotope_Standard == 0)
x<-data.frame(SampleName=x$SampleName, rt=x$rt, Delta_value=x$Delta_value, Peak_Area=x$Peak_Area, FA=x$FA, AcquisitionDate=x$AcquisitionDate, stringsAsFactors = FALSE)
x<-x[is.na(x$FA)==FALSE,]
x$Quantitative_Standard <- NULL

# Normalize Data if Desired (first check if all samples have FA present necessary for normalization AND print histogram)

writeLines("\n\nSelect which standard you have used for an internal standard.")
switch(menu(c(paste(standard_FA[1]),paste(standard_FA[2]),paste(standard_FA[3]))), opera<-(paste(standard_FA[1])), opera<-(paste(standard_FA[2])), opera<-(paste(standard_FA[3])))

opera<-gsub("_", ":", opera)

writeLines("List of unique samples found in your dataset:\n")
writeLines(paste(unique(x$SampleName)))

writeLines(paste("\nNumber of unique samples found in your dataset:\n", length(unique(x$SampleName))))

writeLines(paste("\nNumber of samples containing your internal standard:\n", length(grep(paste("^",opera,"$",sep=""),x$FA)),"\n"))

c16_plot<-data.frame(Peak_Area=x[grep(paste("^",opera,"$",sep=""),x$FA), match("Peak_Area", colnames(x))])
graphics.off()
p <- ggplot(c16_plot, aes(x=Peak_Area)) + geom_histogram(binwidth = mean(c16_plot$Peak_Area)/8)
p
ggsave((filename="Output\\c16_histogram.pdf"), plot = last_plot())

shap<-shapiro.test(c16_plot$Peak_Area)

writeLines(paste("\nIf these numbers do not match, accounting for samples you didn't run with an internal standard, you may not be able to normalize your data. To see whether your data is roughly normal without any transformation, a histogram showing the distribution of Peak_Area from FA c16:0, a commonly occuring FA, has been created in the \\Output\\ folder.\n\nAccording to the Wilks-Shapiro test, you should reject the assumption c16:0 is normally distributed if your p-value is greater than 0.05. Your p-value is:\n\n",paste(shap$p.value),"\n\n c16:0 is found in the following number of samples:\n\n",length(grep("^c16:0$",x$FA))))

writeLines("\nWould you like to normalize your samples to your internal standard?\n\nNote: Only samples containing the internal standard will be normalized, likely biasing your dataset.")
switch(menu(c("Yes", "No")), op<-"Yes", op<-"No")

if (op == "Yes") {
  
  norm<-data.frame(cbind(SampleName=x[grep(paste("^",opera,"$", sep=""),x$FA), match("SampleName", colnames(x))], Norm=(x[grep(paste("^",opera,"$", sep=""),x$FA), match("Peak_Area", colnames(x))])), stringsAsFactors = FALSE)
  norm$Norm<-as.numeric(norm$Norm)
  norm$Norm<-norm$Norm/max(norm$Norm)
  dummyappend<-data.frame(cbind(SampleName=unique(x$SampleName), Norm=NA), stringsAsFactors = FALSE)
  norm<-rbind(norm,dummyappend)
  norm<-norm[-which(duplicated(norm$SampleName)),]

    for (YOLO in 1:nrow(x)) {
      foo<-norm[grep(x$SampleName[[YOLO]],norm$SampleName),]
      foo$Norm<-as.numeric(foo$Norm)
      
        if (is.na(foo$Norm) == TRUE) {
          x$Peak_Area[[YOLO]]<-x$Peak_Area[[YOLO]]
        } else {
            x$Peak_Area[[YOLO]]<-x$Peak_Area[[YOLO]]*foo$Norm
      }
    }
  } else {
      writeLines("\nYou have opted not to normalize peak area.\n\n")
}

#Remove internal standard from downstream analyses
x<-x[-grep(paste("^",opera,"$", sep=""), x$FA),]
standard_FA<-gsub(":", "_", standard_FA)

#Quantitate Peaks

if (force == "Yes") {
  
if (operator == "Average") {
  
  for(n in 1:length(standard_FA)) {
    foo1 <- get(standard_FA[n])
    molC <- foo1[1,] * x$Peak_Area
    assign(standard_FA[n], molC)
    }
  
  foox<-do.call(cbind, lapply(standard_FA, as.symbol))
  x<-cbind(x,foox)
  
  #Average all of the columns produced for each FA standard
  goo <- list((ncol(x)-(length(standard_FA)-1)):ncol(x))
  x$molC_per_uL <- do.call(cbind, lapply(goo, function(i) rowMeans(x[, i])))
  
}  

if (operator == paste(standard_FA[1])) {
  
  foo1 <- get(standard_FA[1])
  x$molC_per_uL <- foo1[1,] * x$Peak_Area
  
}

if (operator == paste(standard_FA[2])) {
  
  foo1 <- get(standard_FA[2])
  x$molC_per_uL <- foo1[1,] * x$Peak_Area
  
}

if (operator == paste(standard_FA[3])) {
  
  foo1 <- get(standard_FA[3])
  x$molC_per_uL <- foo1[1,] * x$Peak_Area
  
}
} else {
  
  if (operator == "Average") {
    
    for(n in 1:length(standard_FA)) {
      foo1 <- get(standard_FA[n])
      molC <- foo1[2,] * x$Peak_Area + foo1[1,]
      assign(standard_FA[n], molC)
    }
    
    foox<-do.call(cbind, lapply(standard_FA, as.symbol))
    x<-cbind(x,foox)
    
    #Average all of the columns produced for each FA standard
    goo <- list((ncol(x)-(length(standard_FA)-1)):ncol(x))
    x$molC_per_uL <- do.call(cbind, lapply(goo, function(i) rowMeans(x[, i])))
    
  }  
  
  if (operator == paste(standard_FA[1])) {
    
    foo1 <- get(standard_FA[1])
    x$molC_per_uL <- foo1[2,] * x$Peak_Area + foo1[1,]
    
  }
  
  if (operator == paste(standard_FA[2])) {
    
    foo1 <- get(standard_FA[2])
    x$molC_per_uL <- foo1[2,] * x$Peak_Area + foo1[1,]
    
  }
  
  if (operator == paste(standard_FA[3])) {
    
    foo1 <- get(standard_FA[3])
    x$molC_per_uL <- foo1[2,] * x$Peak_Area + foo1[1,]
    
  } 
}
#Get user to input volume of hexane in which the sample FAs were suspended during 

cat("In your final step of preparing your PLFA sample(s), you suspended the dried FAs in hexane. Please enter the volume, in microlitres, (uL) which you re-suspended your sample in:","\n") # prompt
vol_hexane<-scan(n=1)
x$molC <- x$molC_per_uL * vol_hexane

##Account for soil dry weights
writeLines("\n\nA spreadsheet entitled *factor_list.csv* has been created in the output directory. \n\nPlease fill in the dry weight of each sample AND include any additional factors you are interested in analyzing (go ahead and rename the column header *Factor* to better describe your factor. \n\nYou can rename your samples using the *Sample_Name* column provided [do not alter the *Sample_List* column]. \n\nIf unenriched control samples are present, please mark these with a zero in the column designated *enriched*. If unenriched controls are absent the final datasheet will contain a column *C13excess_umols_per_drywt* which is empty.")

Sample_list<-unique(x$SampleName)
factor_list<-data.frame(Sample_List=Sample_list, Sample_Name = NA, Dry_wt_Sample_in_grams= 0, Enriched = 1, Factor = NA)
write.csv(factor_list, file="Output\\factor_list.csv")


writeLines("\nPlease copy the *factor_list.csv* to the Input folder and hit Proceed.")
switch(menu(c("Proceed")), operator5<-("Proceed"))

if(operator5 == "Proceed") {
  factor_list<-read.csv("Input\\factor_list.csv", stringsAsFactors = FALSE)
}

#Assigning soil wt to your dataset
for(factor_list_names in colnames(factor_list)) {
  Soil_column<- data.frame(1:nrow(x))
  Soil_column[,]<- NA
  colnames(Soil_column)<-paste(factor_list_names)
  assign(factor_list_names, Soil_column)
}

foox<-do.call(cbind, lapply(colnames(factor_list), as.symbol))
x<-cbind(x,foox)

#Assigning factors to your dataset
for(ID in Sample_list) {
  for(Column_names in colnames(factor_list)) {
    x[grep(paste(ID), x$SampleName), match(Column_names, names(x))] <- factor_list[grep(paste(ID), factor_list$Sample_List), match(Column_names, names(factor_list))]  
  }
}

for (i in 1:nrow(x)) {
  num<-grep(x[i, match("SampleName", colnames(x))], factor_list$Sample_List)
    if (is.na(factor_list[num, match("Sample_Name", colnames(factor_list))]) == FALSE) {
    x[i, match("SampleName", colnames(x))] <- factor_list[num, match("Sample_Name", colnames(factor_list))]
    } 
}

x$Sample_Name<-NULL
x$Sample_List<- NULL
x$X<- NULL

#Rename SampleName according to user inputted "Sample_Name"

#Create function to convert delta_value into mass fraction
percent_13C <- function(x) {
  calc1<-((x/1000+1)*0.0112372)
  calc2<-(calc1/(calc1+1))
  return(calc2)
}

##Convert to nmol 13C Carbon
x$umols13C_per_drywt<-as.numeric(((x$molC*percent_13C(x$Delta_value))*1000000)/x$Dry_wt_Sample_in_grams)
x$umols12C_per_drywt<-as.numeric(((x$molC*(1-percent_13C(x$Delta_value)))*1000000)/x$Dry_wt_Sample_in_grams)

## Calculate the 13C excess - try to match FAs to FAs in control, otherwise take average of controls and use that.
unlabeled_controls<-subset(x, x$Enriched == 0)
mean_unlabeled_controls<-ddply(unlabeled_controls, ~ FA, summarize, mean_delta=mean(Delta_value))
FA_list<-unique(x$FA)
FA_list_controls<-unique(mean_unlabeled_controls$FA)
FA_available<-rep(0, length(FA_list))

count=1
  for (u in FA_list) {
    if (any(grepl(paste("^",u,"$",sep=""), FA_list_controls)) == TRUE) {
      FA_available[[count]]<-1
      count=count+1
    } else {
      FA_available[[count]]<-0
      count=count+1
    }
}

FA_list<-as.data.frame(cbind(FA_list, FA_available), stringsAsFactors = FALSE)
FA_list$FA_available<-as.numeric(FA_list$FA_available)
FA_list$mean_delta<- NA

for (count in 1:nrow(FA_list)) {
  if (FA_list$FA_available[[count]] == 1) {
    FA_list$mean_delta[[count]]<- mean_unlabeled_controls[grep(paste("^",FA_list$FA_list[[count]],"$",sep=""), mean_unlabeled_controls$FA), 2]
    } else {
    FA_list$mean_delta[[count]]<- mean(unlabeled_controls$Delta_value)
    }
}

#Loop taking the mean delta-value for every case where FA_list=1 and ELSE taking the mean of all FAs
x$i13excess<-NA

for (i in 1:nrow(x)) {
  for (FA in FA_list$FA_list) {
    if(any(grep(paste(FA), x[i,match("FA", colnames(x))])) == TRUE) {
      
      x$i13excess[i]<-x$umols13C_per_drywt[i]-(as.numeric((x$molC[i]*percent_13C(FA_list[grep(paste(FA), FA_list$FA_list), match("mean_delta", colnames(FA_list))])*1000000)/x$Dry_wt_Sample_in_grams[i]))        
    
    } 
  }
}

writeLines("\n\nHave you used a glassware control in preparation of your samples? And, would you like to *blank* your samples against it?")
switch(menu(c("Yes", "No")), op<-"Yes", op<-"No")

if (op == "Yes") {
  writeLines("\n\nWhat have you named your glassware control sample in your datasheet?")
  glass_inp<-scan(n=1, what = character())
  glass<-data.frame(cbind(FA=x[grep(paste(glass_inp), x$SampleName), match("FA", colnames(x))], umols13C_per_drywt=x[grep(paste(glass_inp), x$SampleName), match("umols13C_per_drywt", colnames(x))], umols12C_per_drywt=x[grep(paste(glass_inp), x$SampleName), match("umols12C_per_drywt", colnames(x))], i13excess=x[grep(paste(glass_inp), x$SampleName), match("i13excess", colnames(x))]), stringsAsFactors = FALSE)
  glass$umols13C_per_drywt<-as.numeric(glass$umols13C_per_drywt)
  glass$umols12C_per_drywt<-as.numeric(glass$umols12C_per_drywt)
  glass$i13excess<-as.numeric(glass$i13excess)
  
  for (FA in glass$FA) {
    for (count in 1:nrow(x)) {
      if (any(grep(paste("^",FA,"$", sep=""), x[count,match("FA", colnames(x))])) == TRUE) {
        x$umols13C_per_drywt[[count]]<- x$umols13C_per_drywt[[count]]-glass[grep(paste("^",FA,"$", sep=""), glass$FA), match("umols13C_per_drywt", colnames(glass))]
        x$umols12C_per_drywt[[count]]<- x$umols12C_per_drywt[[count]]-glass[grep(paste("^",FA,"$", sep=""), glass$FA), match("umols12C_per_drywt", colnames(glass))]
        x$i13excess[[count]]<- x$i13excess[[count]]-glass[grep(paste("^",FA,"$", sep=""), glass$FA), match("i13excess", colnames(glass))]
      }
    }
  }
  x<-x[-grep(paste(glass_inp), x$SampleName),]
}

x[which(x[, "i13excess"] < 0), match("i13excess", colnames(x))] <- 0

colnames(x)[match("i13excess", colnames(x))]<- "C13excess_umols_per_drywt"

writeLines("\n\nAll negative 13C Excess values have been set to zero, which assumes zero enrichment above background.")
  
write.csv(x, file="Output\\Final_datasheet.csv")
writeLines("\nA spreadsheet containing quantitated FA peaks has been created in the //Output// directory.")
writeLines("\n\nIt is recommended to copy and paste the history of your R session, in order to have the exact details of each session.")
