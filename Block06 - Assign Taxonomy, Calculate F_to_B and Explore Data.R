x<-read.csv("Input\\Final_datasheet.csv", header=TRUE, stringsAsFactors=FALSE)

##Assign taxonomy and Analyze Taxonomic Breakdown
taxonomy<-read.csv("Accessory_Files\\taxonomy.csv", header=TRUE, stringsAsFactors=FALSE)
x$Taxonomy<- NA

##Separate Unclassified FAs from classified
UF<-x[grep("UF", x$FA),]
x<-x[-grep("UF", x$FA),]

for(tax in unique(x$FA)) {
  x[grep(paste(tax), x$FA), match("Taxonomy", names(x))] <- taxonomy[grep(paste("^",tax,"$",sep=""), taxonomy$FA), match("Taxonomy", names(taxonomy))]  
}

x<-rbind(x,UF)
x[which(is.na(x$Taxonomy)), match("Taxonomy", colnames(x))] <- "Unclassified"

writeLines("\nA taxonomical classification has now been assigned to each known FA.\n")

#Subsetting only enriched samples
x_e<-subset(x, Enriched == 1)
x_e$umols13C_per_drywt[which(x_e$umols13C_per_drywt < 0)] <- 0
x_e<-x_e[order(x_e$rt),]
FA_Levels<-na.omit(unique(x_e$FA))
x_e$FA<-factor(x_e$FA, levels = FA_Levels)

x_e$Taxonomy<-factor(x_e$Taxonomy, levels = rev(c("Total Biomass","Fungi","Plant","General Eukaryotes","Gram-positive Bacteria","Gram-negative Bacteria","Unclassified")))

writeLines("\nWould you like to plot a histogram of all FAs and set a threshold for the minimum number of occurrences across samples in order to be inclused in analyses?\n")
switch(menu(c("Yes","No")), proceed<-c("proceed"), proceed<-c("skip"))

if (proceed == "proceed") {
  
  #Plot distribution of FA counts. histogram of FA counts
  x_e$count<- 1
  FA_count<-ddply(x_e, ~ FA, summarize, count=sum(count))
  graphics.off()
  qplot(count, data=FA_count, geom="histogram", binwidth = max(FA_count$count)/5)
  ggsave(filename="Output\\FA_histogram.pdf", plot = last_plot())
  write.csv(FA_count, file="Output\\FA_count.csv")
  x_e$count<- NULL

  writeLines("\nA histogram and spreadsheet have been created in the \\Output\\ folder.\n")

#Remove FAs with a total count less than user input value

  cat("Please enter a minimum number of occurrences for a FA to be included in further analyses")
  FA_remove<-scan(n=1, what = numeric())
  FA_to_Remove<-FA_count[as.numeric(which(FA_count$count < FA_remove[1])),1]

  for(yip in FA_to_Remove) {
    x_e<-x_e[-grep(paste("^",yip,"$",sep=""), x_e$FA),]
  }
}

if (any(is.na(x_e$C13excess_umols_per_drywt)) == TRUE) {
  
  writeLines("It appears that you did not include any unenriched controls. If this is not the case, please select \"Use Excess 13C Anyways\" (Option 2)")
  switch(menu(c("Yes, that is correct", "Use Excess 13C Anyways")), nextup<-c("Correct"), dummy<-c("whatever"))
  
  if (nextup == "Correct") {
    x_e$C13excess_umols_per_drywt<-NULL
  
    ##Calculate Fungal:Bacterial Ratio Modified from (Hogberg 2013: 18:2ω6,9 : (16:1ω7+cy17:0+17:0+18:1ω7+cy19 + i15:0+a15:0+15:0+i16:0+a17:0)). I included additional peaks for fungi and other peaks for bacteria (which I have never seen annotated, so are unlikely to factor)
    cat("We will now calculate the Fungi:Bacteria using a modified version of the equation used by Hogberg et al. 2013 (DOI 10.1007/s11104-013-1742-9). You may calculate the F:B ratio ratio for any one of your factors or for every sample (recommended). To calculate it for all samples, ENTER: SampleName, to use a factor, ENTER: the exact name of your factor, ensuring it matches the exact column heading you created in your \"factor_list.csv\". If you would like to do more than one factor, you'll have to re-run the script. Lastly, if your F:B ratio is returned as \"NoBacteria\" this is because one of your samples within your factor had no identified bacteria, not necessarily all of your samples.")
    FB_factor<-scan(n=1, what = character())

    FB<-ddply(x_e, ~ get(FB_factor) + FA + Taxonomy + Enriched, summarize, umols13C_per_drywt=umols13C_per_drywt, umols12C_per_drywt=umols12C_per_drywt)
    colnames(FB)[1]<-c("Factor")
    Fact_list<-as.vector(na.omit(unique(FB$Factor)))

    #13C enriched samples F:B Ratio
    FB_enriched<-subset(FB, Enriched == 1)

    FB_allC<-rep(1, length(Fact_list))
    FB_13C<-rep(1, length(Fact_list))
    FB_12C<-rep(1, length(Fact_list))
    
    count = 1
    for(fact in Fact_list) {
      gargle<-subset(FB_enriched, Factor == paste(fact))
      gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols13C_per_drywt)+sum(umols12C_per_drywt))
      oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
      if(any(oh) == TRUE) {
        if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
          if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
            
          } else {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

          }

        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        }
      } else {
        gargle<-c("No Bacteria")
      }
      FB_allC[[count]]<-gargle
      count = count + 1
    }

    count = 1
    for(fact in Fact_list) {
      gargle<-subset(FB_enriched, Factor == paste(fact))
      gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols13C_per_drywt))
      oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
      if(any(oh) == TRUE) {
        if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
          if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
            
          } else {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

          }

        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        }
      } else {
        gargle<-c("No Bacteria")
      }
      FB_13C[[count]]<-gargle
      count = count + 1
    }

    count = 1
    for(fact in Fact_list) {
      gargle<-subset(FB_enriched, Factor == paste(fact))
      gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols12C_per_drywt))
      oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
      if(any(oh) == TRUE) {
        if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
          if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
            
          } else {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

          }

        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        }
      } else {
        gargle<-c("No Bacteria")
      }
      FB_12C[[count]]<-gargle
      count = count + 1
      }

    FB_ratio<-data.frame(cbind(Factor=Fact_list,Total_FA=FB_allC, C12_only=FB_12C, C13_only=FB_13C))
    
    if (FB_factor == "SampleName") {
      
      #Add factors to FB ratio
      factors<-colnames(factor_list)[match("Enriched", colnames(factor_list))+1:ncol(factor_list)]
      factors<-factors[which(is.na(factors) == FALSE)]
      
      for (ip in factors) {
        
        fooz <- data.frame(x_e[,match(paste(ip), colnames(x_e))])
        assign(paste(ip),fooz)
      }
      
      
      fooz<-do.call(cbind, lapply(factors, as.symbol))
      fooz<-cbind(fooz, SampleName=x_e[,match("SampleName", colnames(x_e))])
      
      count=1
      for (count in 1:ncol(fooz)-1) {
        
        colnames(fooz)[count] <- factors[count]
        
      }
      
      #Remove duplicates
      fooz<-fooz[-which(duplicated(fooz$SampleName)),]
      
      #Ensure same order
      fooz<-fooz[order(fooz$SampleName),]
      FB_ratio<-FB_ratio[order(FB_ratio$Factor),]
      
      #Bind dataframes
      FB_ratio<-cbind(FB_ratio, fooz)
      
    }
    
    #Create directory according to factor selected
    Workingdirectory<-getwd()
    dir.create(paste(Workingdirectory,"/Output/","/",FB_factor,"/",sep=""))
    write.csv(FB_ratio, file=paste("Output\\",paste(FB_factor),"\\FB_ratios.csv",sep=""))
    writeLines("\n\nA spreadsheet entitled \"FB_ratios.csv\" has been created in the \\Output\\ folder under a folder title corresponding to the factor your selected above. This may simply be \"SampleName\".\n\n")

    ##Now analysis will be done according to a chosen factor
    writeLines("\nWe will now create some additional graphs according to a factor you select.\n")
    cat("Please enter a factor by which you would like to separate your samples. \n\nIf none, ENTER: SampleName")
    Factor<-scan(n=1, what = character())

    #Create empty object with name of factor; this is necessary for calling the factor in subsequent functions like ddply
    empty<-NA
    assign(paste(Factor),empty)
    Factor_data<-x_e[,grep(paste(Factor), colnames(x_e))]

    #Setting the order of variables
    writeLines("\nIf there is a particular order to the variables withing your factor please input the order now. Enter each variable in a separate line according to the order you wish them placed. For example, \nfactor1\nfactor2\n...\nfactor(n). If none, please Enter: NONE. ALSO!! If any of your names contain a space, be sure to use quotation marks around the entire phrase.")
    factor_order<-scan(n=length(unique(Factor_data)), what = character())
    
    if (any(factor_order == "None") == TRUE) {
    print("No Factor Order Input")
    } else {
    x_e<-cbind(x_e, Factor=Factor_data)
    x_e$Factor_data<-factor(Factor_data, levels = factor_order)
  
    }
    ##Show crazy patterns of 13C Enrichment across your factor according to FAs
    FA_13Cenrich<- ddply(x_e, ~  Factor_data + FA, summarize, umols13C_per_drywt=sum(umols13C_per_drywt))
    
    #Create directory according to factor selected
    Workingdirectory<-getwd()
    dir.create(paste(Workingdirectory,"/Output/","/",Factor,"/",sep=""))
    
    graphics.off()
    p<-ggplot(FA_13Cenrich, aes(Factor_data, umols13C_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enriched13C_FAs.pdf",sep=""), height=12, width=10, plot = last_plot())
  
    graphics.off()
    p<-ggplot(FA_13Cenrich, aes(Factor_data, umols13C_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enriched13C_FAs.pdf",sep=""), height=12, width=10, plot = last_plot())
    
    
    #Same graph in black and white
    myColors_FA_BnW<-rev(gray.colors(length(levels(FA_13Cenrich$FA)), start = 0, end = 1, gamma = 2.2, alpha = NULL))
    names(myColors_FA_BnW) <- levels(FA_13Cenrich$FA)
    colScale_FA_BnW <- scale_fill_manual(name = "FA",values = myColors_FA_BnW)
  
    graphics.off()
    p<-ggplot(FA_13Cenrich, aes(Factor_data, umols13C_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
    p + colScale_FA_BnW 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enriched13C_FAs_BnW.pdf",sep=""), height=12, width=10, plot = last_plot())
      
    #Incorrect graph, but with the legend you want.
    #graphics.off()
    #p<-ggplot(FA_13Cenrich, aes(Factor_data, umols13C_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1)) + geom_bar(colour = "black", show_guide=FALSE)
    #p
    #ggsave(filename="Output\\Enriched13C_Proportional_FAs.pdf", height=12, width=10, plot = last_plot())
  
    #writeLines("\n Two general graphs have been printed in the \Output\ folder. Pleas take a look at them. These maybe useful, but are mostly to inform whether you'd like colour graphs or black and white. Please consult them now.\n")
    
    #writeLines("\nPlease select a colour scheme.\n")
    #switch(menu(c("Colourful","Black and White", "Colour-blindness Aware")), colz<-c("Colourful"), colz<-c("BnW"), colz<-c("CBA"))
  
    ## Show Breakdown of enrichment across Taxonomy in Bar charts according to a factor
    Taxo_13Cenrich<- ddply(x_e, ~  Factor_data + Taxonomy, summarize, umols13C_per_drywt=sum(umols13C_per_drywt))
            
    #if (colz == "Colourful") {
    #  ##Create a custom color scale for taxonomy
    #  myColors_tax <- NA
    #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
    
    #}
  
    #if (colz == "BnW") {
    #  #Black and White Scale
    #  myColors_tax<-rev(gray.colors(6, start = 0.1, end = 1, gamma = 2.2, alpha = NULL))
    #  names(myColors_tax) <- levels(Taxo_13Cenrich$Taxonomy)
    #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
    
    #}
  
    #if (colz == "CBA") {
    #  #Color-blind palette
    #  myColors_tax<-rev(c("#331A00","#996136","#D9AF98","#99F8FF","#33E4FF","#007A99"))
    #  names(myColors_tax) <- levels(Taxo_13Cenrich$Taxonomy)
    #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
    
    #}
  
    ##TEST YOUR PALETTE FOR COLOR BLINDNESS
    #myColors_tax<-dichromat(myColors_tax)
    #colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
  
    #Bar plots showing distribution of Enrich 13C in PLFAs With fill (i.e. by percent of total)
    graphics.off()
    p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass"), aes(Factor_data, umols13C_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Taxonomic Distribution of Enriched Carbon in PLFA", axis.title=element_text(hjust=0.1))
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enrich13C_Taxa.pdf",sep=""), plot = last_plot())
  
    #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
    graphics.off()
    p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass"), aes(Factor_data, umols13C_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Taxonomic Distribution of Enriched Carbon in PLFA", axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enrich13C_Taxa.pdf",sep=""), plot = last_plot())
  
    #Bar plots showing distribution of Enrich 13C in PLFAs With fill W/O Fungi (i.e. by percent of total)
    graphics.off()
    p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass" & Taxonomy != "Fungi"), aes(Factor_data, umols13C_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = expression(paste("Taxonomic Distribution of Enriched Carbon in PLFA", italic("(Excluding Fungal Markers)"), sep="")), axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enrich13C_wo_Fungi_Taxa.pdf",sep=""), plot = last_plot())
  
    #Bar plots showing distribution of Enrich 13C in PLFAs With fill W/O Fungi (i.e. by percent of total)
    graphics.off()
    p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass" & Taxonomy != "Fungi"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = expression(paste("Taxonomic Distribution of Enriched Carbon in PLFA", italic("(Excluding Fungal Markers)"), sep="")), axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enrich13C_wo_Fungi_Taxa.pdf",sep=""), plot = last_plot())
    
    
    # Look at the levels of enrichment in the UF Fatty acid assignments to see if we're missing anything interesting
    writeLines("\nWe will now separately plot only \"Unidentified Fatty Acids\" (UF). We will create two graphs: i) one with all UF present and ii) one with all UF found in your control samples removed. Please enter a string which identifies your controls.\n\n Enter: NONE, if no controls present (in this case you will be prompted to chose an average delta-value threshold).\n\n[note: rationale for this step is: for any unidentified legitimate FA, it should be enriched and, therefore controlling for the presence of UFs in unenriched control is critical]")
    control_ID<-scan(n=length(unique(Factor_data)), what = character())
    
    if (control_ID != "NONE") {
      
      writeLines("\n If controls are specific to a factor, please enter the name of the factor now OR Enter: NONE: \n")
      factor_ID<-scan(n=1, what = character())
      
      empty<-paste(factor_ID)
      assign(paste(factor_ID),empty)
      UF<-x_e[grep("UF", x_e$FA),]
      UF<-ddply(UF, ~ Factor_data + FA + get(factor_ID), summarize, sum_enrich13C=sum(umols13C_per_drywt))
      colnames(UF)[3]<- paste(factor_ID)
      
      dir.create(paste(Workingdirectory,"/Output/","/",Factor,"/",sep=""))
      
      graphics.off()
      p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
      p 
      ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
      
      #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
      graphics.off()
      p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
      p  
      ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
      
      
      if (any(grepl(paste(control_ID), UF$Factor_data)) == TRUE) {
        
        removables<-UF[(grep(paste(control_ID), UF$Factor_data)),]
        removables<-ddply(removables, ~ get(factor_ID) + Factor_data, summarize, FA = FA )  
        cownter<-unique(UF[,3])
        
        for (fazz in cownter) {
          
          foo<-subset(removables, removables[,1] == paste(fazz))
          foo<-unique(foo$FA)
          
          for (ubu in foo) {
            
            UF[grep(paste("^",ubu,'$', sep=""), UF$FA),] <- NA
            
          }
        }
        
        #remove NA from UF
        UF<-UF[-which(is.na(UF)),]
        
        graphics.off()
        p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
        p 
        ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_adjusted",sep=""), plot = last_plot())
        
        #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
        graphics.off()
        p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
        p  
        ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_adjusted.pdf",sep=""), plot = last_plot())
        
    } 
      
    n_count<-ddply(x_e, ~ Factor_data + SampleName, summarize, count=1)
    colnames(n_count)[1]<-paste(Factor)
    empty2<-NA
    assign(paste(Factor), empty2)
    n_count<-ddply(n_count, ~ get(paste(Factor)), summarize, count=sum(count))
    colnames(n_count)[1]<-paste(Factor)
    write.csv(n_count, file=paste("Output\\",paste(Factor),"\\Sample Size (n) of Selected Factor.csv",sep=""))
    writeLines("\n\nSeven or eight new plots, depending on whether you had unenriched controls, have been created in the //Output// folder. Examine the title and axis to understand what they are showing.\n\nAlso, a new spreadsheet has been created showing the sample sizes across your factor")
    
    writeLines("\n\nAlso the final datasheet entitle \"Final with Taxonomy\" has been created in the \"Output\" directory.")
    write.csv(x_e, file="Output\\Final with Taxonomy.csv")
      
} else {
    
  writeLines("\n\nPlease input the delta-value cut-off for unidentified FAs. Any UF with a mean delta-value below this threshold will be removed. For SIP data substrate ammendments, a cut off of ~ 50 is suggested, while experiments with CO2 labeling or other marginal labeling concentrations are recommended between -10 and 0.")
  UF_thresh<-scan(n=1, what = numeric())
  
    UF<-x_e[grep("UF", x_e$FA),]
    
    #Remove all UF with an average delta-value < UF_thresh[1]
    FA_exclude<-ddply(UF, ~ FA, summarize, delta=mean(Delta_value))
    FA_exclude_names<-FA_exclude[which(FA_exclude$delta < UF_thresh[1]), 1]
    FA_exclude<-data.frame(Excluded_FA=FA_exclude_names, mean_delta=FA_exclude[which(FA_exclude$delta < UF_thresh[1]), 2])
    
    FA_remove<-data.frame(FA_to_Remove)
    dummy<-data.frame(rep("Found in fewer than 3 Samples", length(FA_to_Remove)))
    FA_remove<-droplevels(data.frame(FA_remove,mean_delta=dummy))
    colnames(FA_remove)[2]<-"mean_delta"
    colnames(FA_remove)[1]<-"Excluded_FA"
    FA_exclude<-rbind(FA_exclude, FA_to_Remove)
  
    write.csv(FA_exclude, file=paste("Output\\",paste(Factor),"\\FA Removed from Analysis.csv",sep=""))
    
    for (poi in FA_exclude_names) {
      
      UF<-UF[-grep(paste("^",poi,"$",sep=""), UF$FA),]
      
    }
    
    #p<-ggplot(UF, aes(FA, Delta_value)) + geom_point()
    #p
    
    UF<-ddply(UF, ~ Factor_data + FA, summarize, sum_enrich13C=sum(umols13C_per_drywt))
    
    dir.create(paste(Workingdirectory,"/Output/","/",Factor,"/",sep=""))
    
    graphics.off()
    p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
    
    #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
    graphics.off()
    p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
    
    n_count<-ddply(x_e, ~ Factor_data + SampleName, summarize, count=1)
    colnames(n_count)[1]<-paste(Factor)
    empty2<-NA
    assign(paste(Factor), empty2)
    n_count<-ddply(n_count, ~ get(paste(Factor)), summarize, count=sum(count))
    colnames(n_count)[1]<-paste(Factor)
    write.csv(n_count, file=paste("Output\\",paste(Factor),"\\Sample Size (n) of Selected Factor.csv",sep=""))
    writeLines("\n\nSeven or eight new plots, depending on whether you had unenriched controls, have been created in the //Output// folder. Examine the title and axis to understand what they are showing.\n\nAlso, a new spreadsheet has been created showing the sample sizes across your factor")
    
    FA_count<-unique(x_e$SampleName)
        
    for (yip in FA_count) {
      
      len<-length(grep(paste(yip), x_e$SampleName))
      assign(yip, len)
      
    }
    
    foox<-do.call(rbind, lapply(FA_count, as.symbol))
    foox<-data.frame(foox)
    foox$SampleName<-row.names(foox)
    foox<-foox[order(foox$SampleName),]
    foox$SampleName<-NULL
    foox<-cbind(foox,fooz)
    colnames(foox)[match("foox",colnames(foox))] <- "FAs"
    
    write.csv(foox, file="Output\\Table of FA Counts Across All Factors.csv")
    writeLines("\n\nA Table of Counts of All FAs in each Sample and corresponding factor variables has been created in the \"Output\" folder. Please consult, since repeated samples will be treated as unique samples unless removed from original batch file.")
    
    writeLines("\n\nAlso the final datasheet entitle \"Final with Taxonomy\" has been created in the \"Output\" directory.")
    write.csv(x_e, file="Output\\Final with Taxonomy.csv")
  } 
options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
} else {
    
    ##Calculate Fungal:Bacterial Ratio Modified from (Hogberg 2013: 18:2ω6,9 : (16:1ω7+cy17:0+17:0+18:1ω7+cy19 + i15:0+a15:0+15:0+i16:0+a17:0)). I included additional peaks for fungi and other peaks for bacteria (which I have never seen annotated, so are unlikely to factor)
    cat("We will now calculate the Fungi:Bacteria using a modified version of the equation used by Hogberg et al. 2013 (DOI 10.1007/s11104-013-1742-9). You may calculate the F:B ratio ratio for any one of your factors or for every sample (recommended). To calculate it for all samples, ENTER: SampleName, to use a factor, ENTER: the exact name of your factor, ensuring it matches the exact column heading you created in your \"factor_list.csv\". If you would like to do more than one factor, you'll have to re-run the script. Lastly, if your F:B ratio is returned as \"NoBacteria\" this is because one of your samples within your factor had no identified bacteria, not necessarily all of your samples.")
    FB_factor<-scan(n=1, what = character())
    
    FB<-ddply(x_e, ~ get(FB_factor) + FA + Taxonomy + Enriched, summarize, umols13C_per_drywt=umols13C_per_drywt, umols12C_per_drywt=umols12C_per_drywt, C13excess_umols_per_drywt=C13excess_umols_per_drywt)
    colnames(FB)[1]<-c("Factor")
    Fact_list<-as.vector(na.omit(unique(FB$Factor)))
    
    #F:B Ratio for only enriched samples (i.e. no controls)
    FB_enriched<-subset(FB, Enriched == 1)
    
    #Setting up the list for each combination
    FB_allC<-rep(1, length(Fact_list))
    FB_13C<-rep(1, length(Fact_list))
    FB_12C<-rep(1, length(Fact_list))
    FB_13Cexcess<-rep(1, length(Fact_list))
    
    count = 1
    for(fact in Fact_list) {
      gargle<-subset(FB_enriched, Factor == paste(fact))
      gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols13C_per_drywt)+sum(umols12C_per_drywt))
      oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
      if(any(oh) == TRUE) {
        if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
          if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
            
          } else {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

          }

        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        }
      } else {
        gargle<-c("No Bacteria")
      }
      FB_allC[[count]]<-gargle
      count = count + 1
    }
    
    count = 1
    for(fact in Fact_list) {
      gargle<-subset(FB_enriched, Factor == paste(fact))
      gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols13C_per_drywt))
      oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
      if(any(oh) == TRUE) {
        if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
          if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
            
          } else {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

          }

        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        }
      } else {
        gargle<-c("No Bacteria")
      }
      FB_13C[[count]]<-gargle
      count = count + 1
    }
    
    count = 1
    for(fact in Fact_list) {
      gargle<-subset(FB_enriched, Factor == paste(fact))
      gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols12C_per_drywt))
      oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
      if(any(oh) == TRUE) {
        if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
          if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
            
          } else {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

          }

        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        }
      } else {
        gargle<-c("No Bacteria")
      }
      FB_12C[[count]]<-gargle
      count = count + 1
    }
    
    count = 1
    for(fact in Fact_list) {
      gargle<-subset(FB_enriched, Factor == paste(fact))
      gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(C13excess_umols_per_drywt))
      oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
      if(any(oh) == TRUE) {
        if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
          if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
            
          } else {
            
            gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

          }

        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        }
      } else {
        gargle<-c("No Bacteria")
      }
      FB_13Cexcess[[count]]<-gargle
      count = count + 1
    }
    
    FB_ratio<-data.frame(cbind(Factor=Fact_list,Total_FA=FB_allC, C12_only=FB_12C, C13_only=FB_13C, Excess.13C=FB_13Cexcess))
    
    if (FB_factor == "SampleName") {
      
      #Add factors to FB ratio
      factors<-colnames(factor_list)[match("Enriched", colnames(factor_list))+1:ncol(factor_list)]
      factors<-factors[which(is.na(factors) == FALSE)]
      
      for (ip in factors) {
        
        fooz <- data.frame(x_e[,match(paste(ip), colnames(x_e))])
        assign(paste(ip),fooz)
      }
      
      
      fooz<-do.call(cbind, lapply(factors, as.symbol))
      fooz<-cbind(fooz, SampleName=x_e[,match("SampleName", colnames(x_e))])
      
      count=1
      for (count in 1:ncol(fooz)-1) {
        
        colnames(fooz)[count] <- factors[count]
        
      }
      
      #Remove duplicates
      fooz<-fooz[-which(duplicated(fooz$SampleName)),]
      
      #Ensure same order
      fooz<-fooz[order(fooz$SampleName),]
      FB_ratio<-FB_ratio[order(FB_ratio$Factor),]
      
      #Bind dataframes
      FB_ratio<-cbind(FB_ratio, fooz)
      
    }
    
    Workingdirectory<-getwd()
    dir.create(paste(Workingdirectory,"/Output/","/",FB_factor,"/",sep=""))
    write.csv(FB_ratio, file=paste("Output\\",paste(FB_factor),"\\FB_ratios.csv",sep=""))
    writeLines("\n\nA spreadsheet entitled \"FB_ratios.csv\" has been created in the \\Output\\ folder under a folder title corresponding to the factor your selected above. This may simply be \"SampleName\".\n\n")
    
    ##Now analysis will be done according to a chosen factor
    writeLines("\nWe will now create some additional graphs according to a factor you select.\n")
    cat("Please enter a factor by which you would like to separate your samples. \n\nIf none, ENTER: SampleName")
    Factor<-scan(n=1, what = character())
    
    #Create empty object with name of factor; this is necessary for calling the factor in subsequent functions like ddply
    empty<-NA
    assign(paste(Factor),empty)
    Factor_data<-x_e[,grep(paste(Factor), colnames(x_e))]
    
    #Setting the order of variables
    writeLines("\nIf there is a particular order to the variables withing your factor please input the order now. Enter each variable in a separate line according to the order you wish them placed. For example, \nfactor1\nfactor2\n...\nfactor(n). If none, please Enter: NONE. ALSO!! If any of your names contain a space, be sure to use quotation marks around the entire phrase.")
    factor_order<-scan(n=length(unique(Factor_data)), what = character())
    
    if (any(factor_order == "None") == TRUE) {
      print("No Factor Order Input")
    } else {
      x_e<-cbind(x_e, Factor=Factor_data)
      x_e$Factor_data<-factor(Factor_data, levels = factor_order)
 
    }
    ##Show crazy patterns of 13C Enrichment across your factor according to FAs
    FA_13Cenrich<- ddply(x_e, ~  Factor_data + FA, summarize, C13excess_umols_per_drywt=sum(C13excess_umols_per_drywt))
    
    graphics.off()
    p<-ggplot(FA_13Cenrich, aes(Factor_data, C13excess_umols_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enriched13C_FAs.pdf",sep=""), height=12, width=10, plot = last_plot())
    
    graphics.off()
    p<-ggplot(FA_13Cenrich, aes(Factor_data, C13excess_umols_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enriched13C_FAs.pdf",sep=""), height=12, width=10, plot = last_plot())
    
    
    #Same graph in black and white
    myColors_FA_BnW<-rev(gray.colors(length(levels(FA_13Cenrich$FA)), start = 0, end = 1, gamma = 2.2, alpha = NULL))
    names(myColors_FA_BnW) <- levels(FA_13Cenrich$FA)
    colScale_FA_BnW <- scale_fill_manual(name = "FA",values = myColors_FA_BnW)
    
    graphics.off()
    p<-ggplot(FA_13Cenrich, aes(Factor_data, C13excess_umols_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
    p + colScale_FA_BnW 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enriched13C_FAs_BnW.pdf",sep=""), height=12, width=10, plot = last_plot())
        
    #Incorrect graph, but with the legend you want.
    #graphics.off()
    #p<-ggplot(FA_13Cenrich, aes(Factor_data, C13excess_umols_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1)) + geom_bar(colour = "black", show_guide=FALSE)
    #p
    #ggsave(filename="Output\\Enriched13C_Proportional_FAs.pdf", height=12, width=10, plot = last_plot())
    
    #writeLines("\n Two general graphs have been printed in the \Output\ folder. Pleas take a look at them. These maybe useful, but are mostly to inform whether you'd like colour graphs or black and white. Please consult them now.\n")
    
    #writeLines("\nPlease select a colour scheme.\n")
    #switch(menu(c("Colourful","Black and White", "Colour-blindness Aware")), colz<-c("Colourful"), colz<-c("BnW"), colz<-c("CBA"))
    
    ## Show Breakdown of enrichment across Taxonomy in Bar charts according to a factor
    Taxo_13Cenrich<- ddply(x_e, ~  Factor_data + Taxonomy, summarize, C13excess_umols_per_drywt=sum(C13excess_umols_per_drywt))
    
    #if (colz == "Colourful") {
    #  ##Create a custom color scale for taxonomy
    #  myColors_tax <- NA
    #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
    
    #}
    
    #if (colz == "BnW") {
    #  #Black and White Scale
    #  myColors_tax<-rev(gray.colors(6, start = 0.1, end = 1, gamma = 2.2, alpha = NULL))
    #  names(myColors_tax) <- levels(Taxo_13Cenrich$Taxonomy)
    #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
    
    #}
    
    #if (colz == "CBA") {
    #  #Color-blind palette
    #  myColors_tax<-rev(c("#331A00","#996136","#D9AF98","#99F8FF","#33E4FF","#007A99"))
    #  names(myColors_tax) <- levels(Taxo_13Cenrich$Taxonomy)
    #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
    
    #}
    
    ##TEST YOUR PALETTE FOR COLOR BLINDNESS
    #myColors_tax<-dichromat(myColors_tax)
    #colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
    
    #Bar plots showing distribution of Enrich 13C in PLFAs With fill (i.e. by percent of total)
    graphics.off()
    p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Taxonomic Distribution of Enriched Carbon in PLFA", axis.title=element_text(hjust=0.1))
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enrich13C_Taxa.pdf",sep=""), plot = last_plot())
    
    #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
    graphics.off()
    p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Taxonomic Distribution of Enriched Carbon in PLFA", axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enrich13C_Taxa.pdf",sep=""), plot = last_plot())
    
    #Bar plots showing distribution of Enrich 13C in PLFAs With fill W/O Fungi (i.e. by percent of total)
    graphics.off()
    p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass" & Taxonomy != "Fungi"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = expression(paste("Taxonomic Distribution of Enriched Carbon in PLFA", italic("(Excluding Fungal Markers)"), sep="")), axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enrich13C_wo_Fungi_Taxa.pdf",sep=""), plot = last_plot())
    
    #Bar plots showing distribution of Enrich 13C in PLFAs With fill W/O Fungi (i.e. by percent of total)
    graphics.off()
    p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass" & Taxonomy != "Fungi"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = expression(paste("Taxonomic Distribution of Enriched Carbon in PLFA", italic("(Excluding Fungal Markers)"), sep="")), axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enrich13C_wo_Fungi_Taxa.pdf",sep=""), plot = last_plot())
        
    # Look at the levels of enrichment in the UF Fatty acid assignments to see if we're missing anything interesting
    writeLines("\nWe will now separately plot only \"Unidentified Fatty Acids\" (UF). We will create two graphs: i) one with all UF present and ii) one with all UF found in your control samples removed. Please enter a string which identifies your controls.\n\n Enter: NONE, if no controls present (in this case all UF with an average delta-value < UF_thresh[1] will be excluded).\n\n[note: rationale for this step is: for any unidentified legitimate FA, it should be enriched and, therefore controlling for the presence of UFs in unenriched control is critical]")
    control_ID<-scan(n=length(unique(Factor_data)), what = character())
    
    if (control_ID != "NONE") {
      
      writeLines("\n If controls are specific to a factor, please enter the name of the factor now OR Enter: NONE: \n")
      factor_ID<-scan(n=1, what = character())
      
      empty<-paste(factor_ID)
      assign(paste(factor_ID),empty)
      UF<-x_e[grep("UF", x_e$FA),]
      UF<-ddply(UF, ~ Factor_data + FA + get(factor_ID), summarize, sum_enrich13C=sum(C13excess_umols_per_drywt))
      colnames(UF)[3]<- paste(factor_ID)
      
      dir.create(paste(Workingdirectory,"/Output/","/",Factor,"/",sep=""))
      
      graphics.off()
      p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
      p 
      ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
      
      #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
      graphics.off()
      p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
      p  
      ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
      
      
      if (any(grepl(paste(control_ID), UF$Factor_data)) == TRUE) {
        
        removables<-UF[(grep(paste(control_ID), UF$Factor_data)),]
        removables<-ddply(removables, ~ get(factor_ID) + Factor_data, summarize, FA = FA )  
        cownter<-unique(UF[,3])
        
        for (fazz in cownter) {
          
          foo<-subset(removables, removables[,1] == paste(fazz))
          foo<-unique(foo$FA)
          
          for (ubu in foo) {
            
            UF[grep(paste("^",ubu,'$', sep=""), UF$FA),] <- NA
            
          }
        }
        
        #remove NA from UF
        UF<-UF[-which(is.na(UF)),]
        
        graphics.off()
        p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
        p 
        ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_adjusted.pdf",sep=""), plot = last_plot())
        
        #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
        graphics.off()
        p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
        p  
        ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_adjusted.pdf",sep=""), plot = last_plot())
        
      } 
    
    n_count<-ddply(x_e, ~ Factor_data + SampleName, summarize, count=1)
    colnames(n_count)[1]<-paste(Factor)
    empty2<-NA
    assign(paste(Factor), empty2)
    n_count<-ddply(n_count, ~ get(paste(Factor)), summarize, count=sum(count))
    colnames(n_count)[1]<-paste(Factor)
    write.csv(n_count, file=paste("Output\\",paste(Factor),"\\Sample Size (n) of Selected Factor.csv",sep=""))
    writeLines("\n\nFive or six new plots, depending on whether you had unenriched controls, have been created in the //Output// folder. Examine the title and axis to understand what they are showing.\n\nAlso, a new spreadsheet has been created showing the sample sizes across your factor")
    
    } else {
      
      UF<-x_e[grep("UF", x_e$FA),]
      UF<-ddply(UF, ~ Factor_data + FA, summarize, sum_enrich13C=sum(C13excess_umols_per_drywt))
      
      dir.create(paste(Workingdirectory,"/Output/","/",Factor,"/",sep=""))
      
      graphics.off()
      p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
      p 
      ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
      
      #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
      graphics.off()
      p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
      p  
      ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
      
      n_count<-ddply(x_e, ~ Factor_data + SampleName, summarize, count=1)
      colnames(n_count)[1]<-paste(Factor)
      empty2<-NA
      assign(paste(Factor), empty2)
      n_count<-ddply(n_count, ~ get(paste(Factor)), summarize, count=sum(count))
      colnames(n_count)[1]<-paste(Factor)
      write.csv(n_count, file=paste("Output\\",paste(Factor),"\\Sample Size (n) of Selected Factor.csv",sep=""))
      writeLines("\n\nSeven or eight new plots, depending on whether you had unenriched controls, have been created in the //Output// folder. Examine the title and axis to understand what they are showing.\n\nAlso, a new spreadsheet has been created showing the sample sizes across your factor")
      
    } 
    
    writeLines("\n\nPlease input the delta-value cut-off for unidentified FAs. Any UF with a mean delta-value below this threshold will be removed. For SIP data substrate ammendments, a cut off of ~ 50 is suggested, while experiments with CO2 labeling or other marginal labeling concentrations are recommended between -10 and 0.")
    UF_thresh<-scan(n=1, what = numeric())
    
    UF<-x_e[grep("UF", x_e$FA),]
    
    #Remove all UF with an average delta-value < UF_thresh[1]
    FA_exclude<-ddply(UF, ~ FA, summarize, delta=mean(Delta_value))
    FA_exclude_names<-FA_exclude[which(FA_exclude$delta < UF_thresh[1]), 1]
    FA_exclude<-data.frame(Excluded_FA=FA_exclude_names, mean_delta=FA_exclude[which(FA_exclude$delta < UF_thresh[1]), 2])
    
    FA_remove<-data.frame(FA_to_Remove)
    dummy<-data.frame(rep("Found in fewer than 3 Samples", length(FA_to_Remove)))
    FA_remove<-droplevels(data.frame(FA_remove,mean_delta=dummy))
    colnames(FA_remove)[2]<-"mean_delta"
    colnames(FA_remove)[1]<-"Excluded_FA"
    FA_exclude<-rbind(FA_exclude, FA_to_Remove)
    
    write.csv(FA_exclude, file=paste("Output\\",paste(Factor),"\\FA Removed from Analysis.csv",sep=""))
    
    for (poi in FA_exclude_names) {
      
      UF<-UF[-grep(paste("^",poi,"$",sep=""), UF$FA),]
      
    }
    
    #p<-ggplot(UF, aes(FA, Delta_value)) + geom_point()
    #p
    
    FA_count<-unique(x_e$SampleName)
    
    for (yip in FA_count) {
      
      len<-length(grep(paste(yip), x_e$SampleName))
      assign(yip, len)
      
    }
    
    foox<-do.call(rbind, lapply(FA_count, as.symbol))
    foox<-data.frame(foox)
    foox$SampleName<-row.names(foox)
    foox<-foox[order(foox$SampleName),]
    foox$SampleName<-NULL
    foox<-cbind(foox,fooz)
    colnames(foox)[match("foox",colnames(foox))] <- "FAs"
    
    write.csv(foox, file="Output\\Table of FA Counts Across All Factors.csv")
    writeLines("\n\nA Table of Counts of All FAs in each Sample and corresponding factor variables has been created in the \"Output\" folder. Please consult, since repeated samples will be treated as unique samples unless removed from original batch file.")
    
    writeLines("\n\nAlso the final datasheet entitle \"Final with Taxonomy\" has been created in the \"Output\" directory.")
    write.csv(x_e, file="Output\\Final with Taxonomy.csv")
    
}
options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
  } else {
    
  ##Calculate Fungal:Bacterial Ratio Modified from (Hogberg 2013: 18:2ω6,9 : (16:1ω7+cy17:0+17:0+18:1ω7+cy19 + i15:0+a15:0+15:0+i16:0+a17:0)). I included additional peaks for fungi and other peaks for bacteria (which I have never seen annotated, so are unlikely to factor)
  cat("We will now calculate the Fungi:Bacteria using a modified version of the equation used by Hogberg et al. 2013 (DOI 10.1007/s11104-013-1742-9). You may calculate the F:B ratio ratio for any one of your factors or for every sample (recommended). To calculate it for all samples, ENTER: SampleName, to use a factor, ENTER: the exact name of your factor, ensuring it matches the exact column heading you created in your \"factor_list.csv\". If you would like to do more than one factor, you'll have to re-run the script. Lastly, if your F:B ratio is returned as \"NoBacteria\" this is because one of your samples within your factor had no identified bacteria, not necessarily all of your samples.")
  FB_factor<-scan(n=1, what = character())
  
  FB<-ddply(x_e, ~ get(FB_factor) + FA + Taxonomy + Enriched, summarize, umols13C_per_drywt=umols13C_per_drywt, umols12C_per_drywt=umols12C_per_drywt, C13excess_umols_per_drywt=C13excess_umols_per_drywt)
  colnames(FB)[1]<-c("Factor")
  Fact_list<-as.vector(na.omit(unique(FB$Factor)))
  
  #13C enriched samples F:B Ratio
  FB_enriched<-subset(FB, Enriched == 1)
  
  FB_allC<-rep(1, length(Fact_list))
  FB_13C<-rep(1, length(Fact_list))
  FB_12C<-rep(1, length(Fact_list))
  FB_13Cexcess<-rep(1, length(Fact_list))
  
  count = 1
  for(fact in Fact_list) {
    gargle<-subset(FB_enriched, Factor == paste(fact))
    gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols13C_per_drywt)+sum(umols12C_per_drywt))
    oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
    if(any(oh) == TRUE) {
      if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
        if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])
          
        }
        
      } else {
        
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
                
      }
          } else {
      gargle<-c("No Bacteria")
    }
    FB_allC[[count]]<-gargle
    count = count + 1
  }
  
  count = 1
  for(fact in Fact_list) {
    gargle<-subset(FB_enriched, Factor == paste(fact))
    gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols13C_per_drywt))
    oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
    if(any(oh) == TRUE) {
      if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
        if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

        }

      } else {
        
        gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
        
      }
    } else {
      gargle<-c("No Bacteria")
    }
    FB_13C[[count]]<-gargle
    count = count + 1
  }
  
  count = 1
  for(fact in Fact_list) {
    gargle<-subset(FB_enriched, Factor == paste(fact))
    gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(umols12C_per_drywt))
    oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
    if(any(oh) == TRUE) {
      if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
        if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

        }

      } else {
        
        gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
        
      }
    } else {
      gargle<-c("No Bacteria")
    }
    FB_12C[[count]]<-gargle
    count = count + 1
  }
  
  count = 1
  for(fact in Fact_list) {
    gargle<-subset(FB_enriched, Factor == paste(fact))
    gargle<-ddply(gargle, ~Taxonomy, summarize, umols=sum(C13excess_umols_per_drywt))
    oh<-grepl("Gram-negative Bacteria|Gram-positive Bacteria", gargle$Taxonomy)
    if(any(oh) == TRUE) {
      if (any(grepl("Gram-negative Bacteria", gargle$Taxonomy)) == TRUE) {
        if (any(grepl("Gram-positive Bacteria", gargle$Taxonomy)) == TRUE) {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2]+gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
          
        } else {
          
          gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-negative Bacteria", gargle$Taxonomy), 2])

        }

      } else {
        
        gargle<-gargle[grep("Fungi", gargle$Taxonomy), 2]/(gargle[grep("Gram-positive Bacteria", gargle$Taxonomy), 2])
        
      }
    } else {
      gargle<-c("No Bacteria")
    }
    FB_13Cexcess[[count]]<-gargle
    count = count + 1
  }
  
  FB_ratio<-data.frame(cbind(Factor=Fact_list,Total_FA=FB_allC, C12_only=FB_12C, C13_only=FB_13C, Excess.13C=FB_13Cexcess))
  
  if (FB_factor == "SampleName") {
    
    #Add factors to FB ratio
    factors<-colnames(factor_list)[match("Enriched", colnames(factor_list))+1:ncol(factor_list)]
    factors<-factors[which(is.na(factors) == FALSE)]
    
    for (ip in factors) {
      
      fooz <- data.frame(x_e[,match(paste(ip), colnames(x_e))])
      assign(paste(ip),fooz)
    }
    
    
    fooz<-do.call(cbind, lapply(factors, as.symbol))
    fooz<-cbind(fooz, SampleName=x_e[,match("SampleName", colnames(x_e))])
    
    count=1
    for (count in 1:ncol(fooz)-1) {
      
      colnames(fooz)[count] <- factors[count]
      
    }
    
    #Remove duplicates
    fooz<-fooz[-which(duplicated(fooz$SampleName)),]
    
    #Ensure same order
    fooz<-fooz[order(fooz$SampleName),]
    FB_ratio<-FB_ratio[order(FB_ratio$Factor),]
    
    #Bind dataframes
    FB_ratio<-cbind(FB_ratio, fooz)
    
  }
    
  Workingdirectory<-getwd()
  dir.create(paste(Workingdirectory,"/Output/","/",FB_factor,"/",sep=""))
  write.csv(FB_ratio, file=paste("Output\\",paste(FB_factor),"\\FB_ratios.csv",sep=""))
  writeLines("\n\nA spreadsheet entitled \"FB_ratios.csv\" has been created in the \\Output\\ folder under a folder title corresponding to the factor your selected above. This may simply be \"SampleName\".\n\n")
  
  ##Now analysis will be done according to a chosen factor
  writeLines("\nWe will now create some additional graphs according to a factor you select.\n")
  cat("Please enter a factor by which you would like to separate your samples. \n\nIf none, ENTER: SampleName")
  Factor<-scan(n=1, what = character())
  
  #Create empty object with name of factor; this is necessary for calling the factor in subsequent functions like ddply
  empty<-NA
  assign(paste(Factor),empty)
  Factor_data<-x_e[,grep(paste(Factor), colnames(x_e))]
  
  #Setting the order of variables
  writeLines("\nIf there is a particular order to the variables withing your factor please input the order now. Enter each variable in a separate line according to the order you wish them placed. For example, \nfactor1\nfactor2\n...\nfactor(n). If none, please Enter: NONE. ALSO!! If any of your names contain a space, be sure to use quotation marks around the entire phrase.")
  factor_order<-scan(n=length(unique(Factor_data)), what = character())
  
  if (any(factor_order == "None") == TRUE) {
    print("No Factor Order Input")
  } else {
    x_e<-cbind(x_e, Factor=Factor_data)
    x_e$Factor_data<-factor(Factor_data, levels = factor_order)
  }
  
  ##Identify whether the groups are equal or not AND explain that unequal sample sizes will be averaged (retaining any unique FA in any of the replicates) for all plots involving comparisons of totals
  n_count<-ddply(x_e, ~ Factor_data + SampleName, summarize, count=1)
  colnames(n_count)[1]<-paste(Factor)
  empty2<-NA
  assign(paste(Factor), empty2)
  n_count<-ddply(n_count, ~ get(paste(Factor)), summarize, count=sum(count))
  colnames(n_count)[1]<-paste(Factor)
  write.csv(n_count, file=paste("Output\\",paste(Factor),"\\Sample Size (n) of Selected Factor (pre-normalization).csv",sep=""))
  
  min_count<-min(n_count[,2])
  equality_test<- min_count == n_count[,2]
  
    if (any(equality_test == FALSE)) {
    
  writeLines("\nIt appears as though the factor you selected has unequal group sizes. This will lead to unbalanced comparisons in many of the graphs which will now be plotted. This script will continue by averaging any replicates in excess of the minimum equal number, BUT you may make adjustments to the data yourself in the \"Final datasheet.csv\" and re-run this script.\n\n[Note: this averaging method will only be used for graphs depicting Total biomass, not proportional biomass.]")
  switch(menu(c("Okay")), dummy<-c("okay"))
      
  group_ex<-n_count[which(n_count[,2] > min_count), 1]  
  
  for (ex in group_ex) {
  
    inexcess<-n_count[grep(ex, n_count[,1]), 2] - min_count
    seq<-seq(from = min_count, to = min_count+inexcess, by = 1)
    master_ID<-unique(x_e[grep(ex, x_e[,match(paste(Factor), colnames(x_e))]), match("SampleName", colnames(x_e))])
    master_ID<-master_ID[seq]
    
      for (inex in master_ID) {
        
      aloha<-subset(x_e, x_e$SampleName == inex)
      assign(paste(inex),aloha)
      
      }
    
    master_compiled<-do.call(rbind, lapply(master_ID, as.symbol))
    amalgam<-ddply(master_compiled, ~FA, summarize, SampleName=master_ID[1], Delta_value=mean(Delta_value), Peak_Area=mean(Peak_Area), molC_per_uL=mean(molC_per_uL), molC=mean(molC), umols13C_per_drywt=mean(umols13C_per_drywt), umols12C_per_drywt=mean(umols12C_per_drywt), C13excess_umols_per_drywt=mean(C13excess_umols_per_drywt))
    
    #For some reaons ddply wont summarize character vectors... so getting errors when taxonomy overlap
    master_compiled<-master_compiled[order(master_compiled$FA),]
    master_compiled<-master_compiled[-which(duplicated(master_compiled$FA)),]
    amalgam<-cbind(amalgam, master_compiled$Taxonomy)    
    colnames(amalgam)[ncol(amalgam)]<-"Taxonomy"
    
    #Add all column from one of samples to new hybrid
    new_ID<-master_ID[1]
    new_ID<-grep(new_ID, x_e$SampleName)
    new_ID<-x_e[new_ID[1],]
    
    #keep the ID of the one used for future use cbind(fooz, foox) down below
    safe_master<-master_ID[-1]
    
    #Remove all traces of replicates from main dataframe and fooz
    for (eliminate in master_ID) {
      
      x_e<-x_e[-grep(eliminate, x_e$SampleName),]    
          
    }
    
    for (eliminate in safe_master) {
      
      fooz<-fooz[-grep(eliminate, fooz$SampleName),]   
      
    }
        
    #Introduce new amalgamated averaged sample
    d <- data.frame(x = new_ID)
    colnames(d) <- colnames(x_e)
    n <- nrow(amalgam)
    d<-do.call("rbind", replicate(n, d, simplify = FALSE))
    d$FA<-amalgam$FA    
    d$umols13C_per_drywt<-amalgam$umols13C_per_drywt
    d$umols12C_per_drywt<-amalgam$umols12C_per_drywt
    d$C13excess_umols_per_drywt<-amalgam$C13excess_umols_per_drywt
    d$Taxonomy<-amalgam$Taxonomy
    d$Delta_value<-amalgam$Delta_value
    d$Peak_Area<-amalgam$Peak_Area
    d$molC_per_uL<-amalgam$molC_per_uL
    d$molC<-amalgam$molC
    
    x_e<-rbind(x_e,d)
    
  }
    
}
    
  ##Show crazy patterns of 13C Enrichment across your factor according to FAs
  FA_13Cenrich<- ddply(x_e, ~  Factor_data + FA, summarize, C13excess_umols_per_drywt=sum(C13excess_umols_per_drywt))
  
  graphics.off()
  p<-ggplot(FA_13Cenrich, aes(Factor_data, C13excess_umols_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
  p 
  ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enriched13C_FAs.pdf",sep=""), height=12, width=10, plot = last_plot())
  
  graphics.off()
  p<-ggplot(FA_13Cenrich, aes(Factor_data, C13excess_umols_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
  p 
  ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enriched13C_FAs.pdf",sep=""), height=12, width=10, plot = last_plot())
  
  #Same graph in black and white
  myColors_FA_BnW<-rev(gray.colors(length(levels(FA_13Cenrich$FA)), start = 0, end = 1, gamma = 2.2, alpha = NULL))
  names(myColors_FA_BnW) <- levels(FA_13Cenrich$FA)
  colScale_FA_BnW <- scale_fill_manual(name = "FA",values = myColors_FA_BnW)
  
  graphics.off()
  p<-ggplot(FA_13Cenrich, aes(Factor_data, C13excess_umols_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1))
  p + colScale_FA_BnW 
  ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enriched13C_FAs_BnW.pdf",sep=""), height=12, width=10, plot = last_plot())
    
  #Incorrect graph, but with the legend you want.
  #graphics.off()
  #p<-ggplot(FA_13Cenrich, aes(Factor_data, C13excess_umols_per_drywt, fill = FA, order = desc(FA))) + geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Proportions of Enriched FAs Along a Gradient According To FA Chain Length", axis.title=element_text(hjust=0.1)) + geom_bar(colour = "black", show_guide=FALSE)
  #p
  #ggsave(filename="Output\\Enriched13C_Proportional_FAs.pdf", height=12, width=10, plot = last_plot())
  
  #writeLines("\n Two general graphs have been printed in the \Output\ folder. Pleas take a look at them. These maybe useful, but are mostly to inform whether you'd like colour graphs or black and white. Please consult them now.\n")
  
  #writeLines("\nPlease select a colour scheme.\n")
  #switch(menu(c("Colourful","Black and White", "Colour-blindness Aware")), colz<-c("Colourful"), colz<-c("BnW"), colz<-c("CBA"))
  
  ## Show Breakdown of enrichment across Taxonomy in Bar charts according to a factor
  Taxo_13Cenrich<- ddply(x_e, ~  Factor_data + Taxonomy, summarize, C13excess_umols_per_drywt=sum(C13excess_umols_per_drywt))
  
  #if (colz == "Colourful") {
  #  ##Create a custom color scale for taxonomy
  #  myColors_tax <- NA
  #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
  
  #}
  
  #if (colz == "BnW") {
  #  #Black and White Scale
  #  myColors_tax<-rev(gray.colors(6, start = 0.1, end = 1, gamma = 2.2, alpha = NULL))
  #  names(myColors_tax) <- levels(Taxo_13Cenrich$Taxonomy)
  #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
  
  #}
  
  #if (colz == "CBA") {
  #  #Color-blind palette
  #  myColors_tax<-rev(c("#331A00","#996136","#D9AF98","#99F8FF","#33E4FF","#007A99"))
  #  names(myColors_tax) <- levels(Taxo_13Cenrich$Taxonomy)
  #  colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
  
  #}
  
  ##TEST YOUR PALETTE FOR COLOR BLINDNESS
  #myColors_tax<-dichromat(myColors_tax)
  #colScale_tax <- scale_fill_manual(name = "Taxonomy",values = myColors_tax)
  
  #Bar plots showing distribution of Enrich 13C in PLFAs With fill (i.e. by percent of total)
  graphics.off()
  p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Taxonomic Distribution of Enriched Carbon in PLFA", axis.title=element_text(hjust=0.1))
  p 
  ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enrich13C_Taxa.pdf",sep=""), plot = last_plot())
  
  #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
  graphics.off()
  p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Taxonomic Distribution of Enriched Carbon in PLFA", axis.title=element_text(hjust=0.1)) 
  p  
  ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enrich13C_Taxa.pdf",sep=""), plot = last_plot())
  
  #Bar plots showing distribution of Enrich 13C in PLFAs With fill W/O Fungi (i.e. by percent of total)
  graphics.off()
  p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass" & Taxonomy != "Fungi"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = expression(paste("Taxonomic Distribution of Enriched Carbon in PLFA", italic("(Excluding Fungal Markers)"), sep="")), axis.title=element_text(hjust=0.1)) 
  p  
  ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_Enrich13C_wo_Fungi_Taxa.pdf",sep=""), plot = last_plot())
  
  #Bar plots showing distribution of Enrich 13C in PLFAs With fill W/O Fungi (i.e. by percent of total)
  graphics.off()
  p<-ggplot(subset(Taxo_13Cenrich, Taxonomy != "Total Biomass" & Taxonomy != "Fungi"), aes(Factor_data, C13excess_umols_per_drywt, fill = Taxonomy, order = desc(Taxonomy))) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = expression(paste("Taxonomic Distribution of Enriched Carbon in PLFA", italic("(Excluding Fungal Markers)"), sep="")), axis.title=element_text(hjust=0.1)) 
  p  
  ggsave(filename=paste("Output\\",paste(Factor),"\\Total_Enrich13C_wo_Fungi_Taxa.pdf",sep=""), plot = last_plot())
    
  # Look at the levels of enrichment in the UF Fatty acid assignments to see if we're missing anything interesting
  writeLines("\nWe will now separately plot only \"Unidentified Fatty Acids\" (UF). We will create two graphs: i) one with all UF present and ii) one with all UF found in your control samples removed. Please enter a string which identifies your controls.\n\n Enter: NONE, if no controls present (in this case you will be prompted to chose an average delta-value threshold).\n\n[note: rationale for this step is: for any unidentified legitimate FA, it should be enriched and, therefore controlling for the presence of UFs in unenriched control is critical]")
  control_ID<-scan(n=length(unique(Factor_data)), what = character())
  
  if (control_ID != "NONE") {
    
    writeLines("\n If controls are specific to a factor, please enter the name of the factor now OR Enter: NONE: \n")
    factor_ID<-scan(n=1, what = character())
    
    empty<-paste(factor_ID)
    assign(paste(factor_ID),empty)
    UF<-x_e[grep("UF", x_e$FA),]
    UF<-ddply(UF, ~ Factor_data + FA + get(factor_ID), summarize, sum_enrich13C=sum(C13excess_umols_per_drywt))
    colnames(UF)[3]<- paste(factor_ID)
    
    dir.create(paste(Workingdirectory,"/Output/","/",Factor,"/",sep=""))
    
    graphics.off()
    p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
    
    #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
    graphics.off()
    p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
    
    
    if (any(grepl(paste(control_ID), UF$Factor_data)) == TRUE) {
      
      removables<-UF[(grep(paste(control_ID), UF$Factor_data)),]
      removables<-ddply(removables, ~ get(factor_ID) + Factor_data, summarize, FA = FA )  
      cownter<-unique(UF[,3])
      
      for (fazz in cownter) {
        
        foo<-subset(removables, removables[,1] == paste(fazz))
        foo<-unique(foo$FA)
        
        for (ubu in foo) {
          
          UF[grep(paste("^",ubu,'$', sep=""), UF$FA),] <- NA
          
        }
      }
      
      #remove NA from UF
      UF<-UF[-which(is.na(UF)),]
      
      graphics.off()
      p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
      p 
      ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_adjusted.pdf",sep=""), plot = last_plot())
      
      #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
      graphics.off()
      p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
      p  
      ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_adjusted.pdf",sep=""), plot = last_plot())
      
    } 
  
  writeLines("\n\nFive or six new plots, depending on whether you had unenriched controls, have been created in the //Output// folder. Examine the title and axis to understand what they are showing.\n\nAlso, a new spreadsheet has been created showing the sample sizes across your factor")
    
  } else {
    
    writeLines("\n\nPlease input the delta-value cut-off for unidentified FAs. Any UF with a mean delta-value below this threshold will be removed. For SIP data substrate ammendments, a cut off of ~ 50 is suggested, while experiments with CO2 labeling or other marginal labeling concentrations are recommended between -10 and 0.")
    UF_thresh<-scan(n=1, what = numeric())
    
    UF<-x_e[grep("UF", x_e$FA),]
    
    #Remove all UF with an average delta-value < 0
    FA_exclude<-ddply(UF, ~ FA, summarize, delta=mean(Delta_value))
    FA_exclude_names<-FA_exclude[which(FA_exclude$delta < UF_thresh[1]), 1]
    FA_exclude<-data.frame(Excluded_FA=FA_exclude_names, mean_delta=FA_exclude[which(FA_exclude$delta < UF_thresh[1]), 2])
    
    FA_remove<-data.frame(FA_to_Remove)
    dummy<-data.frame(rep("Found in fewer than 3 Samples", length(FA_to_Remove)))
    FA_remove<-droplevels(data.frame(Excluded_FA=FA_remove,mean_delta=dummy))
    colnames(FA_remove)[2]<-"mean_delta"
    colnames(FA_remove)[1]<-"Excluded_FA"
    FA_exclude<-rbind(FA_exclude, FA_remove)
    
    write.csv(FA_exclude, file=paste("Output\\",paste(Factor),"\\FA Removed from Analysis.csv",sep=""))
    
    for (poi in FA_exclude_names) {
      
      UF<-UF[-grep(paste("^",poi,"$",sep=""), UF$FA),]
      
    }
    
    p<-ggplot(UF, aes(FA, Delta_value)) + geom_point()
    p
    
    FA_count<-unique(x_e$SampleName)
    
    for (yip in FA_count) {
      
      len<-length(grep(paste(yip), x_e$SampleName))
      assign(yip, len)
      
    }
    
    UF<-ddply(UF, ~ Factor_data + FA, summarize, sum_enrich13C=sum(C13excess_umols_per_drywt))
    
    dir.create(paste(Workingdirectory,"/Output/","/",Factor,"/",sep=""))
    
    graphics.off()
    p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="fill", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab("Percent of Total 13C Enrichment in PLFA") + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
    p 
    ggsave(filename=paste("Output\\",paste(Factor),"\\Proportion_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
    
    #Bar plots showing distribution of Enrich 13C in PLFAs With stacking (i.e. not percent, but total)
    graphics.off()
    p<-ggplot(UF, aes(Factor_data, sum_enrich13C, fill = FA)) + geom_bar(position="stack", stat="identity",colour = "black") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y=element_text(vjust=0.1)) + xlab("") + ylab(expression(paste("Total Excess 13C Enrichment", "(µmols • ", g^-1, "dry wt sample)", sep=""))) + labs(title = "Distribution of Unidentified Enriched PLFAs", axis.title=element_text(hjust=0.1)) 
    p  
    ggsave(filename=paste("Output\\",paste(Factor),"\\Total_enrich13C_UFA_unadjusted.pdf",sep=""), plot = last_plot())
    
    n_count<-ddply(x_e, ~ Factor_data + SampleName, summarize, count=1)
    colnames(n_count)[1]<-paste(Factor)
    empty2<-NA
    assign(paste(Factor), empty2)
    n_count<-ddply(n_count, ~ get(paste(Factor)), summarize, count=sum(count))
    colnames(n_count)[1]<-paste(Factor)
    write.csv(n_count, file=paste("Output\\",paste(Factor),"\\Sample Size (n) of Selected Factor.csv",sep=""))
    writeLines("\n\nSeven or eight new plots, depending on whether you had unenriched controls, have been created in the //Output// folder. Examine the title and axis to understand what they are showing.\n\nAlso, a new spreadsheet has been created showing the sample sizes across your factor")
    
  }
  
  foox<-do.call(rbind, lapply(FA_count, as.symbol))
  foox<-data.frame(foox)
  foox$SampleName<-row.names(foox)
  foox<-foox[order(foox$SampleName),]
  foox$SampleName<-NULL
  foox<-cbind(foox,fooz)
  colnames(foox)[match("foox",colnames(foox))] <- "FAs"
  
  write.csv(foox, file="Output\\Table of FA Counts Across All Factors.csv")
  writeLines("\n\nA Table of Counts of All FAs in each Sample and corresponding factor variables has been created in the \"Output\" folder. Please consult, since repeated samples will be treated as unique samples unless removed from original batch file.")
  
  writeLines("\n\nAlso the final datasheet entitle \"Final with Taxonomy\" has been created in the \"Output\" directory.")
  write.csv(x_e, file="Output\\Final with Taxonomy.csv")
  
  }

writeLines("\n\nIt is recommended to copy and paste the history of your R session, in order to have the exact details of each session.")

options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
options(warn=-1)
