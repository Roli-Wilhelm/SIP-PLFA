SIP-PLFA
========

R script which cleans and processes fatty acid methylester (FAME) data from the Stable Isotope Facility at the University of British Columbia (http://isotopes.forestry.ubc.ca/)

This script is geared towards processing data from stable isotope probing (C13) studies, but can be used (up to a point) for processing unenriched, natural abundance data. It takes raw data and (if followed to the furtherst step) outputs ready-made graphs, Fungal:Bacterial ratios and other goodies which help analyze teh data. It maybe a bit demanding on the front-end (i.e. many files, a number of manual interventions), but the results are worth it.

IMPORTANT: This script is highly interactive and is most effective if you have included i) an internal standard, ii) quantitative standards, iii) isotopic standards AND you have labeled them according to the sample submission form included in this repository. Also, please note that this script makes use of NO sophisticated peak identification algorithms, instead preparing the data in a form which can be manually curated. To that end, the peak identification file (included here0 may become out of date due to changes in the column and gradual shift in retention time). Please see the manager of the facility for an updated peak identification file if this one proves unsatisfactory.

- M.O. -
The pipeline will extract your data from multiple files, which can contain multiple projects/experiments, aid you in identifying peaks with internal standards; calculate excess 13C (i.e. 13C above background levels); and help visualize the distribution of enriched PLFAs (i.e. your enriched community). 

- Critical Information - 
This script is not fool proofed. The potential for weirdness will be vastly minimized if the user follows the naming outlined in the *Sample_Submission.xls* (located in /Accessory_Files/); if the user ensures the raw data is the FULL data (not the "report") and is saved as a *.csv* (NOTE: The .csv format outputted by the version of Excel connected to the instrument is faulty (so try saving again as .csv if you get errors early on in the pipeline). The *Sample_Submission.xls* is the current (Oct. 2013) folder used for submiting samples. It will definitely help to have a basic understanding of the R language. 


Usage:
========
- Directions for Starting the Process -
Input fils must be in the form of a ".csv" and stored in the /Input Folder/. 

Necessary files:

names.csv		modify only if necessary
peakID.csv		modify only if necessary 
taxonomy.csv		modify only if necessary
std_conc.csv		update every time	contains details of the concentrations of standard used
standards.csv		update every time	contains details on which dilutions of standard were used


Created files:

confirm_peaks.csv	List of peaks requiring manual assignment based on rt and BAME, FAME standards + peak assignment data
duplicates_FAs.csv	List of peaks following manual assignment showing whether multiple peaks have been assigned to the same sample
factor_list.csv		List of treatment factors or other designations corresponding to sample names
Final_datasheet.csv	Dataset following cleaning and quantitating
FB_ratio.csv		Fungal:Bacerial ratio 

(additionally: standard curves, histograms etc.) 


Taxonomic Assignment:
========
All taxonomic assignments of FAMEs were based on multiple research articles. The best synthesis is from Ruess 2010 (http://www.sciencedirect.com/science/article/pii/S0038071710002683) . MOST IMPORTANTLY, the assignment is biased for SIP-decomopsition type studies. This enables the assignment of long-chain FAs to fungi as it is unlikely that plants are uptaking carbon via decomposition. The calculation of Fungi:bacteria ratio was done according to Hogberg 2013 (http://link.springer.com/article/10.1007%2Fs11104-013-1742-9).


Caveats and things to read when frustrated
========
If you incorrectly enter information at some point in the pipeline you will not only have to re-run that specific "Block", but also the preceeding one (b/c some of the data may have been erased or renamed). If things appear broken, try re-running from the very beginning again before e-mailing : roliwilhelm@gmail.com.
 
Decided to use manual assignment of peaks b/c i) the ability to identify potential peaks not found in the standard peak ID list based on the fact they are enriched (only useful for SIP studies) and ii) the extent of variability in rt can lead to false-negatives (sadly too common) and false-positives (rare). If you don't believe it, try running your data through an automated pipeline and compare it to one manually curated (the effort, and diligence is worth it for understanding your data). 

Further, I was able to run my samples on a second GC-MS with a spectra library (GC-MS in the Mohn Lab) to confirm the identify of extra peaks BUT, for SIP, it is pretty clear that if you're seeing enrichment of C13 you've got PLFAs.

The way in which quantitative standards have been used yields the umols12C and umols13C, i.e. to the umols of CARBON present in each FA peak. In other words, this has not been normalized to the # of carbon in each FA and is not to be confused with the umols of a specific FA. 

C13excess is calculated by averaging the delta value for each FA present in controls, and using this value to calculate the "natural abundance" umols13C present in the enriched sample. This value is substracted from the total umols13C, yielding the "C13excess." Each FA is paired with its control to account for differences in natural abundance across FAs. However, if a FA that is present in an enriched sample is not present in any control sample, the delta value is averaged across all FAs.  


Notes:

Did not include step for standardizing delta-values to standards. This was rarely applied during my analysis b/c the delta-values were much higher than natural abundances and outside of the range of the standards. If you wish to do this, adjust the delta_value at the same time you manually assign your peaks.

