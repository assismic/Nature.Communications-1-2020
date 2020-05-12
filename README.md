# Nature.Communications-1-2020
Source code of data analysis for Nature Comm manuscript 05/2020

requirements:

  tidyverse 1.3.0;
  
  matrixStats 0.52.0;
  
  preprocessCore (Bioconductor 3.11);
  
  sjmisc 2.8.4;
  
  caret 6-0.86;
  
  ROCR 1.0-11;
  
  pROC 1.16.2;
  
  ggpubr 0.3.0;
  
  ggplot2 3.3.0;
  

The analysis is divided into two scripts:
  1. MakeData.NatComC.R is designed to normalize the input data and generate all data frames necessary for analysis
  2. MainAnalysis_NatCom.R is the main analysis script.
  
Analysis folder must follow the following structure:

ROOT /   
....... InputFolder/

....... Output/

................ rocDFs/

....... Figures/ 

................ MultiAUC/

................ Antigens/

....... Source/

....... SavedSessions/

  
All input files must be copied into the 'InputFolder' folder  
Both script files must be copied into the 'Source' folder

1. Make data:
    setwd() #line 12 must be set to the root folder of analysis.      
    As input, MakeData.NatComC.R takes a single input for the Raw data, a file containing the meta data (stored as NatCom.pheno.csv), a annotation file with descriptions of the features (ann.csv) for full array (all replicates), an annotation file (ann.targets.avr.csv) with the annotations with only unique Antigen names for the averaged data frames and a file with the the preferred order to arrange the features on the data frames (Order.csv).
    
2. MainAnalysis_NatCom:
  setwd() #line 1 must be set to the root folder of analysis. 
  First the formated data created on the "Make Data" is loaded from a saved session (if not ran on the same session).
  top antigens (stored variable top.antigens.vec) were selected manualy.  
  
  
