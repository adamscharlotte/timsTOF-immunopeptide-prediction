#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#-----------------------------------------------------

#Labeling efficiency
  library(readr)
  library(dplyr)
  library(stringr)
  
  #load data
  
    o <- args[1] #"Z:\\internal_projects/active/TOPAS/WP31/Searches/JR/Biology_pilot/Batch2_var/combined/txt"
  o <- paste0(o, "/combined/txt")
  setwd(o)
  evid_var <- read.delim("evidence.txt", stringsAsFactors = FALSE, sep="\t")
  #mod_pep_var <- read.delim("modificationSpecificPeptides.txt", stringsAsFactors = FALSE, sep="\t")
  
  
  #filter data
  #remove contaminants and decoys
  evid_var <- filter(evid_var, evid_var$Reverse!="+"&evid_var$Potential.contaminant!="+") 
  
        # #Counting AS IDs
        # Cys <- sum(str_count(evid_var$Sequence, "C"))
        # Lys <- sum(str_count(evid_var$Sequence, "K"))
        # Met <- sum(str_count(evid_var$Sequence, "M"))
        # TMT_Nterm_free <- count(filter(evid_var, evid_var$Acetyl..Protein.N.term.!="1"))
        # TMT_Nterm_label <- sum(evid_var$TMT10plex.Nter)
        # 
        # #Calculate Labeling efficiencies
        # TMT_Nterm <- as.numeric(TMT_Nterm_label/TMT_Nterm_free)
        # TMT_Lys <- as.numeric(sum(evid_var$TMT10plex.Lys)/Lys)
        # Carba_Cys <- as.numeric(sum(evid_var$Carbamidomethyl..C.)/Cys)
        # Ox_Met <- as.numeric(sum(evid_var$Oxidation..M.)/Met)
 
  #Counting AS Intensities
  evid_var$Cys = str_count(evid_var$Sequence, "C")
  Cys <- sum(evid_var[which(evid_var$Cys!='0'),"Intensity"], na.rm=TRUE)
  
  evid_var$Lys = str_count(evid_var$Sequence, "K")
  Lys <- sum(evid_var[which(evid_var$Lys!='0'),"Intensity"], na.rm=TRUE)
  
  evid_var$Met = str_count(evid_var$Sequence, "M")
  Met <- sum(evid_var[which(evid_var$Met!='0'),"Intensity"], na.rm=TRUE)
  
  TMT_Nterm_free <- sum(evid_var[which(evid_var$Acetyl..Protein.N.term.!='1'),"Intensity"], na.rm=TRUE)
  TMT_Nterm_label <- sum(evid_var[which(evid_var$TMT10plex.Nter!='0'),"Intensity"], na.rm=TRUE)
  
  #Calculate Labeling efficiencies
  TMT_Nterm <- as.numeric(TMT_Nterm_label/TMT_Nterm_free)
  TMT_Lys <- as.numeric(sum(evid_var[which(evid_var$TMT10plex.Lys!="0"), "Intensity"], na.rm=TRUE)/Lys)
  Carba_Cys <- as.numeric(sum(evid_var[which(evid_var$Carbamidomethyl..C.!="0"), "Intensity"], na.rm=TRUE)/Cys)
  Ox_Met <- as.numeric(sum(evid_var[which(evid_var$Oxidation..M.!="0"), "Intensity"], na.rm=TRUE)/Met)

  #ploting Labeling Efficiencies
  setwd(o)
  png("labeling_efficiency.png")
  barplot(c(TMT_Lys, TMT_Nterm, Carba_Cys, Ox_Met), ylab="labeling efficiency",
          names.arg = c("TMT_Lys", "TMT_Nterm", "Carba_Cys", "Ox_Met"), main = c(evid_var[1,"Experiment"]) )
  dev.off()
  
#----------------------------------------------------- 
