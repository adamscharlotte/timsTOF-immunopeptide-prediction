#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(readr)
library(dplyr)
library(stringr)

#-----------------------------------------------------

#load data

m <- args[1] #"Z:\\internal_projects/active/TOPAS/WP31/Searches/JR/Biology_pilot/Batch3_FP/combined/txt"

m <- paste0(m, "/combined/txt")
cat(m, "\n")
setwd(m)
fp_evid <- read.delim("evidence.txt", stringsAsFactors = FALSE, sep="\t")
fp_prot <- read.delim("proteinGroups.txt", stringsAsFactors = FALSE, sep="\t")
fp_sum <- read.delim("summary.txt", stringsAsFactors = FALSE, sep="\t", nrows = 48)
fp_sum2 <- read.delim("summary.txt", stringsAsFactors = FALSE, sep="\t")

#-----------------------------------------------------

#filter data
 #remove contaminants and decoys
fp_evid_filt <- filter(fp_evid, fp_evid$Reverse!="+"&fp_evid$Potential.contaminant!="+")
fp_prot_filt <- filter(fp_prot, fp_prot$Reverse!="+"&fp_prot$Potential.contaminant!="+")


 #keep useful columns
# FP_evid <- data.frame(fp_evid_filt$Sequence, fp_evid_filt$Modified.sequence, fp_evid_filt$Intensity, fp_evid_filt$Missed.cleavages, fp_evid_filt$Proteins,
#                    fp_evid_filt$Raw.file, fp_evid_filt$Fraction, fp_evid_filt$Experiment, fp_evid_filt$PIF, fp_evid_filt$Score, fp_evid_filt$id)
# FP_prot <- data.frame(fp_prot_filt$Protein.IDs)
# 


#-----------------------------------------------------
  
#ID's counting

  #check MS and MS/MS scan numbers
    setwd(m)
    png("MS.png")
    plot(fp_sum$Fraction, fp_sum$MS, ylim = c(0,2500), ylab = "MS scans", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
    
    png("MSMS.png")
    plot(fp_sum$Fraction, fp_sum$MS.MS, ylim = c(0,10000), ylab = "MS/MS scans", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
  
   
  
  #check ID rates per fraction
    setwd(m)
    png("ID_rate.png")
    plot(fp_sum$Fraction, fp_sum$MS.MS.Identified...., ylim = c(0,70), ylab = "ÏD rate [%]", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
 
      
  #check PSMs per fraction
    setwd(m) 
    png("PSMs.png")
    plot(fp_sum$Fraction, fp_sum$Peptide.Sequences.Identified, ylim = c(0,7000), ylab = "PSMs", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
    
   
    
  #check Number of protein groups
    
    #filter
    Num_proteins_FP <- length(fp_prot_filt$Protein.IDs[which(fp_prot_filt$Intensity!="0")])
    
  #check number of peptide sequences
    
    #filter
    Num_peptide_FP <- length(unique(fp_evid_filt$Sequence))
   
    Num_modpeptide_FP <- length(unique(fp_evid_filt$Modified.sequence))
    
  #check distribution of missed cleavages
    
    #filter and plot
    fp_miscleav <- count(fp_evid_filt, fp_evid_filt$Missed.cleavages)
     
    setwd(m)
    png("Miscleavages.png")
    barplot(fp_miscleav$n,  ylab = "Counts", xlab = "miscleavages", main = c(fp_evid_filt[1,"Experiment"]), 
            names.arg =fp_miscleav$`fp_evid_filt$Missed.cleavages` )
    dev.off()
    
  #check number peptide sequences per fraction
    
    #filter
    Num_peptide_FP_per_frac <- aggregate(Sequence ~ Fraction, fp_evid_filt, function(x) length(unique(x)))  #Alternativ dt <-  fp_evid_filt %>% group_by(Fraction) %>% summarise(n_distinct(Sequence))

    Num_modpeptide_FP_per_frac <- aggregate(Modified.sequence ~ Fraction, fp_evid_filt, function(x) length(unique(x)))
   
    #plot
    setwd(m)
    png("Peptides.png")
    barplot(Num_peptide_FP_per_frac$Sequence, ylab = "peptides", xlab = "Fraction", 
            main = c(fp_evid_filt[1,"Experiment"]), names.arg = Num_peptide_FP_per_frac$Fraction, las=3)
    dev.off()
    
  
  #check summed peptide intensity per fraction
    
    #filter
    Int_peptide_FP_per_frac <- fp_evid_filt %>% group_by(Fraction) %>% summarize(sum(Intensity, na.rm = TRUE))
    
    #plot
    setwd(m)
    png("Peptide_intensity.png")
    barplot(Int_peptide_FP_per_frac$`sum(Intensity, na.rm = TRUE)`, ylab = "peptide intensity", xlab = "Fraction", 
            main = c(fp_evid_filt[1,"Experiment"]), names.arg = Int_peptide_FP_per_frac$Fraction, las=3)
    dev.off()
    
    setwd(m)
    png("Peptide_intensity2.png")
    boxplot(log(fp_evid_filt$Intensity) ~ fp_evid_filt$Fraction, ylab = "peptide intensity", xlab = "Fraction", main = c(fp_evid_filt[1,"Experiment"]))
    dev.off()
    
    setwd(m)
    repm <- data.frame(fp_evid_filt$Reporter.intensity.1,fp_evid_filt$Reporter.intensity.2,fp_evid_filt$Reporter.intensity.3,fp_evid_filt$Reporter.intensity.4,fp_evid_filt$Reporter.intensity.5,fp_evid_filt$Reporter.intensity.6,fp_evid_filt$Reporter.intensity.7,fp_evid_filt$Reporter.intensity.8,fp_evid_filt$Reporter.intensity.9,fp_evid_filt$Reporter.intensity.10,fp_evid_filt$Reporter.intensity.11)
    logrepm <- log(repm)
    png("Peptide_intensity3.png")
    boxplot(logrepm, ylab = "peptide intensity", xlab = "reporter ion", main = c(fp_evid_filt[1,"Experiment"]))
    dev.off()    
    
 
#-----------------------------------------------------

#Write table  
  
 
  x <- rep(NA, 48)
  a <- c(Num_proteins_FP, x)
  b <- c(Num_peptide_FP, Num_peptide_FP_per_frac$Sequence)
  c <- c(Num_modpeptide_FP, Num_modpeptide_FP_per_frac$Modified.sequence)
  d <- c(fp_sum2[50,"Peptide.Sequences.Identified"], fp_sum$Peptide.Sequences.Identified)
  e <- c(fp_sum2[50,"MS.MS.Identified...."], fp_sum$MS.MS.Identified....)
  f <- c(NA, Int_peptide_FP_per_frac$`sum(Intensity, na.rm = TRUE)`)
  table <- data.frame(a,b,c,d,e,f)
  row.names(table) = c("total", 1:48)
  colnames(table) = c("Protein groups","Peptides", "Mod_peptides", "PSMs", "ID rate", "peptide intensity")
  setwd(m)
  write.csv(table, file = "ID_table_FP.csv")

