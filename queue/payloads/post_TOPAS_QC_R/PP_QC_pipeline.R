#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


library(readr)
library(dplyr)
library(stringr)

#-----------------------------------------------------

#load data


n <- args[1] #"Z:\\internal_projects/active/TOPAS/WP31/Searches/JR/Biology_pilot/Batch3_PP/combined/txt"

n <- paste0(n, "/combined/txt")
cat(n, "\n")
setwd(n)
pp_evid <- read.delim("evidence.txt", stringsAsFactors = FALSE, sep="\t")
pp_site <- read.delim("Phospho (STY)Sites.txt", stringsAsFactors = FALSE, sep="\t")
pp_prot <- read.delim("proteinGroups.txt", stringsAsFactors = FALSE, sep="\t")
pp_sum <- read.delim("summary.txt", stringsAsFactors = FALSE, sep="\t", nrows = 12)
pp_sum2 <- read.delim("summary.txt", stringsAsFactors = FALSE, sep="\t")

#-----------------------------------------------------

#filter data
 #remove contaminants and decoys

pp_evid_filt <- filter(pp_evid, pp_evid$Reverse!="+"&pp_evid$Potential.contaminant!="+")
pp_site_filt <- filter(pp_site, pp_site$Reverse!="+"&pp_site$Potential.contaminant!="+")
pp_prot_filt <- filter(pp_prot, pp_prot$Reverse!="+"&pp_prot$Potential.contaminant!="+")


 #keep useful columns

# PP_evid <- data.frame(pp_evid_filt$Sequence)
# PP_site <- data.frame(pp_site_filt$Proteins)

#-----------------------------------------------------
  
#ID's counting

  #check MS and MS/MS scan numbers
   
    setwd(n)
    png("MS.png")
    plot(pp_sum$Fraction, pp_sum$MS, ylim = c(0,5000), ylab = "MS scans", xlab = "Fraction", pch=16, cex=2, main = c(pp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(pp_sum$Fraction), at = c(pp_sum$Fraction), las=3)
    dev.off()
    
    png("MSMS.png")
    plot(pp_sum$Fraction, pp_sum$MS.MS, ylim = c(0,20000), ylab = "MS/MS scans", xlab = "Fraction", pch=16, cex=2, main = c(pp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(pp_sum$Fraction), at = c(pp_sum$Fraction), las=3)
    dev.off()
  
  #check ID rates per fraction
   
    setwd(n)
    png("ID_rate.png")
    plot(pp_sum$Fraction, pp_sum$MS.MS.Identified...., ylim = c(0,70), ylab = "ÏD rate [%]", xlab = "Fraction", pch=16, cex=2, main = c(pp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(pp_sum$Fraction), at = c(pp_sum$Fraction), las=3)
    dev.off()
      
  #check PSMs per fraction
   
    setwd(n)
    png("PSMs.png")
    plot(pp_sum$Fraction, pp_sum$Peptide.Sequences.Identified, ylim = c(0,7000), ylab = "PSMs", xlab = "Fraction", pch=16, cex=2, main = c(pp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(pp_sum$Fraction), at = c(pp_sum$Fraction), las=3)
    dev.off()
    
  #check Number of protein groups
    
    #filter
    Num_proteins_PP <- length(pp_prot_filt$Protein.IDs[which(pp_prot_filt$Intensity!="0")])
    Num_phosphoproteins_PP <- length(pp_prot_filt$Protein.IDs[which(pp_prot_filt$Phospho..STY..site.positions!="" & pp_prot_filt$Intensity!="0")])
    
  #check number of peptide sequences
    
    #filter
    Num_peptide_PP <- length(unique(pp_evid_filt$Sequence))
    Num_phosphopeptide_PP <- length(unique(pp_evid_filt$Sequence[which(pp_evid_filt$Phospho..STY.!="0")]))
    
    Num_modpeptide_PP <- length(unique(pp_evid_filt$Modified.sequence))
    Num_modphosphopeptide_PP <- length(unique(pp_evid_filt$Modified.sequence[which(pp_evid_filt$Phospho..STY.!="0")]))
    
  #check distribution of missed cleavages
    
    #filter and plot
    pp_miscleav <- count(pp_evid_filt, pp_evid_filt$Missed.cleavages)
    
   
    setwd(n)
    png("Miscleavages.png")
    barplot(pp_miscleav$n,  ylab = "Counts", xlab = "miscleavages", main = c(pp_evid_filt[1,"Experiment"]), 
            names.arg =pp_miscleav$`pp_evid_filt$Missed.cleavages` )
    dev.off()
    
   #check number peptide sequences per fraction
    
    #filter

    Num_peptide_PP_per_frac <- aggregate(Sequence ~ Fraction, pp_evid_filt, function(x) length(unique(x)))
    Num_modpeptide_PP_per_frac <- aggregate(Modified.sequence ~ Fraction, pp_evid_filt, function(x) length(unique(x)))
    Num_phosphopeptide_PP_per_frac <- aggregate(Sequence ~ Fraction, filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0"), function(x) length(unique(x)))
    Num_modphosphopeptide_PP_per_frac <- aggregate(Modified.sequence ~ Fraction, filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0"), function(x) length(unique(x)))
   
    #plot
    
    setwd(n)
    png("Peptides.png")
    barplot(Num_peptide_PP_per_frac$Sequence, ylab = "peptides", xlab = "Fraction", 
            main = c(pp_evid_filt[1,"Experiment"]), names.arg = Num_peptide_PP_per_frac$Fraction, las=3)
    dev.off()
    
    png("Peptides.png")
    barplot(Num_phosphopeptide_PP_per_frac$Sequence, ylab = "peptides", xlab = "Fraction", 
            main = c(pp_evid_filt[1,"Experiment"]), names.arg = Num_phosphopeptide_PP_per_frac$Fraction, las=3)
    dev.off()
    
  #check summed peptide intensity per fraction
    
    #filter

    Int_peptide_PP_per_frac <- pp_evid_filt %>% group_by(Fraction) %>% summarize(sum(Intensity, na.rm = TRUE))
    Int_phosphopeptide_PP_per_frac <- filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0") %>% group_by(Fraction) %>% summarize(sum(Intensity, na.rm = TRUE))
    
    #plot
   
    setwd(n)
    png("Peptide_intensity.png")
    barplot(Int_peptide_PP_per_frac$`sum(Intensity, na.rm = TRUE)`, ylab = "peptide intensity", xlab = "Fraction", 
            main = c(pp_evid_filt[1,"Experiment"]), names.arg = Int_peptide_PP_per_frac$Fraction, las=3)
    dev.off()
    
    png("Ppeptide_intensity.png")
    barplot(Int_phosphopeptide_PP_per_frac$`sum(Intensity, na.rm = TRUE)`, ylab = "phosphopeptide intensity", xlab = "Fraction", 
            main = c(pp_evid_filt[1,"Experiment"]), names.arg = Int_phosphopeptide_PP_per_frac$Fraction, las=3)
    dev.off()
    
    png("Peptide_intensity2.png")
    boxplot(log(pp_evid_filt$Intensity) ~ pp_evid_filt$Fraction, ylab = "peptide intensity", xlab = "Fraction", main = c(pp_evid_filt[1,"Experiment"]))
    dev.off()
    
    repm <- data.frame(pp_evid_filt$Reporter.intensity.1,pp_evid_filt$Reporter.intensity.2,pp_evid_filt$Reporter.intensity.3,pp_evid_filt$Reporter.intensity.4,pp_evid_filt$Reporter.intensity.5,pp_evid_filt$Reporter.intensity.6,pp_evid_filt$Reporter.intensity.7,pp_evid_filt$Reporter.intensity.8,pp_evid_filt$Reporter.intensity.9,pp_evid_filt$Reporter.intensity.10,pp_evid_filt$Reporter.intensity.11)
    logrepm <- log(repm)
    png("Peptide_intensity3.png")
    boxplot(logrepm, ylab = "peptide intensity", xlab = "reporter ion", main = c(pp_evid_filt[1,"Experiment"]))
    dev.off()
    
    pp_evid_filt_filt <- filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0")
    repm <- data.frame(pp_evid_filt_filt$Reporter.intensity.1,pp_evid_filt_filt$Reporter.intensity.2,pp_evid_filt_filt$Reporter.intensity.3,pp_evid_filt_filt$Reporter.intensity.4,pp_evid_filt_filt$Reporter.intensity.5,pp_evid_filt_filt$Reporter.intensity.6,pp_evid_filt_filt$Reporter.intensity.7,pp_evid_filt_filt$Reporter.intensity.8,pp_evid_filt_filt$Reporter.intensity.9,pp_evid_filt_filt$Reporter.intensity.10,pp_evid_filt_filt$Reporter.intensity.11)
    logrepm <- log(repm)
    png("Ppeptide_intensity3.png")
    boxplot(logrepm, ylab = "phosphopeptide intensity", xlab = "reporter ion", main = c(pp_evid_filt[1,"Experiment"]))
    dev.off() 
    
#-----------------------------------------------------

#IMAC selectivity (ID and intensity based)

        # #ID based
        #   Num_peptide_PP_per_frac <- aggregate(Sequence ~ Fraction, pp_evid_filt, function(x) length(unique(x)))
        #   Num_phosphopeptide_PP_per_frac <- aggregate(Sequence ~ Fraction, filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0"), function(x) length(unique(x)))
        #   IMAC_sel <- Num_phosphopeptide_PP_per_frac/Num_peptide_PP_per_frac

  #Intensity based
  Int_peptide_PP_per_frac <- pp_evid_filt %>% group_by(Fraction) %>% summarize(sum(Intensity, na.rm = TRUE))
  Int_phosphopeptide_PP_per_frac <- filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0") %>% group_by(Fraction) %>% summarize(sum(Intensity, na.rm = TRUE))
  IMAC_sel <- data.frame(Int_peptide_PP_per_frac$Fraction,Int_phosphopeptide_PP_per_frac$`sum(Intensity, na.rm = TRUE)`/Int_peptide_PP_per_frac$`sum(Intensity, na.rm = TRUE)`)
        
  #plot
  setwd(n)
  png("IMAC_selectivity.png")
  barplot(IMAC_sel$Int_phosphopeptide_PP_per_frac..sum.Intensity..na.rm...TRUE...Int_peptide_PP_per_frac..sum.Intensity..na.rm...TRUE.., ylab ='Selectivity IMAC', xlab ='Fraction', 
          main = c(pp_evid_filt[1,"Experiment"]), names.arg = IMAC_sel$Int_peptide_PP_per_frac.Fraction)
  dev.off()

#-----------------------------------------------------
  
#Write table  
  
  y <- rep(NA, 12)
  a <- c(Num_proteins_PP, y)
  b <- c(Num_phosphoproteins_PP, y)
  c <- c(Num_peptide_PP, Num_peptide_PP_per_frac$Sequence)
  d <- c(Num_phosphopeptide_PP, Num_phosphopeptide_PP_per_frac$Sequence)
  e <- c(Num_modpeptide_PP, Num_modpeptide_PP_per_frac$Modified.sequence)
  f <- c(Num_modphosphopeptide_PP, Num_modphosphopeptide_PP_per_frac$Modified.sequence)
  g <- c(pp_sum2[14,"Peptide.Sequences.Identified"], pp_sum$Peptide.Sequences.Identified)
  h <- c(pp_sum2[14,"MS.MS.Identified...."], pp_sum$MS.MS.Identified....)
  i <- c(NA, Int_peptide_PP_per_frac$`sum(Intensity, na.rm = TRUE)`)
  j <- c(NA, Int_phosphopeptide_PP_per_frac$`sum(Intensity, na.rm = TRUE)`)
  k <- c(NA, IMAC_sel$Int_phosphopeptide_PP_per_frac..sum.Intensity..na.rm...TRUE...Int_peptide_PP_per_frac..sum.Intensity..na.rm...TRUE..)
  table <- data.frame(a,b,c,d,e,f,g,h,i,j,k)
  row.names(table) = c("total", 1:12)
  colnames(table) = c("Protein groups","Phosphoprotein groups", "Peptides", "Phosphopeptides", "Mod_peptides",
                      "Mod_phosphopeptides", "PSMs", "ID rate", "peptide intensity", "phosphopeptide intensity", "IMAC selectivity")
  setwd(n)
  write.csv(table, file = "ID_table_PP.csv")
  
#-----------------------------------------------------

