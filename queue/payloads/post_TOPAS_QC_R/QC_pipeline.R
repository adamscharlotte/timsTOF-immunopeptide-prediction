library(readr)
library(dplyr)
library(stringr)

#-----------------------------------------------------

#load data

m <- "Z:\\internal_projects/active/TOPAS/WP31/Searches/JR/Biology_pilot/Batch3_FP/combined/txt"
n <- "Z:\\internal_projects/active/TOPAS/WP31/Searches/JR/Biology_pilot/Batch3_PP/combined/txt"

setwd(m)
fp_evid <- read.delim("evidence.txt", stringsAsFactors = FALSE, sep="\t")
fp_prot <- read.delim("proteinGroups.txt", stringsAsFactors = FALSE, sep="\t")
fp_sum <- read.delim("summary.txt", stringsAsFactors = FALSE, sep="\t", nrows = 48)
fp_sum2 <- read.delim("summary.txt", stringsAsFactors = FALSE, sep="\t")

setwd(n)
pp_evid <- read.delim("evidence.txt", stringsAsFactors = FALSE, sep="\t")
pp_site <- read.delim("Phospho (STY)Sites.txt", stringsAsFactors = FALSE, sep="\t")
pp_prot <- read.delim("proteinGroups.txt", stringsAsFactors = FALSE, sep="\t")
pp_sum <- read.delim("summary.txt", stringsAsFactors = FALSE, sep="\t", nrows = 12)
pp_sum2 <- read.delim("summary.txt", stringsAsFactors = FALSE, sep="\t")

#-----------------------------------------------------

#filter data
 #remove contaminants and decoys
fp_evid_filt <- filter(fp_evid, fp_evid$Reverse!="+"&fp_evid$Potential.contaminant!="+")
fp_prot_filt <- filter(fp_prot, fp_prot$Reverse!="+"&fp_prot$Potential.contaminant!="+")

pp_evid_filt <- filter(pp_evid, pp_evid$Reverse!="+"&pp_evid$Potential.contaminant!="+")
pp_site_filt <- filter(pp_site, pp_site$Reverse!="+"&pp_site$Potential.contaminant!="+")
pp_prot_filt <- filter(pp_prot, pp_prot$Reverse!="+"&pp_prot$Potential.contaminant!="+")


 #keep useful columns
# FP_evid <- data.frame(fp_evid_filt$Sequence, fp_evid_filt$Modified.sequence, fp_evid_filt$Intensity, fp_evid_filt$Missed.cleavages, fp_evid_filt$Proteins,
#                    fp_evid_filt$Raw.file, fp_evid_filt$Fraction, fp_evid_filt$Experiment, fp_evid_filt$PIF, fp_evid_filt$Score, fp_evid_filt$id)
# FP_prot <- data.frame(fp_prot_filt$Protein.IDs)
# 
# PP_evid <- data.frame(pp_evid_filt$Sequence)
# PP_site <- data.frame(pp_site_filt$Proteins)

#-----------------------------------------------------
  
#ID's counting

  #check MS and MS/MS scan numbers
    setwd(m)
    png("MS.png")
    plot(fp_sum$Fraction, fp_sum$MS, ylim = c(0,70), ylab = "MS scans", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
    
    png("MSMS.png")
    plot(fp_sum$Fraction, fp_sum$MS.MS, ylim = c(0,70), ylab = "MS/MS scans", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
  
    setwd(n)
    png("MS.png")
    plot(pp_sum$Fraction, pp_sum$MS, ylim = c(0,70), ylab = "MS scans", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
    
    png("MSMS.png")
    plot(pp_sum$Fraction, pp_sum$MS.MS, ylim = c(0,70), ylab = "MS/MS scans", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
  
  #check ID rates per fraction
    setwd(m)
    png("ID_rate.png")
    plot(fp_sum$Fraction, fp_sum$MS.MS.Identified...., ylim = c(0,70), ylab = "ÏD rate [%]", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
    
    setwd(n)
    png("ID_rate.png")
    plot(pp_sum$Fraction, pp_sum$MS.MS.Identified...., ylim = c(0,70), ylab = "ÏD rate [%]", xlab = "Fraction", pch=16, cex=2, main = c(pp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(pp_sum$Fraction), at = c(pp_sum$Fraction), las=3)
    dev.off()
      
  #check PSMs per fraction
    setwd(m) 
    png("PSMs.png")
    plot(fp_sum$Fraction, fp_sum$Peptide.Sequences.Identified, ylim = c(0,7000), ylab = "PSMs", xlab = "Fraction", pch=16, cex=2, main = c(fp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(fp_sum$Fraction), at = c(fp_sum$Fraction), las=3)
    dev.off()
    
    setwd(n)
    png("PSMs.png")
    plot(pp_sum$Fraction, pp_sum$Peptide.Sequences.Identified, ylim = c(0,7000), ylab = "PSMs", xlab = "Fraction", pch=16, cex=2, main = c(pp_sum[1,"Experiment"]), xaxt='n')
    axis(1, labels = c(pp_sum$Fraction), at = c(pp_sum$Fraction), las=3)
    dev.off()
    
  #check Number of protein groups
    
    #filter
    Num_proteins_FP <- length(fp_prot_filt$Protein.IDs[which(fp_prot_filt$Intensity!="0")])
    Num_proteins_PP <- length(pp_prot_filt$Protein.IDs[which(pp_prot_filt$Intensity!="0")])
    Num_phosphoproteins_PP <- length(pp_prot_filt$Protein.IDs[which(pp_prot_filt$Phospho..STY..site.positions!="" & pp_prot_filt$Intensity!="0")])
    
  #check number of peptide sequences
    
    #filter
    Num_peptide_FP <- length(unique(fp_evid_filt$Sequence))
    Num_peptide_PP <- length(unique(pp_evid_filt$Sequence))
    Num_phosphopeptide_PP <- length(unique(pp_evid_filt$Sequence[which(pp_evid_filt$Phospho..STY.!="0")]))
    
    Num_modpeptide_FP <- length(unique(fp_evid_filt$Modified.sequence))
    Num_modpeptide_PP <- length(unique(pp_evid_filt$Modified.sequence))
    Num_modphosphopeptide_PP <- length(unique(pp_evid_filt$Modified.sequence[which(pp_evid_filt$Phospho..STY.!="0")]))
    
  #check distribution of missed cleavages
    
    #filter and plot
    fp_miscleav <- count(fp_evid_filt, fp_evid_filt$Missed.cleavages)
    pp_miscleav <- count(pp_evid_filt, pp_evid_filt$Missed.cleavages)
    
    setwd(m)
    png("Miscleavages.png")
    barplot(fp_miscleav$n,  ylab = "Counts", xlab = "miscleavages", main = c(fp_evid_filt[1,"Experiment"]), 
            names.arg =fp_miscleav$`fp_evid_filt$Missed.cleavages` )
    dev.off()
    
    setwd(n)
    png("Miscleavages.png")
    barplot(pp_miscleav$n,  ylab = "Counts", xlab = "miscleavages", main = c(pp_evid_filt[1,"Experiment"]), 
            names.arg =pp_miscleav$`pp_evid_filt$Missed.cleavages` )
    dev.off()
    
   #check number peptide sequences per fraction
    
    #filter
    Num_peptide_FP_per_frac <- aggregate(Sequence ~ Fraction, fp_evid_filt, function(x) length(unique(x)))  #Alternativ dt <-  fp_evid_filt %>% group_by(Fraction) %>% summarise(n_distinct(Sequence))
    Num_peptide_PP_per_frac <- aggregate(Sequence ~ Fraction, pp_evid_filt, function(x) length(unique(x)))
    Num_modpeptide_FP_per_frac <- aggregate(Modified.sequence ~ Fraction, fp_evid_filt, function(x) length(unique(x)))
    Num_modpeptide_PP_per_frac <- aggregate(Modified.sequence ~ Fraction, pp_evid_filt, function(x) length(unique(x)))
    Num_phosphopeptide_PP_per_frac <- aggregate(Sequence ~ Fraction, filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0"), function(x) length(unique(x)))
    Num_modphosphopeptide_PP_per_frac <- aggregate(Modified.sequence ~ Fraction, filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0"), function(x) length(unique(x)))
   
    #plot
    setwd(m)
    png("Peptides.png")
    barplot(Num_peptide_FP_per_frac$Sequence, ylab = "peptides", xlab = "Fraction", 
            main = c(fp_evid_filt[1,"Experiment"]), names.arg = Num_peptide_FP_per_frac$Fraction, las=3)
    dev.off()
    
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
    Int_peptide_FP_per_frac <- fp_evid_filt %>% group_by(Fraction) %>% summarize(sum(Intensity, na.rm = TRUE))
    Int_peptide_PP_per_frac <- pp_evid_filt %>% group_by(Fraction) %>% summarize(sum(Intensity, na.rm = TRUE))
    Int_phosphopeptide_PP_per_frac <- filter(pp_evid_filt, pp_evid_filt$Phospho..STY.!="0") %>% group_by(Fraction) %>% summarize(sum(Intensity, na.rm = TRUE))
    
    #plot
    setwd(m)
    png("Peptide_intensity.png")
    barplot(Int_peptide_FP_per_frac$`sum(Intensity, na.rm = TRUE)`, ylab = "peptide intensity", xlab = "Fraction", 
            main = c(fp_evid_filt[1,"Experiment"]), names.arg = Int_peptide_FP_per_frac$Fraction, las=3)
    dev.off()
    
    setwd(n)
    png("Peptide_intensity.png")
    barplot(Int_peptide_PP_per_frac$`sum(Intensity, na.rm = TRUE)`, ylab = "peptide intensity", xlab = "Fraction", 
            main = c(pp_evid_filt[1,"Experiment"]), names.arg = Int_peptide_PP_per_frac$Fraction, las=3)
    dev.off()
    
    png("Peptide_intensity.png")
    barplot(Int_phosphopeptide_PP_per_frac$`sum(Intensity, na.rm = TRUE)`, ylab = "phosphopeptide intensity", xlab = "Fraction", 
            main = c(pp_evid_filt[1,"Experiment"]), names.arg = Int_phosphopeptide_PP_per_frac$Fraction, las=3)
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

#Labeling efficiency
  library(readr)
  library(dplyr)
  library(stringr)
  
  #load data
  
    o <- "Z:\\internal_projects/active/TOPAS/WP31/Searches/JR/Biology_pilot/Batch2_var/combined/txt"
  
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