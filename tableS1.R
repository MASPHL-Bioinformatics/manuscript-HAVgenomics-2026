readRenviron("_environment")

OUTPUT_DIR=Sys.getenv("OUTPUT_DIR")
TABLE_S1=Sys.getenv("TABLE_S1")
CLEANING_SCRIPT=Sys.getenv("CLEANING_SCRIPT")
source(CLEANING_SCRIPT)

s1_metadata <- metadata %>%
  subset(sequenced=="Yes", select=-county_of_residence)

s1_metadata$year <- format(s1_metadata$Specimen_Collection_Date, "%Y")
# one sample is missing collection date, but has other dates showing it was a 2024 sample (e.g., event date)
s1_metadata$year[is.na(s1_metadata$year)] <- "2024"

# combine genotype columns
s1_metadata$final_subgenotype <- ifelse(s1_metadata$Subgenotype==s1_metadata$GHOST.Result, s1_metadata$Subgenotype,
                                 ifelse(s1_metadata$Subgenotype=="" & s1_metadata$GHOST.Result!="", s1_metadata$GHOST.Result,
                                 ifelse(s1_metadata$Subgenotype!="" & s1_metadata$GHOST.Result=="", s1_metadata$Subgenotype,
                                 ifelse(s1_metadata$Subgenotype=="" & s1_metadata$GHOST.Result=="", "Missing", "Check"))))

#table(s1_metadata$final_subgenotype, s1_metadata$Subgenotype)

s1_metadata <- subset(s1_metadata, select=-c(GHOST.Result, Subgenotype))

s1_table <- cbind.data.frame("Sample ID"=s1_metadata$MASPHL_Sequencing_ID,
                             "BioSample Accession"=s1_metadata$Biosample_ID,
                             "Sample from Same Patient"=s1_metadata$Sample.from.Same.Patient,
                             "State of Residence"=s1_metadata$State.of.Residence,
                             "Year"=s1_metadata$year,
                             "MA Outbreak?"=s1_metadata$outbreak,
                             "Total Reads *"=s1_metadata$Total.Reads...,
                             "Aligned Reads"=s1_metadata$Aligned.Reads,
                             "Assembly Length"=s1_metadata$Assembly.Length,
                             "Unambig Bases"=s1_metadata$Unambig.Bases,
                             "Mean Read Depth"=s1_metadata$Mean.Read.Depth,
                             "Included?"=s1_metadata$Sample_Passed_QC,
                             "Subgenotype"=s1_metadata$final_subgenotype,
                             "Metagenomic"=s1_metadata$Metagenomic,
                             "Metagenomic Assembly Length"=s1_metadata$Metagenomic.Assembly.Length,
                             "Metagenomic Unambig Bases"=s1_metadata$Metagenomic.Unambig.Bases,
                             "Percent of Unclassified Reads Metagenomic"=s1_metadata$Percent.of.Unclassified.Reads.Metagenomic,
                             "Percent of HAV Reads Metagenomic"=s1_metadata$Percent.of.HAV.Reads.Metagenomic,
                             "Percent of Other Viral Reads Metagenomic"=s1_metadata$Percent.of.Other.Viral.Reads.Metagenomic,
                             "Twist Capture"=s1_metadata$Twist.Capture,
                             "Twist Assembly Length"=s1_metadata$Twist.Assembly.Length,
                             "Twist Unambig Bases"=s1_metadata$Twist.Unambig.Bases,
                             "Percent of Unclassified Reads Twist"=s1_metadata$Percent.of.Unclassified.Reads.Twist,
                             "Percent of HAV Reads Twist"=s1_metadata$Percent.of.HAV.Reads.Twist,	
                             "Percent of Other Viral Reads Twist"=s1_metadata$Percent.of.Other.Viral.Reads.Twist,
                             "HAV Panel Capture"=s1_metadata$HAV.Panel.Capture,
                             "HAV Panel Assembly Length"=s1_metadata$HAV.Panel.Assembly.Length,
                             "HAV Panel Unambig Bases"=s1_metadata$HAV.Panel.Unambig.Bases,
                             "Percent of Unclassified Reads HAV Panel"=s1_metadata$Percent.of.Unclassified.Reads.HAV.Panel,
                             "Percent of HAV Reads HAV Panel"=s1_metadata$Percent.of.HAV.Reads.HAV.Panel,
                             "Percent of Other Viral Reads HAV Panel"=s1_metadata$Percent.of.Other.Viral.Reads.HAV.Panel)

s1_table[is.na(s1_table)] <- ""

s1_table

write.csv(s1_table, TABLE_S1, row.names=FALSE)