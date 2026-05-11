########################
##### Folder paths #####
########################

readRenviron("_environment")

OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR")
METADATA_FILE <- Sys.getenv("METADATA_FILE")
TOWN_SHAPE <- Sys.getenv("TOWN_SHAPE")
STATE_POP <- Sys.getenv("STATE_POP")
TOWN_POP <- Sys.getenv("TOWN_POP")
TOWN_COUNTRY <- Sys.getenv("TOWN_COUNTRY")
MMWR_KEY <- Sys.getenv("MMWR_KEY")
IMMUNE <- Sys.getenv("IMMUNE")
NH_SAMPLES <- Sys.getenv("NH_SAMPLES")
NH_METADATA <- Sys.getenv("NH_METADATA")


#############################
##### Global parameters #####
#############################
knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)

#####################
##### Libraries #####
#####################

# install required packages with Bioconductor
pkglist = c("dplyr", "stringi", "sf", "readxl","geodist")

lapply(pkglist, require, character.only = TRUE)

############################
##### Custom functions #####
############################
##### Rolling average
# calculates a rolling average; supply the vector x and the length of the rolling average (in this example, it is 5 units)
ma <- function(x, n = 5){
  # sides=1 is saying go from current row backwards to calculate rolling average
  stats::filter(x, rep(1 / n, n), sides = 1)
}


##### Metadata summary table by vaccination status
# calculates summary data frames for vaccination data in both wide and long format
immune_sum <- function(category="Sex") {
  # get summary table of vaccines administered by selected variable's categories, week, and year
  df1 <- subset(immune_long, BIN==category) %>%
    group_by(week, year, wkYr, binCat) %>%
    summarize(wkImmuneCt = sum(vaccine_count))
  
  # merge in actual dates (just have MMWR weeks at this point)
  df1 <- merge(df1,
               subset(mmwr_key, select=c(wkYr, endDt)),
               by.all="wkYr",
               all.x=T)
  
  # merge in state population by year so we can do vaccine rate calculations
  df1 <- merge(df1,
               state_pop,
               by.x="year",
               by.y="Year",
               all.x=TRUE) %>%
    arrange(endDt)
  
  # go from long to wide format
  df2 <- subset(df1, select=c(wkYr, endDt, wkImmuneCt, binCat)) %>%
    tidyr::pivot_wider(id_cols=c(wkYr, endDt),
                       names_from=binCat,
                       values_from=wkImmuneCt) %>%
    arrange(endDt)
  
  # rename variables
  df2 <- df2 %>%
    rename_with(.fn=~paste0("wkImmuneCt_", category, .), .cols=all_of(colnames(df2[,3:ncol(df2)])))
  
  # output both long and wide formats (made both so figured it might be nice to have both potentially)
  out <- list(df1, df2)
  
  return(out)
}


##### Table 1 variable summary
# gives summary/descriptive stats for metadata variables in the format needed for Table 1
summary_table_char <- function(data, var_int, name, columnnames) {
  # count number of times each category within the chosen variable appears, not dropping the 0 count variables
  df_1 <- data %>%
    count({{var_int}}, .drop=FALSE)
  
  # get the % of each category
  df_1$`%` <- round(df_1$n/sum(df_1$n)*100, 1)
  
  # creating a label column for the percentage and sample size (i.e., X% (N=X)) 
  df_1$label <- paste0(df_1$`%`, "% (", df_1$n, ")")
  
  # selecting the desired columns (just need the names of the categories and label to display)
  df_1 <- subset(df_1, select=c(1, 4))
  
  # rename columns
  colnames(df_1) <- c("Variable", columnnames)
  
  # fix levels from original assignment in the metadata
  df_1$Variable <- as.character(df_1$Variable)
  
  # adding a variable header to the category %'s and counts
  df_1 <- rbind.data.frame(c(name, ""), df_1)
  
  return(df_1)
}


##############################################
##### Shapefile for MA towns and NH & CT #####
##############################################
gis_ma_town_shape <- read_sf(TOWN_SHAPE)

# converting cities and counties to have capital letters for the start of each word
gis_ma_town_shape$town <- stri_trans_totitle(gis_ma_town_shape$TOWN)

# rename geometry variable to something a little more user-friendly (and allows for merging multiple geometries into one dataframe later)
names(gis_ma_town_shape)[names(gis_ma_town_shape)=="geometry"] <- "townGeo"



################################
##### Donahue denominators #####
################################
# state population
state_pop <- read.csv(STATE_POP)

# change Population column name to statePop
names(state_pop)[names(state_pop)=="Population"] <- "statePop"


# town population
town_pop <- read.csv(TOWN_POP)

# I am trying to match town naming conventions across all of the datasets, which requires a little manual processing
town_pop$town[town_pop$town=="Manchester"] <- "Manchester-By-The-Sea"

# merge in town shape with the town shape data
town_pop <- merge(town_pop,
                  gis_ma_town_shape,
                  by.all="town",
                  all=T)

# change Population column name to townPop
names(town_pop)[names(town_pop)=="Population"] <- "townPop"




###############################
##### Town and county key #####
###############################
town_county <- read.csv(TOWN_COUNTRY)



####################
##### Week key #####
####################
mmwr_key <- read_xlsx(MMWR_KEY)

# properly formatting the start and end dates for each MMWR week
mmwr_key$startDt <- as.Date(mmwr_key$WeekStartDate, format="%Y-%m-%d")
mmwr_key$endDt <- as.Date(mmwr_key$WeekEndDate...5, format="%Y-%m-%d")

# creating a combined wk_yr column
mmwr_key$wkYr <- paste0(mmwr_key$`MMWR Week`, "_", mmwr_key$Year)


####################
##### Metadata #####
####################
# MA MAVEN data
metadata <- read.csv(METADATA_FILE)

# conversion to proper date format
metadata$event_date <- as.Date(metadata$event_date, format="%m/%d/%Y")
metadata$Specimen_Collection_Date <- as.Date(metadata$Specimen_Collection_Date, format="%m/%d/%Y")
metadata$diagnostic_specimen_date <- as.Date(metadata$diagnostic_specimen_date, format="%m/%d/%Y")
metadata$symptom_onset_date <- as.Date(metadata$symptom_onset_date, format="%m/%d/%Y")
metadata$First_Positive_Specimen_Date <- as.Date(metadata$First_Positive_Specimen_Date, format="%m/%d/%Y")

# convert age to numeric from integer
metadata$age <- as.numeric(metadata$age)

# format caseid
metadata$caseid <- as.character(metadata$caseid)

# now drop HepA/USA/MA/0046/2018
metadata <- metadata[metadata$MASPHL_Sequencing_ID != "HepA/USA/MA/0046/2018", ]


####################
#### MIIS data #####
####################
immune <- read_xlsx(IMMUNE, sheet="Output 1-2-2025")


###################
##### NH data #####
###################
# NH sequences sent for sequencing...the sequenced cases can be found in qc_samples
nh_sequenced <- read.table(NH_SAMPLES, header=TRUE, sep ="\t") 

# all NH outbreak cases
nh_outbreak <- read_xlsx(NH_METADATA, sheet="Epi Curve Dates")

# NH metadata (summary, not row-level)
nh_metadata <- read_xlsx(NH_METADATA, sheet="HAV Data")



## Get vector of dates for each day within each MMWR week of interest
# creating a list of dates for every day within each MMWR week
full_week_wide <- lubridate::Date()

for (i in 1:nrow(mmwr_key)) {
  # this creates a vector of dates from the start date to the end date of each MMWR week
  date_vec <- seq(mmwr_key$startDt[i], mmwr_key$endDt[i], by=1)
  # bind by column to create a dataframe
  full_week_wide <- cbind(full_week_wide, date_vec)
}

# I am sure I could put this step into the loop, but it doesn't work because the empty dataframe in the very first loop step. It became easier to just pull that step out of the loop.
full_week_wide <- as.data.frame(full_week_wide)

# get the correct date formats for each row
for (i in 1:ncol(full_week_wide)) {
  full_week_wide[,i] <- as.Date(full_week_wide[,i], format="%Y-%m-%d")
}

colnames(full_week_wide) <- paste0(mmwr_key$`MMWR Week`, "_", mmwr_key$Year)

# creates a long version of the full_week_wide with just wkYr and date
full_week_long <- full_week_wide %>%
  tidyr::pivot_longer(cols=1:ncol(full_week_wide),
                      names_to="wkYr",
                      values_to="dates")




##### Gender #####
# gender
#unique(metadata$gender)
metadata$genderGrp3 <- ifelse(metadata$gender=="", "Unknown", metadata$gender)
#unique(metadata$genderGrp3)
metadata$genderGrp3 <- factor(metadata$genderGrp3, levels=c("Female", "Male", "Unknown"))
#table(metadata$genderGrp3, metadata$gender)


##### Age #####
#summary(metadata$age)
metadata$ageGrp5 <- ifelse(metadata$age<20 & !(is.na(metadata$age)), "<20 years old",
                           ifelse(metadata$age>=20 & metadata$age<=39 & !(is.na(metadata$age)), "20-39 years old",
                                  ifelse(metadata$age>=40 & metadata$age<=59 & !(is.na(metadata$age)), "40-59 years old",
                                         ifelse(metadata$age>=60 & !(is.na(metadata$age)), ">=60 years old",
                                                # I put in a Check grouping category so that any mistakes in assignment here are glaringly obvious when I check values       
                                                ifelse(is.na(metadata$age), "Unknown", "Check grouping")))))
#unique(metadata$ageGrp5)
metadata$ageGrp5 <- factor(metadata$ageGrp5, levels=c("<20 years old", "20-39 years old", "40-59 years old", ">=60 years old", "Unknown"))
#table(metadata$ageGrp5, metadata$age)




##### County #####
# county groups - making sure NA is showing as Unknown/missing
#unique(metadata$county_of_residence)
metadata$countyGrp <- ifelse(metadata$county_of_residence=="", "Unknown", metadata$county_of_residence)
#unique(metadata$countyGrp)
#table(metadata$countyGrp, metadata$county_of_residence)




##### Town locations #####
### residence, work, and clinician cities
## residence
# this function makes sure that the start of every word is capitalized
#unique(metadata$official_city)
metadata$cityResidence <- stri_trans_totitle(ifelse(is.na(metadata$official_city), "Unknown", 
                                                    ifelse(metadata$official_city=="", "Unknown", 
                                                           ifelse(metadata$official_city=="N/A", "Unknown", metadata$official_city))))
#table(metadata$cityResidence, metadata$official_city)




##### Food handler #####
#unique(metadata$IS_CASE_A_FOODHANDLER)
metadata$foodhandler <- ifelse(metadata$IS_CASE_A_FOODHANDLER=="", "Unknown", metadata$IS_CASE_A_FOODHANDLER)
#table(metadata$foodhandler, metadata$IS_CASE_A_FOODHANDLER)




##### Symptoms #####
# symptoms present
#unique(metadata$symptoms)
metadata$symptomsGrp3 <- ifelse(metadata$symptoms=="UNKNOWN", "Unknown",
                                ifelse(metadata$symptoms=="", "Unknown",metadata$symptoms))
#table(metadata$symptomsGrp3, metadata$symptoms)


# nausea
#unique(metadata$nausea)
metadata$nauseaGrp3 <- ifelse(metadata$nausea=="", "Unknown", metadata$nausea)
#table(metadata$nauseaGrp3, metadata$nausea)


# vomiting
#unique(metadata$vomiting)
metadata$vomitGrp3 <- ifelse(metadata$vomiting=="", "Unknown", metadata$vomiting)
#table(metadata$vomitGrp3, metadata$vomiting)


# diarrhea
#unique(metadata$diarrhea)
metadata$diarrheaGrp3 <- ifelse(metadata$diarrhea=="", "Unknown", metadata$diarrhea)
#table(metadata$diarrheaGrp3, metadata$diarrhea)


# jaundice
#unique(metadata$jaundice)
metadata$jaundiceGrp3 <- ifelse(metadata$jaundice=="", "Unknown", metadata$jaundice)
#table(metadata$jaundiceGrp3, metadata$jaundice)




##### Drug use #####
#unique(metadata$recent_illicit)
metadata$recentIllicitGrp3 <- ifelse(metadata$recent_illicit=="", "Unknown", metadata$recent_illicit)
#unique(metadata$recentIllicitGrp3)
metadata$recentIllicitGrp3 <- factor(metadata$recentIllicitGrp3, levels=c("Yes", "No", "Unknown"))
#table(metadata$recentIllicitGrp3, metadata$recent_illicit)


#unique(metadata$illicit_ever)
metadata$illicitEverGrp3 <- ifelse(metadata$illicit_ever=="", "Unknown", metadata$illicit_ever)
#unique(metadata$illicitEverGrp3)
metadata$illicitEverGrp3 <- factor(metadata$illicitEverGrp3, levels=c("Yes", "No", "Unknown"))
#table(metadata$illicitEverGrp3, metadata$illicit_ever)


#unique(metadata$recent_injection)
metadata$recentInjectGrp3 <- ifelse(metadata$recent_injection=="", "Unknown", metadata$recent_injection)
#unique(metadata$recentInjectGrp3)
metadata$recentInjectGrp3 <- factor(metadata$recentInjectGrp3, levels=c("Yes", "No", "Unknown"))
#table(metadata$recentInjectGrp3, metadata$recent_injection)

#unique(metadata$injection_ever)
metadata$injectEverGrp3 <- ifelse(metadata$injection_ever=="", "Unknown", metadata$injection_ever)
#unique(metadata$injectEverGrp3)
metadata$injectEverGrp3 <- factor(metadata$injectEverGrp3, levels=c("Yes", "No", "Unknown"))
#table(metadata$injectEverGrp3, metadata$injection_ever)




##### Vaccines #####
#unique(metadata$vaccine)
metadata$vacGrp2 <- ifelse(metadata$vaccine=="yes", "Yes",
                           ifelse(metadata$vaccine=="no/unknown", "No/unknown", 
                                  ifelse(metadata$vaccine=="", "No/unknown", "Check grouping")))
#unique(metadata$vacGrp2)
metadata$vacGrp2 <- factor(metadata$vacGrp2, levels=c("Yes", "No/unknown"))
#table(metadata$vacGrp2, metadata$vaccine)




##### Housing status #####
# creating a flag that says if homeless and/or unstably housed, than they are homeless/unstably housed
#unique(metadata$unstablehousing)
metadata$unstableHousingGrp3 <- ifelse(metadata$unstablehousing=="", "Unknown", metadata$unstablehousing)
#table(metadata$unstableHousingGrp3, metadata$unstablehousing)


#unique(metadata$homelessness)
metadata$homelessGrp3 <- ifelse(metadata$homelessness=="", "Unknown", metadata$homelessness)
#unique(metadata$homelessGrp3)
metadata$homelessGrp3 <- factor(metadata$homelessGrp3, levels=c("Yes", "No", "Unknown"))
#table(metadata$homelessGrp3, metadata$homelessness)


metadata$homeless_unstable <- ifelse(metadata$homelessGrp3=="Yes" | metadata$unstableHousingGrp3=="Yes", "Yes", 
                                     ifelse(metadata$homelessGrp3=="Unknown" & metadata$unstableHousingGrp3=="Unknown", "Unknown", 
                                            ifelse(metadata$homelessGrp3=="No" & metadata$unstableHousingGrp3=="Unknown", "No",
                                                   ifelse(metadata$homelessGrp3=="Unknown" & metadata$unstableHousingGrp3=="No", "No", 
                                                          ifelse(metadata$homelessGrp3=="No" & metadata$unstableHousingGrp3=="No", "No", "Check grouping")))))
#unique(metadata$homeless_unstable)
metadata$homeless_unstable <- factor(metadata$homeless_unstable, levels=c("Yes", "No", "Unknown"))
#table(metadata$homeless_unstable, metadata$homelessGrp3)
#table(metadata$homeless_unstable, metadata$unstableHousingGrp3)




##### Create race/ethnicity variable #####
#unique(metadata$hispanic)
metadata$hispanicGrp3 <- ifelse(metadata$hispanic=="", "Unknown", metadata$hispanic)
#table(metadata$hispanicGrp3, metadata$hispanic)


#unique(metadata$race)
metadata$raceGrp7 <- ifelse(metadata$race=="", "Unknown",
                            ifelse(metadata$race %in% c("Native Hawaiian/Pacific Islander", "American Indian/Alaska Native"), "Other", metadata$race))
#table(metadata$raceGrp7, metadata$race)


metadata$raceEth <- ifelse(metadata$hispanicGrp3=="Unknown", metadata$raceGrp7,
                           ifelse(metadata$hispanicGrp3=="No", metadata$raceGrp7,
                                  ifelse(metadata$hispanicGrp3=="Yes", "Hispanic/Latinx", "Check grouping")))
#unique(metadata$raceEth)
metadata$raceEth <- factor(metadata$raceEth, levels=c("White", "Hispanic/Latinx", "Black/African American", "Asian", "Multirace", "Other", "Unknown"))
#table(metadata$raceEth, metadata$hispanicGrp3)
#table(metadata$raceEth, metadata$raceGrp7)




##### MSM #####
#unique(metadata$MSM)
metadata$msmGrp3 <- ifelse(metadata$MSM=="YES", "Yes",
                           ifelse(metadata$MSM=="No", "No",
                                  ifelse(metadata$MSM %in% c("Unknown", ""), "Unknown", "Check grouping")))
#table(metadata$msmGrp3, metadata$MSM)




##### Sequenced & outbreak status #####
# checking to make sure sequenced and outbreak are not missing or anything weird
# the only missing values are those for NH, so filling those in the best I can
#unique(metadata$sequenced)
#nrow(subset(metadata, sequenced=="" & State.of.Residence=="MA"))
#nrow(subset(metadata, sequenced=="" & State.of.Residence=="NH"))
#nrow(subset(metadata, sequenced==""))
# All of the samples missing whether or not they were sequenced are from NH. All of these samples were sent for sequencing by the very fact that they are here in MA. So changing to a yes.
metadata$sequenced[metadata$sequenced==""] <- "Yes"
#unique(metadata$sequenced)
metadata$sequenced <- factor(metadata$sequenced, levels=c("Yes", "No"))

#unique(metadata$Sample_Passed_QC)
#nrow(subset(metadata, Sample_Passed_QC=="" & State.of.Residence=="MA"))
#nrow(subset(metadata, sequenced=="Yes" & Sample_Passed_QC=="No" & State.of.Residence=="MA"))
#nrow(subset(metadata, sequenced=="Yes" & Sample_Passed_QC=="" & State.of.Residence=="MA"))
#nrow(subset(metadata, Sample_Passed_QC=="" & State.of.Residence=="NH"))
#nrow(subset(metadata, Sample_Passed_QC==""))


# all cases missing QC are not sequenced, except for 8 cases which were sent for sequencing but failed sequencing.
metadata$Sample_Passed_QC[metadata$Sample_Passed_QC=="" & metadata$sequenced=="Yes"] <- "No"
#nrow(subset(metadata, sequenced=="Yes" & Sample_Passed_QC=="No" & State.of.Residence=="MA"))
#unique(metadata$Sample_Passed_QC)
metadata$Sample_Passed_QC[metadata$Sample_Passed_QC==""] <- "Not sequenced"
#unique(metadata$Sample_Passed_QC)
metadata$Sample_Passed_QC <- factor(metadata$Sample_Passed_QC, levels=c("Yes", "No", "Not sequenced"))

#unique(metadata$outbreak)
#nrow(subset(metadata, outbreak=="" & State.of.Residence=="MA"))
#nrow(subset(metadata, outbreak=="" & State.of.Residence=="NH"))
#nrow(subset(metadata, outbreak==""))
# all the cases missing outbreak information are NH cases; we can assume they were outbreak cases, but to avoid assumptions changing to label them as unknown
metadata$outbreak[metadata$outbreak==""] <- "Unknown - NH Case"
#unique(metadata$outbreak)
metadata$outbreak <- factor(metadata$outbreak, levels=c("Yes", "No", "Unknown - NH Case"))


##### Remove duplicates from metadata summaries #####
# manually removing the duplicates we don't want by setting a dropped duplicates flag
metadata$droppedDups <- ifelse(!(metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0265/2018"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0305/2018"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0128/2018"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0142/2018"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0290/2019"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0318/2019"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0324/2019"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0328/2019"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0331/2019"|
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0333/2019"|
                                   # added on 7/18/2025 as it has come along with the non-outbreak cases. 
                                   metadata$MASPHL_Sequencing_ID=="HepA/USA/MA/0371/2024"), 0, 1)


# # double checking we have no duplicates
# metadata %>%
#  subset(droppedDups==0) %>%
#  subset(!(is.na(caseid))) %>%
#  subset(duplicated(caseid)) %>%
#  nrow()


##### Create outbreak/sequenced concat. variable #####
# create a group that gives outbreak and sequenced status
metadata$outSeq <- paste0(metadata$outbreak, "-", metadata$Sample_Passed_QC)
#unique(metadata$outSeq)
metadata$outSeq[metadata$outSeq=="Yes-Not sequenced"] <- "Yes-No"
metadata$outSeq[metadata$outSeq=="No-Not sequenced"] <- "No-No"
metadata$outSeq[metadata$outSeq=="Unknown - NH Case-No"] <- "Likely-No"
metadata$outSeq[metadata$outSeq=="Unknown - NH Case-Yes"] <- "Likely-Yes"
#unique(metadata$outSeq)
metadata$outSeq <- factor(metadata$outSeq, levels=c("Yes-Yes", "Yes-No", "No-Yes", "No-No", "Likely-No", "Likely-Yes"))



##### Factoring necessary variables #####
# for plotting and keeping things in good order, making factors of a few variables
metadata$above18 <- ifelse(metadata$age >= 18, 1, 0)




##### Analysis date #####
# setting the analysis date to be the specimen collection date except if the case is from NH, in which case we have to sue the event date
metadata$analysis_date <- as.Date(ifelse(metadata$State.of.Residence != "NH", metadata$event_date, metadata$Specimen_Collection_Date))




##### Merging in MMWR week with metadata #####
metadata <- merge(metadata,
                  full_week_long,
                  by.x="analysis_date",
                  by.y="dates",
                  all.x=TRUE) %>%
  arrange(analysis_date)

# get week and year variables
metadata$week <- as.numeric(ifelse(nchar(metadata$wkYr)==6, 
                                   substr(metadata$wkYr, 1, 1),
                                   substr(metadata$wkYr, 1, 2)))
metadata$year <- as.numeric(ifelse(nchar(metadata$wkYr)==6, 
                                   substr(metadata$wkYr, 3, 6),
                                   substr(metadata$wkYr, 4, 7)))




##### Merge in town shape and population #####
metadata <- merge(metadata,
                  subset(town_pop, select=-Year),
                  by.x="cityResidence",
                  by.y="town",
                  all.x=T)

# column for organizing table 1
metadata$tbl1_outbreak_yr <- ifelse(metadata$analysis_date < "2021-01-01" & 
                                 !(is.na(metadata$analysis_date)) & 
                                 metadata$outbreak=="Yes", "2018-2020",
                               ifelse(metadata$analysis_date < "2021-01-01" & 
                                        !(is.na(metadata$analysis_date)) & 
                                        metadata$outbreak=="No" &
                                        metadata$Sample_Passed_QC=="Yes", "2018-2020", 
                                      ifelse("2023-01-01" <= metadata$analysis_date & 
                                               metadata$analysis_date <= "2024-12-31"  & 
                                               !(is.na(metadata$analysis_date)) & 
                                               metadata$outbreak=="Yes", "2023-2024",
                                             ifelse("2023-01-01" <= metadata$analysis_date & 
                                                      metadata$analysis_date <= "2024-12-31"  & 
                                                      !(is.na(metadata$analysis_date)) & 
                                                      metadata$outbreak=="No" &
                                                      metadata$Sample_Passed_QC=="Yes", "2023-2024",
                                                    ifelse(is.na(metadata$analysis_date), "Missing", 
                                                           ifelse(metadata$outbreak=="No", "Non-outbreak", "Check grouping"))))))




##### Weekly NH sequenced cases #####
nhSeqCases_week <- metadata %>%
  subset(droppedDups==0 & Sample_Passed_QC=="Yes" & State.of.Residence=="NH") %>%
  group_by(wkYr) %>%
  summarize(nhSeqWkCaseCt = n())


##### Weekly MA case counts (no outbreak or sequenced status included) #####
cases_week_all <- metadata %>%
  subset(droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(wkYr) %>%
  summarize(totalWkCaseCt = n())

##### Weekly MA cases by outbreak #####
cases_week_out_long <- metadata %>%
  subset(droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(wkYr, outbreak) %>%
  summarize(outWkCaseCt = n()) %>%
  tidyr::complete(outbreak, fill = list(outWkCaseCt = 0))

# go from long to wide
cases_week_out_wide <- cases_week_out_long %>%
  subset(outbreak != "") %>%
  tidyr::pivot_wider(names_from=outbreak,
                     values_from=outWkCaseCt) %>%
  rename("outbreakYes"="Yes",
         "outbreakNo"="No")


##### Weekly MA cases by sequencing #####
cases_week_seq_long <- metadata %>%
  subset(droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(wkYr, Sample_Passed_QC) %>%
  summarize(seqWkCaseCt = n())%>%
  tidyr::complete(Sample_Passed_QC, fill = list(seqWkCaseCt = 0))

# go from long to wide
cases_week_seq_wide <- cases_week_seq_long %>%
  tidyr::pivot_wider(names_from=Sample_Passed_QC,
                     values_from=seqWkCaseCt) %>%
  rename("sequenceYes"="Yes",
         "sequenceNo"="No")


##### Weekly MA cases by outbreak-sequenced combinations #####
cases_week_out_seq_long <- metadata %>%
  subset(droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(wkYr, outSeq) %>%
  summarize(outSeqWkCaseCt = n()) %>%
  tidyr::complete(outSeq, fill = list(outSeqWkCaseCt = 0))

# go from long to wide
cases_week_out_seq_wide <- cases_week_out_seq_long %>%
  tidyr::pivot_wider(names_from=outSeq,
                     values_from=outSeqWkCaseCt) %>%
  rename("outYesseqYes"="Yes-Yes",
         "outYesseqNo"="Yes-No",
         "outNoseqYes"="No-Yes",
         "outNoseqNo"="No-No")


##### Merging in weekly case summaries #####
cases_week_all <- merge(cases_week_all,
                        cases_week_out_wide,
                        by.all="wkYr",
                        all=T)
cases_week_all <- merge(cases_week_all,
                        cases_week_seq_wide,
                        by.all="wkYr",
                        all=T)
cases_week_all <- merge(cases_week_all,
                        nhSeqCases_week,
                        by.all="wkYr",
                        all=T)
cases_week_all <- merge(cases_week_all,
                        cases_week_out_seq_wide,
                        by.all="wkYr",
                        all=T)

##### Merge in MMWR weeks #####
cases_week_all <- merge(cases_week_all,
                        subset(mmwr_key, select=c(wkYr, endDt)),
                        by.all="wkYr",
                        all=T)

# if a case count is missing (this happens when merging if a week did not have an outbreak or sequenced case, but did have a case that week), fill in with a 0
cases_week_all[is.na(cases_week_all)] <- 0


range(cases_week_all$wkYr, na.rm=TRUE)
# [1] "0"      "9_2025"


# remove first junk row
nh_outbreak <- subset(nh_outbreak, `...2`!="Count")

# rename columns
colnames(nh_outbreak) <- c("collectionDt", "nhCaseCount")

# format variables
nh_outbreak$collectionDt <- as.numeric(nh_outbreak$collectionDt)
nh_outbreak$collectionDt <- as.Date(nh_outbreak$collectionDt, format="%Y-%m-%d", origin = "1899-12-30")
nh_outbreak$nhCaseCount <- as.numeric(nh_outbreak$nhCaseCount)

# merge in MMWR weeks
nh_outbreak <- merge(nh_outbreak,
                     full_week_long,
                     by.x="collectionDt",
                     by.y="dates",
                     all.x=T)

# get weekly counts of total NH detected outbreak cases
nh_outbreak_summ <- nh_outbreak %>%
  group_by(wkYr) %>%
  summarize(nhTotWkCaseCt = sum(nhCaseCount))

# merge in with total weekly count dataframe
cases_week_all <- merge(cases_week_all,
                        nh_outbreak_summ,
                        by="wkYr",
                        all=T)

# add 0's for NH data with missing data post-merge (i.e., weeks where no cases were detected)
cases_week_all$nhTotWkCaseCt[is.na(cases_week_all$nhTotWkCaseCt)] <- 0




##### Subset the data to proper date range #####
# This step subsets the immune data to just the dates where we have case data. The vaccine data starts on 1_2018 and the case data starts on 3_2018, so I kept the min immune date so that when we calculate 4-week rolling averages, we have a value by 4_2018.
immune <- immune %>%
  select(BIN:Shots_2024_47)




##### Rename oddly formatted variables
colnames(immune)[colnames(immune)=="Shots_2022_52...211"] <- "Shots_2021_52"
colnames(immune)[colnames(immune)=="Shots_2022_52...263"] <- "Shots_2022_52"





##### Go to long format #####
immune_long <- immune %>%
  tidyr::pivot_longer(cols=3:ncol(immune),
                      names_to="mmwrweek",
                      values_to="vaccine_count")




##### Formatting variables #####
# Get week & year numbers
# if the string is 12 characters long, it's a single digit month; if it's greater than 12, than it is a double digit month
immune_long$week <- as.numeric(ifelse(nchar(immune_long$mmwrweek)==12, 
                                      substr(immune_long$mmwrweek, 12, 12), 
                                      substr(immune_long$mmwrweek, 12, 13)))

# get year variable from old column names
immune_long$year <- as.numeric(substr(immune_long$mmwrweek, 7, 10))


# the data pull from MIIS names all categories of interest in a variable called CITY, so renaming so it makes sense
immune_long$binCat <- immune_long$CITY


# formatting a wkYr variable
immune_long$wkYr <- paste0(as.character(immune_long$week), "_", as.character(immune_long$year))




##### Merge in week date data #####
immune_long <- merge(immune_long,
                     subset(mmwr_key, select=c(endDt, `MMWR Week`, Year)),
                     by.x=c("week", "year"),
                     by.y=c("MMWR Week", "Year"),
                     all.x=TRUE)




##### Weekly vaccines administered #####
# I don't have a total for each week, just totals for each of the categories for each week (e.g., residence, provider location, age group, race/ethnicity, etc.). So I get the totals in each category and use the max value in each week across all categories as the total vaccines given. Provider location appears to have some data missing (e.g., a few missing data), but all of the other categories appear to have the same value. The max value for any week would represent that week's total vaccines administered count.
immune_weekly_total <- immune_long %>%
  group_by(year, week, BIN) %>%
  summarize(vax_weekly_count = sum(vaccine_count))


# now get the max value from the different vaccine descriptor categories
immune_weekly_total <- immune_weekly_total %>%
  group_by(year, week) %>%
  summarize(vax_weekly_count = max(vax_weekly_count))


# merge in MMWR weeks
immune_weekly_total <- merge(immune_weekly_total,
                             subset(mmwr_key, select=c(wkYr, `MMWR Week`, Year, startDt, endDt)),
                             by.x=c("year", "week"),
                             by.y=c("Year", "MMWR Week"),
                             all.x=TRUE)


# merge in population totals
immune_weekly_total <- merge(immune_weekly_total,
                             state_pop,
                             by.x="year",
                             by.y="Year",
                             all.x=TRUE) %>%
  arrange(endDt)


# get vaccination rates per 100,000
immune_weekly_total$vax_per_100000 <- immune_weekly_total$vax_weekly_count/(immune_weekly_total$statePop/100000)


# calculate a 4-week rolling vaccination rate average
immune_weekly_total$rollingVaxPer100000 <- ma(immune_weekly_total$vax_per_100000, n=4)




##### Weekly vaccines administered by provider #####
# get summary stats for immunizations by provider type (see immune_sum function at start of document, but just gives summary stat counts by week by the vaccine descriptor of interest)
immune_weekly_provType <- immune_sum(category="provider_type")


# choosing a the long format output; immune_sum() has the option of seeing the data in long format, df[[1]], or wide format, df[[2]]
immune_weekly_provType_long <- immune_weekly_provType[[1]]


# renaming provider type categories for display - creating two, as not fully sure what grouping we want yet, The first (here) is 9 categories for provider type (the next will be 5 groups)
#unique(immune_weekly_provType_long$binCat)
immune_weekly_provType_long$provTypeGrp9 <- ifelse(immune_weekly_provType_long$binCat %in% c("Assisted Living/Adult Day Care",
                                                                                             "LTCF/Nursing Home/Rest Home", 
                                                                                             "VNA", 
                                                                                             "LTC Public"), 
                                                   "LTCF", 
                                                   ifelse(immune_weekly_provType_long$binCat %in% c("Board of Health/Health Dept", 
                                                                                                    "Dept of Mental Health", 
                                                                                                    "State Agency", 
                                                                                                    "Council on Aging"), 
                                                          "State Agency",
                                                          ifelse(immune_weekly_provType_long$binCat %in% c("Community Health Center", 
                                                                                                           "Free Clinic",
                                                                                                           "STD walk in clinic"), 
                                                                 "Community Health Center",
                                                                 ifelse(immune_weekly_provType_long$binCat %in% c("Other (Private)",
                                                                                                                  "Other (Public)", 
                                                                                                                  "Out Of State IIS (RI/NY)", 
                                                                                                                  "Home Health Agency", 
                                                                                                                  "Employee Health"), 
                                                                        "Other",
                                                                        ifelse(immune_weekly_provType_long$binCat %in% c("College (Private)", 
                                                                                                                         "College (Public)"), 
                                                                               "College",
                                                                               ifelse(immune_weekly_provType_long$binCat %in% c("Family Planning", 
                                                                                                                                "Family Practice",
                                                                                                                                "Hospital (Private)", 
                                                                                                                                "Hospital(Public)", 
                                                                                                                                "Internal Medicine", 
                                                                                                                                "Multi-Specialty Center", 
                                                                                                                                "OB/GYN", 
                                                                                                                                "Specialty Practice", 
                                                                                                                                "Urgent Care Centers"), 
                                                                                      "Healthcare Center",
                                                                                      ifelse(immune_weekly_provType_long$binCat %in% c("Pediatric Practice", 
                                                                                                                                       "School (Public)", 
                                                                                                                                       "School (Special Education)", 
                                                                                                                                       "School Based Health Center"), 
                                                                                             "Pediatric",
                                                                                             ifelse(immune_weekly_provType_long$binCat=="Correctional Facility", 
                                                                                                    "Correctional Facility",
                                                                                                    ifelse(immune_weekly_provType_long$binCat=="Commercial Pharmacy", 
                                                                                                           "Commercial Pharmacy", "Check")))))))))


# now the 5 category provider type variable
immune_weekly_provType_long$provTypeGrp5 <- ifelse(immune_weekly_provType_long$binCat %in% c("Assisted Living/Adult Day Care", 
                                                                                             "LTCF/Nursing Home/Rest Home", 
                                                                                             "VNA", 
                                                                                             "LTC Public",
                                                                                             "Board of Health/Health Dept",
                                                                                             "Dept of Mental Health", 
                                                                                             "State Agency", 
                                                                                             "Council on Aging",
                                                                                             "Other (Private)",
                                                                                             "Other (Public)", 
                                                                                             "Out Of State IIS (RI/NY)",
                                                                                             "Home Health Agency", 
                                                                                             "Employee Health",
                                                                                             "College (Private)", 
                                                                                             "College (Public)",
                                                                                             "Correctional Facility"), 
                                                   "Other",
                                                   ifelse(immune_weekly_provType_long$binCat %in% c("Community Health Center", 
                                                                                                    "Free Clinic",
                                                                                                    "STD walk in clinic"), 
                                                          "Community Health Center",
                                                          ifelse(immune_weekly_provType_long$binCat %in% c("Family Planning", 
                                                                                                           "Family Practice", 
                                                                                                           "Hospital (Private)", 
                                                                                                           "Hospital(Public)",
                                                                                                           "Internal Medicine",
                                                                                                           "Multi-Specialty Center",
                                                                                                           "OB/GYN", 
                                                                                                           "Specialty Practice", 
                                                                                                           "Urgent Care Centers"), 
                                                                 "Healthcare Center",
                                                                 ifelse(immune_weekly_provType_long$binCat %in% c("Pediatric Practice", 
                                                                                                                  "School (Public)", 
                                                                                                                  "School (Special Education)", 
                                                                                                                  "School Based Health Center"), 
                                                                        "Pediatric",
                                                                        ifelse(immune_weekly_provType_long$binCat=="Commercial Pharmacy", 
                                                                               "Commercial Pharmacy", "Check")))))


# Now resummarizing to the 9 provider categories
immune_weekly_provType_long9 <- immune_weekly_provType_long %>%
  group_by(year, week, provTypeGrp9) %>%
  summarize(vax_weekly_provType_count = sum(wkImmuneCt))


# merge in MMWR weeks
immune_weekly_provType_long9 <- merge(immune_weekly_provType_long9,
                                      subset(mmwr_key, select=c(wkYr, `MMWR Week`, Year, startDt, endDt)),
                                      by.x=c("year", "week"),
                                      by.y=c("Year", "MMWR Week"),
                                      all.x=TRUE)


# go wide now so that we can unify with other datasets
immune_weekly_provType_wide9 <- subset(immune_weekly_provType_long9, select=c(wkYr, endDt, vax_weekly_provType_count, provTypeGrp9)) %>%
  tidyr::pivot_wider(id_cols=c(wkYr, endDt),
                     names_from=provTypeGrp9,
                     values_from=vax_weekly_provType_count) %>%
  arrange(endDt)


# rename columns, starting at the 3rd column, and just add a prefix so we know what these counts are
immune_weekly_provType_wide9 <- immune_weekly_provType_wide9 %>%
  rename_with(.fn=~paste0("wkImmuneCt_prov", .), .cols=all_of(colnames(immune_weekly_provType_wide9[,3:ncol(immune_weekly_provType_wide9)])))


# get 4-week rolling averages for all of the vaccine counts for each provider type
immune_weekly_provType_wide9$rollingWkImmCt_provCollege <- ma(immune_weekly_provType_wide9$wkImmuneCt_provCollege, n=4)
immune_weekly_provType_wide9$rollingWkImmCt_provCommPharm <- ma(immune_weekly_provType_wide9$`wkImmuneCt_provCommercial Pharmacy`, n=4)
immune_weekly_provType_wide9$rollingWkImmCt_provCommHealthCent <- ma(immune_weekly_provType_wide9$`wkImmuneCt_provCommunity Health Center`, n=4)
immune_weekly_provType_wide9$rollingWkImmCt_provCorrFac <- ma(immune_weekly_provType_wide9$`wkImmuneCt_provCorrectional Facility`, n=4)
immune_weekly_provType_wide9$rollingWkImmCt_provHealtCent <- ma(immune_weekly_provType_wide9$`wkImmuneCt_provHealthcare Center`, n=4)
immune_weekly_provType_wide9$rollingWkImmCt_provLTCF <- ma(immune_weekly_provType_wide9$wkImmuneCt_provLTCF, n=4)
immune_weekly_provType_wide9$rollingWkImmCt_provOther <- ma(immune_weekly_provType_wide9$wkImmuneCt_provOther, n=4)
immune_weekly_provType_wide9$rollingWkImmCt_provPediatric <- ma(immune_weekly_provType_wide9$wkImmuneCt_provPediatric, n=4)
immune_weekly_provType_wide9$rollingWkImmCt_provStateAgency <- ma(immune_weekly_provType_wide9$`wkImmuneCt_provState Agency`, n=4)


# repeating everything above but with the 5 categories of provider
immune_weekly_provType_long5 <- immune_weekly_provType_long %>%
  group_by(year, week, provTypeGrp5) %>%
  summarize(vax_weekly_provType2_count = sum(wkImmuneCt))


# merge in MMWR weeks
immune_weekly_provType_long5 <- merge(immune_weekly_provType_long5,
                                      subset(mmwr_key, select=c(wkYr, `MMWR Week`, Year, startDt, endDt)),
                                      by.x=c("year", "week"),
                                      by.y=c("Year", "MMWR Week"),
                                      all.x=TRUE)


# go wide now so that we can unify with other datasets
immune_weekly_provType_wide5 <- subset(immune_weekly_provType_long5, select=c(wkYr, endDt, vax_weekly_provType2_count, provTypeGrp5)) %>%
  tidyr::pivot_wider(id_cols=c(wkYr, endDt),
                     names_from=provTypeGrp5,
                     values_from=vax_weekly_provType2_count) %>%
  arrange(endDt)


# rename columns
immune_weekly_provType_wide5 <- immune_weekly_provType_wide5 %>%
  # didn't choose a great prefix for these but don't want to change it now and have to alter the analysis code potentially.
  rename_with(.fn=~paste0("wkImmuneCt_prov2", .), .cols=all_of(colnames(immune_weekly_provType_wide5[,3:ncol(immune_weekly_provType_wide5)])))


# get 4-week rolling averages for all of the vaccine counts for each provider type
immune_weekly_provType_wide5$rollingWkImmCt_prov2Pediatric <- ma(immune_weekly_provType_wide5$wkImmuneCt_prov2Pediatric, n=4)
immune_weekly_provType_wide5$rollingWkImmCt_prov2CommPharm <- ma(immune_weekly_provType_wide5$`wkImmuneCt_prov2Commercial Pharmacy`, n=4)
immune_weekly_provType_wide5$rollingWkImmCt_prov2CommHealthCent <- ma(immune_weekly_provType_wide5$`wkImmuneCt_prov2Community Health Center`, n=4)
immune_weekly_provType_wide5$rollingWkImmCt_prov2Other <- ma(immune_weekly_provType_wide5$wkImmuneCt_prov2Other, n=4)
immune_weekly_provType_wide5$rollingWkImmCt_prov2HealtCent <- ma(immune_weekly_provType_wide5$`wkImmuneCt_prov2Healthcare Center`, n=4)




##### Weekly vaccines administered by vaccine type #####
immune_weekly_vax <- immune_sum(category="cvx_code")


# choosing the wide format this time, as opposed to the long format chosen for the provider type vax rate summary tables. That's because I had to recombine provider categories before getting the final rates
immune_weekly_vax_wide <- immune_weekly_vax[[2]]


# get 4-week rolling count averages
immune_weekly_vax_wide$rollingWkImmCt_hepAB <- ma(immune_weekly_vax_wide$`wkImmuneCt_cvx_code104 HepA-HepB`, n=4)
immune_weekly_vax_wide$rollingWkImmCt_hepAPed <- ma(immune_weekly_vax_wide$`wkImmuneCt_cvx_code31 & 83 - HepA  Pediatric doses`, n=4)
immune_weekly_vax_wide$rollingWkImmCt_hepA2 <- ma(immune_weekly_vax_wide$`wkImmuneCt_cvx_code52 HepA 2-dose`, n=4)
immune_weekly_vax_wide$rollingWkImmCt_hepAUnk <- ma(immune_weekly_vax_wide$`wkImmuneCt_cvx_code85 HepA Unspecified`, n=4)


# making another version of the vaccine categories where I combine pediatric and unknown doses
# choosing an immune sum output (long format summary stats for weekly vax rate by provider)
immune_weekly_vax_long <- immune_weekly_vax[[1]]


# renaming provider type categories for display
immune_weekly_vax_long$vaxTypeGrp3 <- ifelse(immune_weekly_vax_long$binCat %in% c("31 & 83 - HepA  Pediatric doses", 
                                                                                  "85 HepA Unspecified"), 
                                             "Other", immune_weekly_vax_long$binCat)


# getting counts again based on 3 category vaccine grouping
immune_weekly_vax_long3 <- immune_weekly_vax_long %>%
  group_by(year, week, vaxTypeGrp3) %>%
  summarize(vax_weekly_vax_count = sum(wkImmuneCt))


# merge in MMWR weeks
immune_weekly_vax_long3 <- merge(immune_weekly_vax_long3,
                                 subset(mmwr_key, select=c(wkYr, `MMWR Week`, Year, startDt, endDt)),
                                 by.x=c("year", "week"),
                                 by.y=c("Year", "MMWR Week"),
                                 all.x=TRUE)


# go wide now so that we can unify with other datasets
immune_weekly_vax_wide3 <- subset(immune_weekly_vax_long3, select=c(wkYr, endDt, vax_weekly_vax_count, vaxTypeGrp3)) %>%
  tidyr::pivot_wider(id_cols=c(wkYr, endDt),
                     names_from=vaxTypeGrp3,
                     values_from=vax_weekly_vax_count) %>%
  arrange(endDt)


# rename columns
immune_weekly_vax_wide3 <- immune_weekly_vax_wide3 %>%
  # again, better names should have been chosen. Original vaccine type categories are the cvx prefixes, and the vx prefixes are for the newly made vaccine categories.
  rename_with(.fn=~paste0("wkImmuneCt_vax", .), .cols=all_of(colnames(immune_weekly_vax_wide3[,3:ncol(immune_weekly_vax_wide3)])))


# get 4-week rolling averages for all of the vaccine counts for each provider type
immune_weekly_vax_wide3$rollingWkImmCt_vaxAB <- ma(immune_weekly_vax_wide3$`wkImmuneCt_vax104 HepA-HepB`, n=4)
immune_weekly_vax_wide3$rollingWkImmCt_vax2Dose <- ma(immune_weekly_vax_wide3$`wkImmuneCt_vax52 HepA 2-dose`, n=4)
immune_weekly_vax_wide3$rollingWkImmCt_vaxOther <- ma(immune_weekly_vax_wide3$wkImmuneCt_vaxOther, n=4)




##### Merge everything all together, including cases #####
# create a list of dataframes we want to combine
df_list <- list(cases_week_all,
                immune_weekly_provType_wide9,
                immune_weekly_provType_wide5,
                immune_weekly_vax_wide,
                immune_weekly_vax_wide3,
                subset(immune_weekly_total, select=c(endDt, wkYr, rollingVaxPer100000)))


# now join all of the dataframes from the list
cases_week_all <- df_list %>%
  purrr::reduce(full_join, by=c("endDt", "wkYr"))


# arrange by date
cases_week_all <- cases_week_all %>%
  arrange(endDt)


# vaccine data extends past our study date, so selecting the last date we want included so we don't just plot vaccine data with no case data
cases_week_all <- subset(cases_week_all, endDt < "2024-11-30")
range(cases_week_all$wkYr, na.rm=TRUE)
# [1] "0"      "9_2025"

# "0" will get inappropriately coerced to the Unix epoch date "1970-01-01"
# need to remove the wkYr=0
cases_week_all <- cases_week_all %>%
  filter(wkYr != "0")

range(cases_week_all$wkYr, na.rm=TRUE)
#[1] "1_2018" "9_2024"


##########################################
##### Generate categorical summaries #####
##########################################
##### Gender #####
##### 2018-2020
# 2018-2020 outbreak cases
gender_1820_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                      genderGrp3, 
                                      "Gender", 
                                      columnnames="Outbreak cases: % (N=X)")

# 2018-2020 sequenced cases
gender_1820_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                       genderGrp3, 
                                       "Gender", 
                                       columnnames="Sequenced cases: % (N=X)")


##### 2023-2024
# 2023-2024 outbreak cases
gender_2324_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                      genderGrp3, 
                                      "Gender", 
                                      columnnames="Outbreak cases: % (N=X)")

# 2023-2024 sequenced cases
gender_2324_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                       genderGrp3, 
                                       "Gender", 
                                       columnnames="Sequenced cases: % (N=X)")


##### All
gender_all <- summary_table_char(subset(metadata, droppedDups==0 & State.of.Residence=="MA"),
                                 genderGrp3, 
                                 "Gender", 
                                 columnnames="All cases: % (N=X)")



##### Age #####
# Note - not using the summary_table_char() function for this one because I want to create some unique rows for age (e.g., median and IQR)

##### 2018-2020
# 2018-2020 outbreak cases
age_1820_out <- metadata %>%
  subset(tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(ageGrp5, .drop=FALSE) %>%
  count()
# calculate the percentage for each age category
age_1820_out$`%` <- round(age_1820_out$n/sum(age_1820_out$n)*100, 1)
# get the proper label (% (N=X))
age_1820_out$`Outbreak cases: % (N=X)` <- paste0(age_1820_out$`%`, "% (", age_1820_out$n, ")")
# select the desired columns
age_1820_out <- subset(age_1820_out, select=c(ageGrp5, `Outbreak cases: % (N=X)`))
# rename columns
colnames(age_1820_out) <- c("Variable", "Outbreak cases: % (N=X)")
# in order to build up the table to display the way I want, I have to convert this variable to a character from a factor so that I can create a header column (the header isn't apart of the factor levels, and it was too complicated to re-level so this was a simpler fix for me).
age_1820_out$Variable <- as.character(age_1820_out$Variable)
# Now put together the header columns and median
age_1820_out <- rbind.data.frame(c("Age group (yrs)", ""), c("Median (IQR):", paste0(median(subset(metadata, tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T), " (", quantile(subset(metadata, tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[2], "-", quantile(subset(metadata, tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[4], ")")), age_1820_out)

# 2018-2020 sequenced cases
age_1820_seq2 <- metadata %>%
  subset(tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(ageGrp5, .drop=FALSE) %>%
  count()
age_1820_seq2$`%` <- round(age_1820_seq2$n/sum(age_1820_seq2$n)*100, 1)
age_1820_seq2$`Sequenced cases: % (N=X)` <- paste0(age_1820_seq2$`%`, "% (", age_1820_seq2$n, ")")
age_1820_seq2 <- subset(age_1820_seq2, select=c(ageGrp5, `Sequenced cases: % (N=X)`))
colnames(age_1820_seq2) <- c("Variable", "Sequenced cases: % (N=X)")
age_1820_seq2$Variable <- as.character(age_1820_seq2$Variable)
age_1820_seq2 <- rbind.data.frame(cbind.data.frame("Variable"="Age group (yrs)", `Sequenced cases: % (N=X)`=""), cbind.data.frame("Variable"="Median (IQR):", `Sequenced cases: % (N=X)`=paste0(median(subset(metadata, tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T), " (", quantile(subset(metadata, tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[2], "-", quantile(subset(metadata, tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[4], ")")), age_1820_seq2)


##### 2023-2024
# 2023-2024 outbreak cases
age_2324_out <- metadata %>%
  subset(tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(ageGrp5, .drop=FALSE) %>%
  count()
age_2324_out$`%` <- round(age_2324_out$n/sum(age_2324_out$n)*100, 1)
age_2324_out$`Outbreak cases: % (N=X)` <- paste0(age_2324_out$`%`, "% (", age_2324_out$n, ")")
age_2324_out <- subset(age_2324_out, select=c(ageGrp5, `Outbreak cases: % (N=X)`))
colnames(age_2324_out) <- c("Variable", "Outbreak cases: % (N=X)")
age_2324_out$Variable <- as.character(age_2324_out$Variable)
age_2324_out <- rbind.data.frame(cbind("Variable"="Age group (yrs)", `Outbreak cases: % (N=X)`=""), cbind("Variable"="Median (IQR):", `Outbreak cases: % (N=X)`=paste0(median(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T), " (", quantile(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[2], "-", quantile(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[4], ")")), age_2324_out)

# 2023-2024 sequenced cases
age_2324_seq2 <- metadata %>%
  subset(tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(ageGrp5, .drop=FALSE) %>%
  count()
age_2324_seq2$`%` <- round(age_2324_seq2$n/sum(age_2324_seq2$n)*100, 1)
age_2324_seq2$`Sequenced cases: % (N=X)` <- paste0(age_2324_seq2$`%`, "% (", age_2324_seq2$n, ")")
age_2324_seq2 <- subset(age_2324_seq2, select=c(ageGrp5, `Sequenced cases: % (N=X)`))
colnames(age_2324_seq2) <- c("Variable", "Sequenced cases: % (N=X)")
age_2324_seq2$Variable <- as.character(age_2324_seq2$Variable)
age_2324_seq2 <- rbind.data.frame(cbind.data.frame("Variable"="Age group (yrs)", `Sequenced cases: % (N=X)`=""), cbind.data.frame("Variable"="Median (IQR):", `Sequenced cases: % (N=X)`=paste0(median(subset(metadata, tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T), " (", quantile(subset(metadata, tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[2], "-", quantile(subset(metadata, tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[4], ")")), age_2324_seq2)


##### All
age_all <- metadata %>%
  subset(droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(ageGrp5, .drop=FALSE) %>%
  count()
age_all$`%` <- round(age_all$n/sum(age_all$n)*100, 1)
age_all$`All cases: % (N=X)` <- paste0(age_all$`%`, "% (", age_all$n, ")")
age_all <- subset(age_all, select=c(ageGrp5, `All cases: % (N=X)`))
colnames(age_all) <- c("Variable", "All cases: % (N=X)")
age_all$Variable <- as.character(age_all$Variable)
age_all <- rbind.data.frame(cbind.data.frame("Variable"="Age group (yrs)", `All cases: % (N=X)`=""), cbind.data.frame("Variable"="Median (IQR):", `All cases: % (N=X)`=paste0(median(subset(metadata, droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T), " (", quantile(subset(metadata, droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[2], "-", quantile(subset(metadata, droppedDups==0 & State.of.Residence=="MA", select=age)[,1], na.rm=T)[4], ")")), age_all)




##### Race/ethnicity #####
##### 2018-2020
# 2018-2020 outbreak cases
raceEth_1820_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                       raceEth, 
                                       "Race/ethnicity combined", 
                                       columnnames="Outbreak cases: % (N=X)")

# 2018-2020 sequenced cases
raceEth_1820_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                        raceEth, 
                                        "Race/ethnicity combined", 
                                        columnnames="Sequenced cases: % (N=X)")


##### 2023-2024
# 2023-2024 outbreak cases
raceEth_2324_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                       raceEth, 
                                       "Race/ethnicity combined", 
                                       columnnames="Outbreak cases: % (N=X)")

# 2023-2024 sequenced cases
raceEth_2324_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                        raceEth, 
                                        "Race/ethnicity combined", 
                                        columnnames="Sequenced cases: % (N=X)")


##### All
raceEth_all <- summary_table_char(subset(metadata, droppedDups==0 & State.of.Residence=="MA"), 
                                  raceEth, 
                                  "Race/ethnicity combined", 
                                  columnnames="All cases: % (N=X)")





##### Illicit drug use (ever) #####
##### 2018-2020
# 2018-2020 outbreak cases
illicitEver_1820_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                           illicitEverGrp3, 
                                           "History of drug use", 
                                           columnnames="Outbreak cases: % (N=X)")


# 2018-2020 sequenced cases
illicitEver_1820_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"),
                                            illicitEverGrp3, 
                                            "History of drug use", 
                                            columnnames="Sequenced cases: % (N=X)")


##### 2023-2024
# 2023-2024 outbreak cases
illicitEver_2324_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                           illicitEverGrp3, 
                                           "History of drug use", 
                                           columnnames="Outbreak cases: % (N=X)")

# 2023-2024 sequenced cases
illicitEver_2324_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                            illicitEverGrp3, 
                                            "History of drug use", 
                                            columnnames="Sequenced cases: % (N=X)")


##### All
illicitEver_all <- summary_table_char(subset(metadata, droppedDups==0 & State.of.Residence=="MA"), 
                                      illicitEverGrp3, 
                                      "History of drug use", 
                                      columnnames="All cases: % (N=X)")




##### Unhoused #####
##### 2018-2020
# 2018-2020 outbreak cases
unhoused_final_1820_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                              homelessGrp3, 
                                              "Experiencing homelessness", 
                                              columnnames="Outbreak cases: % (N=X)")


# 2018-2020 sequenced cases
unhoused_final_1820_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                               homelessGrp3, 
                                               "Experiencing homelessness", 
                                               columnnames="Sequenced cases: % (N=X)")


##### 2023-2024
# 2023-2024 outbreak cases
unhoused_final_2324_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                              homelessGrp3, 
                                              "Experiencing homelessness", 
                                              columnnames="Outbreak cases: % (N=X)")

# 2023-2024 sequenced cases
unhoused_final_2324_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                               homelessGrp3, 
                                               "Experiencing homelessness", 
                                               columnnames="Sequenced cases: % (N=X)")


##### All
unhoused_final_all <- summary_table_char(subset(metadata, droppedDups==0 & State.of.Residence=="MA"),
                                         homelessGrp3, 
                                         "Experiencing homelessness", 
                                         columnnames="All cases: % (N=X)")




##### Vaccination status #####
##### 2018-2020
# 2018-2020 outbreak cases
vax_1820_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                   vacGrp2, 
                                   "Vaccinated", 
                                   columnnames="Outbreak cases: % (N=X)")


# 2018-2020 sequenced cases
vax_1820_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                    vacGrp2, 
                                    "Vaccinated", 
                                    columnnames="Sequenced cases: % (N=X)")


##### 2023-2024
# 2023-2024 outbreak cases
vax_2324_out <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                   vacGrp2, 
                                   "Vaccinated", 
                                   columnnames="Outbreak cases: % (N=X)")

# 2023-2024 sequenced cases
vax_2324_seq2 <- summary_table_char(subset(metadata, tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA"), 
                                    vacGrp2, 
                                    "Vaccinated", 
                                    columnnames="Sequenced cases: % (N=X)")


##### All
vax_all <- summary_table_char(subset(metadata, droppedDups==0 & State.of.Residence=="MA"),
                              vacGrp2, 
                              "Vaccinated", 
                              columnnames="All cases: % (N=X)")




#########################################
##### Combine categorical summaries #####
#########################################
##### 2018-2020 #####
##### Outbreak
# combine all of the summary tables into the desired order. Starting with 2018-2020, combine the outbreak and sequenced columns. Then do the same for the 2023-2024 outbreak. Then combine both of those outputs to get the final summary table.
summary_1820_out <- rbind.data.frame(c("Total cases", 
                                       paste0("N=", 
                                              nrow(subset(metadata, 
                                                          tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA")))), 
                                     gender_1820_out, 
                                     age_1820_out, 
                                     raceEth_1820_out, 
                                     illicitEver_1820_out,
                                     unhoused_final_1820_out, 
                                     vax_1820_out)


##### Sequenced
summary_1820_seq2 <- rbind.data.frame(c("Total cases", 
                                        paste0("N=",
                                               nrow(subset(metadata, 
                                                           tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA")))), 
                                      gender_1820_seq2, 
                                      age_1820_seq2, 
                                      raceEth_1820_seq2, 
                                      illicitEver_1820_seq2,
                                      unhoused_final_1820_seq2, 
                                      vax_1820_seq2)


##### combine 2018-2020 summary for sequenced and outbreaks
final_summary_1820_2 <- cbind(summary_1820_out, 
                              subset(summary_1820_seq2, 
                                     select=`Sequenced cases: % (N=X)`), 
                              c(paste0("N=", 
                                       nrow(subset(metadata, Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="NH"))), 
                                rep("", 
                                    times=nrow(summary_1820_out)-1)))

# correct column names
colnames(final_summary_1820_2) <- c("Variable",
                                    "Outbreak-associated MA cases: % (N=X) - 2018-2020",
                                    "Successfully sequenced MA cases: % (N=X) - 2018-2020",
                                    "Outbreak-associated NH cases: % (N=X)")

# create an rn variable to make sure we can merge to the 2023-2024 outbreak. The rn variable stands for row number and is just the row number for the desired order of the final summary table
final_summary_1820_2$rn <- seq(1, nrow(final_summary_1820_2), by=1)




##### 2023-2024 #####
##### Outbreak
summary_2324_out <- rbind.data.frame(c("Total cases", 
                                       paste0("N=",
                                              nrow(subset(metadata, tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA")))), 
                                     gender_2324_out, 
                                     age_2324_out, 
                                     raceEth_2324_out,  
                                     illicitEver_2324_out,
                                     unhoused_final_2324_out, 
                                     vax_2324_out)


##### Sequenced
summary_2324_seq2 <- rbind.data.frame(c("Total cases", 
                                        paste0("N=",
                                               nrow(subset(metadata, 
                                                           tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & droppedDups==0 & State.of.Residence=="MA")))), 
                                      gender_2324_seq2, 
                                      age_2324_seq2, 
                                      raceEth_2324_seq2,  
                                      illicitEver_2324_seq2,
                                      unhoused_final_2324_seq2, 
                                      vax_2324_seq2)


##### Combine summaries for sequenced and outbreak cases
final_summary_2324_2 <- cbind(summary_2324_out, subset(summary_2324_seq2, select=`Sequenced cases: % (N=X)`))

# rename variables
colnames(final_summary_2324_2) <- c("Variable",
                                    "Outbreak-associated MA cases: % (N=X) - 2023-2024",
                                    "Successfully sequenced MA cases: % (N=X) - 2023-2024")

# provide row number for linking
final_summary_2324_2$rn <- seq(1, nrow(final_summary_2324_2), by=1)




##### All #####
summary_all <- rbind.data.frame(c("Total cases", 
                                  paste0("N=", nrow(subset(metadata, droppedDups==0 & State.of.Residence=="MA")))), 
                                gender_all,
                                age_all, 
                                raceEth_all, 
                                illicitEver_all,
                                unhoused_final_all, 
                                vax_all)

# rename columns
colnames(summary_all) <- c("Variable",
                           "All MA cases: % (N=X)*")

# provide row number for linking
summary_all$rn <- seq(1, nrow(summary_all), by=1)




##### Merge for final tables #####
##### Outbreak and sequenced (all years)
# merge the 2018-2020 and 2023-2024 summaries
final_summary_2 <- merge(final_summary_1820_2,
                         final_summary_2324_2,
                         by.all=c("Variable", "rn"))

# now add the total summary
final_summary_2 <- merge(final_summary_2,
                         summary_all,
                         by.all=c("Variable", "rn"))

# make sure the summary table is in the proper order and then drop the rn variable
final_summary_2 <- final_summary_2 %>%
  arrange(rn) %>%
  subset(select=-rn)




##### Counts #####
correct_nh_ages <- cbind.data.frame(c("<20 years old", "20-39 years old", "40-59 years old"),
                                    c(as.numeric(subset(as.data.frame(nh_metadata), `Hepatitis A Outbreak, NH, 2018-2020`=="0-9", select=...2)[,1]) + as.numeric(subset(as.data.frame(nh_metadata), `Hepatitis A Outbreak, NH, 2018-2020`=="10-19", select=...2)[,1]),
                                      as.numeric(subset(as.data.frame(nh_metadata), `Hepatitis A Outbreak, NH, 2018-2020`=="20-29", select=...2)[,1]) + as.numeric(subset(as.data.frame(nh_metadata), `Hepatitis A Outbreak, NH, 2018-2020`=="30-39", select=...2)[,1]),
                                      as.numeric(subset(as.data.frame(nh_metadata), `Hepatitis A Outbreak, NH, 2018-2020`=="40-49", select=...2)[,1]) + as.numeric(subset(as.data.frame(nh_metadata), `Hepatitis A Outbreak, NH, 2018-2020`=="50-59", select=...2)[,1])))

colnames(correct_nh_ages) <- colnames(nh_metadata)

nh_metadata <- subset(nh_metadata, !(`Hepatitis A Outbreak, NH, 2018-2020` %in% c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "American Indian/Alaska Native", "Native Hawaiian/Pacific Islander")))

nh_metadata <- rbind(subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020` %in% c("Total Outbreak Cases", "Total Cases", "Age")), correct_nh_ages, subset(nh_metadata, !(`Hepatitis A Outbreak, NH, 2018-2020` %in% c("Total Outbreak Cases", "Total Cases", "Age"))))

nh_metadata_test <- rbind(subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Gender"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Female"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Male"),
                          c(`Hepatitis A Outbreak, NH, 2018-2020`="Unknown", ...2="0"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Age"),
                          c(`Hepatitis A Outbreak, NH, 2018-2020`="Median (IQR):", ...2=""),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="<20 years old"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="20-39 years old"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="40-59 years old"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="≥60"),
                          nh_metadata[8,],
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Race/Ethnicity"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="White"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Hispanic/Latinx"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Black/African American"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Asian"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Multirace"),
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Other"),
                          nh_metadata[19,],
                          subset(nh_metadata, `Hepatitis A Outbreak, NH, 2018-2020`=="Lifetime Illicit Drug Use"),
                          nh_metadata[22,],
                          nh_metadata[21,],
                          nh_metadata[c(23:24),],
                          nh_metadata[26,],
                          nh_metadata[25,],
                          nh_metadata[c(27:28),],
                          nh_metadata[30,],
                          c(`Hepatitis A Outbreak, NH, 2018-2020`="No/unknown", ...2=as.numeric(nh_metadata[29,2])+as.numeric(nh_metadata[31,2])))

# replace the count label with NAs
nh_metadata_test$`...2`[nh_metadata_test$`...2`=="Count"] <- ""

# get percentages and labels together
nh_metadata_test$percent <- ""
for (i in 1:nrow(nh_metadata_test)) {
  nh_metadata_test$percent[i] <- paste0(round((as.numeric(nh_metadata_test$`...2`[i])/338)*100, 1), "% (", nh_metadata_test$`...2`[i], ")")
}

# remove the NA
nh_metadata_test$percent[nh_metadata_test$percent=="NA% ()"] <- ""

nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`[nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`=="Age"] <- "Age group (yrs)"
nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`[nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`=="\u226560"] <- ">=60 years old"
nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`[nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`=="Race/Ethnicity"] <- "Race/ethnicity combined"
nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`[nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`=="Lifetime Illicit Drug Use"] <- "History of drug use"
nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`[nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`=="Unhoused"] <- "Experiencing homelessness"
nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`[nh_metadata_test$`Hepatitis A Outbreak, NH, 2018-2020`=="Ever Vaccinated (≥1 Dose)"] <- "Yes"

# subset to the final tally column: now we have a row with counts and percentages, in the desired format and order to merge with Table 1
nh_metadata_test <- subset(nh_metadata_test, select=-...2)

# match column names between the NH data and Table 1
colnames(nh_metadata_test) <- colnames(subset(final_summary_2, select=c(Variable, `Outbreak-associated NH cases: % (N=X)`)))

# correcting the first row to show the data we want (total outbreak and total sequenced)
nh_metadata_test <- rbind(c(Variable="Total cases", `Outbreak-associated NH cases: % (N=X)`=paste0("N=338 (", nrow(subset(metadata, Sample_Passed_QC=="Yes" & State.of.Residence=="NH")), " sequenced)")), nh_metadata_test)

# now put it all together (MA and NH cases) into a final Table 1.
final_summary_2 <- cbind.data.frame(final_summary_2[,c(1:3)],
                                    nh_metadata_test[,2],
                                    final_summary_2[,c(5:7)])

# change column names and we're done!
colnames(final_summary_2)[colnames(final_summary_2) == "nh_table1_final[, 2]"] <- "Outbreak-associated NH cases: % (N=X)"




##############################
##### Case total by town #####
##############################
metadata_town <- metadata %>%
  subset(droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(overallTownCt = n())

# merge in town population, as well as getting all MA towns represented as metadata_town up to this point only shows the towns that had cases.
metadata_town <- merge(subset(metadata_town, select=-townGeo), # dropping townGeo because those are in the town_pop dataframe
                        town_pop,
                        by.x="cityResidence",
                        by.y="town",
                        all=T)

# change NA to 0 as these are the towns with 0 counts that were merged in from town_pop
metadata_town$overallTownCt[is.na(metadata_town$overallTownCt)] <- 0

# get a state case rate per 100000 people
metadata_town$overallTownRate <- metadata_town$overallTownCt/(metadata_town$townPop/100000)

# create a case rate grouping and factoring
metadata_town$overallRtGrp6 <- factor(ifelse(metadata_town$overallTownRate == 0, "0",
                                              ifelse(metadata_town$overallTownRate > 0 & metadata_town$overallTownRate <=5, ">0-5",
                                                     ifelse(metadata_town$overallTownRate > 5 & metadata_town$overallTownRate <= 10, ">5-10",
                                                            ifelse(metadata_town$overallTownRate > 10 & metadata_town$overallTownRate <= 30, ">10-30",
                                                                   ifelse(metadata_town$overallTownRate > 30 & metadata_town$overallTownRate <= 60, ">30-60",
                                                                          ifelse(metadata_town$overallTownRate > 60, ">60", "Check")))))), 
                                       levels=c("0", ">0-5", ">5-10", ">10-30", ">30-60", ">60"))

# do the same grouping and factoring for the raw count as was done for the rate
metadata_town$overallCtGrp6 <- factor(ifelse(metadata_town$overallTownCt == 0, "0",
                                              ifelse(metadata_town$overallTownCt > 0 & metadata_town$overallTownCt <=5, "1-5",
                                                     ifelse(metadata_town$overallTownCt > 5 & metadata_town$overallTownCt <= 10, "6-10",
                                                            ifelse(metadata_town$overallTownCt > 10 & metadata_town$overallTownCt <= 30, "11-30",
                                                                   ifelse(metadata_town$overallTownCt > 30 & metadata_town$overallTownCt <= 60, "31-60",
                                                                          ifelse(metadata_town$overallTownCt > 60, ">60", "Check")))))), 
                                       levels=c("0", "1-5", "6-10", "11-30", "31-60", ">60"))




########################################
##### Sequenced case total by town #####
########################################
# rinse and repeat the above process
sequenced_cases_town <- metadata %>%
  subset(Sample_Passed_QC=="Yes" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(sequencedTownCt = n())

sequenced_cases_town <- merge(subset(sequenced_cases_town, select=-townGeo),
                              town_pop,
                              by.x="cityResidence",
                              by.y="town",
                              all=T)

sequenced_cases_town$sequencedTownCt[is.na(sequenced_cases_town$sequencedTownCt)] <- 0

sequenced_cases_town$sequencedTownRate <- sequenced_cases_town$sequencedTownCt/(sequenced_cases_town$townPop/100000)

sequenced_cases_town$sequencedRtGrp6 <- factor(ifelse(sequenced_cases_town$sequencedTownRate == 0, "0",
                                                      ifelse(sequenced_cases_town$sequencedTownRate > 0 & 
                                                               sequenced_cases_town$sequencedTownRate <=5, ">0-5",
                                                             ifelse(sequenced_cases_town$sequencedTownRate > 5 & 
                                                                      sequenced_cases_town$sequencedTownRate <= 10, ">5-10",
                                                                    ifelse(sequenced_cases_town$sequencedTownRate > 10 & 
                                                                             sequenced_cases_town$sequencedTownRate <= 30, ">10-30",
                                                                           ifelse(sequenced_cases_town$sequencedTownRate > 30 & 
                                                                                    sequenced_cases_town$sequencedTownRate <= 60, ">30-60",
                                                                                  ifelse(sequenced_cases_town$sequencedTownRate > 60, ">60", "Check")))))), 
                                               levels=c("0", ">0-5", ">5-10", ">10-30", ">30-60", ">60"))

sequenced_cases_town$sequencedCtGrp6 <- factor(ifelse(sequenced_cases_town$sequencedTownCt == 0, "0",
                                                      ifelse(sequenced_cases_town$sequencedTownCt > 0 & 
                                                               sequenced_cases_town$sequencedTownCt <=5, "1-5",
                                                             ifelse(sequenced_cases_town$sequencedTownCt > 5 & 
                                                                      sequenced_cases_town$sequencedTownCt <= 10, "6-10",
                                                                    ifelse(sequenced_cases_town$sequencedTownCt > 10 & 
                                                                             sequenced_cases_town$sequencedTownCt <= 30, "11-30",
                                                                           ifelse(sequenced_cases_town$sequencedTownCt > 30 & 
                                                                                    sequenced_cases_town$sequencedTownCt <= 60, "31-60",
                                                                                  ifelse(sequenced_cases_town$sequencedTownCt > 60, ">60", "Check")))))),
                                               levels=c("0", "1-5", "6-10", "11-30", "31-60", ">60"))




##################################
##### Outbreak cases by town #####
##################################

outbreak_cases_town <- metadata %>%
  subset(outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(outbreakTownCt = n())

outbreak_cases_town <- merge(subset(outbreak_cases_town, select=-townGeo),
                             town_pop,
                             by.x="cityResidence",
                             by.y="town",
                             all=T)

outbreak_cases_town$outbreakTownCt[is.na(outbreak_cases_town$outbreakTownCt)] <- 0

outbreak_cases_town$outbreakTownRate <- outbreak_cases_town$outbreakTownCt/(outbreak_cases_town$townPop/100000)

outbreak_cases_town$outbreakRtGrp6 <- factor(ifelse(outbreak_cases_town$outbreakTownRate == 0, "0",
                                                    ifelse(outbreak_cases_town$outbreakTownRate > 0 & 
                                                             outbreak_cases_town$outbreakTownRate <=5, ">0-5",
                                                           ifelse(outbreak_cases_town$outbreakTownRate > 5 & 
                                                                    outbreak_cases_town$outbreakTownRate <= 10, ">5-10",
                                                                  ifelse(outbreak_cases_town$outbreakTownRate > 10 & 
                                                                           outbreak_cases_town$outbreakTownRate <= 30, ">10-30",
                                                                         ifelse(outbreak_cases_town$outbreakTownRate > 30 & 
                                                                                  outbreak_cases_town$outbreakTownRate <= 60, ">30-60",
                                                                                ifelse(outbreak_cases_town$outbreakTownRate > 60, ">60", "Check")))))), 
                                             levels=c("0", ">0-5", ">5-10", ">10-30", ">30-60", ">60"))

outbreak_cases_town$outbreakCtGrp6 <- factor(ifelse(outbreak_cases_town$outbreakTownCt == 0, "0",
                                                    ifelse(outbreak_cases_town$outbreakTownCt > 0 & 
                                                             outbreak_cases_town$outbreakTownCt <=5, "1-5",
                                                           ifelse(outbreak_cases_town$outbreakTownCt > 5 & 
                                                                    outbreak_cases_town$outbreakTownCt <= 10, "6-10",
                                                                  ifelse(outbreak_cases_town$outbreakTownCt > 10 & 
                                                                           outbreak_cases_town$outbreakTownCt <= 30, "11-30",
                                                                         ifelse(outbreak_cases_town$outbreakTownCt > 30 & 
                                                                                  outbreak_cases_town$outbreakTownCt <= 60, "31-60",
                                                                                ifelse(outbreak_cases_town$outbreakTownCt > 60, ">60", "Check")))))), 
                                             levels=c("0", "1-5", "6-10", "11-30", "31-60", ">60"))




#################################################
##### 2018-2020 outbreak period total cases #####
#################################################

outbreak1_cases_town <- metadata %>%
  subset(tbl1_outbreak_yr=="2018-2020" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(outbreak1OverallTownCt = n())

outbreak1_cases_town <- merge(subset(outbreak1_cases_town, select=-townGeo),
                              town_pop,
                              by.x="cityResidence",
                              by.y="town",
                              all=T)

outbreak1_cases_town$outbreak1OverallTownCt[is.na(outbreak1_cases_town$outbreak1OverallTownCt)] <- 0

outbreak1_cases_town$outbreak1OverallTownRate <- outbreak1_cases_town$outbreak1OverallTownCt/(outbreak1_cases_town$townPop/100000)

outbreak1_cases_town$outbreak1OverallRtGrp6 <- factor(ifelse(outbreak1_cases_town$outbreak1OverallTownRate == 0, "0",
                                                             ifelse(outbreak1_cases_town$outbreak1OverallTownRate > 0 & 
                                                                      outbreak1_cases_town$outbreak1OverallTownRate <=5, ">0-5",
                                                                    ifelse(outbreak1_cases_town$outbreak1OverallTownRate > 5 & 
                                                                             outbreak1_cases_town$outbreak1OverallTownRate <= 10, ">5-10",
                                                                           ifelse(outbreak1_cases_town$outbreak1OverallTownRate > 10 & 
                                                                                    outbreak1_cases_town$outbreak1OverallTownRate <= 30, ">10-30",
                                                                                  ifelse(outbreak1_cases_town$outbreak1OverallTownRate > 30 & 
                                                                                           outbreak1_cases_town$outbreak1OverallTownRate <= 60, ">30-60",
                                                                                         ifelse(outbreak1_cases_town$outbreak1OverallTownRate > 60, ">60", "Check")))))),
                                                      levels=c("0", ">0-5", ">5-10", ">10-30", ">30-60", ">60"))

outbreak1_cases_town$outbreak1OverallCtGrp6 <- factor(ifelse(outbreak1_cases_town$outbreak1OverallTownCt == 0, "0",
                                                             ifelse(outbreak1_cases_town$outbreak1OverallTownCt > 0 & 
                                                                      outbreak1_cases_town$outbreak1OverallTownCt <=5, "1-5",
                                                                    ifelse(outbreak1_cases_town$outbreak1OverallTownCt > 5 & 
                                                                             outbreak1_cases_town$outbreak1OverallTownCt <= 10, "6-10",
                                                                           ifelse(outbreak1_cases_town$outbreak1OverallTownCt > 10 & 
                                                                                    outbreak1_cases_town$outbreak1OverallTownCt <= 30, "11-30",
                                                                                  ifelse(outbreak1_cases_town$outbreak1OverallTownCt > 30 & 
                                                                                           outbreak1_cases_town$outbreak1OverallTownCt <= 60, "31-60",
                                                                                         ifelse(outbreak1_cases_town$outbreak1OverallTownCt > 60, ">60", "Check")))))), 
                                                      levels=c("0", "1-5", "6-10", "11-30", "31-60", ">60"))




##############################################
##### 2018-2020 sequenced outbreak cases #####
##############################################

outbreak1_seq_cases_town <- metadata %>%
  subset(tbl1_outbreak_yr=="2018-2020" & Sample_Passed_QC=="Yes" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(outbreak1SequencedTownCt = n())

outbreak1_seq_cases_town <- merge(subset(outbreak1_seq_cases_town, select=-townGeo),
                                  town_pop,
                                  by.x="cityResidence",
                                  by.y="town",
                                  all=T)

outbreak1_seq_cases_town$outbreak1SequencedTownCt[is.na(outbreak1_seq_cases_town$outbreak1SequencedTownCt)] <- 0

outbreak1_seq_cases_town$outbreak1SequencedTownRate <- outbreak1_seq_cases_town$outbreak1SequencedTownCt/(outbreak1_seq_cases_town$townPop/100000)

outbreak1_seq_cases_town$outbreak1SequencedRtGrp6 <- factor(ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownRate == 0, "0",
                                                                   ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownRate > 0 & 
                                                                            outbreak1_seq_cases_town$outbreak1SequencedTownRate <=5, ">0-5",
                                                                          ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownRate > 5 & 
                                                                                   outbreak1_seq_cases_town$outbreak1SequencedTownRate <= 10, ">5-10",
                                                                                 ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownRate > 10 & 
                                                                                          outbreak1_seq_cases_town$outbreak1SequencedTownRate <= 30, ">10-30",
                                                                                        ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownRate > 30 & 
                                                                                                 outbreak1_seq_cases_town$outbreak1SequencedTownRate <= 60, ">30-60",
                                                                                               ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownRate > 60, ">60", "Check")))))), 
                                                            levels=c("0", ">0-5", ">5-10", ">10-30", ">30-60", ">60"))

outbreak1_seq_cases_town$outbreak1SequencedCtGrp6 <- factor(ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownCt == 0, "0",
                                                                   ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownCt > 0 & 
                                                                            outbreak1_seq_cases_town$outbreak1SequencedTownCt <=5, "1-5",
                                                                          ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownCt > 5 & 
                                                                                   outbreak1_seq_cases_town$outbreak1SequencedTownCt <= 10, "6-10",
                                                                                 ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownCt > 10 & 
                                                                                          outbreak1_seq_cases_town$outbreak1SequencedTownCt <= 30, "11-30",
                                                                                        ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownCt > 30 & 
                                                                                                 outbreak1_seq_cases_town$outbreak1SequencedTownCt <= 60, "31-60",
                                                                                               ifelse(outbreak1_seq_cases_town$outbreak1SequencedTownCt > 60, ">60", "Check")))))), 
                                                            levels=c("0", "1-5", "6-10", "11-30", "31-60", ">60"))




####################################
##### 2018-2020 outbreak cases #####
####################################

outbreak1_out_cases_town <- metadata %>%
  subset(tbl1_outbreak_yr=="2018-2020" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(outbreak1OutbreakTownCt = n())

outbreak1_out_cases_town <- merge(subset(outbreak1_out_cases_town, select=-townGeo),
                                  town_pop,
                                  by.x="cityResidence",
                                  by.y="town",
                                  all=T)

outbreak1_out_cases_town$outbreak1OutbreakTownCt[is.na(outbreak1_out_cases_town$outbreak1OutbreakTownCt)] <- 0

outbreak1_out_cases_town$outbreak1OutbreakTownRate <- outbreak1_out_cases_town$outbreak1OutbreakTownCt/(outbreak1_out_cases_town$townPop/100000)

outbreak1_out_cases_town$outbreak1OutbreakRtGrp6 <- factor(ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownRate == 0, "0",
                                                                  ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownRate > 0 & 
                                                                           outbreak1_out_cases_town$outbreak1OutbreakTownRate <=5, ">0-5",
                                                                         ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownRate > 5 & 
                                                                                  outbreak1_out_cases_town$outbreak1OutbreakTownRate <= 10, ">5-10",
                                                                                ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownRate > 10 & 
                                                                                         outbreak1_out_cases_town$outbreak1OutbreakTownRate <= 30, ">10-30",
                                                                                       ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownRate > 30 & 
                                                                                                outbreak1_out_cases_town$outbreak1OutbreakTownRate <= 60, ">30-60",
                                                                                              ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownRate > 60, ">60", "Check")))))), 
                                                           levels=c("0", ">0-5", ">5-10", ">10-30", ">30-60", ">60"))

outbreak1_out_cases_town$outbreak1OutbreakCtGrp6 <- factor(ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownCt == 0, "0",
                                                                  ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownCt > 0 & 
                                                                           outbreak1_out_cases_town$outbreak1OutbreakTownCt <=5, "1-5",
                                                                         ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownCt > 5 & 
                                                                                  outbreak1_out_cases_town$outbreak1OutbreakTownCt <= 10, "6-10",
                                                                                ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownCt > 10 & 
                                                                                         outbreak1_out_cases_town$outbreak1OutbreakTownCt <= 30, "11-30",
                                                                                       ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownCt > 30 & 
                                                                                                outbreak1_out_cases_town$outbreak1OutbreakTownCt <= 60, "31-60",
                                                                                              ifelse(outbreak1_out_cases_town$outbreak1OutbreakTownCt > 60, ">60", "Check")))))), 
                                                           levels=c("0", "1-5", "6-10", "11-30", "31-60", ">60"))




#################################################
##### 2023-2024 outbreak period total cases #####
#################################################

outbreak2_cases_town <- metadata %>%
  subset(tbl1_outbreak_yr=="2023-2024" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(outbreak2OverallTownCt = n())

outbreak2_cases_town <- merge(subset(outbreak2_cases_town, select=-townGeo),
                              town_pop,
                              by.x="cityResidence",
                              by.y="town",
                              all=T)

outbreak2_cases_town$outbreak2OverallTownCt[is.na(outbreak2_cases_town$outbreak2OverallTownCt)] <- 0

outbreak2_cases_town$outbreak2OverallTownRate <- outbreak2_cases_town$outbreak2OverallTownCt/(outbreak2_cases_town$townPop/100000)

outbreak2_cases_town$outbreak2OverallRtGrp3 <- factor(ifelse(outbreak2_cases_town$outbreak2OverallTownRate == 0, "0",
                                                             ifelse(outbreak2_cases_town$outbreak2OverallTownRate > 0 & 
                                                                      outbreak2_cases_town$outbreak2OverallTownRate <=5, ">0-5",
                                                                    ifelse(outbreak2_cases_town$outbreak2OverallTownRate > 5, ">5", "Check"))), 
                                                      levels=c("0", ">0-5", ">5"))

outbreak2_cases_town$outbreak2OverallCtGrp3 <- factor(ifelse(outbreak2_cases_town$outbreak2OverallTownCt == 0, "0",
                                                             ifelse(outbreak2_cases_town$outbreak2OverallTownCt > 0 & 
                                                                      outbreak2_cases_town$outbreak2OverallTownCt <=5, "1-5",
                                                                    ifelse(outbreak2_cases_town$outbreak2OverallTownCt > 5, ">5", "Check"))), 
                                                      levels=c("0", "1-5", ">5"))




##############################################
##### 2023-2024 sequenced outbreak cases #####
##############################################

outbreak2_seq_cases_town <- metadata %>%
  subset(tbl1_outbreak_yr=="2023-2024" & Sample_Passed_QC=="Yes" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(outbreak2SequencedTownCt = n())

outbreak2_seq_cases_town <- merge(subset(outbreak2_seq_cases_town, select=-townGeo),
                                  town_pop,
                                  by.x="cityResidence",
                                  by.y="town",
                                  all=T)

outbreak2_seq_cases_town$outbreak2SequencedTownCt[is.na(outbreak2_seq_cases_town$outbreak2SequencedTownCt)] <- 0

outbreak2_seq_cases_town$outbreak2SequencedTownRate <- outbreak2_seq_cases_town$outbreak2SequencedTownCt/(outbreak2_seq_cases_town$townPop/100000)

outbreak2_seq_cases_town$outbreak2SequencedRtGrp3 <- factor(ifelse(outbreak2_seq_cases_town$outbreak2SequencedTownRate == 0, "0",
                                                                   ifelse(outbreak2_seq_cases_town$outbreak2SequencedTownRate > 0 & 
                                                                            outbreak2_seq_cases_town$outbreak2SequencedTownRate <=5, ">0-5",
                                                                          ifelse(outbreak2_seq_cases_town$outbreak2SequencedTownRate > 5, ">5", "Check"))), 
                                                            levels=c("0", ">0-5", ">5"))

outbreak2_seq_cases_town$outbreak2SequencedCtGrp3 <- factor(ifelse(outbreak2_seq_cases_town$outbreak2SequencedTownCt == 0, "0",
                                                                   ifelse(outbreak2_seq_cases_town$outbreak2SequencedTownCt > 0 & 
                                                                            outbreak2_seq_cases_town$outbreak2SequencedTownCt <=5, "1-5",
                                                                          ifelse(outbreak2_seq_cases_town$outbreak2SequencedTownCt > 5, ">5", "Check"))), 
                                                            levels=c("0", "1-5", ">5"))




####################################
##### 2023-2024 outbreak cases #####
####################################

outbreak2_out_cases_town <- metadata %>%
  subset(tbl1_outbreak_yr=="2023-2024" & outbreak=="Yes" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(outbreak2OutbreakTownCt = n())

outbreak2_out_cases_town <- merge(subset(outbreak2_out_cases_town, select=-townGeo),
                                  town_pop,
                                  by.x="cityResidence",
                                  by.y="town",
                                  all=T)

outbreak2_out_cases_town$outbreak2OutbreakTownCt[is.na(outbreak2_out_cases_town$outbreak2OutbreakTownCt)] <- 0

outbreak2_out_cases_town$outbreak2OutbreakTownRate <- outbreak2_out_cases_town$outbreak2OutbreakTownCt/(outbreak2_out_cases_town$townPop/100000)

outbreak2_out_cases_town$outbreak2OutbreakRtGrp3 <- factor(ifelse(outbreak2_out_cases_town$outbreak2OutbreakTownRate == 0, "0",
                                                                  ifelse(outbreak2_out_cases_town$outbreak2OutbreakTownRate > 0 & 
                                                                           outbreak2_out_cases_town$outbreak2OutbreakTownRate <=5, ">0-5",
                                                                         ifelse(outbreak2_out_cases_town$outbreak2OutbreakTownRate > 5, ">5", "Check"))), 
                                                           levels=c("0", ">0-5", ">5"))

outbreak2_out_cases_town$outbreak2OutbreakCtGrp3 <- factor(ifelse(outbreak2_out_cases_town$outbreak2OutbreakTownCt == 0, "0",
                                                                  ifelse(outbreak2_out_cases_town$outbreak2OutbreakTownCt > 0 & 
                                                                           outbreak2_out_cases_town$outbreak2OutbreakTownCt <=5, "1-5",
                                                                         ifelse(outbreak2_out_cases_town$outbreak2OutbreakTownCt > 5, ">5", "Check"))),
                                                           levels=c("0", "1-5", ">5"))




######################
##### IIIA cases #####
######################

iiia_cases_town <- metadata %>%
  subset(Subgenotype=="IIIA" & droppedDups==0 & State.of.Residence=="MA") %>%
  group_by(cityResidence, townGeo) %>%
  summarise(iiiaTownCt = n())

iiia_cases_town <- merge(subset(iiia_cases_town, select=-townGeo),
                         town_pop,
                         by.x="cityResidence",
                         by.y="town",
                         all=T)

iiia_cases_town$iiiaTownCt[is.na(iiia_cases_town$iiiaTownCt)] <- 0

iiia_cases_town$iiiaTownRate <- iiia_cases_town$iiiaTownCt/(iiia_cases_town$townPop/100000)

iiia_cases_town$iiiaRtGrp3 <- factor(ifelse(iiia_cases_town$iiiaTownRate == 0, "0",
                                            ifelse(iiia_cases_town$iiiaTownRate > 0 & iiia_cases_town$iiiaTownRate <=5, ">0-5",
                                                   ifelse(iiia_cases_town$iiiaTownRate > 5 & iiia_cases_town$iiiaTownRate <=10, ">5-10",
                                                          ifelse(iiia_cases_town$iiiaTownRate > 10 & iiia_cases_town$iiiaTownRate <=20, ">10-20",
                                                                 ifelse(iiia_cases_town$iiiaTownRate > 20 & iiia_cases_town$iiiaTownRate <=30, ">20-30",
                                                                        ifelse(iiia_cases_town$iiiaTownRate > 30 & iiia_cases_town$iiiaTownRate <=50, ">30-50",
                                                                               ifelse(iiia_cases_town$iiiaTownRate > 50, ">=50", "Check"))))))), 
                                     levels=c("0", ">0-5", ">5-10", ">10-20", ">20-30", ">30-50", ">=50"))

iiia_cases_town$iiiaCtGrp3 <- factor(ifelse(iiia_cases_town$iiiaTownCt == 0, "0",
                                            ifelse(iiia_cases_town$iiiaTownCt > 0 & iiia_cases_town$iiiaTownCt <=5, "1-5",
                                                   ifelse(iiia_cases_town$iiiaTownCt > 5 & iiia_cases_town$iiiaTownCt <=10, "6-10",
                                                          ifelse(iiia_cases_town$iiiaTownCt > 10 & iiia_cases_town$iiiaTownCt <=20, "11-20",
                                                                 ifelse(iiia_cases_town$iiiaTownCt > 20 & iiia_cases_town$iiiaTownCt <=30, "21-30",
                                                                        ifelse(iiia_cases_town$iiiaTownCt > 30, ">=31", "Check")))))), 
                                     levels=c("0", "1-5", "6-10", "11-20", "21-30", ">=31"))




##############################################
##### Combine all geographical summaries #####
##############################################

# now combine all of the geographical summaries into one dataset, much like how the cases_week_all dataset has all of the metadata variables. First, creating a list of all the dataframes I want to combine so that I can combine them in one step instead of multiple single join steps.
# I am not including townGeo as when exporting into excel from R, that column goes blank. I am thus merging in the town shapes in the analysis program as I can do it straight from the gis file.
df_list <- list(subset(metadata_town, select=c(cityResidence, overallTownCt, overallTownRate, overallRtGrp6, overallCtGrp6)),
                subset(outbreak_cases_town, select=c(cityResidence, outbreakTownCt, outbreakTownRate, outbreakRtGrp6, outbreakCtGrp6)),
                subset(sequenced_cases_town, select=c(cityResidence, sequencedTownCt, sequencedTownRate, sequencedRtGrp6, sequencedCtGrp6)),
                subset(outbreak1_cases_town, select=c(cityResidence, outbreak1OverallTownCt, outbreak1OverallTownRate, outbreak1OverallRtGrp6, outbreak1OverallCtGrp6)),
                subset(outbreak2_cases_town, select=c(cityResidence, outbreak2OverallTownCt, outbreak2OverallTownRate, outbreak2OverallRtGrp3, outbreak2OverallCtGrp3)),
                subset(outbreak1_seq_cases_town, select=c(cityResidence, outbreak1SequencedTownCt, outbreak1SequencedTownRate, outbreak1SequencedRtGrp6, outbreak1SequencedCtGrp6)),
                subset(outbreak2_seq_cases_town, select=c(cityResidence, outbreak2SequencedTownCt, outbreak2SequencedTownRate, outbreak2SequencedRtGrp3, outbreak2SequencedCtGrp3)),
                subset(outbreak1_out_cases_town, select=c(cityResidence, outbreak1OutbreakTownCt, outbreak1OutbreakTownRate, outbreak1OutbreakRtGrp6, outbreak1OutbreakCtGrp6)),
                subset(outbreak2_out_cases_town, select=c(cityResidence, outbreak2OutbreakTownCt, outbreak2OutbreakTownRate, outbreak2OutbreakRtGrp3, outbreak2OutbreakCtGrp3)),
                subset(iiia_cases_town, select=c(cityResidence, iiiaTownCt, iiiaTownRate, iiiaRtGrp3, iiiaCtGrp3)))

# now join all of the dataframes from the list
metadata_town_clean <- df_list %>%
  purrr::reduce(full_join, by="cityResidence")

# merge in county
metadata_town_clean <- merge(metadata_town_clean,
                              town_county,
                              by.x="cityResidence",
                              by.y="town",
                              all.x=T)
