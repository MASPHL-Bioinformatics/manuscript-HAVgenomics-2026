# create direct/known epi and all-vs-all epi links derived from linkage groups
# load juniper and output direct transmission to file
# create pairwise epi-and-juniper links dataframe object

readRenviron("_environment")
METADATA_FILE <- Sys.getenv("METADATA_FILE")
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR")
JUNIPER_RES <- Sys.getenv("JUNIPER_RES")
JUNIPER_OUT <- Sys.getenv("JUNIPER_OUT")
COMBINED_OUT <- Sys.getenv("COMBINED_OUT")

###################################################################
###################################################################
###################################################################

pkglist = c("dplyr", "igraph", "ggplot2", "stringr", "splitstackshape", "data.table", "tidyr")
lapply(pkglist, require, character.only = TRUE)

devtools::install_github('xavierdidelot/TransPhylo')
devtools::install_github("broadinstitute/juniper0", dependencies = FALSE) #because t kept failing on TransPhylo install

library(juniper0)

###############################################################
###############################################################
###############################################################=
# get list of all direct/known epi links from metadata
metadata <- as.data.frame(read.csv(METADATA_FILE))
metadata$caseid <- as.character(metadata$caseid) # format caseid
# apply trimws() to strip leading/trailing whitespace for all character columns - doesn't affect numeric/date columns
metadata <- metadata %>%
  mutate(across(where(is.character), trimws))

# not necessary, but for good practice, drop HepA/USA/MA/0046/2018
metadata <- metadata[metadata$MASPHL_Sequencing_ID != "HepA/USA/MA/0046/2018", ]

# create a table with all epi links. n = 44
sum(ifelse(metadata$individual_links=="", 0, stringr::str_count(metadata$individual_links, pattern="\\,")+1))/2 # total number of pairs. Giving a 0 for a link count if there is no individual link, and then giving a 1 plus the number of commas for the number of links for those not having a null individual link value.
every_link <- subset(metadata, individual_links!="", select=c(caseid, individual_links))

# create a long format count of each epi link, both directions, so it has double the rows as the total epi link count
every_link <- every_link %>%
  splitstackshape::cSplit("individual_links", sep=",", direction="wide") %>%
  tidyr::pivot_longer(cols=c(2:5), names_to="linkNum", values_to="linked_caseid") %>%
  subset(!is.na(linked_caseid))

every_link$caseid <- as.character(every_link$caseid)
every_link$linked_caseid <- as.character(every_link$linked_caseid)
pair_link <- every_link
nrow(pair_link) #88


# remove redundant pairs in 88-row set using vectorized operations - can scale to larger data
pair_link_unique <- pair_link %>%
  mutate(
    id1 = pmin(caseid, linked_caseid), # parallel min/max are vectorized operations. compare values across two vectors element-by-element
    id2 = pmax(caseid, linked_caseid) # regardless of order, pmin and pmax will always put the smaller/lower in caseid/linked_caseid respectively
  ) %>%
  distinct(id1, id2, .keep_all = TRUE) %>% # distinct() drops duplicates and leaves unique pairs
  select(-id1, -id2)

# filter df to only contain pairs where BOTH caseids have passed QC
qc_passed <- metadata %>%
  filter(Sample_Passed_QC == "Yes") %>%
  pull(caseid)

pair_link_filtered <- pair_link_unique %>%
  filter(
    caseid %in% qc_passed,
    linked_caseid %in% qc_passed
  )

pair_link_filtered #n = 16 

# what is the distribution of subgenotypes among these pairs?
pair_link_filtered_sb <- pair_link_filtered |>
  left_join(metadata |> select(caseid, sb = Subgenotype), by = "caseid") |>
  left_join(metadata |> select(caseid, sb_linked = Subgenotype), by = c("linked_caseid" = "caseid"))

print(sum(pair_link_filtered_sb$sb == "IIIA" & pair_link_filtered_sb$sb_linked == "IIIA", na.rm = TRUE))
# 16

##########################
# # Concordance test - Another route for me to explode all caseid/individual_links, regardless of sequencing status
# pair_link <- metadata %>%
#   filter(!is.na(individual_links) & individual_links != "") %>%
#   separate_rows(individual_links, sep = ",\\s*") %>%
#   select(caseid, individual_links) %>%
#   dplyr::rename(linked_caseid = individual_links)
# nrow(pair_link)
# # 88 - dedup and filtering gives me 16 again

###############################################################
###############################################################
###############################################################

# Build undirected graph from first two columns
g <- graph_from_data_frame(
  pair_link_filtered[, c("caseid", "linked_caseid")],
  directed = FALSE
)

# Find connected components
comps <- components(g)

# Map component number → letter label (A, B, C, ...)
comp_labels <- setNames(
  LETTERS[comps$membership],
  names(comps$membership)
)

# Build output dataframe
linkage_groups <- data.frame(
  strain         = names(comp_labels),
  linkage_group  = unname(comp_labels),
  stringsAsFactors = FALSE
) |>
  arrange(linkage_group, strain)

linkage_groups # n=24

# figure out what letter label the added jail group will have
jail_label <- LETTERS[which(LETTERS == tail(linkage_groups$linkage_group, 1)) + 1]

# # figure out what letter label the added restaurant group will have
restaurant_label <- LETTERS[which(LETTERS==jail_label) + 1]

# get the Cluster 2 (jail) cases that passed QC - n=5
jail_cluster <- metadata[(metadata$cluster_links == "Cluster 2 (jail)" & metadata$Sample_Passed_QC=="Yes"), "caseid"]
new_rows <- data.frame(strain = jail_cluster, linkage_group = jail_label)
new_rows # 5 samples

linkage_groups <- rbind(linkage_groups, new_rows) # 24 lines (A-J) without jail cluster; adding those gives 29 lines (A-K)
# for comparison, Lydia's original contains 27, and she had 25 in the juniper-epi link combined set. I suspect that Lydia just started with a smaller set than I did, possibly because of those epi links that cut off in earlier versions.
linkage_groups

# # get the Cluster 1 (restaurant) cases that passed QC
# rest_cluster <- metadata[(metadata$cluster_links == "Cluster 1 (restaurant)" & metadata$Sample_Passed_QC=="Yes"), "caseid"]
# new_rest_rows <- data.frame(strain = rest_cluster, linkage_group = restaurant_label)
# new_rest_rows
# linkage_groups <- rbind(linkage_groups, new_rest_rows)
# linkage_groups
# This only adds a single case passing QC; this is excluded from the linkage group definition criteria (must contain 2+ cases passing QC)

# Check that I only have IIIA isolates in this set - all IIIA
metadata$Subgenotype[match(trimws(linkage_groups$strain), trimws(metadata$caseid))]

# swap out caseid for MASPHL_Sequencing_ID
linkage_groups <- linkage_groups %>%
  left_join(metadata %>% select(caseid, MASPHL_Sequencing_ID),
            by = c("strain" = "caseid")) %>%
  mutate(strain = coalesce(MASPHL_Sequencing_ID, strain)) %>% # used coalesce() because it will default back to original strain value if there's no match in the metadata caseid column
  select(-MASPHL_Sequencing_ID)

linkage_groups

###############################################################
###############################################################
###############################################################
# now explode all pairwise combinations within linkage groups
# "linkage pairs"
setDT(linkage_groups)

# unique pairs in a group of n items = n(n-1)/2. We should have 31 unique pairs
linkage_pairs <- linkage_groups[, {
  ids <- strain
  if (length(ids) > 1) {
    pairs <- combn(ids, 2) #combn(ids, 2) generates all non-redundant pairs within each group automatically
    .(case1 = pairs[1,], case2 = pairs[2,])
  } else {
    .(case1 = character(0), case2 = character(0))
  }
}, by = linkage_group]
linkage_pairs #n=31

# # concordance test
# pairs <- linkage_groups[
#   , as.data.table(t(combn(strain, 2))), 
#   by = linkage_group
# ]
# setnames(pairs, c("linkage_group","strain1","strain2"))
# pairs # n=31

##################################################################
##################################################################
##################################################################
# load juniper RData for "pi 0.4" and write direct transmissions to file

load(JUNIPER_RES)
if(exists("my_res")) {
  assign(sub("_downsampled\\.RData$", "",substring(JUNIPER_RES, regexpr("clade", JUNIPER_RES), nchar(JUNIPER_RES))), my_res)
  rm(my_res)
} else if(exists("small_res")) {
  assign(sub("_downsampled\\.RData$", "",substring(JUNIPER_RES, regexpr("clade", JUNIPER_RES), nchar(JUNIPER_RES))), small_res)
  rm(small_res)
} else {
  print(paste("ERROR: res not read in for file ", JUNIPER_RES, sep=""))
}

# summarize with juniper default function
clade3a_juniper_pi_0.4_summary = juniper0::summarize(clade3a_juniper_pi_0.4_juniper,0)

clade3a_juniper_pi_0.4_tmrca = c(quantile(clade3a_juniper_pi_0.4_summary$time_of_MRCA, 0.025, type=1),
                      median(clade3a_juniper_pi_0.4_summary$time_of_MRCA),
                      quantile(clade3a_juniper_pi_0.4_summary$time_of_MRCA, 0.975, type=1))

# check out results
clade3a_juniper_pi_0.4_transmissions = as.data.table(as.table(clade3a_juniper_pi_0.4_summary$direct_transmissions))[N > 0.5,]
epi_pairs <- linkage_pairs
clade3a_juniper_pi_0.4_epi_transmissions = rbind(clade3a_juniper_pi_0.4_transmissions[epi_pairs, on=.(V1=case1, V2=case2)],
                                      clade3a_juniper_pi_0.4_transmissions[epi_pairs, on=.(V1=case2, V2=case1)])
clade3a_juniper_pi_0.4_epi_transmissions = clade3a_juniper_pi_0.4_epi_transmissions[!is.na(N),]


clade3a_juniper_pi_0.4_direct_transmissions <- as.data.table(as.table(clade3a_juniper_pi_0.4_summary$direct_transmissions))[N>0,]
colnames(clade3a_juniper_pi_0.4_direct_transmissions) = c("from","to","weight")

write.table(clade3a_juniper_pi_0.4_direct_transmissions, JUNIPER_OUT, sep = "\t", row.names = FALSE, quote = FALSE)

nrow(clade3a_juniper_pi_0.4_direct_transmissions)
#n = 1963

##################################################################
##################################################################
##################################################################
# read in juniper direct transmissions
junipert <- read.table(JUNIPER_OUT, sep="\t", header=TRUE)
nrow(junipert) #1963
junipert_0.5 <- junipert[junipert$weight >= 0.5,] # filter by > 0.5 
nrow(junipert_0.5) #n=107

# double-check - are these redundant? 
junipert_0.5[junipert_0.5$from == "HepA/USA/MA/0078/2018" | junipert_0.5$to   == "HepA/USA/MA/0078/2018",]
# no, not redundant

# remove any inferred transmissions including HepA/USA/MA/0046/2018
junipert_0.5 <- junipert_0.5[junipert_0.5$from != "HepA/USA/MA/0046/2018", ]
junipert_0.5 <- junipert_0.5[junipert_0.5$to != "HepA/USA/MA/0046/2018", ]
nrow(junipert_0.5) # n=78

#need to combine with the direct epi links
# sanity checkL recreate keys for each source separately, check overlap between them
juniper_keys <- paste(pmin(junipert_0.5$from, junipert_0.5$to), 
                      pmax(junipert_0.5$from, junipert_0.5$to), sep = "|") #n=107

epi_keys <- paste(pmin(linkage_pairs$case1, linkage_pairs$case2), 
                  pmax(linkage_pairs$case1, linkage_pairs$case2), sep = "|") #n=31

intersect(juniper_keys, epi_keys) #6

# Standardize column names
juniper_dt <- as.data.table(junipert_0.5)[, .(case1 = from, case2 = to, support = weight)]
juniper_dt[, type := "juniper"]
juniper_dt

epi_dt <- copy(linkage_pairs)[, support := NA_real_]
epi_dt[, type := "epi"]

# Combine both link sets
combined <- rbind(epi_dt, juniper_dt, fill = TRUE)
# Create a sorted key to identify duplicates regardless of direction
combined[, key := paste(pmin(case1, case2), pmax(case1, case2), sep = "|")]
# if the same pair appears in both, collapse to "epi and juniper"
combined <- combined[, .(
  case1   = first(case1),
  case2   = first(case2),
  support = max(support, na.rm = TRUE),
  type    = ifelse(.N > 1, "epi and juniper", first(type))
), by = key]


# Clean up
combined[support == -Inf, support := NA_real_]  # max of all-NA returns -Inf
combined[, key := NULL]

setnames(combined, c("case1", "case2"), c("from", "to"))
nrow(combined) # n = 103

# double-check
combined[combined$type %in% c("epi", "epi and juniper")] #n=31
combined[combined$type %in% c("epi and juniper")] #n=6


write.table(combined, COMBINED_OUT, sep = "\t", row.names = FALSE, quote = FALSE)


# double-check to see if any juniper links include samples that didn't pass QC (just in case)
passed <- metadata$MASPHL_Sequencing_ID[metadata$Sample_Passed_QC == "Yes"]
combined[, QC_check := ifelse(from %in% passed & to %in% passed, "both passed", "one or both failed")]
nrow(combined[combined$QC_check=="one or both failed", ]) #0
combined[, QC_check := NULL]
