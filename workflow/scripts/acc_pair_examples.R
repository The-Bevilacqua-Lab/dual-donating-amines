library(dplyr)

#### SET INITIAL VARIABLES ####

# Creates data frames from the combined data.
combined_df <- read.csv(snakemake@input[["combined"]], header = TRUE, na.strings = "NaN", comment.char = "#",
                        colClasses = c("PDB" = "character"))

# Specify the columns to write to the output CSV file.
output_columns <- c('don_index', 'don_name', 'don_resn', 'don_resi', 'don_chain', 'eta', 'theta', 'chi', 'type',
                    'acc_pair_1_name', 'acc_pair_1_resi', 'acc_pair_1_chain',
                    'acc_pair_2_name', 'acc_pair_2_resi', 'acc_pair_2_chain',
                    'model', 'PDB', 'eq_class_member')

#### 5A - SHEARED GA BASE PAIR ####

# Extract the dual donor A(N6)'s that interact with the N3 and O2' atoms of guanines.
n6a_df <- combined_df %>%
  filter(don_resn == "A", type == 2, acc_pair_combined == "G(N3), N(O2')", same_resi == "True") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Which guanines form a G(N2)-to-A(N7) H-bond?
n2g_n7a_h_bond_df <- merge(combined_df,
                           n6a_df[c('don_resi', 'don_chain', 'acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')],
                           by.x = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB'),
                           by.y = c('don_resi', 'don_chain', 'acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_name == "N7") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Backtrack to identify the A(N6)'s involved in the sheared GA base pair.
n6a_sheared_ga_df <- merge(n6a_df,
                           n2g_n7a_h_bond_df[c('don_resi', 'don_chain', 'model', 'PDB')],
                           by.x = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB'),
                           by.y = c('don_resi', 'don_chain', 'model', 'PDB'))

# Write the data to a csv.
write.csv(n6a_sheared_ga_df[output_columns], snakemake@output[["supplemental_data_S1"]], na = "NaN")

# Which guanines form a G(N2)-to-A(N7) H-bond and bear a dual-donating N2?
n2g_n7a_h_bond_dual_df <- merge(combined_df,
                                n6a_df[c('don_resi', 'don_chain', 'acc_pair_1_resi', 'acc_pair_1_chain',
                                         'model', 'PDB')],
                                by.x = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB'),
                                by.y = c('don_resi', 'don_chain', 'acc_pair_1_resi', 'acc_pair_1_chain',
                                         'model', 'PDB')) %>%
  filter(type == 2, h_bond == 1, acc_name == "N7") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(n2g_n7a_h_bond_dual_df[output_columns], snakemake@output[["supplemental_data_S2"]], na = "NaN")

#### 5B - WC/H A-MINOR MOTIF ####

# Which guanines form a G(N2)-to-C(O2) H-bond?
n2g_o2c_h_bond_df <- merge(combined_df, n6a_df[c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')],
                           by.x = c('don_resi', 'don_chain', 'model', 'PDB'),
                           by.y = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1 & acc_resn == "C" & acc_name == "O2") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Regarding the cytosines accepting an H-bond via their O2, which ones form a C(N4)-to-G(O6) H-bond, involving the same
# guanine that donates to the O2 of the cytosine? The G(N2)-to-C(O2) and C(N4)-to-G(O6) H-bonds suggests that the C and
# G are forming a canonical base pair. Meanwhile, the G is engaged with a dual-donating A(N6) via its sugar edge. The
# entire unit is a WC/H A-minor motif.
n4c_o6g_h_bond_df <- merge(combined_df,
                           n2g_o2c_h_bond_df[c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB')],
                           by.x = c('don_resi', 'don_chain', 'acc_resi', 'acc_chain', 'model', 'PDB'),
                           by.y = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_name == "O6") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Backtrack to identify the A's involved in the WC/H A-minor motifs.
a_minor_df <- merge(n6a_df, n4c_o6g_h_bond_df[c('acc_resi', 'acc_chain', 'model', 'PDB')],
                    by.x = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB'),
                    by.y = c('acc_resi', 'acc_chain', 'model', 'PDB'))

# Write the data to a csv.
write.csv(a_minor_df[output_columns], snakemake@output[["supplemental_data_S3"]], na = "NaN")

# Which guanines form a G(N2)-to-A(N1) H-bond, where the A(N1) is from the adenine in the A-minor motif?
n2g_n1a_h_bond_df <- merge(combined_df, a_minor_df[c('acc_pair_1_resi', 'acc_pair_1_chain', 'don_resi', 'don_chain',
                                                     'model', 'PDB')],
                           by.x = c('don_resi', 'don_chain', 'acc_resi', 'acc_chain', 'model', 'PDB'),
                           by.y = c('acc_pair_1_resi', 'acc_pair_1_chain', 'don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_name == "N1") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Make the adenines the focus of the data frame (e.g., don_resi refers to the A, not the G).
adenine_df <- merge(a_minor_df, n2g_n1a_h_bond_df[c('acc_resi', 'acc_chain', 'model', 'PDB')],
                    by.x = c('don_resi', 'don_chain', 'model', 'PDB'),
                    by.y = c('acc_resi', 'acc_chain', 'model', 'PDB'))

# Write the data to a csv.
write.csv(adenine_df[output_columns], snakemake@output[["supplemental_data_S5"]], na = "NaN")

# Which guanines form a G(N2)-to-A(N7) H-bond, where the A(N7) is from the adenine in the A-minor motif?
n2g_n7a_h_bond_df <- merge(combined_df, a_minor_df[c('acc_pair_1_resi', 'acc_pair_1_chain', 'don_resi', 'don_chain',
                                                     'model', 'PDB')],
                           by.x = c('don_resi', 'don_chain', 'acc_resi', 'acc_chain', 'model', 'PDB'),
                           by.y = c('acc_pair_1_resi', 'acc_pair_1_chain', 'don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_name == "N7") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Make the adenines the focus of the data frame (e.g., don_resi refers to the A, not the G).
adenine_df <- merge(a_minor_df, n2g_n7a_h_bond_df[c('acc_resi', 'acc_chain', 'model', 'PDB')],
                    by.x = c('don_resi', 'don_chain', 'model', 'PDB'),
                    by.y = c('acc_resi', 'acc_chain', 'model', 'PDB'))

# Write the data to a csv.
write.csv(adenine_df[output_columns], snakemake@output[["supplemental_data_S6"]], na = "NaN")

# Which of the guanines also bear a dual-donating G(N2)?
dual_g_df <- merge(combined_df, n4c_o6g_h_bond_df[c('acc_resi', 'acc_chain', 'model', 'PDB')],
                   by.x = c('don_resi', 'don_chain', 'model', 'PDB'), 
                   by.y = c('acc_resi', 'acc_chain', 'model', 'PDB')) %>%
  filter(type == 2) %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(dual_g_df[output_columns], snakemake@output[["supplemental_data_S4"]], na = "NaN")

#### 5C - A DUAL-DONATING AMINE TO A HOOGSTEEN ####

# Extract the adenines amines that dual donate to the N7 and NPO of another A.
a_hoogsteen_df <- combined_df %>%
  filter(don_resn == "A", type == 2, acc_pair_combined == "A(N7), N(NPO)", same_resi == "True") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Which adenines bearing the acceptor pair form a A(N6)-to-A(N7) H-bond?
n6a_n7a_h_bond_df <- merge(combined_df, a_hoogsteen_df[c('don_resi', 'don_chain',
                                                         'acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')],
                           by.x = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB'),
                           by.y = c('don_resi', 'don_chain', 'acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_name == "N7") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Backtrack to identify the dual-donating A(N6)'s involved in the tHH AA base pair.
tHH_aa_df <- merge(a_hoogsteen_df, n6a_n7a_h_bond_df[c('don_resi', 'don_chain', 'model', 'PDB')],
                   by.x = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB'),
                   by.y = c('don_resi', 'don_chain', 'model', 'PDB'))

# Write the data to a csv.
write.csv(tHH_aa_df[output_columns], snakemake@output[["supplemental_data_S7"]], na = "NaN")

# Now, find which adenines bearing the acceptor pair form an A(N6)-to-A(N7) H-bond and bear a dual-donating N6.
n6a_n7a_h_bond_dual_df <- merge(combined_df, a_hoogsteen_df[c('don_resi', 'don_chain',
                                                              'acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')],
                                by.x = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB'),
                                by.y = c('don_resi', 'don_chain', 'acc_pair_1_resi', 'acc_pair_1_chain',
                                         'model', 'PDB')) %>%
  filter(type == 2, h_bond == 1, acc_name == "N7") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(n6a_n7a_h_bond_dual_df[output_columns], snakemake@output[["supplemental_data_S8"]], na = "NaN")

#### 5D - C AMINE DUAL DONATION AND GC BASE PAIR ####

# Extract the cytosines that interact with G(O6) and N(NPO).
cytosines_df <- combined_df %>% filter(don_resn == "C", type == 2, acc_pair_combined == "G(O6), N(NPO)") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Create new columns that specify the G(O6) acceptor.
cytosines_df[c('guanine_resi', 'guanine_chain')] <- NA
cytosines_df[cytosines_df['acc_pair_1_name'] == "O6", c('guanine_resi', 'guanine_chain')] <- 
  cytosines_df[cytosines_df['acc_pair_1_name'] == "O6", c('acc_pair_1_resi', 'acc_pair_1_chain')]
cytosines_df[cytosines_df['acc_pair_1_name'] != "O6", c('guanine_resi', 'guanine_chain')] <- 
  cytosines_df[cytosines_df['acc_pair_1_name'] != "O6", c('acc_pair_2_resi', 'acc_pair_2_chain')]

# Which guanines form a G(N2)-to-C(O2) H-bond?
n2g_o2c_df <- merge(combined_df,
                    cytosines_df[c('don_resi', 'don_chain', 'guanine_resi', 'guanine_chain', 'model', 'PDB')],
                    by.x = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB'),
                    by.y = c('don_resi', 'don_chain', 'guanine_resi', 'guanine_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_name == "O2") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Backtrack to identify the C's involved in GC base pairs.
base_pairing_cytosines_df <- merge(cytosines_df,
                                   n2g_o2c_df[c('don_resi', 'don_chain', 'model', 'PDB')],
                                   by.x = c('guanine_resi', 'guanine_chain', 'model', 'PDB'),
                                   by.y = c('don_resi', 'don_chain', 'model', 'PDB'))

# Write the data to a csv.
write.csv(base_pairing_cytosines_df[output_columns], snakemake@output[["supplemental_data_S9"]], na = "NaN")

#### 5E - G AMINE DUAL DONATION AND GC BASE PAIR ####

# Extract the guanines that interact with A(N3) and C(O2).
guanines_df <- combined_df %>% filter(don_resn == "G", type == 2, acc_pair_combined == "A(N3), C(O2)") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Create new columns that specify the C(O2) acceptor.
guanines_df[c('cytosine_resi', 'cytosine_chain')] <- NA
guanines_df[guanines_df['acc_pair_1_name'] == "O2", c('cytosine_resi', 'cytosine_chain')] <- 
  guanines_df[guanines_df['acc_pair_1_name'] == "O2", c('acc_pair_1_resi', 'acc_pair_1_chain')]
guanines_df[guanines_df['acc_pair_1_name'] != "O2", c('cytosine_resi', 'cytosine_chain')] <- 
  guanines_df[guanines_df['acc_pair_1_name'] != "O2", c('acc_pair_2_resi', 'acc_pair_2_chain')]

# Which cytosines form a C(N4)-to-G(O6) H-bond?
n4c_o6g_df <- merge(combined_df,
                    guanines_df[c('cytosine_resi', 'cytosine_chain', 'don_resi', 'don_chain', 'model', 'PDB')],
                    by.x = c('don_resi', 'don_chain', 'acc_resi', 'acc_chain', 'model', 'PDB'),
                    by.y = c('cytosine_resi', 'cytosine_chain', 'don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_name == "O6") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Backtrack to identify the G's involved in GC base pairs.
base_pairing_guanines_df <- merge(combined_df,
                                  n4c_o6g_df[c('acc_resi', 'acc_chain', 'model', 'PDB')],
                                  by.x = c('don_resi', 'don_chain', 'model', 'PDB'),
                                  by.y = c('acc_resi', 'acc_chain', 'model', 'PDB')) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(base_pairing_guanines_df[output_columns], snakemake@output[["supplemental_data_S10"]], na = "NaN")

# Create new columns that specify the A(N3) acceptor.
base_pairing_guanines_df[c('adenine_resi', 'adenine_chain')] <- NA
base_pairing_guanines_df[base_pairing_guanines_df['acc_pair_1_name'] == "N3", c('adenine_resi', 'adenine_chain')] <-
  base_pairing_guanines_df[base_pairing_guanines_df['acc_pair_1_name'] == "N3",
                           c('acc_pair_1_resi', 'acc_pair_1_chain')]
base_pairing_guanines_df[base_pairing_guanines_df['acc_pair_1_name'] != "N3", c('adenine_resi', 'adenine_chain')] <-
  base_pairing_guanines_df[base_pairing_guanines_df['acc_pair_1_name'] != "N3",
                           c('acc_pair_2_resi', 'acc_pair_2_chain')]

# Identify the adenines that accept an H-bond from the dual-donating G(N2) and that also bear a dual-donating N6.
n6a_dual_df <- merge(combined_df, base_pairing_guanines_df[c('adenine_resi', 'adenine_chain', 'model', 'PDB')],
                     by.x = c('don_resi', 'don_chain', 'model', 'PDB'),
                     by.y = c('adenine_resi', 'adenine_chain', 'model', 'PDB')) %>%
  filter(type == 2) %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(n6a_dual_df[output_columns], snakemake@output[["supplemental_data_S11"]], na = "NaN")

#### 5F - G DUAL-DONATING AMINE TO A(N7) AND NPO ####

# Extract the guanine amines that dual donate to A(N7) and an NPO. The acceptors belong to the same residue.
g_dual_n7a_npo_same_df <- combined_df %>%
  filter(don_resn == "G", type == 2, acc_pair_combined == "A(N7), N(NPO)", same_resi == "True") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(g_dual_n7a_npo_same_df[output_columns], snakemake@output[["supplemental_data_S12"]], na = "NaN")

# Create new columns that specify the A(N7) acceptor.
g_dual_n7a_npo_same_df[c('adenine_resi', 'adenine_chain')] <- NA
g_dual_n7a_npo_same_df[g_dual_n7a_npo_same_df['acc_pair_1_name'] == "N7", c('adenine_resi', 'adenine_chain')] <-
  g_dual_n7a_npo_same_df[g_dual_n7a_npo_same_df['acc_pair_1_name'] == "N7",
                           c('acc_pair_1_resi', 'acc_pair_1_chain')]
g_dual_n7a_npo_same_df[g_dual_n7a_npo_same_df['acc_pair_1_name'] != "N7", c('adenine_resi', 'adenine_chain')] <-
  g_dual_n7a_npo_same_df[g_dual_n7a_npo_same_df['acc_pair_1_name'] != "N7",
                           c('acc_pair_2_resi', 'acc_pair_2_chain')]

# Identify the adenines that accept an H-bond from the dual-donating G(N2) and that also bear a dual-donating N6.
n6a_dual_df <- merge(combined_df, g_dual_n7a_npo_same_df[c('adenine_resi', 'adenine_chain', 'model', 'PDB')],
                     by.x = c('don_resi', 'don_chain', 'model', 'PDB'),
                     by.y = c('adenine_resi', 'adenine_chain', 'model', 'PDB')) %>%
  filter(type == 2) %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(n6a_dual_df[output_columns], snakemake@output[["supplemental_data_S13"]], na = "NaN")

#### S1A - C DUAL-DONATING AMINE TO AA(O) AND G(O6) ####

# Extract the cytosine amines that dual donate to AA(O) and G(O6).
c_aa_g_df <- combined_df %>% filter(don_resn == "C", type == 2, acc_pair_combined == "AA(O), G(O6)") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(c_aa_g_df[output_columns], snakemake@output[["supplemental_data_S14"]], na = "NaN")

#### S1B - G DUAL-DONATING AMINE TO AA(O) AND C(O2) ####

# Extract the guanine amines that dual donate to AA(O) and C(O2).
g_aa_c_df <- combined_df %>% filter(don_resn == "G", type == 2, acc_pair_combined == "AA(O), C(O2)") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(g_aa_c_df[output_columns], snakemake@output[["supplemental_data_S15"]], na = "NaN")

#### PSEUDO-TORSION LOCATION 1, LOCATION 2, AND NEIGHBORS ####

# Extract the rows from the combined data frame.
pair_1_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G") & type == 2) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Translate the pseudo-torsions to range from 0 to 360 degrees.
pair_1_df[which(pair_1_df$eta >= 0), "eta_translated"] <- pair_1_df[which(pair_1_df$eta >= 0), "eta"]
pair_1_df[which(pair_1_df$eta < 0), "eta_translated"] <- pair_1_df[which(pair_1_df$eta < 0), "eta"] + 360
pair_1_df[which(pair_1_df$theta >= 0), "theta_translated"] <- pair_1_df[which(pair_1_df$theta >= 0), "theta"]
pair_1_df[which(pair_1_df$theta < 0), "theta_translated"] <- pair_1_df[which(pair_1_df$theta < 0), "theta"] + 360

# Filter for location 1.
pair_1_df <- pair_1_df %>% filter(eta_translated >= 43.2 & eta_translated < 72 &
                                  theta_translated >= 151.2 & theta_translated < 180)

# Create a column that contains the donor atom name along with the donor residue name.
pair_1_df <- pair_1_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Extract the rows from the combined data frame.
pair_2_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G") & type == 2) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Translate the pseudo-torsions to range from 0 to 360 degrees.
pair_2_df[which(pair_2_df$eta >= 0), "eta_translated"] <- pair_2_df[which(pair_2_df$eta >= 0), "eta"]
pair_2_df[which(pair_2_df$eta < 0), "eta_translated"] <- pair_2_df[which(pair_2_df$eta < 0), "eta"] + 360
pair_2_df[which(pair_2_df$theta >= 0), "theta_translated"] <- pair_2_df[which(pair_2_df$theta >= 0), "theta"]
pair_2_df[which(pair_2_df$theta < 0), "theta_translated"] <- pair_2_df[which(pair_2_df$theta < 0), "theta"] + 360

# Filter for location 2.
pair_2_df <- pair_2_df %>% filter(eta_translated >= 295.2 & eta_translated < 324 &
                                  theta_translated >= 14.4 & theta_translated < 43.2)

# Create a column that contains the donor atom name along with the donor residue name.
pair_2_df <- pair_2_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Create two new columns.
pair_1_df <- pair_1_df %>% mutate(acc_pair_combined_reformat = gsub(", ", "/", acc_pair_combined, fixed = TRUE))
pair_2_df <- pair_2_df %>% mutate(acc_pair_combined_reformat = gsub(", ", "/", acc_pair_combined, fixed = TRUE))
pair_1_df["location"] <- "Location 1"
pair_2_df["location"] <- "Location 2"

# Create csv files for the residues within Location 1 and Location 2.
write.csv(pair_1_df[output_columns], snakemake@output[["supplemental_data_S16"]], quote = FALSE, na = "NaN")
write.csv(pair_2_df[output_columns], snakemake@output[["supplemental_data_S17"]], quote = FALSE, na = "NaN")

# Create a csv file with residues from Location 1 that are connected to the 5'-end of any residue from Location 2. This
# approach may not work if an insertion code is used for this or the downstream residue. If this is the case, print an
# error message, and create an empty csv file.
if (any(grepl("[a-zA-Z]", pair_1_df$don_resi)) | any(grepl("[a-zA-Z]", pair_2_df$don_resi))) {
  cat("Error: At least one residue within one of the pseudo-torsion locations contains an insertion code. The file
      named Supplemental_Data_S18.csv will be empty.")
  write.csv(data.frame(), snakemake@output[["supplemental_data_S18"]], quote = FALSE, na = "NaN")
} else {
  neighbors_df <- merge(pair_1_df %>% mutate(downstream_resi = as.character(as.numeric(don_resi) + 1)), pair_2_df,
                        by.x = c("downstream_resi", "don_chain", "eq_class_member"),
                        by.y = c("don_resi", "don_chain", "eq_class_member"), suffixes = c("",".y"))
  write.csv(neighbors_df[output_columns], snakemake@output[["supplemental_data_S18"]], quote = FALSE, na = "NaN")
}
