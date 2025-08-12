library(dplyr)

#### SET INITIAL VARIABLES ####

# Creates data frames from the combined data.
combined_df <- read.csv(snakemake@input[["combined"]], header = TRUE, na.strings = "NaN", comment.char = "#",
                        colClasses = c("PDB" = "character"))

#### 4A - SHEARED GA BASE PAIR ####

# Extract the guanines that interact with a dual donating A via their N3 and O2' atoms.
guanines_to_consider_df <- combined_df %>% filter(don_resn == "A", type == 2, 
                                                  acc_pair_combined == "G(N3), N(O2')", same_resi == "True") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Only keep relevant columns.
guanines_to_consider_df <- guanines_to_consider_df[c('don_resi', 'don_chain', 
                                                     'acc_pair_1_name', 'acc_pair_2_name', 'acc_pair_1_resn', 
                                                     'acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')]

# Rename select columns.
guanines_to_consider_df <- guanines_to_consider_df %>% rename(adenine_resi=don_resi, adenine_chain=don_chain)

# Which guanines form a G(N2)-to-A(N7) H-bond?
g_a_n7_h_bonding_df <- merge(guanines_to_consider_df, combined_df, 
                             by.x = c('adenine_resi', 'adenine_chain', 'acc_pair_1_resi', 'acc_pair_1_chain', 
                                      'model', 'PDB'),
                             by.y = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, don_name == "N2", acc_name == "N7") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Backtrack to identify the A's involved in the sheared GA base pair.
sheared_ga_df <- merge(g_a_n7_h_bonding_df, combined_df, 
                       by.x = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB'),
                       by.y = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')) %>% 
  filter(don_resn.y == "A", type.y == 2, acc_pair_combined.y == "G(N3), N(O2')", same_resi.y == "True") %>%
  distinct(don_index.y, eq_class_member.y, .keep_all = TRUE)

# Rename and remove select columns.
sheared_ga_df <- sheared_ga_df[c('don_index.y', 'don_name.y', 'don_resn.y', 'don_resi', 'don_chain', 'eta.y', 'theta.y',
                                 'acc_pair_1_name', 'acc_pair_1_resi', 'acc_pair_1_chain', 'acc_pair_2_name', 
                                 'acc_pair_2_resi.y', 'acc_pair_2_chain.y', 'model', 'PDB', 'eq_class_member.y', 
                                 'type.y')] %>% 
  rename(don_index=don_index.y, don_name=don_name.y, don_resn=don_resn.y, eta=eta.y, theta=theta.y,
         acc_pair_2_resi=acc_pair_2_resi.y, acc_pair_2_chain=acc_pair_2_chain.y, eq_class_member=eq_class_member.y,
         type=type.y)

# Write the data to a csv.
write.csv(sheared_ga_df, snakemake@output[["supplemental_data_S1"]], na = "NaN")

# Now, find which guanines form a G(N2)-to-A(N7) H-bond and bear a dual-donating N2.
g_a_n7_h_bonding_dual_df <- merge(guanines_to_consider_df, combined_df,
                                  by.x = c('adenine_resi', 'adenine_chain', 'acc_pair_1_resi', 'acc_pair_1_chain',
                                           'model', 'PDB'),
                                  by.y = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(type == 2, h_bond == 1, don_name == "N2", acc_name == "N7") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Rename and remove select columns.
g_a_n7_h_bonding_dual_df <- g_a_n7_h_bonding_dual_df[c('don_index', 'don_name', 'don_resn', 'acc_pair_1_resi',
                                                       'acc_pair_1_chain', 'eta', 'theta', 'acc_pair_1_name.y',
                                                       'acc_pair_1_resi.y', 'acc_pair_1_chain.y', 'acc_pair_2_name.y',
                                                       'acc_pair_2_resi', 'acc_pair_2_chain', 'model', 'PDB',
                                                       'eq_class_member', 'type')] %>%
  rename(don_resi=acc_pair_1_resi, don_chain=acc_pair_1_chain, acc_pair_1_name=acc_pair_1_name.y,
         acc_pair_1_resi=acc_pair_1_resi.y, acc_pair_1_chain=acc_pair_1_chain.y, acc_pair_2_name=acc_pair_2_name.y)

# Write the data to a csv.
write.csv(g_a_n7_h_bonding_dual_df, snakemake@output[["supplemental_data_S2"]], na = "NaN")

#### 4B - WC/H A-MINOR MOTIF ####

# Which guanines form a G(N2)-to-C(O2) H-bond?
g_c_h_bonding_df <- merge(guanines_to_consider_df, combined_df, 
                          by.x = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB'),
                          by.y = c('don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1 & acc_resn == "C" & acc_name == "O2") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Rename and only keep relevant columns.
g_c_h_bonding_df <- g_c_h_bonding_df %>% rename(guanine_resi=acc_pair_1_resi, guanine_chain=acc_pair_1_chain, 
                                                cytosine_resi=acc_resi, cytosine_chain=acc_chain)
g_c_h_bonding_df <- 
  g_c_h_bonding_df[c('guanine_resi', 'guanine_chain', 'cytosine_resi', 'cytosine_chain', 'model', 'PDB')]

# Which cytosines form a C(N4)-to-G(O6) H-bond? These cytosines also receive an H-bond at their O2 from the N2 of the 
# same guanine. This suggests that the C and G are forming a canonical base pair. Meanwhile, the G is engaged with a 
# dual donating A(N6) via its sugar edge. The entire unit is a WC/H A-minor motif.
c_g_h_bonding_df <- merge(g_c_h_bonding_df, combined_df, 
                          by.x = c('cytosine_resi', 'cytosine_chain', 'model', 'PDB'),
                          by.y = c('don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_resn == "G", acc_name == "O6") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Only keep relevant columns.
c_g_h_bonding_df <- 
  c_g_h_bonding_df[c('cytosine_resi', 'cytosine_chain', 'guanine_resi', 'guanine_chain', 'model', 'PDB')]

# Backtrack to identify the A's involved in the WC/H A-minor motifs.
a_minor_df <- merge(c_g_h_bonding_df, combined_df, 
                    by.x = c('guanine_resi', 'guanine_chain', 'model', 'PDB'),
                    by.y = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')) %>% 
  filter(don_resn == "A", type == 2, acc_pair_combined == "G(N3), N(O2')", same_resi == "True") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Rename select columns.
a_minor_df <- a_minor_df %>% rename(acc_pair_1_resi=guanine_resi, acc_pair_1_chain=guanine_chain)

# Write the data to a csv.
write.csv(a_minor_df, snakemake@output[["supplemental_data_S3"]], na = "NaN")

# Which guanines form a G(N2)-to-A(N1) H-bond?
g_n2_a_n1_h_bonding_df <- merge(a_minor_df, combined_df, 
                                by.x = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB'),
                                by.y = c('don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond.y == 1, acc_resn.y == "A", acc_name.y == "N1") %>% 
  distinct(don_index.y, eq_class_member.y, .keep_all = TRUE)

# Remove and rename select columns.
g_n2_a_n1_h_bonding_df <- g_n2_a_n1_h_bonding_df[c('acc_resi.y', 'acc_chain.y', 'model', 'PDB')] %>% 
  rename(adenine_resi=acc_resi.y, adenine_chain=acc_chain.y)

# Determine which of the adenines are involved in a WC/H A-minor motif.
g_n2_a_n1_h_bonding_filtered_df <- merge(a_minor_df, g_n2_a_n1_h_bonding_df,
                                         by.x = c('don_resi', 'don_chain', 'model', 'PDB'), 
                                         by.y = c('adenine_resi', 'adenine_chain', 'model', 'PDB')) %>% 
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(g_n2_a_n1_h_bonding_filtered_df, snakemake@output[["supplemental_data_S4"]], na = "NaN")

# Which guanines form a G(N2)-to-A(N7) H-bond?
g_n2_a_n7_h_bonding_df <- merge(a_minor_df, combined_df, 
                          by.x = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB'),
                          by.y = c('don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond.y == 1, acc_resn.y == "A", acc_name.y == "N7") %>% 
  distinct(don_index.y, eq_class_member.y, .keep_all = TRUE)

# Remove and rename select columns.
g_n2_a_n7_h_bonding_df <- g_n2_a_n7_h_bonding_df[c('acc_resi.y', 'acc_chain.y', 'model', 'PDB')] %>% 
  rename(adenine_resi=acc_resi.y, adenine_chain=acc_chain.y)

# Determine which of the adenines are involved in a WC/H A-minor motif.
g_n2_a_n7_h_bonding_filtered_df <- merge(a_minor_df, g_n2_a_n7_h_bonding_df,
                                         by.x = c('don_resi', 'don_chain', 'model', 'PDB'), 
                                         by.y = c('adenine_resi', 'adenine_chain', 'model', 'PDB')) %>% 
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(g_n2_a_n7_h_bonding_filtered_df, snakemake@output[["supplemental_data_S5"]], na = "NaN")

# Which of the guanines also bear a dual donating G(N2)?
dual_g_df <- merge(combined_df, c_g_h_bonding_df,
                   by.x = c('don_resi', 'don_chain', 'model', 'PDB'), 
                   by.y = c('guanine_resi', 'guanine_chain', 'model', 'PDB')) %>% 
  filter(don_resn == "G", type == 2) %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(dual_g_df, snakemake@output[["supplemental_data_S6"]], na = "NaN")

#### 4C - A DUAL DONATING AMINE TO A HOOGSTEEN ####

# Extract the adenines amines that dual donate to the N7 and NPO of another A.
a_hoogsteen_df <- combined_df %>% filter(don_resn == "A", type == 2, 
                                         acc_pair_combined == "A(N7), N(NPO)", same_resi == "True") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Only keep relevant columns.
a_hoogsteen_df <- a_hoogsteen_df[c('don_resi', 'don_chain', 'acc_pair_1_name', 'acc_pair_2_name', 'acc_pair_1_resn',
                                   'acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')]

# Rename select columns.
a_hoogsteen_df <- a_hoogsteen_df %>% rename(dual_don_A_resi=don_resi, dual_don_A_chain=don_chain)

# Which adenines bearing the acceptor pair form a A(N6)-to-A(N7) H-bond?
a_n6_a_n7_h_bonding_df <- merge(a_hoogsteen_df, combined_df,
                                by.x = c('dual_don_A_resi', 'dual_don_A_chain', 'acc_pair_1_resi', 'acc_pair_1_chain',
                                         'model', 'PDB'),
                                by.y = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, don_name == "N6", acc_name == "N7") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Backtrack to identify the A's involved in the tHH AA base pair.
tHH_aa_df <- merge(a_n6_a_n7_h_bonding_df, combined_df,
                   by.x = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB'),
                   by.y = c('acc_pair_1_resi', 'acc_pair_1_chain', 'model', 'PDB')) %>%
  filter(don_resn.y == "A", type.y == 2, acc_pair_combined.y == "A(N7), N(NPO)", same_resi.y == "True") %>%
  distinct(don_index.y, eq_class_member.y, .keep_all = TRUE)

# Rename and remove select columns.
tHH_aa_df <- tHH_aa_df[c('don_index.y', 'don_name.y', 'don_resi', 'don_chain', 'eta.y', 'theta.y', 'acc_pair_1_name',
                         'acc_pair_1_resi', 'acc_pair_1_chain', 'acc_pair_2_name', 'acc_pair_2_resi.y',
                         'acc_pair_2_chain.y', 'model', 'PDB', 'eq_class_member.y', 'type.y')] %>%
  rename(don_index=don_index.y, don_name=don_name.y, eta=eta.y, theta=theta.y, acc_pair_2_resi=acc_pair_2_resi.y,
         acc_pair_2_chain=acc_pair_2_chain.y, eq_class_member=eq_class_member.y, type=type.y)

# Write the data to a csv.
write.csv(tHH_aa_df, snakemake@output[["supplemental_data_S7"]], na = "NaN")

# Now, find which adenines form an A(N6)-to-A(N7) H-bond and bear a dual-donating N6.
a_n6_a_n7_h_bonding_dual_df <- merge(a_hoogsteen_df, combined_df,
                                     by.x = c('dual_don_A_resi', 'dual_don_A_chain', 'acc_pair_1_resi',
                                              'acc_pair_1_chain', 'model', 'PDB'),
                                     by.y = c('acc_resi', 'acc_chain', 'don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(type == 2, h_bond == 1, don_name == "N6", acc_name == "N7") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Rename and remove select columns.
a_n6_a_n7_h_bonding_dual_df <- a_n6_a_n7_h_bonding_dual_df[c('don_index', 'don_name', 'don_resn', 'acc_pair_1_resi',
                                                             'acc_pair_1_chain', 'eta', 'theta', 'acc_pair_1_name.y',
                                                             'acc_pair_1_resi.y', 'acc_pair_1_chain.y',
                                                             'acc_pair_2_name.y', 'acc_pair_2_resi',
                                                             'acc_pair_2_chain', 'model', 'PDB', 'eq_class_member',
                                                             'type')] %>%
  rename(don_resi=acc_pair_1_resi, don_chain=acc_pair_1_chain, acc_pair_1_name=acc_pair_1_name.y,
         acc_pair_1_resi=acc_pair_1_resi.y, acc_pair_1_chain=acc_pair_1_chain.y, acc_pair_2_name=acc_pair_2_name.y)

# Write the data to a csv.
write.csv(a_n6_a_n7_h_bonding_dual_df, snakemake@output[["supplemental_data_S8"]], na = "NaN")

#### 4D - C AMINE DUAL DONATION AND GC BASE PAIR ####

# Extract the cytosines that interact with G(O6) and N(NPO).
cytosines_df <- combined_df %>% filter(don_resn == "C", type == 2, acc_pair_combined == "G(O6), N(NPO)") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Create new columns that specify the G(O6) acceptor.
cytosines_df[c('guanine_resi', 'guanine_chain')] <- NA
cytosines_df[cytosines_df['acc_pair_1_name'] == "O6", c('guanine_resi', 'guanine_chain')] <- 
  cytosines_df[cytosines_df['acc_pair_1_name'] == "O6", c('acc_pair_1_resi', 'acc_pair_1_chain')]
cytosines_df[cytosines_df['acc_pair_1_name'] != "O6", c('guanine_resi', 'guanine_chain')] <- 
  cytosines_df[cytosines_df['acc_pair_1_name'] != "O6", c('acc_pair_2_resi', 'acc_pair_2_chain')]

# Only keep relevant columns.
cytosines_df <- cytosines_df[c('don_resi', 'don_chain', 'guanine_resi', 'guanine_chain', 'model', 'PDB')]

# Which guanines form a G(N2)-to-C(O2) H-bond?
g_n2_c_o2_df <- merge(cytosines_df, combined_df, 
                     by.x = c('guanine_resi', 'guanine_chain', 'model', 'PDB'),
                     by.y = c('don_resi', 'don_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_resn == "C", acc_name == "O2") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Backtrack to identify the C's involved in GC base pairs.
base_pairing_cytosines_df <- merge(g_n2_c_o2_df, cytosines_df, 
                                   by = c('guanine_resi', 'guanine_chain', 'model', 'PDB'))

# Make cytosine the focus of the data frame (e.g., the type column for each row should be associated with the cytosine, 
# not the guanine).
base_pairing_cytosines_df <- merge(combined_df, 
                                   base_pairing_cytosines_df[c('don_resi.y', 'don_chain.y', 'model', 'PDB')], 
                                   by.x = c('don_resi', 'don_chain', 'model', 'PDB'),
                                   by.y = c('don_resi.y', 'don_chain.y', 'model', 'PDB')) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(base_pairing_cytosines_df, snakemake@output[["supplemental_data_S9"]], na = "NaN")

#### 4E - G AMINE DUAL DONATION AND GC BASE PAIR ####

# Extract the guanines that interact with A(N3) and C(O2).
guanines_df <- combined_df %>% filter(don_resn == "G", type == 2, acc_pair_combined == "A(N3), C(O2)") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Create new columns that specify the C(O2) acceptor.
guanines_df[c('cytosine_resi', 'cytosine_chain')] <- NA
guanines_df[guanines_df['acc_pair_1_name'] == "O2", c('cytosine_resi', 'cytosine_chain')] <- 
  guanines_df[guanines_df['acc_pair_1_name'] == "O2", c('acc_pair_1_resi', 'acc_pair_1_chain')]
guanines_df[guanines_df['acc_pair_1_name'] != "O2", c('cytosine_resi', 'cytosine_chain')] <- 
  guanines_df[guanines_df['acc_pair_1_name'] != "O2", c('acc_pair_2_resi', 'acc_pair_2_chain')]

# Only keep relevant columns.
guanines_df <- guanines_df[c('don_resi', 'don_chain', 'cytosine_resi', 'cytosine_chain', 'model', 'PDB')] %>%
  rename(guanine_resi=don_resi, guanine_chain=don_chain)

# Which cytosines form a C(N4)-to-G(O6) H-bond?
c_n4_g_o6_df <- merge(guanines_df, combined_df, 
                      by.x = c('cytosine_resi', 'cytosine_chain', 'guanine_resi', 'guanine_chain', 'model', 'PDB'),
                      by.y = c('don_resi', 'don_chain', 'acc_resi', 'acc_chain', 'model', 'PDB')) %>%
  filter(h_bond == 1, acc_resn == "G", acc_name == "O6") %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Make guanines the focus of the data frame (e.g., the type column for each row should be associated with the guanine, 
# not the cytosine).
base_pairing_guanines_df <- merge(c_n4_g_o6_df[c('guanine_resi', 'guanine_chain', 'model', 'PDB')], combined_df,
                                   by.x = c('guanine_resi', 'guanine_chain', 'model', 'PDB'),
                                   by.y = c('don_resi', 'don_chain', 'model', 'PDB')) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE) %>% rename(don_resi=guanine_resi, don_chain=guanine_chain)

# Specify the columns to write to the output CSV file.
output_columns <- c('don_index', 'don_name', 'don_resn', 'don_resi', 'don_chain', 'eta', 'theta', 'chi', 'type',
                    'acc_pair_1_name', 'acc_pair_1_resi', 'acc_pair_1_chain',
                    'acc_pair_2_name', 'acc_pair_2_resi', 'acc_pair_2_chain',
                    'model', 'PDB', 'eq_class_member')

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
a_n6_dual_df <- merge(combined_df, base_pairing_guanines_df[c('adenine_resi', 'adenine_chain', 'model', 'PDB')],
                      by.x = c('don_resi', 'don_chain', 'model', 'PDB'),
                      by.y = c('adenine_resi', 'adenine_chain', 'model', 'PDB')) %>%
  filter(type == 2) %>% distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(a_n6_dual_df[output_columns], snakemake@output[["supplemental_data_S11"]], na = "NaN")

#### 4F - G DUAL DONATING AMINE TO A(N7) AND NPO ####

# Extract the guanine amines that dual donate to A(N7) and an NPO. The acceptors belong to the same residue.
g_a_n7_npo_same_df <- combined_df %>% filter(don_resn == "G", type == 2, 
                                             acc_pair_combined == "A(N7), N(NPO)", same_resi == "True") %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(g_a_n7_npo_same_df, snakemake@output[["supplemental_data_S12"]], na = "NaN")

#### S1A - C DUAL DONATING AMINE TO AA(O) AND G(O6) ####

# Extract the cytosine amines that dual donate to AA(O) and G(O6).
c_g_aa_df <- combined_df %>% filter(don_resn == "C", type == 2, acc_pair_combined == "AA(O), G(O6)") %>% 
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(c_g_aa_df, snakemake@output[["supplemental_data_S13"]], na = "NaN")

#### S1B - G DUAL DONATING AMINE TO AA(O) AND C(O2) ####

# Extract the guanine amines that dual donate to AA(O) and C(O2).
g_c_aa_df <- combined_df %>% filter(don_resn == "G", type == 2, acc_pair_combined == "AA(O), C(O2)") %>% 
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Write the data to a csv.
write.csv(g_c_aa_df, snakemake@output[["supplemental_data_S14"]], na = "NaN")

#### PSEUDO-TORSION LOCATION 1, LOCATION 2, AND NEIGHBORS ####

# Extract the rows from the combined data frame relevant for these plots.
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

# Extract the rows from the combined data frame relevant for these plots.
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
write.csv(pair_1_df[,c("don_index", "don_name", "don_resn", "don_resi", "don_chain",
                       "acc_pair_1_name", "acc_pair_1_resi", "acc_pair_1_chain",
                       "acc_pair_2_name", "acc_pair_2_resi", "acc_pair_2_chain",
                       "don_label", "acc_pair_combined_reformat", "same_resi", "type",
                       "PDB", "model", "eq_class_member",
                       "eta", "theta", "eta_translated", "theta_translated", "chi")],
          snakemake@output[["supplemental_data_S15"]], quote = FALSE, na = "NaN", row.names = FALSE)
write.csv(pair_2_df[,c("don_index", "don_name", "don_resn", "don_resi", "don_chain",
                       "acc_pair_1_name", "acc_pair_1_resi", "acc_pair_1_chain",
                       "acc_pair_2_name", "acc_pair_2_resi", "acc_pair_2_chain",
                       "don_label", "acc_pair_combined_reformat", "same_resi", "type",
                       "PDB", "model", "eq_class_member",
                       "eta", "theta", "eta_translated", "theta_translated", "chi")],
          snakemake@output[["supplemental_data_S16"]], quote = FALSE, na = "NaN", row.names = FALSE)

# Create a csv file with residues from Location 1 that are connected to the 5'-end of any residue from Location 2. This
# approach may not work if an insertion code is used for this or the downstream residue. If this is the case, print an
# error message, and create an empty csv file.
if (any(grepl("[a-zA-Z]", pair_1_df$don_resi)) | any(grepl("[a-zA-Z]", pair_2_df$don_resi))) {
  cat("Error: At least one residue within one of the pseudo-torsion locations contains an insertion code. The file
      named Supplemental_Data_S14.csv will be empty.")
  write.csv(data.frame(), snakemake@output[["supplemental_data_S17"]], quote = FALSE, na = "NaN", row.names = FALSE)
} else {
  neighbors_df <- merge(pair_1_df %>% mutate(downstream_resi = as.character(as.numeric(don_resi) + 1)), pair_2_df,
                        by.x = c("downstream_resi", "don_chain", "eq_class_member"),
                        by.y = c("don_resi", "don_chain", "eq_class_member"), suffixes = c("",".y"))
  neighbors_df <- neighbors_df[,c("don_index", "don_name", "don_resn", "don_resi", "don_chain",
                                  "acc_pair_1_name", "acc_pair_1_resi", "acc_pair_1_chain",
                                  "acc_pair_2_name", "acc_pair_2_resi", "acc_pair_2_chain",
                                  "don_label", "acc_pair_combined_reformat", "same_resi", "type",
                                  "PDB", "model", "eq_class_member",
                                  "eta", "theta", "eta_translated", "theta_translated", "chi")]
  write.csv(neighbors_df, snakemake@output[["supplemental_data_S17"]], quote = FALSE, na = "NaN", row.names = FALSE)
}
