#!/usr/bin/env Rscript

#Load required libraries.
print('Loading required libraries...')
library(SingleCellExperiment)
library(dreamlet)
library(ggplot2)

#Specify settings.
num_cores <- as.integer(Sys.getenv("NSLOTS", "1")) #script will get the number of cores to use from NSLOTS variable if you are submitting as a job or will use 1 if not
dir_path <- "/wynton/group/pollen/reedmcmullen/projects/NEMP25/dreamlet_NEMP25/"
fig_path <- "/wynton/group/pollen/reedmcmullen/projects/NEMP25/dreamlet_NEMP25/figures/"

#Load voom results.
print('Loading voom results...')
voom_rds_file <- paste0(dir_path, "NEMP25_Tel_res_voom.rds")
res_voom <- readRDS(voom_rds_file)

#Contrasts for the effect of a perturbation for a species at a time point.
#c('FGF8b', 'EGF', 'BMP7', 'IGF1LR3', 'FGF2', 'NRG1', 'SHH', 'JAG1', 'DLL1')
#c('LY411575', 'HX531', 'Verteporfin', 'LDN193189', 'Rapamycin', 'ATRA', 'BGJ398', 'TRULI', 'WNTC59', 'AZD8931', 'CHIR99021', 'SANT1')

perty_by_timepoint_by_species_cons <- c(
# Vehicle perturbations
TreMan_6hr_Human = '(groupHuman_6hr_TreMan - groupHuman_6hr_Untreated)',
TreMan_54hr_Human = '(groupHuman_54hr_TreMan - groupHuman_54hr_Untreated)',
TreMan_7day_Human = '(groupHuman_7day_TreMan - groupHuman_7day_Untreated)',
DMSO_6hr_Human = '(groupHuman_6hr_DMSO - groupHuman_6hr_Untreated)',
DMSO_54hr_Human = '(groupHuman_54hr_DMSO - groupHuman_54hr_Untreated)',
DMSO_7day_Human = '(groupHuman_7day_DMSO - groupHuman_7day_Untreated)',
TreMan_6hr_Chimp = '(groupChimp_6hr_TreMan - groupChimp_6hr_Untreated)',
TreMan_54hr_Chimp = '(groupChimp_54hr_TreMan - groupChimp_54hr_Untreated)',
TreMan_7day_Chimp = '(groupChimp_7day_TreMan - groupChimp_7day_Untreated)',
DMSO_6hr_Chimp = '(groupChimp_6hr_DMSO - groupChimp_6hr_Untreated)',
DMSO_54hr_Chimp = '(groupChimp_54hr_DMSO - groupChimp_54hr_Untreated)',
DMSO_7day_Chimp = '(groupChimp_7day_DMSO - groupChimp_7day_Untreated)',

# TreMan controlled perturbations
FGF8b_6hr_Human = '(groupHuman_6hr_FGF8b - groupHuman_6hr_TreMan)',
FGF8b_54hr_Human = '(groupHuman_54hr_FGF8b - groupHuman_54hr_TreMan)',
FGF8b_7day_Human = '(groupHuman_7day_FGF8b - groupHuman_7day_TreMan)',
EGF_6hr_Human = '(groupHuman_6hr_EGF - groupHuman_6hr_TreMan)',
EGF_54hr_Human = '(groupHuman_54hr_EGF - groupHuman_54hr_TreMan)',
EGF_7day_Human = '(groupHuman_7day_EGF - groupHuman_7day_TreMan)',
BMP7_6hr_Human = '(groupHuman_6hr_BMP7 - groupHuman_6hr_TreMan)',
BMP7_54hr_Human = '(groupHuman_54hr_BMP7 - groupHuman_54hr_TreMan)',
BMP7_7day_Human = '(groupHuman_7day_BMP7 - groupHuman_7day_TreMan)',
IGF1LR3_6hr_Human = '(groupHuman_6hr_IGF1LR3 - groupHuman_6hr_TreMan)',
IGF1LR3_54hr_Human = '(groupHuman_54hr_IGF1LR3 - groupHuman_54hr_TreMan)',
IGF1LR3_7day_Human = '(groupHuman_7day_IGF1 - groupHuman_7day_TreMan)',
FGF2_6hr_Human = '(groupHuman_6hr_FGF2 - groupHuman_6hr_TreMan)',
FGF2_54hr_Human = '(groupHuman_54hr_FGF2 - groupHuman_54hr_TreMan)',
FGF2_7day_Human = '(groupHuman_7day_FGF2 - groupHuman_7day_TreMan)',
NRG1_6hr_Human = '(groupHuman_6hr_NRG1 - groupHuman_6hr_TreMan)',
NRG1_54hr_Human = '(groupHuman_54hr_NRG1 - groupHuman_54hr_TreMan)',
NRG1_7day_Human = '(groupHuman_7day_NRG1 - groupHuman_7day_TreMan)',
SHH_6hr_Human = '(groupHuman_6hr_SHH - groupHuman_6hr_TreMan)',
SHH_54hr_Human = '(groupHuman_54hr_SHH - groupHuman_54hr_TreMan)',
SHH_7day_Human = '(groupHuman_7day_SHH - groupHuman_7day_TreMan)',
JAG1_6hr_Human = '(groupHuman_6hr_JAG1 - groupHuman_6hr_TreMan)',
JAG1_54hr_Human = '(groupHuman_54hr_JAG1 - groupHuman_54hr_TreMan)',
JAG1_7day_Human = '(groupHuman_7day_JAG1 - groupHuman_7day_TreMan)',
DLL1_6hr_Human = '(groupHuman_6hr_DLL1 - groupHuman_6hr_TreMan)',
DLL1_54hr_Human = '(groupHuman_54hr_DLL1 - groupHuman_54hr_TreMan)',
DLL1_7day_Human = '(groupHuman_7day_DLL1 - groupHuman_7day_TreMan)',

FGF8b_6hr_Chimp = '(groupChimp_6hr_FGF8b - groupChimp_6hr_TreMan)',
FGF8b_54hr_Chimp = '(groupChimp_54hr_FGF8b - groupChimp_54hr_TreMan)',
FGF8b_7day_Chimp = '(groupChimp_7day_FGF8b - groupChimp_7day_TreMan)',
EGF_6hr_Chimp = '(groupChimp_6hr_EGF - groupChimp_6hr_TreMan)',
EGF_54hr_Chimp = '(groupChimp_54hr_EGF - groupChimp_54hr_TreMan)',
EGF_7day_Chimp = '(groupChimp_7day_EGF - groupChimp_7day_TreMan)',
BMP7_6hr_Chimp = '(groupChimp_6hr_BMP7 - groupChimp_6hr_TreMan)',
BMP7_54hr_Chimp = '(groupChimp_54hr_BMP7 - groupChimp_54hr_TreMan)',
BMP7_7day_Chimp = '(groupChimp_7day_BMP7 - groupChimp_7day_TreMan)',
IGF1LR3_6hr_Chimp = '(groupChimp_6hr_IGF1LR3 - groupChimp_6hr_TreMan)',
IGF1LR3_54hr_Chimp = '(groupChimp_54hr_IGF1LR3 - groupChimp_54hr_TreMan)',
IGF1LR3_7day_Chimp = '(groupChimp_7day_IGF1 - groupChimp_7day_TreMan)',
FGF2_6hr_Chimp = '(groupChimp_6hr_FGF2 - groupChimp_6hr_TreMan)',
FGF2_54hr_Chimp = '(groupChimp_54hr_FGF2 - groupChimp_54hr_TreMan)',
FGF2_7day_Chimp = '(groupChimp_7day_FGF2 - groupChimp_7day_TreMan)',
NRG1_6hr_Chimp = '(groupChimp_6hr_NRG1 - groupChimp_6hr_TreMan)',
NRG1_54hr_Chimp = '(groupChimp_54hr_NRG1 - groupChimp_54hr_TreMan)',
NRG1_7day_Chimp = '(groupChimp_7day_NRG1 - groupChimp_7day_TreMan)',
SHH_6hr_Chimp = '(groupChimp_6hr_SHH - groupChimp_6hr_TreMan)',
SHH_54hr_Chimp = '(groupChimp_54hr_SHH - groupChimp_54hr_TreMan)',
SHH_7day_Chimp = '(groupChimp_7day_SHH - groupChimp_7day_TreMan)',
JAG1_6hr_Chimp = '(groupChimp_6hr_JAG1 - groupChimp_6hr_TreMan)',
JAG1_54hr_Chimp = '(groupChimp_54hr_JAG1 - groupChimp_54hr_TreMan)',
JAG1_7day_Chimp = '(groupChimp_7day_JAG1 - groupChimp_7day_TreMan)',
DLL1_6hr_Chimp = '(groupChimp_6hr_DLL1 - groupChimp_6hr_TreMan)',
DLL1_54hr_Chimp = '(groupChimp_54hr_DLL1 - groupChimp_54hr_TreMan)',
DLL1_7day_Chimp = '(groupChimp_7day_DLL1 - groupChimp_7day_TreMan)',

# DMSO controlled perturbations
LY411575_6hr_Human = '(groupHuman_6hr_LY411575 - groupHuman_6hr_DMSO)',
LY411575_54hr_Human = '(groupHuman_54hr_LY411575 - groupHuman_54hr_DMSO)',
LY411575_7day_Human = '(groupHuman_7day_LY411575 - groupHuman_7day_DMSO)',
HX531_6hr_Human = '(groupHuman_6hr_HX531 - groupHuman_6hr_DMSO)',
HX531_54hr_Human = '(groupHuman_54hr_HX531 - groupHuman_54hr_DMSO)',
HX531_7day_Human = '(groupHuman_7day_HX531 - groupHuman_7day_DMSO)',
Verteporfin_6hr_Human = '(groupHuman_6hr_Verteporfin - groupHuman_6hr_DMSO)',
Verteporfin_54hr_Human = '(groupHuman_54hr_Verteporfin - groupHuman_54hr_DMSO)',
Verteporfin_7day_Human = '(groupHuman_7day_Verteporfin - groupHuman_7day_DMSO)',
LDN193189_6hr_Human = '(groupHuman_6hr_LDN193189 - groupHuman_6hr_DMSO)',
LDN193189_54hr_Human = '(groupHuman_54hr_LDN193189 - groupHuman_54hr_DMSO)',
LDN193189_7day_Human = '(groupHuman_7day_LDN193189 - groupHuman_7day_DMSO)',
Rapamycin_6hr_Human = '(groupHuman_6hr_Rapamycin - groupHuman_6hr_DMSO)',
Rapamycin_54hr_Human = '(groupHuman_54hr_Rapamycin - groupHuman_54hr_DMSO)',
Rapamycin_7day_Human = '(groupHuman_7day_Rapamycin - groupHuman_7day_DMSO)',
ATRA_6hr_Human = '(groupHuman_6hr_ATRA - groupHuman_6hr_DMSO)',
ATRA_54hr_Human = '(groupHuman_54hr_ATRA - groupHuman_54hr_DMSO)',
ATRA_7day_Human = '(groupHuman_7day_ATRA - groupHuman_7day_DMSO)',
BGJ398_6hr_Human = '(groupHuman_6hr_BGJ398 - groupHuman_6hr_DMSO)',
BGJ398_54hr_Human = '(groupHuman_54hr_BGJ398 - groupHuman_54hr_DMSO)',
BGJ398_7day_Human = '(groupHuman_7day_BGJ398 - groupHuman_7day_DMSO)',
TRULI_6hr_Human = '(groupHuman_6hr_TRULI - groupHuman_6hr_DMSO)',
TRULI_54hr_Human = '(groupHuman_54hr_TRULI - groupHuman_54hr_DMSO)',
TRULI_7day_Human = '(groupHuman_7day_TRULI - groupHuman_7day_DMSO)',
WNTC59_6hr_Human = '(groupHuman_6hr_WNTC59 - groupHuman_6hr_DMSO)',
WNTC59_54hr_Human = '(groupHuman_54hr_WNTC59 - groupHuman_54hr_DMSO)',
WNTC59_7day_Human = '(groupHuman_7day_WNTC59 - groupHuman_7day_DMSO)',
AZD8931_6hr_Human = '(groupHuman_6hr_AZD8931 - groupHuman_6hr_DMSO)',
AZD8931_54hr_Human = '(groupHuman_54hr_AZD8931 - groupHuman_54hr_DMSO)',
AZD8931_7day_Human = '(groupHuman_7day_AZD8931 - groupHuman_7day_DMSO)',
CHIR99021_6hr_Human = '(groupHuman_6hr_CHIR99021 - groupHuman_6hr_DMSO)',
CHIR99021_54hr_Human = '(groupHuman_54hr_CHIR99021 - groupHuman_54hr_DMSO)',
CHIR99021_7day_Human = '(groupHuman_7day_CHIR99021 - groupHuman_7day_DMSO)',
SANT1_6hr_Human = '(groupHuman_6hr_SANT1 - groupHuman_6hr_DMSO)',
SANT1_54hr_Human = '(groupHuman_54hr_SANT1 - groupHuman_54hr_DMSO)',
SANT1_7day_Human = '(groupHuman_7day_SANT1 - groupHuman_7day_DMSO)',

LY411575_6hr_Chimp = '(groupChimp_6hr_LY411575 - groupChimp_6hr_DMSO)',
LY411575_54hr_Chimp = '(groupChimp_54hr_LY411575 - groupChimp_54hr_DMSO)',
LY411575_7day_Chimp = '(groupChimp_7day_LY411575 - groupChimp_7day_DMSO)',
HX531_6hr_Chimp = '(groupChimp_6hr_HX531 - groupChimp_6hr_DMSO)',
HX531_54hr_Chimp = '(groupChimp_54hr_HX531 - groupChimp_54hr_DMSO)',
HX531_7day_Chimp = '(groupChimp_7day_HX531 - groupChimp_7day_DMSO)',
Verteporfin_6hr_Chimp = '(groupChimp_6hr_Verteporfin - groupChimp_6hr_DMSO)',
Verteporfin_54hr_Chimp = '(groupChimp_54hr_Verteporfin - groupChimp_54hr_DMSO)',
Verteporfin_7day_Chimp = '(groupChimp_7day_Verteporfin - groupChimp_7day_DMSO)',
LDN193189_6hr_Chimp = '(groupChimp_6hr_LDN193189 - groupChimp_6hr_DMSO)',
LDN193189_54hr_Chimp = '(groupChimp_54hr_LDN193189 - groupChimp_54hr_DMSO)',
LDN193189_7day_Chimp = '(groupChimp_7day_LDN193189 - groupChimp_7day_DMSO)',
Rapamycin_6hr_Chimp = '(groupChimp_6hr_Rapamycin - groupChimp_6hr_DMSO)',
Rapamycin_54hr_Chimp = '(groupChimp_54hr_Rapamycin - groupChimp_54hr_DMSO)',
Rapamycin_7day_Chimp = '(groupChimp_7day_Rapamycin - groupChimp_7day_DMSO)',
ATRA_6hr_Chimp = '(groupChimp_6hr_ATRA - groupChimp_6hr_DMSO)',
ATRA_54hr_Chimp = '(groupChimp_54hr_ATRA - groupChimp_54hr_DMSO)',
ATRA_7day_Chimp = '(groupChimp_7day_ATRA - groupChimp_7day_DMSO)',
BGJ398_6hr_Chimp = '(groupChimp_6hr_BGJ398 - groupChimp_6hr_DMSO)',
BGJ398_54hr_Chimp = '(groupChimp_54hr_BGJ398 - groupChimp_54hr_DMSO)',
BGJ398_7day_Chimp = '(groupChimp_7day_BGJ398 - groupChimp_7day_DMSO)',
TRULI_6hr_Chimp = '(groupChimp_6hr_TRULI - groupChimp_6hr_DMSO)',
TRULI_54hr_Chimp = '(groupChimp_54hr_TRULI - groupChimp_54hr_DMSO)',
TRULI_7day_Chimp = '(groupChimp_7day_TRULI - groupChimp_7day_DMSO)',
WNTC59_6hr_Chimp = '(groupChimp_6hr_WNTC59 - groupChimp_6hr_DMSO)',
WNTC59_54hr_Chimp = '(groupChimp_54hr_WNTC59 - groupChimp_54hr_DMSO)',
WNTC59_7day_Chimp = '(groupChimp_7day_WNTC59 - groupChimp_7day_DMSO)',
AZD8931_6hr_Chimp = '(groupChimp_6hr_AZD8931 - groupChimp_6hr_DMSO)',
AZD8931_54hr_Chimp = '(groupChimp_54hr_AZD8931 - groupChimp_54hr_DMSO)',
AZD8931_7day_Chimp = '(groupChimp_7day_AZD8931 - groupChimp_7day_DMSO)',
CHIR99021_6hr_Chimp = '(groupChimp_6hr_CHIR99021 - groupChimp_6hr_DMSO)',
CHIR99021_54hr_Chimp = '(groupChimp_54hr_CHIR99021 - groupChimp_54hr_DMSO)',
CHIR99021_7day_Chimp = '(groupChimp_7day_CHIR99021 - groupChimp_7day_DMSO)',
SANT1_6hr_Chimp = '(groupChimp_6hr_SANT1 - groupChimp_6hr_DMSO)',
SANT1_54hr_Chimp = '(groupChimp_54hr_SANT1 - groupChimp_54hr_DMSO)',
SANT1_7day_Chimp = '(groupChimp_7day_SANT1 - groupChimp_7day_DMSO)'
)

#Contrasts for the species differences in the effect of a perturbation at a time point.
perturbations <- c('TreMan', 'DMSO', 'FGF8b', 'EGF', 'BMP7', 'IGF1LR3', 'FGF2', 'NRG1', 'SHH', 'JAG1', 'DLL1', 'LY411575', 'HX531', 'Verteporfin', 'LDN193189', 'Rapamycin', 'ATRA', 'BGJ398', 'TRULI', 'WNTC59', 'AZD8931', 'CHIR99021', 'SANT1')
pert_by_timepoint_speciesdiff_cons <- c() #Initialize empty vector.
for (perturbation in perturbations) { #loop through the perturbations
  # Create the line based on the compound
  con1 <- paste0(perturbation, "_6hr_SpeciesDiff = paste0('(', ", perturbation, "_6hr_Human, '-', ", perturbation, "_6hr_Chimp, ')')")
  con2 <- paste0(perturbation, "_54hr_SpeciesDiff = paste0('(', ", perturbation, "_54hr_Human, '-', ", perturbation, "_54hr_Chimp, ')')")
  con3 <- paste0(perturbation, "_7day_SpeciesDiff = paste0('(', ", perturbation, "_7day_Human, '-', ", perturbation, "_7day_Chimp, ')')")
  pert_by_timepoint_speciesdiff_cons <- c(pert_by_timepoint_speciesdiff_cons, con1, con2, con3)
}

# Combine the contrasts.
contrasts <- c(perty_by_timepoint_by_species_cons, pert_by_timepoint_speciesdiff_cons)





