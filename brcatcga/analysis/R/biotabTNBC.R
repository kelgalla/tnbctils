
#---------------------------------------------------------------------
# FILE     : biotabTNBC.R
# AUTHOR   : Kelly E. Craven
# CREATED  : 2017-05-20
# COMMENTS : reads the biotab file and filters it down to only TNBC 
#            cases
#---------------------------------------------------------------------
library(plyr)

# legacy
dir <- "~/brca/brcatcga/dataLegacyReorg/"
biotabDir <- paste0(dir, "biotab/")
clinptFile <- "nationwidechildrens.org_clinical_patient_brca.735bc5ff-86d1-421a-8693-6e6f92055563.txt"
clinpt <- read.delim(paste0(biotabDir, clinptFile), header=F, stringsAsFactors=F)

colnames(clinpt) <- clinpt[2,]
clinpt <- clinpt[c(4:nrow(clinpt)),]
clinpt <- clinpt[,c("bcr_patient_uuid", "bcr_patient_barcode",
                    "breast_carcinoma_estrogen_receptor_status", "er_level_cell_percentage_category",
                    "breast_carcinoma_progesterone_receptor_status",
                    "progesterone_receptor_level_cell_percent_category",
                    "lab_proc_her2_neu_immunohistochemistry_receptor_status", 
                    "her2_erbb_pos_finding_cell_percent_category", "her2_immunohistochemistry_level_result",
                    "lab_procedure_her2_neu_in_situ_hybrid_outcome_type")]

# new
dir <- "~/brca/brcatcga/analysis/R/"
biotabDir <- "~/brca/brcatcga/analysis/bin/"
clinptFile <- "biotab.txt"
clinpt <- read.delim(paste0(biotabDir, clinptFile), header=T, stringsAsFactors=F)

# figure out ER status
erStatus <- ddply(clinpt, .(clinpt$breast_carcinoma_estrogen_receptor_status, 
                            clinpt$er_level_cell_percentage_category), nrow)
names(erStatus) <- c("breast_carcinoma_estrogen_receptor_status", "er_level_cell_percentage_category", "freq")
#erStatus <- unique(clinpt[,c("breast_carcinoma_estrogen_receptor_status", "er_level_cell_percentage_category")])
write.table(erStatus, paste0(dir,"erStatus.txt"), sep="\t", quote=F, row.names=F)

# identify ER positive 1% cut off
erPos1 <- clinpt[clinpt$breast_carcinoma_estrogen_receptor_status=="Positive" | 
                  (clinpt$breast_carcinoma_estrogen_receptor_status=="Negative" & 
                     clinpt$er_level_cell_percentage_category=="<10%") |
                   (clinpt$breast_carcinoma_estrogen_receptor_status=="Indeterminate" &
                      clinpt$er_level_cell_percentage_category=="20-29%"),]
write.table(erPos1, paste0(dir,"erPos1.txt"), sep="\t", quote=F, row.names=F)

erPos10 <- clinpt[(clinpt$breast_carcinoma_estrogen_receptor_status=="Positive" &
                     !(clinpt$er_level_cell_percentage_category %in% c("<10%"))) |
                    (clinpt$breast_carcinoma_estrogen_receptor_status=="Indeterminate" &
                       clinpt$er_level_cell_percentage_category=="20-29%"),]
write.table(erPos10, paste0(dir,"erPos10.txt"), sep="\t", quote=F, row.names=F)

# figure out PR status
prStatus <- ddply(clinpt, .(clinpt$breast_carcinoma_progesterone_receptor_status, 
                            clinpt$progesterone_receptor_level_cell_percent_category), nrow)
names(prStatus) <- c("breast_carcinoma_progesterone_receptor_status", 
                     "progesterone_receptor_level_cell_percent_category", "freq")
#prStatus <- unique(clinpt[,c("breast_carcinoma_progesterone_receptor_status", "progesterone_receptor_level_cell_percent_category")])
write.table(prStatus, paste0(dir,"prStatus.txt"), sep="\t", quote=F, row.names=F)

# identify PR positive 1% cut off
prPos1 <- clinpt[clinpt$breast_carcinoma_progesterone_receptor_status=="Positive" | 
                   (clinpt$breast_carcinoma_progesterone_receptor_status=="Negative" & 
                      clinpt$progesterone_receptor_level_cell_percent_category=="<10%"),]
write.table(prPos1, paste0(dir,"prPos1.txt"), sep="\t", quote=F, row.names=F)

prPos10 <- clinpt[(clinpt$breast_carcinoma_progesterone_receptor_status=="Positive" &
                     !(clinpt$progesterone_receptor_level_cell_percent_category %in% c("<10%"))),]
write.table(prPos10, paste0(dir,"prPos10.txt"), sep="\t", quote=F, row.names=F)

# er and pr positive 1
erprPos1 <- merge(erPos1, prPos1, by=c("bcr_patient_uuid",
                                       "bcr_patient_barcode",
                                       "breast_carcinoma_estrogen_receptor_status",
                                       "er_level_cell_percentage_category",
                                       "breast_carcinoma_progesterone_receptor_status",
                                       "progesterone_receptor_level_cell_percent_category",
                                       "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                                       "her2_erbb_pos_finding_cell_percent_category",
                                       "her2_immunohistochemistry_level_result",
                                       "lab_procedure_her2_neu_in_situ_hybrid_outcome_type"))
write.table(erprPos1, paste0(dir,"erprPos1.txt"), sep="\t", quote=F, row.names=F)

# er and pr positive 10
erprPos10 <- merge(erPos10, prPos10, by=c("bcr_patient_uuid",
                                       "bcr_patient_barcode",
                                       "breast_carcinoma_estrogen_receptor_status",
                                       "er_level_cell_percentage_category",
                                       "breast_carcinoma_progesterone_receptor_status",
                                       "progesterone_receptor_level_cell_percent_category",
                                       "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                                       "her2_erbb_pos_finding_cell_percent_category",
                                       "her2_immunohistochemistry_level_result",
                                       "lab_procedure_her2_neu_in_situ_hybrid_outcome_type"))
write.table(erprPos10, paste0(dir,"erprPos10.txt"), sep="\t", quote=F, row.names=F)

# find ER negative cases
#erNeg <- clinpt[clinpt$breast_carcinoma_estrogen_receptor_status=="Negative" | 
#                  clinpt$er_level_cell_percentage_category=="<10%",]
erNeg <- clinpt[clinpt$breast_carcinoma_estrogen_receptor_status=="Negative",]

# take ER negative cases and only keep PR negative cases
#erprNeg <- erNeg[erNeg$breast_carcinoma_progesterone_receptor_status=="Negative" |
#                   erNeg$progesterone_receptor_level_cell_percent_category=="<10%",]
erprNeg <- erNeg[erNeg$breast_carcinoma_progesterone_receptor_status=="Negative",]

tnbc <- erprNeg[erprNeg$lab_procedure_her2_neu_in_situ_hybrid_outcome_type %in% 
                  c("Negative", "[Not Evaluated]", "[Not Available]", ""),]
tnbc <- tnbc[(tnbc$lab_procedure_her2_neu_in_situ_hybrid_outcome_type=="Negative" &
                tnbc$lab_proc_her2_neu_immunohistochemistry_receptor_status %in%
                c("Negative", "Indeterminate", "Equivocal", "[Not Evaluated]", "[Not Available]", ""))
             |
               (tnbc$lab_proc_her2_neu_immunohistochemistry_receptor_status=="Negative" &
                  tnbc$her2_immunohistochemistry_level_result %in%
                  c("[Not Available]", "0", "1+", "")),]

HER2status <- unique(tnbc[,c("lab_proc_her2_neu_immunohistochemistry_receptor_status",
                                "her2_immunohistochemistry_level_result",
                                "lab_procedure_her2_neu_in_situ_hybrid_outcome_type")])

write.table(HER2status, paste0(dir,"HER2status.txt"), sep="\t", quote=F, row.names=F)

write.table(tnbc, paste0(dir,"tnbc.txt"), sep="\t", quote=F, row.names=F)


# clinpt <- clinpt[,c("bcr_patient_uuid", "bcr_patient_barcode", "histological_type_other", 
#                     "breast_carcinoma_estrogen_receptor_status", "er_level_cell_percentage_category", 
#                     "breast_carcinoma_immunohistochemistry_er_pos_finding_scale", 
#                     "immunohistochemistry_positive_cell_score", 
#                     "positive_finding_estrogen_receptor_other_measurement_scale_text", "er_detection_method_text",
#                     "breast_carcinoma_progesterone_receptor_status",
#                     "progesterone_receptor_level_cell_percent_category", 
#                     "breast_carcinoma_immunohistochemistry_progesterone_receptor_pos_finding_scale",
#                     "breast_carcinoma_immunohistochemistry_pos_cell_score", 
#                     "pos_finding_progesterone_receptor_other_measurement_scale_text", "pgr_detection_method_text",
#                     "lab_proc_her2_neu_immunohistochemistry_receptor_status", 
#                     "her2_erbb_pos_finding_cell_percent_category", "her2_immunohistochemistry_level_result", 
#                     "pos_finding_her2_erbb2_other_measurement_scale_text", 
#                     "her2_erbb_method_calculation_method_text",
#                     "lab_procedure_her2_neu_in_situ_hybrid_outcome_type", 
#                     "her2_neu_breast_carcinoma_copy_analysis_input_total_number",
#                     "fluorescence_in_situ_hybridization_diagnostic_procedure_chromosome_17_signal_result_range",
#                     "her2_neu_and_centromere_17_copy_number_analysis_input_total_number_count",
#                     "her2_neu_chromosone_17_signal_ratio_value", 
#                     "her2_and_centromere_17_positive_finding_other_measurement_scale_text",
#                     "her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text")]

#   [1] "bcr_patient_uuid"                                                                    
#   [2] "bcr_patient_barcode"                                                                                             
#   [3] "form_completion_date"                                                                                            
#   [4] "tissue_prospective_collection_indicator"                                                                         
#   [5] "tissue_retrospective_collection_indicator"                                                                       
#   [6] "days_to_birth"                                                                                                   
#   [7] "gender"                                                                                                          
#   [8] "menopause_status"                                                                                                
#   [9] "race"                                                                                                            
#   [10] "ethnicity"                                                                                                       
#   [11] "other_dx"                                                                                                        
#   [12] "history_of_neoadjuvant_treatment"                                                                                
#   [13] "person_neoplasm_cancer_status"                                                                                   
#   [14] "vital_status"                                                                                                    
#   [15] "days_to_last_followup"                                                                                           
#   [16] "days_to_death"                                                                                                   
#   [17] "radiation_therapy"                                                                                               
#   [18] "postoperative_rx_tx"                                                                                             
#   [19] "histological_type_other"                                                                                         
#   [20] "year_of_initial_pathologic_diagnosis"                                                                            
#   [21] "age_at_initial_pathologic_diagnosis"                                                                             
#   [22] "initial_pathologic_diagnosis_method"                                                                             
#   [23] "init_pathology_dx_method_other"                                                                                  
#   [24] "breast_carcinoma_surgical_procedure_name"                                                                        
#   [25] "surgical_procedure_purpose_other_text"                                                                           
#   [26] "margin_status"                                                                                                   
#   [27] "breast_carcinoma_primary_surgical_procedure_name"                                                                
#   [28] "breast_neoplasm_other_surgical_procedure_descriptive_text"                                                       
#   [29] "breast_cancer_surgery_margin_status"                                                                             
#   [30] "axillary_lymph_node_stage_method_type"                                                                           
#   [31] "axillary_lymph_node_stage_other_method_descriptive_text"                                                         
#   [32] "cytokeratin_immunohistochemistry_staining_method_micrometastasis_indicator"                                      
#   [33] "primary_lymph_node_presentation_assessment"                                                                      
#   [34] "lymph_node_examined_count"                                                                                       
#   [35] "number_of_lymphnodes_positive_by_he"                                                                             
#   [36] "number_of_lymphnodes_positive_by_ihc"                                                                            
#   [37] "system_version"                                                                                                  
#   [38] "pathologic_T"                                                                                                    
#   [39] "pathologic_N"                                                                                                    
#   [40] "pathologic_M"                                                                                                    
#   [41] "pathologic_stage"                                                                                                
#   [42] "metastatic_site_at_diagnosis"                                                                                    
#   [43] "metastatic_site_at_diagnosis_other"                                                                              
#   [44] "breast_carcinoma_estrogen_receptor_status"                                                                       
#   [45] "er_level_cell_percentage_category"                                                                               
#   [46] "breast_carcinoma_immunohistochemistry_er_pos_finding_scale"                                                      
#   [47] "immunohistochemistry_positive_cell_score"                                                                        
#   [48] "positive_finding_estrogen_receptor_other_measurement_scale_text"                                                 
#   [49] "er_detection_method_text"                                                                                        
#   [50] "breast_carcinoma_progesterone_receptor_status"                                                                   
#   [51] "progesterone_receptor_level_cell_percent_category"                                                               
#   [52] "breast_carcinoma_immunohistochemistry_progesterone_receptor_pos_finding_scale"                                   
#   [53] "breast_carcinoma_immunohistochemistry_pos_cell_score"                                                            
#   [54] "pos_finding_progesterone_receptor_other_measurement_scale_text"                                                  
#   [55] "pgr_detection_method_text"                                                                                       
#   [56] "lab_proc_her2_neu_immunohistochemistry_receptor_status"                                                          
#   [57] "her2_erbb_pos_finding_cell_percent_category"                                                                     
#   [58] "her2_immunohistochemistry_level_result"                                                                          
#   [59] "pos_finding_her2_erbb2_other_measurement_scale_text"                                                             
#   [60] "her2_erbb_method_calculation_method_text"                                                                        
#   [61] "lab_procedure_her2_neu_in_situ_hybrid_outcome_type"                                                              
#   [62] "her2_neu_breast_carcinoma_copy_analysis_input_total_number"                                                      
#   [63] "fluorescence_in_situ_hybridization_diagnostic_procedure_chromosome_17_signal_result_range"                       
#   [64] "her2_neu_and_centromere_17_copy_number_analysis_input_total_number_count"                                        
#   [65] "her2_neu_chromosone_17_signal_ratio_value"                                                                       
#   [66] "her2_and_centromere_17_positive_finding_other_measurement_scale_text"                                            
#   [67] "her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text"                                
#   [68] "new_tumor_event_after_initial_treatment"                                                                         
#   [69] "metastatic_breast_carcinoma_estrogen_receptor_status"                                                            
#   [70] "metastatic_breast_carcinoma_estrogen_receptor_level_cell_percent_category"                                       
#   [71] "metastatic_breast_carcinoma_immunohistochemistry_er_pos_cell_score"                                              
#   [72] "pos_finding_metastatic_breast_carcinoma_estrogen_receptor_other_measuremenet_scale_text"                         
#   [73] "metastatic_breast_carcinoma_estrogen_receptor_detection_method_text"                                             
#   [74] "metastatic_breast_carcinoma_progesterone_receptor_status"                                                        
#   [75] "metastatic_breast_carcinoma_progesterone_receptor_level_cell_percent_category"                                   
#   [76] "metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_score"                                              
#   [77] "metastatic_breast_carcinoma_pos_finding_progesterone_receptor_other_measure_scale_text"                          
#   [78] "metastatic_breast_carcinoma_progesterone_receptor_detection_method_text"                                         
#   [79] "metastatic_breast_carcinoma_lab_proc_her2_neu_immunohistochemistry_receptor_status"                              
#   [80] "metastatic_breast_carcinoma_her2_erbb_pos_finding_cell_percent_category"                                         
#   [81] "metastatic_breast_carcinoma_erbb2_immunohistochemistry_level_result"                                             
#   [82] "metastatic_breast_carcinoma_pos_finding_her2_erbb2_other_measure_scale_text"                                     
#   [83] "metastatic_breast_carcinoma_her2_erbb_method_calculation_method_text"                                            
#   [84] "metastatic_breast_carcinoma_lab_proc_her2_neu_in_situ_hybridization_outcome_type"                                
#   [85] "her2_neu_metastatic_breast_carcinoma_copy_analysis_input_total_number"                                           
#   [86] "metastatic_breast_carcinoma_fluorescence_in_situ_hybridization_diagnostic_proc_centromere_17_signal_result_range"
#   [87] "her2_neu_and_centromere_17_copy_number_metastatic_breast_carcinoma_analysis_input_total_number_count"            
#   [88] "metastatic_breast_carcinoma_her2_neu_chromosone_17_signal_ratio_value"                                           
#   [89] "metastatic_breast_carcinoma_pos_finding_other_scale_measurement_text"                                            
#   [90] "metastatic_breast_carcinoma_her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text"    
#   [91] "anatomic_neoplasm_subdivision"                                                                                   
#   [92] "clinical_M"                                                                                                      
#   [93] "clinical_N"                                                                                                      
#   [94] "clinical_T"                                                                                                      
#   [95] "clinical_stage"                                                                                                  
#   [96] "days_to_initial_pathologic_diagnosis"                                                                            
#   [97] "days_to_patient_progression_free"                                                                                
#   [98] "days_to_tumor_progression"                                                                                       
#   [99] "disease_code"                                                                                                    
#   [100] "extranodal_involvement"                                                                                          
#   [101] "histological_type"                                                                                               
#   [102] "icd_10"                                                                                                          
#   [103] "icd_o_3_histology"                                                                                               
#   [104] "icd_o_3_site"                                                                                                    
#   [105] "informed_consent_verified"                                                                                       
#   [106] "distant_metastasis_present_ind2"                                                                                 
#   [107] "patient_id"                                                                                                      
#   [108] "project_code"                                                                                                    
#   [109] "tumor_tissue_site_other"                                                                                         
#   [110] "stage_other"                                                                                                     
#   [111] "tissue_source_site"                                                                                              
#   [112] "tumor_tissue_site" 