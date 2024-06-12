library(tidyverse)
library(struct)
library(SummarizedExperiment)
library(pmp)
library(structToolbox)
library(limma)
library(lme4)
library(nlme)
library(expss)


########LOADING IN METADATA##################################################

clinical <- read_csv("./ROSMAP_clinical.csv",
                     col_types = list(visit=col_double(),
                                      age_at_visit=col_character(),
                                      cts_mmse30=col_double(),
                                      dcfdx=col_factor(),
                                      apoe_genotype=col_factor(),
                                      race=col_factor(),
                                      braaksc = col_factor(),
                                      ceradsc = col_factor(),
                                      cogdx = col_double(),
                                      ad_reagan = col_factor(),
                                      study = col_factor())
) %>%
  as.data.frame() %>%
  mutate(apoe4 = factor(str_count(apoe_genotype,"4"))) %>%
  mutate(projid = as.character(projid)) %>%
  mutate(apoe4_carrier = ifelse(grepl("4", apoe_genotype),1,0) ) %>%
  mutate(phenotype = factor(paste0(cogdx,apoe4_carrier),
                            levels=c("10","11","20","40","21","41"),
                            labels=c("NCI-","NCI+","MCI-","AD-","MCI+","AD+"))) %>%
  mutate(sub_phenotype = "NCI-/AD-") %>%
  mutate(sub_phenotype = ifelse(phenotype=="NCI+" | phenotype=="MCI-","NCI+/MCI-",sub_phenotype)) %>%
  mutate(sub_phenotype = ifelse(phenotype=="MCI+" | phenotype=="AD+", "MCI+/AD+",sub_phenotype))

for(i in 9:12){
  clinical[,i] <-  ifelse(clinical[,i]=="90+","90",clinical[,i])
  clinical[,i] <- as.numeric(clinical[,i])
}

age_death <- read_csv("./dataset_707_ageofdeath.csv") %>%
  select(projid, age_death) %>%
  rename("age_death" = "age_death_complete") %>%
  mutate(projid = as.character(projid)) %>% drop_na(age_death_complete)

clinical <- clinical %>% 
  left_join(age_death , by ="projid")

long_data <- read_csv("./dataset_707_long_12-10-2020.csv") %>%
  mutate(dcfdx = factor(dcfdx, levels= c(1,2,4), labels = c("NCI","MCI","AD"))) %>%
  select(projid, fu_year,dcfdx,age_at_visit)

cog_data <- read_csv("./cog_patho_data.csv", col_types = list(projid = col_character())) %>%
  select(-study,-educ,-msex,-braaksc,-ceradsc, -pmi,-apoe_genotype, -cogdx)

#### Original Brain Data #############################

raw_brain <- read_csv("brain_lipidomics_experiment_output.csv") %>%
  mutate(projid = ifelse(projid=="5236 study pool","QC",projid)) %>%
  drop_na(projid) %>%
  as.data.frame()
raw_brain <- raw_brain[,!grepl("X",colnames(raw_brain))]
dim(raw_brain)

raw_brain$projid[raw_brain$projid=="QC"] <- paste0("QC",
                                                   1:length(raw_brain$projid[raw_brain$projid=="QC"]))

brain_info <- raw_brain[,1:11] %>%
  left_join(clinical,by="projid") %>%
  mutate(ID = paste0(projid)) %>%
  relocate(ID) %>%
  mutate(QC = ifelse(grepl("QC",projid),"QC","Not QC")) %>%
  mutate(Plate.Bar.Code = factor(Plate.Bar.Code, labels=(1:2))) %>%
  mutate(Well.Position = as.numeric(Plate.Bar.Code)*100 + Well.Position) %>%
  select(-contains("Sample"), -contains("Number"),-contains("mmse30")) %>%
  left_join(cog_data, by = "projid")


brain_matrix <- t(raw_brain[,12:ncol(raw_brain)])
colnames(brain_matrix) <- brain_info$projid

idx_df <- data.frame('met' = rownames(brain_matrix))

brain_se <- SummarizedExperiment(assays=list(counts=brain_matrix),
                                 colData=brain_info,
                                 metadata=list("Description"="brain Metabolite Output"),
                                 rowData = idx_df
)


brain_se <- brain_se[,!is.na(colData(brain_se)$cogdx) & colData(brain_se)$QC !="QC"]
brain_se <- brain_se[,(colData(brain_se)$cogdx %in% c(1,2,4))]
brain_se <- brain_se[,!is.na(colData(brain_se)$age_death_complete)]

######LOADING IN SERUM DATA ############

raw_serum <- read_csv("serum_lipidomics_experiment_output.csv") %>%
  mutate(projid = ifelse(Subaliquot.ID=="4688 study pool","QC",projid)) %>%
  drop_na(projid) %>%
  as.data.frame()
raw_serum <- raw_serum[,!grepl("X",colnames(raw_serum))]


raw_serum$projid[raw_serum$projid=="QC"] <- paste0("QC",
                                                   1:length(raw_serum$projid[raw_serum$projid=="QC"]))

####

serum_info <- raw_serum[,1:13] %>%
  left_join(clinical,by="projid") %>%
  mutate(ID = paste0(projid,"-",Visit.ID)) %>%
  relocate(ID) %>%
  mutate(QC = ifelse(grepl("QC",projid),"QC","Not QC")) %>%
  mutate(Plate.Bar.Code = factor(Plate.Bar.Code, labels=(1:8))) %>%
  mutate(Well.Position = as.numeric(Plate.Bar.Code)*100 + Well.Position) %>%
  select(-contains("Sample"), -contains("Number"),-contains("mmse30")) %>%
  left_join(long_data, by = c("projid", "Visit.ID"= "fu_year")) %>%
  mutate(phenotype.y = factor(ifelse(apoe4_carrier == 1, paste0(dcfdx.y,"+"), paste0(dcfdx.y,"-")),
                              levels=c("NCI-","NCI+","MCI-","AD-","MCI+","AD+"))) %>%
  mutate(sub_phenotype.y = "NCI-/AD-") %>%
  mutate(sub_phenotype.y = ifelse(phenotype.y=="NCI+" | phenotype.y=="MCI-","NCI+/MCI-",sub_phenotype.y)) %>%
  mutate(sub_phenotype.y = ifelse(phenotype.y=="MCI+" | phenotype.y=="AD+", "MCI+/AD+",sub_phenotype.y)) %>%
  mutate(prox_age = age_first_ad_dx - age_at_visit.y) %>%
  group_by(projid) %>%
  mutate(consistent = ifelse(length(unique(phenotype.y))==1, 1, 0)) %>%
  mutate(count = n()) %>%
  mutate(bio_rep = ifelse(count==1, "1",projid)) %>%
  select(-count, -consistent, -contains("dcfdx"), -contains(".ID"),-Species,-Material,-Well.Position) %>%
  ungroup() %>%
  left_join(cog_data, by = "projid")

serum_info$phenotype <- NULL
serum_info$sub_phenotype <- NULL
serum_info <- serum_info %>% rename("phenotype.y"= "phenotype") %>%
  rename("sub_phenotype.y" = "sub_phenotype")
rownames(serum_info) <- serum_info$ID


serum_matrix <- t(raw_serum[,14:ncol(raw_serum)])
colnames(serum_matrix) <- serum_info$ID 
idx_df <- data.frame('met' = rownames(serum_matrix))
serum_se <- SummarizedExperiment(assays=list(counts=serum_matrix),
                                 colData=serum_info,
                                 metadata=list("Description"="serum Metabolite Output"),
                                 rowData = idx_df
)
serum_se <- serum_se[,(!(is.na(serum_se$phenotype)) &
                         colData(serum_se)$QC !="QC")]
serum_se <- serum_se[,!is.na(colData(serum_se)$msex)]
serum_se <- serum_se[,!is.na(colData(serum_se)$age_at_visit.y)]

##########QC and Normalization####################################################################

filters_imputation <- function(se_object, norm_method, imp_method = NULL, missing = 0.17){
  
  print(dim(assay(se_object)))
  high_na_lipids <- rownames(se_object)[apply(assay(se_object),1,function(x) {sum(is.na(x))})/ncol(assay(se_object)) >= missing]
  print(high_na_lipids)
  
  print(sum(is.na(assay(se_object)))/product(dim(assay(se_object))))
  se_filtered <- filter_samples_by_mv(df=se_object, max_perc_mv=0.30)
  print(dim(se_filtered))
  se_filtered <- filter_peaks_by_fraction(df=se_filtered, min_frac=1-missing, 
                                          classes=se_filtered$QC, method="across")
  print(dim(se_filtered))
  print(sum(is.na(assay(se_filtered)))/product(dim(assay(se_filtered))))
  
    assay(se_filtered)[is.na(assay(se_filtered))] <- 0

  if(norm_method == "mole_perc"){
    fd <- assay(se_filtered) %>%  expss::prop_col() %>% t() %>% as.data.frame() %>%
      mutate(projid = se_filtered$projid,
             ID = se_filtered$ID,
             phenotype = se_filtered$phenotype,
             age = log(se_filtered$age_death_complete),
             sex = se_filtered$msex,
             apoe4_carrier = se_filtered$apoe4_carrier, educ = log(se_filtered$educ))  %>%
      pivot_longer(-any_of(c("projid", "ID","phenotype", "age","sex","apoe4_carrier","educ")), names_to = "lipid_name", values_to= "value") %>%
      mutate(idx = rep(rowData(se_filtered)$idx, length(unique(ID))))
    loc <- str_extract(se_object@metadata$Description, "brain|serum")
    fd %>% write_tsv(paste0(loc,"_filtered_mole_perc_unnormalized_raw_data_2024.tsv"))
    assay(se_filtered) <- expss::prop_col(assay(se_filtered))
    assay(se_filtered) <- t(scale(t(assay(se_filtered))))
  }
  
 
  
  print("BEC")
  if(grepl("brain",metadata(se_object)$Description)){
    matrix <- removeBatchEffect(assay(se_filtered), batch=colData(se_filtered)$Plate.Bar.Code
                                , covariates=colData(se_filtered)[,c("msex", "age_death_complete","pmi","educ")]
    )
  }
  else{
    
    matrix <- removeBatchEffect(assay(se_filtered), batch=colData(se_filtered)$Plate.Bar.Code,
                                batch2=colData(se_filtered)$bio_rep,
                                covariates=colData(se_filtered)[,c("msex","age_at_visit.y","educ")] 
    )
  }
  colnames(matrix) <- colnames(assay(se_filtered))
  assay(se_filtered) <- matrix 
  
  print(dim(se_filtered))

  return(se_filtered)
}


imputation_choice <- function(brain_se,serum_se,normalization_method=NULL, imp_method=NULL, missing=.19){
  ########### Brain se
  if(!is.null(brain_se)){
    brain_data <- filters_imputation(brain_se, normalization_method,imp_method, missing)
    
    brain_de <<- DatasetExperiment(data=t(assay(brain_data)),
                                   sample_meta=colData(brain_data),
                                   variable_meta = rowData(brain_data),
                                   name="Brain Metabolites Output",
                                   description=normalization_method)
    print(brain_de)
  }
  ############### Serum se
  if(!is.null(serum_se)){
    
    serum_data <- filters_imputation(serum_se,normalization_method, imp_method)
    
      serum_data <- filters_imputation(serum_se, normalization_method,imp_method)
      
      serum_de <<- DatasetExperiment(data=t(assay(serum_data)),
                                     sample_meta=colData(serum_data),
                                     variable_meta = rowData(serum_data),
                                     name="Serum Metabolites Output",
                                     description=normalization_method)
    
  
    print(serum_de)
  }
}


adj_data <- function(de_object){
  
  if(de_object$description=="mole_perc"){
    
  }
  
  first_lip <- colnames(assay(de_object))[1]
  last_lip <- tail(colnames(assay(de_object)),1)
  
  final_data <- as_tibble(assay(de_object)) %>%
    mutate(projid = rowData(de_object)$projid) %>%
    mutate(ID = rowData(de_object)$ID) %>%
    relocate(projid) %>%
    relocate(ID) %>%
    left_join(as_tibble(rowData(de_object)), by = "ID") %>%
    select(-projid.y) %>%
    rename("projid.x" = "projid") %>%
    mutate(cognitive_status = ifelse(cogdx==1,0,1)) %>%
    pivot_longer(all_of(first_lip):all_of(last_lip),names_to="lipid_name",values_to="value") %>%
    left_join(as_tibble(colData(de_object)), by =c("lipid_name"="met")) %>%
    select(-contains(c("QC","flag","filter","fraction")))
  return(final_data)
  
  
}


###############SAVING OBJECTS####################################################

choice <- "mole_perc"
imputation_choice(brain_se,serum_se,normalization_method = choice)

brain_final_data <- adj_data(brain_de)
serum_final_data <- adj_data(serum_de)

save(brain_de, brain_final_data,
     serum_de,serum_final_data,
     file = "mole_perc_objects.rda")



