setwd("~")
setwd("lipid_workflow/bulk_rna_seq")
library(tidyverse)
library(ggrepel)

genes <- read.csv("pathway_genes.csv", header=T, sep=",")
colnames(genes) <- c("Gene","Pathway")
genes$Gene <- str_to_upper(genes$Gene)

matrix <- read.csv("ROSMAP_bulkrnaseq_updated_residuals.csv")
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]
colnames(matrix) <- str_remove(colnames(matrix), "X")
rownames(matrix) <- str_to_upper(rownames(matrix))
matrix <- matrix[rownames(matrix) %in% genes$Gene,]
old_matrix <- matrix

clinical <- read.csv("/mnt/vast/hpc/homes/jm4454/lipid_workflow/ROSMAP_clinical.csv") %>%
  select(projid, study, msex, apoe_genotype, age_first_ad_dx,age_death, cogdx,educ, pmi) %>%
  mutate(age_death = parse_number(age_death),
         age_first_ad_dx = parse_number(age_first_ad_dx)) 

rid_data <- read_csv("/mnt/vast/hpc/homes/jm4454/lipid_workflow/cog_patho_data.csv") %>%
  select(projid, study, cogn_global_random_slope, RIN, plaq_n_sqrt, tangles_sqrt) %>%
  mutate(RIN = as.numeric(RIN)) %>%
  drop_na(RIN)

age_death <- read_csv("/mnt/vast/hpc/homes/jm4454/lipid_workflow/dataset_707_ageofdeath.csv") %>%
  select(projid, study, age_death) %>%
  dplyr::rename("age_death_complete"= "age_death")

#Age of Onset
ad_complete <- read_csv("/mnt/vast/hpc/homes/jm4454/lipid_workflow/dataset_707_long_12-10-2020.csv") %>%
  select(projid,study, dcfdx,age_at_visit) %>%
  mutate(projid = as.numeric(projid)) %>%
  filter(dcfdx %in% c(4)) %>%
  group_by(projid) %>%
  top_n(1, -age_at_visit) %>%
  dplyr::rename("age_first_ad_dx_complete" = "age_at_visit")

clinical <- clinical %>%
  left_join(age_death, by = c("projid")) %>% #uncensored age of death 
  left_join(ad_complete , by = c("projid")) %>% #uncensored age of AD onset
  group_by(projid) %>%
  mutate(
    age_death_complete = max(age_death,age_death_complete,na.rm=T),
    age_death_complete =ifelse(age_death_complete < 0, NA, age_death_complete)) %>%
  mutate(
    age_first_ad_dx_complete = max(age_first_ad_dx,age_first_ad_dx_complete,na.rm=T),
    age_first_ad_dx_complete =ifelse(age_first_ad_dx_complete < 0, NA, age_first_ad_dx_complete)) %>%
  ungroup()


cogn_measures_df <- read.csv("/mnt/vast/hpc/homes/jm4454/lipid_workflow/dataset_707_basic_04-24-2022.csv") %>%
  select(projid, cogng_demog_slope, cogng_path_slope)


matrix <- old_matrix
meta <- clinical %>% left_join(rid_data , by = c("projid")) %>%
  mutate(APOE4 = ifelse(grepl("4",apoe_genotype),1,0),
         msex = factor(msex)) %>%
  filter(projid %in% colnames(matrix)) %>%
  arrange(projid) %>%
  select(-contains("study")) %>%
  drop_na(msex,RIN, age_death_complete)  %>%
  left_join(cogn_measures_df, by ="projid")

matrix <- t(matrix[,colnames(matrix) %in% meta$projid])
matrix <- matrix[order(as.numeric(rownames(matrix))),]

colnames(matrix) %>% as.matrix(ncol=1) %>% write.table("available_genes_in_bulkrna_matrix.csv",quote=F,row.names=F) 

######


model_fn <- function(gene, cogn=T, tan_plaq=F){
  # print(gene)
  meta_new <- meta %>% filter(!is.na(tangles_sqrt), !is.na(plaq_n_sqrt))
  if(!isTRUE(cogn)){
    meta_new <- meta_new # %>% filter(age_first_ad_dx_complete < 90, age_death_complete < 90)
  }
  matrix_sub <- matrix[rownames(matrix) %in% meta_new$projid,]
  meta_new$gene <- matrix_sub[,gene]
  
  dep <- ifelse(isTRUE(cogn),"cogn_global_random_slope", "age_first_ad_dx_complete")
  base <- "~RIN+age_death_complete+msex+gene+educ+pmi"
  tp <- ifelse(tan_plaq, "+tangles_sqrt+plaq_n_sqrt","")
  
  
  modval <- lm(as.formula(paste0(dep, base,tp)),data= meta_new)
  
  # if(isTRUE(cogn)){
  # modval <- lm(cogn_global_random_slope~as.numeric(RIN)+ age_death_complete+msex+gene,data=meta_new)
  # } else{
  # modval <- lm(age_first_ad_dx~as.numeric(RIN)+age_death_complete+msex+gene,data=meta_new)
  #   
  # }
  #  print(dim(summary(modval)$coefficients))
  gene_results <- summary(modval)$coefficients["gene",]
  result <- data.frame(Gene = gene, Estimate = gene_results[1], P = gene_results[4])
  rownames(result) <- gene
  return(result)
}

cogn_results <- map_dfr(colnames(matrix), function(x) model_fn(x, cogn=T, tan_plaq=F))
onset_results <- map_dfr(colnames(matrix), function(x) model_fn(x, cogn=F,tan_plaq=F))
cogn_tp_results <- map_dfr(colnames(matrix), function(x) model_fn(x, cogn=T, tan_plaq=T))
onset_tp_results <- map_dfr(colnames(matrix), function(x) model_fn(x, cogn=F,tan_plaq=T))

volcano_fn <- function(results1,results2, title){
  
  data1 <- results1 %>%
    mutate(P = p.adjust(P, method = "fdr")) %>%
    mutate(Significant = ifelse(Estimate <= 0, "Yes", "No")) %>%
    filter(P <= 0.05)
  data2 <- results2 %>%
    mutate(P = p.adjust(P, method = "fdr")) %>%
    mutate(Significant = ifelse(Estimate <= 0, "Yes", "No")) %>%
    filter(P <= 0.05)
  data <- rbind(data1 %>% mutate(Model = "standard adjustment"),
                data2 %>% mutate(Model = "adjusted for tangles & plaques")) %>%
    mutate(Model = factor(Model, levels = c("standard adjustment",
                                            "adjusted for tangles & plaques")))
  plt <- data %>%
    ggplot(aes(x=-log10(P),y = Estimate, color = Significant))+
    geom_point()+
    geom_label_repel(data = data, aes(label =Gene))+
    theme_bw()+
    scale_color_manual(values = c("Yes"="blue","No"="red"))+
    guides(color = "none")+ ggtitle(title)+
    ylab("Beta estimate")+
    xlab(expression(-log[10](p-value)))+
    theme(
      axis.text = element_text(size=13, color = "black"),
      axis.title = element_text(size=14, color = "black"),
      strip.text=element_text(size=13))+
    facet_wrap(~Model, nrow=1)
  ggsave(paste0("Volcano_Plots_Grid_Significant_",title,"_Educ_PMI.png"), plt, height = 4, width=9, dpi=2200)
  
  return(plt)
}
library(patchwork)
volc <- volcano_fn(cogn_results, cogn_tp_results, "CGRS")


##########


genes <- manual_genes_input

meta_new <- meta %>% filter(!is.na(tangles_sqrt), !is.na(plaq_n_sqrt))
matrix_sub <- matrix[rownames(matrix) %in% meta_new$projid,genes]


genes_df <- meta_new %>%
  left_join(as.data.frame(matrix_sub) %>% rownames_to_column("projid") %>%
              mutate(projid = as.numeric(projid)),
            by = "projid") %>%
  select(any_of(c("projid","cogn_global_random_slope",genes))) %>%
  pivot_longer(genes, names_to="gene", values_to ="value")

plt <- genes_df %>%
  mutate(gene = factor(gene, levels = genes)) %>%
  ggplot(aes(x=value, y = cogn_global_random_slope, color = gene, group=gene))+
  geom_point()+
  theme_bw()+
  geom_smooth(method='lm', formula= y~x, color = "black", linetype="dashed")+
  #  scale_color_brewer(palette="Set1")+
  scale_color_manual(values = brewer.pal(n=9,"Set1")[-6])+
  xlab("Normalized Gene Expression")+ylab("CGRS")+
  facet_wrap(~gene, ncol=2)+guides(color = "none")+
  theme(strip.text = element_text(size=13),
        axis.text = element_text(size=12, color = "black"),
        axis.title = element_text(size=14, color = "black"))
ggsave("Bulk_Gene_Scatter_Plot_Educ_PMI.png", plt, height = 4, width=7, dpi = 2200)

