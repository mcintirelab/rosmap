packages <- c("tidyverse","SummarizedExperiment",
              "struct","RColorBrewer","readr"
              "umap","ggrepel","gridExtra",
              "WGCNA","gridExtra","ggpubr",
              "rstatix","lme4",
              )


load("DLPFC_Serum_objects.rda"))

####DIFFERENTIAL LIPIDS###

differential_lipids <- function(final_data, max_p = 0.05){
  
  pheno_pair <- function(pheno){
    df <- data.frame()
    sub_data <- final_data[final_data$phenotype %in% c("NCI-",pheno),]
    for(i in unique(final_data$lipid_name) ){
      min_data <- sub_data %>%
        filter(lipid_name == i) %>%
        mutate(phenotype = factor(phenotype, levels=c("NCI-",pheno)))
      t <- t.test(value~phenotype,data=min_data)
      lipid_name <- i
      df <- rbind(df, c(lipid_name,-diff(t$estimate),t$p.value))
    }
    df$pheno <- pheno
    colnames(df)[1:3] <- c("lipid_name","estimate","p.value")
    df <- df %>% mutate(direction = ifelse(estimate>0,"NEG","POS")) 
    return(df)
  }
  
  pair_results <- map_dfr(c("NCI+","MCI-","MCI+","AD-","AD+"), pheno_pair)
  diff_lipids <- pair_results %>%
    group_by(pheno) %>%
    mutate(p_adj = p.adjust(p.value, method = "none")) %>%
    ungroup() %>%
    filter(p.value <= max_p) %>%
    distinct(lipid_name, direction) 
  return(diff_lipids)
  
}

p_val <- 0.05
brain_lipids <- differential_lipids(brain_final_data,max_p=p_val)
serum_lipids <- differential_lipids(serum_final_data, max_p =p_val)


differential_lipids <- function(loc){
  if(loc=="serum"){
    return(unique(serum_lipids$lipid_name))
  } else if(loc=="brain"){
    return(unique(brain_lipids$lipid_name))
  }
}
differential_lipids("brain")
differential_lipids("serum")


####MEDIAN CORRELATION HEATMAP####

mean_correlation <- function(loc, lipids_selection=NULL){
  de_object <- get(paste0(loc,"_de"))
  phenotypes <- c("NCI-","NCI+","MCI-","MCI+","AD-","AD+")
  results <- data.frame()
  if(is.null(lipids_selection)){
    lipids_selection <- colData(de_object)$met
    
  }
  for(pheno in phenotypes){

    for(other_pheno in phenotypes){
      pheno_individuals <- de_object[rowData(de_object)$phenotype == pheno,
                                     colnames(de_object) %in% lipids_selection] %>%
        assay() %>% t()
      not_pheno_individuals <- de_object[rowData(de_object)$phenotype == other_pheno,
                                         colnames(de_object) %in% lipids_selection] %>%
        assay() %>% t()
      
      corr_indiv <- c()
      for(i in seq_along(colnames(pheno_individuals))){
        matrix <- cbind(pheno_individuals[,i], not_pheno_individuals)
        corr <- cor(matrix, method = "spearman")[-1,1]
        corr_indiv <- c(corr_indiv, corr)
      } 
      
      vec <- c(pheno, other_pheno, mean(corr_indiv))
      results <- rbind(results, vec)  
    }
    phenotypes <- phenotypes[phenotypes != pheno]  
  }
  colnames(results) <- c("Group A","Group B", "Correlation")
  results <- results %>% mutate(Correlation = as.numeric(Correlation))
  
  return(results)
}


mean_cor_htp <- function(lipids_selection=NULL){
  phenotypes <- rev(c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))
  mat <- rbind(
    mean_correlation("brain", lipids_selection) %>% mutate(Location = "DLPFC"),
    mean_correlation("serum", lipids_selection) %>% mutate(Location = "Serum")
  )

  tops <- mat %>% filter(`Group A` != `Group B`, Location =="Brain") %>% top_n(3, Correlation) %>%
    mutate(transparency = "black")
  mat <- mat %>%
    as_tibble() %>%
    left_join(tops, by= c("Group A","Group B","Location","Correlation")) %>%
    replace_na(list(transparency=alpha("black", 0))) %>%
    filter(`Group A` != `Group B`) %>%
    arrange(transparency)
  print(head(mat))
  
  htp <- mat %>%
    mutate(`Group A` = factor(as.character(`Group A`), levels = c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))) %>%
    mutate(`Group B` = factor(as.character(`Group B`), levels = c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))) %>%
    ggplot(aes(`Group B`,`Group A`))+
    geom_tile(aes(fill=`Correlation`),size=.6)
  
  htp <- htp + geom_text(data=mat %>% mutate(Correlation = round(Correlation,2)),
                         aes(`Group B`, `Group A`, label = `Correlation`,size=12))+
    scale_fill_gradient(low="white", high = "firebrick1", name = "Mean Pairwise \nCorrelation",
    )+
    scale_color_identity()+
    theme_classic()+
    guides(size="none")+
    theme(
      legend.title = element_text(size=8),
      legend.text = element_text(size= 7),
      axis.text = element_text(size=10, color = "black"),
      strip.text = element_text(size=12),
      axis.title = element_blank())+guides(color = "none")+
    facet_grid(~Location)
  
  return(htp)
}
mean_cor_htp()

#####LINEAR ASSOCIATION MODELS#####


associations <- function(location, save_results = T){
  
  final_data <- get(paste0(location,"_final_data"))
  final_data <- final_data %>%
    as_tibble() %>%
    select(any_of(c("cogn_global_random_slope","age_first_ad_dx_complete","lipid_name",
                    "value","age_death_complete","msex",
                    "plaq_n_sqrt","tangles_sqrt", "idx", "apoe4_carrier", "educ","pmi"))) %>%
    drop_na(plaq_n_sqrt,tangles_sqrt)
  
  part_associations <- function(cogn_var){
    dep_term <- ifelse(isTRUE(cogn_var), "cogn_global_random_slope","age_first_ad_dx_complete")
    dep_title <- ifelse(isTRUE(cogn_var), "Cogn. Global Random Slope","Age of Onset")
    
    terms <- list()
    terms[[1]] <- paste0(dep_term," ~ ","value"," + age_death_complete + msex+educ")
    if(location=="brain"){
      terms[[1]] <- paste0(terms[[1]], "+ pmi")
    }
    terms[[2]] <- paste0(terms[[1]]," + plaq_n_sqrt + tangles_sqrt")
    
    all_models <- data.frame()
    for(i in terms){
      models <- final_data %>%
        group_by(lipid_name) %>%
        do(model = lm(as.formula(i), data=.)) 
      actual_models <- map_dfr(1:nrow(models), function(i){return(tidy(models$model[[i]]) %>%
                                                                    filter(term=="value") %>%                                                              mutate(n = paste0("n = ", nobs(models$model[[i]]))))})
      models <- cbind(models,actual_models) %>%
        select(-term, -model) %>%
        mutate(tangles_plaques = ifelse(grepl("tangles",i), "adjusted for tangles & plaques","standard adjustment")) %>%
        drop_na()
      all_models <- rbind(all_models, models)
      
    }
    return(all_models)
  }
  
  all_models <- rbind(part_associations(T) %>% mutate(cogn_var = "Cogn. Global Random Slope"),
                      part_associations(F) %>% mutate(cogn_var= "Age of Onset"))
  all_lipids <- rbind(brain_final_data %>% select(lipid_name,idx),
                      serum_final_data %>% select(lipid_name,idx)) %>%
    distinct(lipid_name, idx)
  
  all_results <- all_models %>%
    left_join(all_lipids, by = "lipid_name") %>%
    relocate(lipid_name,idx) %>%
    mutate(cogn_var = factor(cogn_var, levels = c("Cogn. Global Random Slope","Age of Onset")))%>%
    group_by(cogn_var) %>% mutate(p_adj = p.adjust(p.value,"BH")) %>% ungroup()
  write_tsv(all_results  %>%
              mutate(tangles_plaques =ifelse(tangles_plaques=="adjusted for tangles & plaques",
                                             "adjusted for tangles & plaques",
                                             "standard adjustment")),
            paste0("Association_results_",str_to_upper(location),".tsv"))
  
  plot_fn <- function(data){
    data <- data %>%
      filter(cogn_var !="Age of Onset")
    plot <- data %>%
      ggplot(aes(x=estimate, y= -log10(p_adj),color=idx))+
      geom_point()+
      xlab("Lipid Regr. Coefficient")+
      ylab(expression(-log[10](p-value)))+
      geom_hline(yintercept = -log10(0.05), linetype=3)+
      theme_bw()+
      geom_text_repel(aes(label=ifelse(p_adj <= 0.05,lipid_name,"")),
                      size=2.2,show.legend = F,max.overlaps=20)+
      labs(color = "Lipid Class")+
      scale_color_brewer(palette="Set1")+
      theme(
        strip.text = element_text(size = 10),
        axis.text = element_text(size=8, color = "black"),
        axis.title=element_text(size=11),
        legend.text= element_text(size=11),
        legend.position = "bottom",
        plot.title = element_text())+
      guides(color = guide_legend(ncol=2))+
      facet_wrap(~tangles_plaques, scales = "free_x",ncol=1)
    return(plot)
  }
  label_input <- ifelse(location =="brain","a","b")
  plot1 <- plot_grid(plot_fn(all_results %>%
                               filter(tangles_plaques =="standard adjustment")),
                     labels = label_input ,label_size=16)
  plot2 <- plot_grid(plot_fn(all_results %>%
                               filter(tangles_plaques !="standard adjustment")),
                     labels = label_input ,label_size=16)
  
  if(isTRUE(save_results)){
    ggsave(paste0(location,"_supp_adj.png"),plot1, width = 4, height = 5)
    ggsave(paste0(location,"_main_adj.png"),plot2, width = 4, height = 5)
  }
  
} 

associations("brain")
associations("serum")


####VOLCANO PLOTS#####


volcanos <- function(initial_pheno, save = F){ 

    ttest <- function(pheno){
    df <- data.frame()
    sub_data <- final_data[final_data$phenotype %in% c(initial_pheno,pheno),] %>%
      select(any_of(c("ID","phenotype","value", "lipid_name", "idx")))
    sub_data$phenotype <- factor(sub_data$phenotype, levels=c(initial_pheno,pheno), labels=c(initial_pheno,pheno))
    
    n_pheno <- sub_data %>% distinct(ID, phenotype) %>% filter(phenotype == pheno) %>% nrow()
    
    for(i in unique(final_data$lipid_name) ){
      min_data <- sub_data[sub_data$lipid_name == i,]
       t <- t.test(min_data$value~min_data$phenotype)
      df <- rbind(df, c(i, diff(t$estimate),t$p.value, n_pheno)) 
    }
    colnames(df) <- c("lipid_name","ES", "p", "n")
    df$Phenotype <- pheno
    return(df)
  }
  
  final_data <- get(paste0("brain_final_data"))
  brain_lipids <- final_data %>% select(lipid_name, idx)
  
  remain <- unique(final_data$phenotype)[unique(final_data$phenotype)!=initial_pheno]
  pair_results <- map_dfr(as.character(remain), ttest) %>%
    mutate(Location = "DLPFC")
  
  final_data <- get(paste0("serum_final_data"))
  serum_lipids <- final_data %>% select(lipid_name,idx)
  pair_results <- rbind(pair_results %>% mutate(p_adj = p.adjust(p, method = "BH")),
                        map_dfr(as.character(remain), ttest) %>%
                          mutate(Location = "Serum") %>% mutate(p_adj = p.adjust(p,method="BH"))
  )
  all_lipids <- rbind(brain_lipids, serum_lipids)
  df <- pair_results %>%
    left_join((all_lipids %>% distinct(lipid_name,.keep_all=T)), by = "lipid_name") %>%
    mutate(idx = factor(idx, levels = c("C","PC.aa","PC.ae","SM","lysoPC"))) %>%
    dplyr::rename("Lipid Class"="idx")
  
  volcano <- as_tibble(df) %>%
    mutate(p = as.numeric(p)) %>%
    group_by(Phenotype) %>%
    ungroup() %>%
    mutate(ES = as.numeric(ES)) %>%
    mutate(p = ifelse(-log10(p)>= 5, 10e-5, p)) %>%
    mutate(beyond_sig = ifelse(p ==10e-5, "p<=1e-6", "p>1e-6")) %>%
    mutate(Phenotype = factor(Phenotype, levels= c("NCI+","MCI-","MCI+","AD-","AD+"))) %>%
    ggplot(aes(x=ES,y=-log10(p),
               color=`Lipid Class`
    ))+
    geom_point(size=.5)+
    geom_text_repel(aes(label=ifelse(p <= 0.05,lipid_name,"")),
                    size=2.5,show.legend = F
                    , max.overlaps = 1000
    )+
    geom_text(aes(label = paste0("n = ",n)), x=0.75 , y = 0.2, size = 3, color = "black") +
    xlab("Mean Difference")+
    ylab(expression(-log[10](p-value)))+
    theme_bw()+

    theme(legend.position = "bottom",
          strip.text = element_text(size=14, color = "black"),
          axis.title = element_text(size=14, color = "black"),
          axis.text = element_text(size=10, color = "black"),
          strip.background = element_rect(fill="white"))+
    guides(color = guide_legend(override.aes = list(size = 4)))+
    scale_color_brewer(palette="Set1", name = "Lipid Class")+
    geom_hline(yintercept=-log10(0.05), linetype=3)+
    facet_grid(Phenotype~Location)
  
  if(isTRUE(save)){
    volcano <- plot_grid(volcano, labels = c('a'), label_size = 16)
    ggsave(paste0(initial_pheno,"_volcano_v8.png"), volcano, width=7,height=13,dpi=2200)
  }
  return()
  
}

require(gtable)
volcanos( "NCI-")

########

library(patchwork)
source("03_knn_WGCNA_pipeline.R")

knn_plots("brain",n_clusters=3,
          kmeans_method = "spearman",save=T,paired=F)
knn_plots("serum",n_clusters=3,
          kmeans_method ="spearman",save=T, paired=F)

#################################################
######## Module analysis

source("03_knn_WGCNA_pipeline.R")

k <- "pearson"  
library(WGCNA)
brain_adj_matrix <<- adjacency_matrix("brain", mthd = k,diff_lipids = differential_lipids("brain"))
serum_adj_matrix <<- adjacency_matrix("serum", mthd = k, diff_lipids = differential_lipids("serum"))

met <- "complete"
i <- 6
module_analysis("brain",
                title = paste0("Brain_",i,"_",k),
                min_clust_size = i, 
                remove_imputed=ifelse(choice =="median_norm",T,F),
                cluster_method =met, 
                print_heatmap = T)

j <- 3
module_analysis("serum"
                ,diff_lipids = differential_lipids("serum"),
                title=paste0("Serum_", j,"_",k,"_0515"),
                min_clust_size = j,
                remove_imputed=F,
                cluster_method =met,
                print_heatmap = T) 
rm(brain_adj_matrix)


#######Comparing education residuals in blue module#########
load("Brain_6_pearson_WGCNA_6.rds")


module_eigen <- function(color, de_object, region){
  set.seed(123)
  print(color)
  de_object <- de_object[,colData(de_object)$met %in% differential_lipids(region)]
  assay <- assay(de_object)[,colorh1 == color]
  mat <- prcomp(assay,scale=T,center=T)
  mat_sum <- summary(mat)
  print("Proportion of Variance Explained -")
  try(print(mat_sum$importance[2,1:2]))
  return(as.data.frame(mat$x)$PC1)
}


datPC <- map_dfc(unique(colorh1),function(x) module_eigen(x, brain_de,"brain"))
colnames(datPC) <- paste0("ME",unique(colorh1))

brain_meta <- cbind(rowData(brain_de), datPC) %>% as.data.frame() %>%
  dplyr::select(MEblue, educ,phenotype,age_death_complete,pmi,msex,cogn_global_random_slope) %>%
  drop_na(cogn_global_random_slope)

res <- residuals(lm(MEblue~educ+age_death_complete+pmi+msex,data = brain_meta))
brain_meta$residuals <- res


floor1 <- quantile(brain_meta$residuals, 1)*1.03
ceil1 <- quantile(brain_meta$residuals, 1)*1.28

res.aov <- brain_meta %>% anova_test(residuals ~ phenotype)

pwc <- brain_meta %>%  tukey_hsd(residuals ~ phenotype)
p_df <- pairwise.wilcox.test(brain_meta$residuals, brain_meta$phenotype, p.adjust.method = "none")$p.value %>%
  as.data.frame() %>%
  rownames_to_column("group2") %>%
  pivot_longer(-group2, names_to="group1",values_to="p") %>%
  drop_na(p) %>%
  arrange(group1) %>%
  filter(grepl("NCI-",group1)) %>% mutate(p = p.adjust(p,"BH")) %>%
  mutate(p.signif = ifelse(p>= 0.1,"ns","*")) %>%
  mutate(p.signif = ifelse(p<= 0.05,paste0(p.signif,"*"), p.signif)) 
pwc <- pwc %>% 
  left_join(p_df, by = c("group1","group2")) %>%
  add_xy_position(x="phenotype") %>%
  select(-p.adj,-p.adj.signif) %>%
  filter(p.signif !="ns") %>%
  arrange(y.position) %>%
  mutate(y.position = seq(floor1,ceil1 ,length.out = dplyr::n()))
print(dim(pwc))


plt1 <- brain_meta %>%
  mutate(phenotype = factor(phenotype, levels = c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))) %>%
  ggplot(aes(x = phenotype, y = residuals, color = phenotype))+
  geom_boxplot()+geom_point()+theme_bw()+
  guides(color = "none")+ylab("Adjusted MEblue")+
  theme(axis.text = element_text(color = "black"),axis.title.x = element_blank())+
  stat_pvalue_manual(pwc, label="p.signif", hide.ns = TRUE, size = 6) +
  ggtitle(paste0("MEblue adj. for CGRS (Anova, p =", round(res.aov$p,2),")"))

ggsave("MEblue_adj_for_cgrs_boxplots_v2.png",plt1, height= 7 ,width=7)

#####Correlation between education and blue module######

brain_meta <- cbind(rowData(brain_de), datPC) %>% as.data.frame() %>%
  dplyr::select(MEblue, educ,phenotype,age_death_complete,pmi,msex,cogn_global_random_slope)
res <- residuals(lm(MEblue~educ+age_death_complete+pmi+msex,data = brain_meta))
brain_meta$residuals <- res

educ_cor <- brain_meta %>%
  mutate(phenotype = factor(phenotype, levels = c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))) %>%
  group_by(phenotype) %>% cor_test(MEblue, educ ) %>%ungroup() %>% select(phenotype, cor, p) 
cgrs_cor <- brain_meta %>% drop_na(cogn_global_random_slope) %>% 
  mutate(phenotype = factor(phenotype, levels = c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))) %>%
  group_by(phenotype) %>% cor_test(MEblue, cogn_global_random_slope ) %>%ungroup()%>% select(phenotype, cor, p) 
cgrs_educ_cor <- brain_meta %>% drop_na(cogn_global_random_slope) %>% 
  mutate(phenotype = factor(phenotype, levels = c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))) %>%
  group_by(phenotype) %>% cor_test(residuals, cogn_global_random_slope ) %>%ungroup()%>% select(phenotype, cor, p) 


cors <-  brain_meta %>% group_by(phenotype) %>% summarize(corr = round(cor(MEblue,cogn_global_random_slope,use = "pairwise.complete.obs"),3))
colnames(cors)[2] <- "Pearson\ncorrelation"

plt1 <- brain_meta %>%
  mutate(phenotype = factor(phenotype, levels = c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))) %>%
  ggplot(aes(y = residuals, x = cogn_global_random_slope, color = phenotype, group = phenotype))+
  geom_smooth(method='lm',formula= y~x, linetype="solid",se=F)+
  geom_point()+theme_bw()+
  theme(axis.text = element_text(color = "black"))+
  xlab("CGRS")+ylab("Blue module lipid residuals")+
  ggtitle("Educ-adjusted MEblue vs cgrs")
ggsave("MEblue_adj_educ_vs_cgrs.png",plt1, height= 6 ,width=6)

vals <- brain_meta %>% group_by(phenotype) %>% do(r2 = summary(lm(MEblue~educ, data = .))$r.squared )
vals <- vals %>% unnest(r2)

####MEDIATION ANALYSIS####

library(psych)

brain_ncis <- brain_meta %>% as.data.frame() %>% filter(phenotype %in% c("NCI-","NCI+")) %>%
  mutate(pheno = as.numeric(factor(phenotype, levels= c("NCI-","NCI+"), labels=  c(0,1)))-1) %>%
  dplyr::select(MEblue, pheno, educ, pmi, msex, age_death_complete)
rownames(brain_ncis) <- NULL
png("NCI_NCI_mediation.png", height=6,width=6, units = "in",res=1000)

med_nci_nci = psych::mediate(
  pheno ~ MEblue + (educ) - pmi - msex- age_death_complete, plot=T,
  data = brain_ncis,  n.iter = 500
)
dev.off()

brain_nci_mci <- brain_meta %>% as.data.frame() %>% filter(phenotype %in% c("NCI-","MCI-")) %>%
  mutate(pheno = as.numeric(factor(phenotype, levels= c("NCI-","MCI-"), labels=  c(0,1)))-1) %>%
  dplyr::select(MEblue, pheno, educ, pmi, msex, age_death_complete)
rownames(brain_nci_mci) <- NULL
png("NCI_MCI_mediation.png", height=8,width=6, units = "in",res=1000)
med_nci_mci = psych::mediate(
  pheno ~ MEblue + (educ) - pmi - msex - age_death_complete, plot = T,
  data = brain_nci_mci,  n.iter = 500
)
dev.off()
##########
