
library(patchwork)
library(rstatix)
library(ggpubr)
library(grid)
library(gridExtra)
library(umap)
library(RColorBrewer)
library(circlize)

knn_plots <- function(location, n_clusters=3, kmeans_method, diff_lipids=NULL, title, save=T, paired =T){
  
  de_object <- get(paste0(location,"_de"))
  
  if(!is.null(diff_lipids)){
    de_object <- de_object[,(colData(de_object)$met %in% diff_lipids)]
  }
  assay <- assay(de_object)
  
  set.seed(123)
  km_cluster <- pam(assay,n_clusters)$clustering %>% unlist()
  color_code= scales::hue_pal()(6)[c(1,5,2,3,4,6)]
  color_code <- c("red","orange","#9590FF","black","skyblue","grey")
  colors_6 <- data.frame(groups=1:6, group =c("red","orange","purple","black","skyblue","grey"))

  set.seed(123)
  
  umap_data <- as.data.frame(umap::umap(assay)$layout) %>%
    mutate(projid = rownames(assay)) %>%
    mutate(phenotype = rowData(de_object)$phenotype) %>%
    mutate(phenotype = factor(phenotype, levels=c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))) %>%
    mutate(groups = km_cluster) %>%
    left_join(colors_6, by = "groups") %>%
    mutate(group = factor(group, levels=colors_6$group[1:length(unique(km_cluster))])) 
  
  pheno_types <- "phenotype"
  
  #  for(pheno_col in pheno_types){
  pheno_col <- pheno_types
  umap_plot <- umap_data %>%
    ggplot(aes(x=V1,y=V2,color=group,shape=umap_data[,pheno_col]))+
    geom_point(size=2)+ # Brain has aes(size=2)
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    theme_bw()+
    scale_color_manual(values=color_code[1:max(umap_data$groups)], name = "Cluster\nAssignment")+
    scale_shape_manual(name = "Diagnosis", values =c(0,15,1,16,2,17))+
    guides(size = "none")+
    theme(axis.text=element_text(size=13),
          axis.title=element_text(size=15))

  ##############################
  
  if(!isTRUE(paired)){
    
    heatmap0 <- umap_data %>%
      dplyr::rename("Phenotype"=all_of(pheno_col)) %>%
      group_by(Phenotype,groups) %>%
      mutate(count = n()) %>%
      ungroup(groups) %>%
      mutate(pheno_proportion = count/n()) %>%
      ungroup() %>%
      distinct(groups,Phenotype,.keep_all=T) %>%
      select(Phenotype, group, count, pheno_proportion)
    phenoTitle <- "Unpaired"
  } else{
    heatmap0 <- umap_data %>%
      dplyr::rename("Phenotype"=all_of(pheno_col)) %>%
      mutate(Phenotype = as.character(Phenotype)) %>%
      mutate(Phenotype = ifelse(Phenotype %in% c("NCI-","AD-"),"NCI-/AD-",Phenotype),
             Phenotype = ifelse(Phenotype %in% c("MCI-","NCI+"),"NCI+/MCI-",Phenotype),
             Phenotype = ifelse(Phenotype %in% c("MCI+","AD+"),"MCI+/AD+",Phenotype),
             Phenotype = factor(Phenotype , levels = c("NCI-/AD-","NCI+/MCI-","MCI+/AD+"))) %>%
      group_by(Phenotype,groups) %>%
      mutate(count = n()) %>%
      ungroup(groups) %>%
      mutate(pheno_proportion = count/n()) %>%
      ungroup() %>%
      distinct(groups,Phenotype,.keep_all=T) %>%
      select(Phenotype, group, count, pheno_proportion)
    
    phenoTitle <- "Paired"
  }
  
  
  
  heatmap <- heatmap0 %>%
    tidyr::expand(Phenotype,group) %>%
    left_join(heatmap0, by=c("Phenotype","group")) %>%
    mutate(pheno_proportion = replace_na(pheno_proportion,0)) %>%
    mutate(count = replace_na(count,0))
  
  if(n_clusters >2){
    matrix_form <- heatmap %>% select(Phenotype, group, count) %>%
      pivot_wider(names_from = group, values_from= count) %>%
      select(-Phenotype) %>%
      as.matrix()
    fisher <- function(matrix){
      p.val <- matrix
      for(i in 1:nrow(matrix)){
        for(j in 1:ncol(matrix)){
          red_matrix <- matrix(c(matrix[i,j],
                                 sum(matrix[-i,j]),
                                 sum(matrix[i,-j]),
                                 sum(matrix[-i,-j])),nrow=2)
          p.val[i,j] <- round(fisher.test(red_matrix)$p.value,3)
        }
      }
      return(p.val)
    }
    
    heatmap <- heatmap %>% arrange(group) %>%
      mutate(count = paste0(matrix_form,"\n (", fisher(matrix_form),")"))
    
  }
  
  if(n_clusters == 2){
    unique_pheno <- c("NCI-","NCI+","MCI-","MCI+","AD-","AD+")
    
    pheno_data <- rowData(de_object) %>%
      fastDummies::dummy_cols(pheno_col) %>%
      select(contains(pheno_col))
    
    unique_pheno_p <- c()
    vec <- paste0("phenotype_",c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"))
    for(i in vec){
      print(i)
      sub_data <- pheno_data %>%
        select(contains(i)) %>%
        pull()
      sum <- fisher.test(sub_data,km_cluster)$p.value
      print(sum)
      unique_pheno_p <- c(unique_pheno_p, round(sum,2))
    }
    
    heatmap <- heatmap %>%
      mutate(Phenotype = factor(Phenotype, levels=c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"),
                                labels=paste0(unique_pheno,"\n(",unique_pheno_p,")")))
  }
  fill_var <- ifelse(n_clusters==2, "pheno_proportion","mlog10p")

  heatmap <- heatmap %>%
    mutate(mlog10p = -log(parse_number(str_remove(count,"^[0-9]+\n ")),base=10)) #%>%
  print(head(heatmap))
  heatmap_plot <- heatmap %>%
    mutate(sig = ifelse(mlog10p >= 1,"Yes","No")) %>%
    mutate(Color = ifelse(sig=="Yes",T,F)) %>%
    mutate(Size = ifelse(sig=="Yes", T, F)) %>% arrange(sig) %>%
    ggplot(aes_string(x="Phenotype", y="group", label = "count"))
  
  if(n_clusters==2){
    heatmap_plot <- heatmap_plot + geom_tile(aes_string(fill=fill_var))
  } else{
    heatmap_plot <- heatmap_plot + geom_tile(aes_string(fill=fill_var,color ="Color", size="Size")) 
    
    
  }
  heatmap_plot <- heatmap_plot + scale_colour_manual(values = rev(c("black", "white"))) + 
    scale_size_manual(values = rev(c(1.2, 0)))+
    geom_text() + theme_bw()+
    xlab("Diagnosis Group")+
    ylab("Cluster Assignment")+
    scale_fill_gradient(name = "Group\nDistribution",low = "grey", high = "red") +
    theme(panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_blank(),
          legend.title = element_text(size=11),
          legend.text = element_text(size=10),
          axis.text = element_text(size=ifelse(n_clusters!=2,10,7.5), color = "black"),
          axis.title=element_text(size=11),
          legend.position = "bottom")+
    guides(fill = "none",size= "none",color = "none")
  
  if(isTRUE(save)){
    ggsave(paste0(title,"_knn_heatmap_",phenoTitle,"_",n_clusters,"_",kmeans_method,"_Heatmap.png"),
           heatmap_plot,height=4,width=5, dpi = 2100)
    ggsave(paste0(title,"_knn_",phenoTitle,"_",n_clusters,"_",kmeans_method,"_UMAP.png"),
           umap_plot, width=5,height = 3,dpi=2100)  
  }
}

##########




adjacency_matrix <- function(location, mthd, diff_lipids=NULL){
  
  de_object <- get(paste0(location,"_de"))
  if(!is.null(diff_lipids)){
    de_object <- de_object[,colData(de_object)$met %in% diff_lipids]
  }
  assay <- assay(de_object)
  
  optimal_power <- function(power){
    ADJ1 <- abs(cor(assay,use="p", method=mthd))^power
    k <- as.vector(apply(ADJ1,2,sum, na.rm=T))
    return(ifelse(scaleFreeFitIndex(k)$slope.SFT>0,0,scaleFreeFitIndex(k)$Rsquared.SFT))
  }
  
  pwr <- map_dbl(2:20,optimal_power)
  print(pwr)
  print(pwr[which.max(pwr)])
  pwr <- c(2:10)[which.max(pwr)]
  print(pwr)
  
  ADJ1=abs(cor(assay,use="p", method=mthd))^pwr
  return(ADJ1)
  
}

####### MDS Plot
mds_plot <- function(v1,v2, legend.display =T, mds_data, var_exp){
  dim_one <- parse_number(as.character(v1))
  dim_two <- parse_number(as.character(v2))
  
  mod_colors <- sort(unique(mds_data$module))
  mod_colors1 <- replace(mod_colors, mod_colors == "yellow", "gold")
  mod_colors1 <- replace(mod_colors1, mod_colors1 == "turquoise",
                         "darkturquoise")
  
  mds_graph <- mds_data %>%
    ggplot(aes(x=.data[[v1]],y=.data[[v2]],color=module,
               shape=class,
               size=significance)) +
    geom_text(aes(label=met),size=4)+
    scale_color_manual(values=mod_colors1,
                       labels=mod_colors, name="Module")+
    scale_shape_manual(values=c(16,17,15,10,7), name = "Class")+
    xlab(paste0("Scaling Dimension ", dim_one ," (",round(var_exp[dim_one],3)*100,"%)" ))+
    ylab(paste0("Scaling Dimension ", dim_two ," (",round(var_exp[dim_two],3)*100,"%)" ))+
    theme_bw()+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"),
          legend.text = element_text(size=14),
          legend.key = element_rect(fill = "white"),
          legend.title = element_text(size=16))
  
  if(legend.display == T){
    mds_graph <- mds_graph +
      guides(size = "none",
             color = guide_legend(override.aes = list(size=5)),
             shape = guide_legend(override.aes = list(size=5)))
  }
  else{
    mds_graph <- mds_graph + guides(size="none",color="none",shape="none")
  }
  
  
  return(mds_graph)
}


tukey_gg <- function(target_module, modules_data){
  modules_data1 <- modules_data %>%
    filter(module==target_module) %>%
    mutate(phenotype = factor(phenotype, levels = c("NCI-","NCI+","MCI-",
                                                    "MCI+","AD-","AD+")))
  
  
  floor1 <- quantile(modules_data1$module_value, 1)*1.03
  ceil1 <- quantile(modules_data1$module_value, 1)*1.28
  
  res.aov <- modules_data1 %>% anova_test(module_value ~ phenotype)
  
  if(res.aov$p <= .1){
    pwc <- modules_data1 %>%  tukey_hsd(module_value ~ phenotype)
    p_df <- pairwise.wilcox.test(modules_data1$module_value, modules_data1$phenotype, p.adjust.method = "none")$p.value %>%
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
  }
  gg <- ggboxplot(modules_data1 %>% filter(module_value<=15), x = "phenotype", y = "module_value", color="phenotype") +
    ggtitle(paste0(target_module, " (Anova, "," p = ",res.aov$p,")"))+
    theme_bw()+
    guides(color = "none")+
    theme(axis.title = element_text(size=12, color = "black"),
          axis.text = element_text(size=11,color = "black"),
          plot.title = element_text(size=14),
          axis.title.x= element_blank())  +
    ylab("Eigenlipid Value")
  
  if(res.aov$p <= .1){
    my_comparisons <- data.frame(group1 = c("NCI-","NCI+","MCI-","MCI+","AD-","AD+"),group2= c("NCI-","NCI+","MCI-","MCI+","AD-","AD+")) %>% tidyr::expand(group1,group2)
    my_comparisons <- lapply(1:nrow(my_comparisons), function(i) return(unname(unlist(my_comparisons[i,]))))
    
    gg <- gg +
      stat_compare_means(aes(label = ..p.signif..),
                         comparisons = my_comparisons,# label.y = max_jp
                         size=3, test = "t.test",hide.ns=T)
    stat_pvalue_manual(pwc, label="p.signif", hide.ns = TRUE, size = 6) 
  }
  
  
  return(gg)
  
}

##############

module_analysis <- function(location, title, min_clust_size,diff_lipids = NULL,
                            remove_imputed=F, cluster_method="complete",
                            print_heatmap=T){
  de_object <- get(paste0(location,"_de"))
  adj_matrix <- get(paste0(location,"_adj_matrix"))

    if(!is.null(diff_lipids)){
    de_object <- de_object[,colData(de_object)$met %in% diff_lipids]
  } else{
    diff_lipids <- colnames(de_object)
  }
  assay <- assay(de_object)
  
  lipid_classes <- as_tibble(colData(de_object)) %>%
    select(any_of(c("met","idx"))) %>%
    distinct(.keep_all=T)
  
  
  dissTOM=TOMdist(adj_matrix)
  collectGarbage()
  
  cmd1=cmdscale(as.dist(dissTOM),3, eig=T)
  var_exp <- (cmd1$eig[1:10] / sum(cmd1$eig))
  
  hierTOM = hclust(as.dist(dissTOM),method=cluster_method)
  colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , 
                                                      cutHeight = 0.995,
                                                      deepSplit=4, minClusterSize =min_clust_size,
                                                      pamRespectsDendro = FALSE))
  color_labels <- data.frame(colorDynamicHybridTOM)
  rownames(color_labels) <- colnames(assay)
  
  dendro <- plotDendroAndColors(hierTOM,
                                colors=color_labels,
                                dendroLabels = NULL,marAll = c(1, 8, 3, 1),
                                main = "Lipid dendrogram and module colors, TOM dissimilarity")

  
  colorh1 <- colorDynamicHybridTOM
  save(dissTOM,cmd1, hierTOM, color_labels, var_exp,colorh1,
       file = paste0(title,"_WGCNA_",min_clust_size,".rds"))
  # load(paste0(title,"_WGCNA_",min_clust_size,".rds"))
  
  mds_data <- cmd1$points %>% as.data.frame() %>%
    mutate(class = lipid_classes$idx) %>%
    mutate(met = lipid_classes$met) %>%
    mutate(module = as.character(colorh1))
  mds_data %>% write_tsv(paste0(title,"_mds_data.tsv"))
  
  print(colorh1)

  datME=moduleEigengenes(assay,colorh1)$eigengenes
  print("Boxplots")
  
  module_eigen <- function(color){
    set.seed(123)
    print(color)
    assay <- assay(de_object)[,colorh1 == color]
    mat <- prcomp(assay,scale=T,center=T)
    mat_sum <- summary(mat)
    print("Proportion of Variance Explained -")
    try(print(mat_sum$importance[2,1:2]))
    return(as.data.frame(mat$x)$PC1)
  }
  
  
  datPC <- map_dfc(unique(colorh1),module_eigen)
  colnames(datPC) <- paste0("ME",unique(colorh1))

  modules_data <- datME %>%
    mutate(phenotype = rowData(de_object)$phenotype)%>%
    mutate(sub_phenotype = "NCI-/AD-") %>%
    mutate(sub_phenotype = ifelse(phenotype=="NCI+" | phenotype=="MCI-","NCI+/MCI-",sub_phenotype)) %>%
    mutate(sub_phenotype = ifelse(phenotype=="MCI+" | phenotype=="AD+", "MCI+/AD+",sub_phenotype)) %>%
    pivot_longer(contains("^ME"),names_to="module",values_to="module_value")
  
  
  
  gg_list <- c()
  i <- 1
  for(j in unique(modules_data$module)){
    print(i)
    gg_list[[i]] <- tukey_gg(j, modules_data)
    i <- i+1
  }
  n_mod <- length(unique(modules_data$module))
  module_grid <- grid.arrange(grobs = gg_list, as.table=T, nrow = ifelse(n_mod<=2,1,2)) 
  ggsave(paste0(title,"_module_boxplots_",cluster_method,"_",min_clust_size,"_table_v2.png"),
         module_grid, width=12,height=8, dpi=1200)
  
  #####
  
  m2 <- dissTOM
  rownames(m2) <- colnames(adj_matrix)
  colnames(m2) <- rownames(m2)
  m2 <- m2[lipid_df$Lipid,lipid_df$Lipid]
  m2[upper.tri(m2, diag=T)] <- NA
  m2 <- m2[,ncol(m2):1]
  
  lipid_df <- lipid_df[nrow(lipid_df):1,]
  lipid_df_v2 <- lipid_df %>% select(Module) %>% as.data.frame()
  rownames(lipid_df_v2) <- lipid_df$Lipid
  colors_vec <- unique(lipid_df$Module)
  names(colors_vec) <- colors_vec
  ha1 = columnAnnotation(
    df=lipid_df_v2,col=list(Module=colors_vec), show_legend=T)
  row_split_labels <- lipid_df$Module[nrow(lipid_df):1] ; names(row_split_labels) <- lipid_df$Lipid[nrow(lipid_df):1]
  col_split_labels <- lipid_df$Module ; names(col_split_labels) <- lipid_df$Lipid
  
  for(i in unique(col_split_labels)){
    names(col_split_labels[col_split_labels==i]) <- rev(names(col_split_labels[col_split_labels==i])) 
  }
  col_split_labels  = factor(col_split_labels, levels= unique(col_split_labels))
 
  m2 <- Heatmap(m2, 
          cluster_rows=F,cluster_columns=F,show_row_names=T,
         na_col="white",show_column_names=F,
          bottom_annotation = ha1,

          row_split = row_split_labels, column_split = col_split_labels,
          name = "Dissimilarity\nCorrelation", 
          column_names_gp = grid::gpar(fontsize = 10),

  )
  
  if(isTRUE(print_heatmap)){
    png(paste0(str_to_title(title),"_heatmap_Complex_Correlation_adj.png"),width=11,height=10,units="in",res=1200)
    print(m2)
    dev.off()
     ggsave(paste0(title,"_module_boxplots_",cluster_method,"_",min_clust_size,"_adj.png"), module_grid,
            width=12,height=8, dpi=350)
  }
}


