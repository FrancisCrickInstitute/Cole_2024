

#### Scripts for Treg paper - Figure 2 & Supplementary figure 2 ####

## Megan Cole (megan.cole@crick.ac.uk)

rm(list=ls())

# Load relevant libraries 
library(ggplot2)
library(dplyr) 
library(RColorBrewer) 
library(ComplexHeatmap) 
library(gplots) 
library(ggpubr) 
library(textshape)
library(tibble)
library(tidyr) 

##########################
#### GLobal variables ####

date = format(Sys.Date(), "%Y%m%d")

path = "/Users/colem/Documents/Projects/Treg_paper/Golden/Data/"

neighb1 = read.csv(paste0(path, "20231002_neighb1.csv"))

out_fig2 = "/Users/colem/Documents/Projects/Treg_paper/Figures/Figure2/"
dir.create(out_fig2)
out_supp2 = "/Users/colem/Documents/Projects/Treg_paper/Figures/Supplementary_Figure2/"
dir.create(out_supp2)

colours1 =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
              "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
              "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
              "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")

treat_col = c("Vehicle" = "#F8766D", "MRTX" = "#00BFC4")
domain_col = c("Normal" = "#F8766D", "Interface" = "#00BA38", "Tumour" = "#619CFF")

## Dataset 2
neighb2 = read.csv(paste0(path, "20231002_neighb2.csv"))

###############################################################################################

###################
#### Functions ####
###################

# Commulative flow diagram
cumulative_flow = function(data,
                           treat,
                           fontsize,
                           out,
                           name) {
  
  # Cumulative flow diagram - MRTX
  p = ggplot(data[which(
    data$scaled_Y <= 0.5 &
      data$scaled_Y >= -0.5 &
      data$domain != "n/a" &
      data$treatment == treat
  ), ],
  aes(x = scaled_X)) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = -1.1,
      ymin = 0,
      ymax = Inf,
      alpha = 0.5,
      fill = "#F3F3F3"
    ) +
    annotate(
      "rect",
      xmin = -1.1,
      xmax = -0.9,
      ymin = 0,
      ymax = Inf,
      alpha = 0.5,
      fill = "#979797"
    ) +
    annotate(
      "rect",
      xmin = -0.9,
      xmax = 0.9,
      ymin = 0,
      ymax = Inf,
      alpha = 0.5,
      fill = "#000000"
    ) +
    annotate(
      "rect",
      xmin = 0.9,
      xmax = 1.1,
      ymin = 0,
      ymax = Inf,
      alpha = 0.5,
      fill = "#979797"
    ) +
    annotate(
      "rect",
      xmin = 1.1,
      xmax = Inf,
      ymin = 0,
      ymax = Inf,
      alpha = 0.5,
      fill = "#F3F3F3"
    ) +
    geom_area(aes(y = ..count.., fill = factor(agglom18_average)), stat =
                'bin') +
    ylim(0, 4500) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 14)
    ) +
    guides(fill = guide_legend(title = "Community")) +
    xlab("Cross section through tissue") +
    ylab("Cell count") +
    scale_fill_manual(values = random_cols)
  p
  ggsave(
    plot = p,
    device = "png",
    width = 7,
    height = 5.5,
    bg = 'white',
    dpi = 300,
    path = out,
    filename = paste0(date, name)
  )
}

# Stacked bar function 
stacked_bar = function(data,
                       X,
                       fill,
                       select,
                       path,
                       ROI_name,
                       position = "fill",
                       colour,
                       width = 10,
                       height = 10,
                       split = 'single') {
  p = ggplot(data, aes(
    x = reorder(get(X), get(fill) == select, FUN = mean),
    fill = as.factor(get(fill))
  )) +
    theme_classic() +
    coord_flip() +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text = element_text(size = 16)
    ) +
    theme(legend.text = element_text(size = 16)) +
    xlab("Community") +
    labs(fill = "")
  if (position == "fill") {
    p = p + geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution")
  } else if (position == "count") {
    p = p + geom_bar(stat = "count") + ylab("Cell count")
  }
  if (split == "multiple"){
    p = p + facet_wrap(treatment~.)
    p = p + theme(strip.text = element_text(size = 16))
    
  }
  p
  ggsave(
    plot = p,
    device = "png",
    width = width,
    height = height,
    dpi = 300,
    path = path,
    filename = paste(ROI_name, ".png", sep = "")
  )
} 


# Community proportions heatmap
community_prop_heatmap = function(data,
                                  ROI,
                                  name,
                                  out_path){
  
  #Generate table for heatmap that contains community size and proportion in treatment groups per ROI
  cd_prop = data.frame()
  for (r in unique(data[,ROI])){
    print(r)
    cd = data[which(data[,ROI] == r),]
    cd_prop[r, "treatment"] = as.character(unique(cd$treatment))
    cd_prop[r, ROI] = r
    total_cc = dim(cd)[1]
    for (cl in as.character(sort(unique(data$agglom18_average)))){
      print(cl)
      cd_prop[r, cl] = dim(cd[which(cd$agglom18_average == cl),])[1]/total_cc
    }
  }
  
  names(cd_prop)
  # Heatmap of cell type proportions
  heatmap_prop = cd_prop[,]
  row.names(heatmap_prop) = cd_prop[,ROI]
  head(heatmap_prop)
  heatmap_prop2 = heatmap_prop[,c(3:ncol(heatmap_prop))]
  heatmap_prop2 = as.matrix(heatmap_prop2)
  col1= colorRampPalette(brewer.pal(8, "Blues"))(25)
  col3 = list(Treatment = c("Vehicle" = "#F8766D", "MRTX" = "#00BFC4"))
  dd = hclust(dist(heatmap_prop2), method = "average")
  ha = HeatmapAnnotation(Treatment = heatmap_prop$treatment,
                         which = "row",
                         col = c(col3))
  

  pdf(file=paste(out_path, date, name, sep = ""), width=9, height=5)
  ht = Heatmap(t(scale(t(heatmap_prop2))), name = "heatmap", col = col1, cluster_rows = dd, row_dend_side = "left", right_annotation = ha)
  draw(ht)
  dev.off()
}

# Create list of cell type pairs 
marker_list = function(celltypes){
  
  # Create df of each unique cell type 
  ct_pairs = data.frame('name1' = select)
  # Repeat each row in df by number of rows 
  ct_pairs = ct_pairs %>% slice(rep(1:n(), each = nrow(ct_pairs)))
  # Create another df on unique cell types and repeat block by number of rows
  pairs2 = data.frame('name2' = select)
  pairs2 = bind_rows(replicate(nrow(pairs2), pairs2, simplify = FALSE))
  # Merge columns into 1 df 
  ct_pairs = cbind(ct_pairs, pairs2)
  
  return(ct_pairs)
}

# Correlation of cell type pairs per community
community_correlations = function(data,
                                  name1,
                                  name2,
                                  df){
  
  # Create data frame with calculated proportions
  # of each cell type in each community
  prop = data %>%
    dplyr::select('cellType', 'agglom18_average') %>%
    dplyr::group_by(agglom18_average) %>%
    dplyr::mutate(community_count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cellType %in% c(name1, name2)) %>%
    dplyr::group_by(agglom18_average) %>%
    dplyr::summarise(community_count = unique(community_count),
                     count1 = sum(cellType == name1),
                     count2 = sum(cellType == name2),
                     prop1 = (count1/community_count)*100,
                     prop2 = (count2/community_count)*100
    )
  
  df$name1[row] = name1
  df$name2[row] = name2
  df$corr[row] = cor(prop$prop1, prop$prop2, method = "pearson")
  df$p_val[row] = cor.test(prop$prop1, prop$prop2, method = "pearson")$p.value
  
  return(df)
}

# Correlation of cell type pairs per ROI 
ct_correlations = function(data,
                           name1,
                           name2,
                           df){
  
  prop = data %>% 
    dplyr::select('cellType', 'ROI_name') %>% 
    dplyr::group_by(ROI_name) %>% 
    dplyr::mutate(roi_count = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(cellType %in% c(name1, name2)) %>% 
    dplyr::group_by(ROI_name) %>% 
    dplyr::summarise(roi_count = unique(roi_count),
                     count1 = sum(cellType == name1),
                     count2 = sum(cellType == name2),
                     prop1 = (count1/roi_count)*100,
                     prop2 = (count2/roi_count)*100
    )
  
  df$name1[row] = name1
  df$name2[row] = name2
  df$corr[row] = cor(prop$prop1, prop$prop2, method = "pearson")
  df$p_val[row] = cor.test(prop$prop1, prop$prop2, method = "pearson")$p.value
  
  return(df)
}

###############################################################################################

###################
#### Figure 2a ####

#### Cumulative flow diagram of top 5 communities for dataset 1 ####
## Likely to go in figure 3 or Supp figure 3 

random_cols = c("1" = "#9c4767",
                "2" = "#5bb84d",
                "3" = "#a35cca",
                "4" = "#b6b236",
                "5" = "#586ccd",
                "6" = "#dc9541",
                "7" = "#4ca2d5",
                "8" = "#cf4e32",
                "9" = "#54be9d",
                "10" = "#c947a0",
                "11" = "#3c8149",
                "12" = "#d84169",
                "13" = "#99ab5b",
                "14" = "#9080c3",
                "15" = "#737129",
                "16" = "#da84bb",
                "17" = "#a16434",
                "18" = "#dc8074")

for (r in unique(neighb1$ROI_name)){
  cd = neighb1[which(neighb1$ROI_name == r),]
  tum_x_min = min(cd[which(cd$domain == "Tumour"), "Location_Center_X"]) 
  tum_y_min = min(cd[which(cd$domain == "Tumour"), "Location_Center_Y"]) 
  tum_x_max = max(cd[which(cd$domain == "Tumour"), "Location_Center_X"]) 
  tum_y_max = max(cd[which(cd$domain == "Tumour"), "Location_Center_Y"]) 
  tum_x_centre = (tum_x_min + tum_x_max)/2
  tum_y_centre = (tum_y_min + tum_y_max)/2
  print(r)
  # set centre to zero, tum_x_min to -1 and tum_x_max to +1
  neighb1[which(neighb1$ROI_name == r) , "scaled_X"] = (cd$Location_Center_X - tum_x_centre)/(tum_x_max- tum_x_centre)
  neighb1[which(neighb1$ROI_name == r) , "scaled_Y"] = (cd$Location_Center_Y - tum_y_centre)/(tum_y_max- tum_y_centre)
}


#### Figure 2a-b & Supp Figure 2a-d ####

# Fig 2a - Vehicle 
cumulative_flow(neighb1, "Vehicle", 5, out_fig2,  "_dataset1_cumulative_flow_18communities_Vehicle_onMRTXscale.png")
# Fig 2b - MRTX
cumulative_flow(neighb1, "MRTX", 4.5, out_fig2, "_dataset1_cumulative_flow_18communities_MRTX_4500threshold.png")

# Supp Fig 2a
com18 = neighb1 %>% filter(agglom18_average == 18)
cumulative_flow(com18, "Vehicle", 4.5, out_supp2, "_dataset1_cumulative_flow_com18_Vehicle_onMRTX.png")
cumulative_flow(com18, "MRTX", 4.5, out_supp2, "_dataset1_cumulative_flow_com18_MRTX_4500threshold.png")

# Supp Fig 2b
com10 = neighb1 %>% filter(agglom18_average == 10)
cumulative_flow(com10, "Vehicle", 4.5, out_supp2, "_dataset1_cumulative_flow_com10_Vehicle_onMRTX.png")
cumulative_flow(com10, "MRTX", 4.5, out_supp2, "_dataset1_cumulative_flow_com10_MRTX_4500threshold.png")

# Supp Fig 2c
com2 = neighb1 %>% filter(agglom18_average == 2)
cumulative_flow(com2, "Vehicle", 4.5, out_supp2, "_dataset1_cumulative_flow_com2_Vehicle_onMRTX.png")
cumulative_flow(com2, "MRTX", 4.5, out_supp2, "_dataset1_cumulative_flow_com2_MRTX_4500threshold.png")

# Supp Fig 2d
com3 = neighb1 %>% filter(agglom18_average == 3)
cumulative_flow(com3, "Vehicle", 4.5, out_supp2, "_dataset1_cumulative_flow_com3_Vehicle_onMRTX.png")
cumulative_flow(com3, "MRTX", 4.5, out_supp2, "_dataset1_cumulative_flow_com3_MRTX_4500threshold.png")


################
#### Fig 2c ####

# Stacked bar
stacked_bar(neighb1, "agglom18_average", "treatment", "Vehicle", out_fig2, paste0(date, "_dataset1_stacked_bar_18communities_proportions_treatment"),
            colour = treat_col, width = 6.5, height = 6.5)

stacked_bar(neighb2, "agglom18_average", "treatment", "Vehicle", out_fig2, paste0(date, "_dataset2_stacked_bar_18communities_proportions_treatment"),
            colour = treat_col, width = 6.5, height = 6.5)

#####################
#### Supp Fig 2e ####

stacked_bar(neighb1[which(neighb1$domain != "n/a"),], "agglom18_average", "domain", "Normal", out_supp2, paste0(date, "_dataset1_stacked_bar_18communities_proportions_domain"),
            position = "count", colour = domain_col, width = 8, height = 6.5, split = "multiple")


#################################
#### Figure 2d / Supp Fig 2f ####

#### Dataset 1 - cell type proportion per community or per ROI ####

# Shorten ROI names 
neighb1 = neighb1 %>% 
  mutate(ROI_short = case_when(ROI_name == "20190913_BRAC3529.2d_ROI1_MRTX"~"BRAC3529.2d_ROI1",
                               ROI_name == "20200130_BRAC4002.3c_ROI2_MRTX_crop1"~"BRAC4002.3c_ROI2_t1",
                               ROI_name == "20200130_BRAC4002.3c_ROI2_MRTX_crop2"~"BRAC4002.3c_ROI2_t2",
                               ROI_name == "20200130_BRAC4002.3c_ROI3_MRTX"~"BRAC4002.3c_ROI3",
                               ROI_name == "20190917_BRAC3495.3f_ROI1_Vehicle_crop1"~"BRAC3495.3f_ROI1_t1",
                               ROI_name == "20190917_BRAC3495.3f_ROI1_Vehicle_crop2"~"BRAC3495.3f_ROI1_t2",
                               ROI_name == "20190927_BRAC3529.2b_ROI1_MRTX_crop2"~"BRAC3529.2b_ROI1",
                               ROI_name == "20191119_BRAC3326.4e_ROI1_Vehicle_crop1"~"BRAC3326.4e_ROI1",
                               ROI_name == "20191121_BRAC3438.6f_ROI1_Vehicle"~"BRAC3438.6f_ROI1",
                               ROI_name == "20191121_BRAC3438.6f_ROI2_Vehicle"~"BRAC3438.6f_ROI2",
                               ROI_name == "20191121_BRAC3438.6f_ROI3_Vehicle"~"BRAC3438.6f_ROI3",
                               ROI_name == "20200130_BRAC4002.3c_ROI1_MRTX"~"BRAC4002.3c_ROI1",
                               TRUE ~ ROI_name))


# Dataset 1 - Fig 2d
community_prop_heatmap(neighb1, "ROI_short", "_neighb1_heatmap_community_proportions_per_ROI.png", out_fig2)
# Dataset 2 - Supp Fig 2f
community_prop_heatmap(neighb2, "ROI_name", "_neighb2_heatmap_community_proportions_per_ROI.pdf", out_supp2)


#################################
#### Figure 2e / Supp Fig 2i ####

neighb2 = neighb2 %>% 
  mutate(ROI_short = case_when(ROI_name == "BRAC3326.4e_ROI1_Vehicle"~"BRAC3326.4e_Vehicle",
                               ROI_name == "BRAC3438.6f_ROI1_Vehicle"~"BRAC3438.6f_t1_Vehicle",
                               ROI_name == "BRAC3438.6f_ROI3_Vehicle"~"BRAC3438.6f_t2_Vehicle",
                               ROI_name == "BRAC3495.3f_ROI1_Vehicle"~"BRAC3495.3f_Vehicle",
                               ROI_name == "BRAC3529.2a_ROI1_MRTX"~"BRAC3529.2a_MRTX",
                               ROI_name == "BRAC3529.2b_ROI1_MRTX"~"BRAC3529.2b_MRTX",
                               ROI_name == "BRAC3529.2d_ROI3_MRTX"~"BRAC3529.2d_MRTX",
                               ROI_name == "BRAC3708.2d_ROI1_Vehicle"~"BRAC3708.2d_Vehicle",
                               ROI_name == "BRAC4002.3c_ROI1_MRTX"~"BRAC4002.3c_t1_MRTX",
                               ROI_name == "BRAC4002.3c_ROI2_MRTX"~"BRAC4002.3c_t2_MRTX",
                               ROI_name == "BRAC4002.3c_ROI3_MRTX"~"BRAC4002.3c_t3_MRTX",
                               TRUE ~ ROI_name))


select = c("B cells", "Dendritic cells", "Dendritic cells CD103", "Fibroblasts", "Macrophages type 1",
           "Macrophages type 2", "T cells CD4", "T cells CD8", "T reg cells", "Tumour")

# Create a list of cell type pairs 
pairs = marker_list(select)

#### Figure 2e ####

# Create new df for saving cor value and p-value
corr_df = data.frame(matrix(nrow = nrow(pairs), ncol = 4))
names(corr_df) = c("name1", "name2", "corr", "p_val")

## Correlation of cell types per community ##
for(row in 1:nrow(pairs)){
  print(paste0(pairs$name1[row], ", ", pairs$name2[row]))
  corr_df = community_correlations(neighb1, pairs$name1[row], pairs$name2[row], corr_df)
}

# Cluster cell types based on similar correlations 
corr_df_wide = corr_df %>% 
  select(name1, name2, corr) %>% 
  pivot_wider(names_from = name2, values_from = corr, values_fill = list(name1 = 0)) %>% 
  column_to_rownames(var = "name1")

corr_df_wide = as.matrix(corr_df_wide)


clust = hclust(dist(corr_df_wide), method = "average")
# Only plot the heatmap to see the order of the cell types following clustering 
heatmap.2(corr_df_wide, cluster_rows = clust,trace = "none", margins = c(10,10))
#dev.off()

# Order based on results from clust (plotted as heatmap for visualisation)
order = c("Dendritic cells", "Dendritic cells CD103", "T cells CD4", "T reg cells", "T cells CD8", 
          "Macrophages type 1", "B cells", "Fibroblasts", "Macrophages type 2", "Tumour")


# Add significance column to corr_df
corr_df$sig = ""
corr_df <- within(corr_df, sig[p_val < 0.05] <- '*')
corr_df <- within(corr_df, sig[p_val < 0.01] <- "**")
corr_df <- within(corr_df, sig[p_val > 0 & p_val < 0.001] <- "***")
corr_df <- within(corr_df, sig[corr == 0 & p_val == 0] <- "" )
# Replace NA values with zeros
corr_df[is.na(corr_df)] <- 0

# Create new DF for correlations with self 
corr_with_self = filter(corr_df, name1 == name2)
corr_df = filter(corr_df, !(name1 == name2))

# Remove repeated correlations
corr_df = corr_df[-c(20,56,57,74,75,79,65,66,70,71,38,39,42:44,1,2,4,6:8,
                      28:31,33:35,46:53,82:90),]
# Plot 
p = ggplot(corr_df,
           aes(x = factor(name1), y = factor(name2), fill = corr)) + 
  geom_tile()+
  theme_classic()+
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1)) + 
  geom_tile(data = corr_with_self, fill = "red") + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(label=sig), size = 5) + 
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0)) +
  scale_x_discrete(limits = order, NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0), limits=rev(order))
p
ggsave(plot = p, device = 'png', height = 5, width = 6.5, path = out_fig2,
       filename = paste0(date, "_dataset1_cellType_correlation_matrix_perCommunity_half_reordered.png"))

################################################

#### Supp Fig 2i ####

## Correlation of cell types per ROI ##

corr_df = data.frame(matrix(nrow = nrow(pairs), ncol = 4))
names(corr_df) = c("name1", "name2", "corr", "p_val")

# Get correlation values for each pair, per community 
for(row in 1:nrow(pairs)){
  print(paste0(pairs$name1[row], ", ", pairs$name2[row]))
  corr_df = ct_correlations(neighb1, pairs$name1[row], pairs$name2[row], corr_df)
}

# Cluster cell types based on similar correlations 
corr_df_wide = corr_df %>% 
  select(name1, name2, corr) %>% 
  pivot_wider(names_from = name2, values_from = corr, values_fill = list(name1 = 0)) %>% 
  column_to_rownames(var = "name1")

corr_df_wide = as.matrix(corr_df_wide)

# # Only plot the heatmap to see the order of the cell types following clustering 
clust = hclust(dist(corr_df_wide), method = "average")
heatmap.2(corr_df_wide, cluster_rows = clust,trace = "none", margins = c(10,10))
dev.off()

# Add significance column to corr_df
corr_df$sig = ""
corr_df <- within(corr_df, sig[p_val < 0.05] <- '*')
corr_df <- within(corr_df, sig[p_val < 0.01] <- "**")
corr_df <- within(corr_df, sig[p_val > 0 & p_val < 0.001] <- "***")
corr_df <- within(corr_df, sig[corr == 0 & p_val == 0] <- "" )
# Replace NA values with zeros
corr_df[is.na(corr_df)] <- 0

# # Create new DF for correlations with self 
corr_with_self = filter(corr_df, name1 == name2)
corr_df = filter(corr_df, !(name1 == name2))

# Based on order of community clustering results
corr_df = corr_df[-c(20,56,57,74,75,79,65,66,70,71,38,39,42:44,1,2,4,6:8,
                     28:31,33:35,46:53,82:90),]

# Plot 
p = ggplot(corr_df,
           aes(x = factor(name1), y = factor(name2), fill = corr)) + 
  geom_tile()+
  theme_classic()+
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1)) + 
  geom_tile(data = corr_with_self, fill = "red") + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(label=sig), size = 5) + 
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0)) +
  scale_x_discrete(limits = order, NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0), limits=rev(order))
p
ggsave(plot = p, device = 'png', height = 5, width = 6.5, path = out_supp2,
       filename = paste(date, "_dataset1_cellType_correlation_matrix_perROI_half_reordered_based_on_community_clustering.png"))


#####################
#### Supp Fig 2g ####

# Showing an example of cell type correlation across the 18 communities 
pearson_corr_community = function(data, name1, name2){
  
  prop = data %>%
    dplyr::select(cellType, agglom18_average) %>%
    dplyr::group_by(agglom18_average) %>%
    dplyr::mutate(community_count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cellType %in% c(name1, name2)) %>%
    group_by(agglom18_average) %>%
    dplyr::summarise(community_count = unique(community_count),
                     count1 = sum(cellType == name1),
                     count2 = sum(cellType == name2),
                     prop1 = (count1/community_count)*100,
                     prop2 = (count2/community_count)*100) 
  
  p = ggplot(prop, aes(x = prop1, y = prop2, label = agglom18_average)) + 
    geom_point(size = 2) + 
    theme_minimal() + 
    geom_smooth(method = "lm") + 
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 20)) +
    stat_cor() + 
    xlab(name1) + 
    ylab(name2)
  p
  ggsave(plot = p, device = 'png', height = 5, width = 5, path = out_supp2, bg = "white",
         filename = paste(date, "_dataset1_cellType_correlation_per_community_", name1, "_", name2, ".png"))
}

# Strong vs no correlation per community 
pearson_corr_community(neighb1, "T cells CD8", "Dendritic cells CD103")


#####################
#### Supp Fig 2h ####

# Showing an example of cell type correlation across the 12 ROIs
pearson_corr_roi = function(data, name1, name2){
  
  prop = data %>%
    dplyr::select(cellType, ROI_name) %>%
    dplyr::group_by(ROI_name) %>%
    dplyr::mutate(community_count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cellType %in% c(name1, name2)) %>%
    group_by(ROI_name) %>%
    dplyr::summarise(community_count = unique(community_count),
                     count1 = sum(cellType == name1),
                     count2 = sum(cellType == name2),
                     prop1 = (count1/community_count)*100,
                     prop2 = (count2/community_count)*100) 
  
  p = ggplot(prop, aes(x = prop1, y = prop2, label = ROI_name)) + 
    geom_point(size = 2) + 
    theme_minimal() + 
    geom_smooth(method = "lm") + 
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 20)) +
    stat_cor() + 
    xlab(name1) + 
    ylab(name2)
  p
  ggsave(plot = p, device = 'png', height = 5, width = 5, path = out_supp2, bg = "white",
         filename = paste(date, "_dataset1_cellType_correlation_per_ROI_", name1, "_", name2, ".png"))
}

pearson_corr_roi(neighb1, "T cells CD8", "Dendritic cells CD103")


#########################################################################################################




