
#### Scripts for Treg paper - Figure 4 & Supplementary figure 4 ####

## Megan Cole (megan.cole@crick.ac.uk)

rm(list=ls())

library(ggplot2)
library(gplots)
library(dplyr)
library(ggridges) 
library(scales) 

date = format(Sys.Date(), "%Y%m%d")

path = ""

neighb1 = read.csv(paste0(path, "neighb1.csv"))
neighb2 = read.csv(paste0(path, "neighb2.csv"))

out_fig4 = paste0(path, "Figure4/")
dir.create(out_fig4)
out_supp4 = paste0(path, "Supplementary_Figure4/")
dir.create(out_supp4)

top5 = c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2", "Other")
top5_col = c("T/NA" = "#F8766D", "T/M1" = "#FCBF49", "T/DC" = "#00BF7D", "T/M2_1" = "#00B0F6", "T/M2_2" = "#E76BF3", "Other" = "#E1E1E1")
treat_col = c("Vehicle" = "#F8766D", "MRTX" = "#00BFC4")

#########################################################################################################

###################
#### Functions ####

# Plotting expression of markers across communities & treatment groups - box plot 
expression_boxplot = function(data, select_ct, X, marker, dir, name, strip_text){
  
  p = ggplot(transform(data[which(data$top5 %in% c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2") & data$cellType %in% select_ct),],
                       treatment = factor(treatment, levels = c("Vehicle", 'MRTX'))),
             aes(x = treatment, y = get(X), fill = top5))  + 
    geom_boxplot() + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          strip.text = element_text(size=strip_text, face = "bold"),
          legend.text = element_text(size = 16),
          legend.title = element_blank(),
          legend.position = "right") +
    ylab(paste0(marker, " expression")) +
    scale_y_continuous(trans='log2', labels = label_number(accuracy = 0.0001)) +
    scale_fill_manual(values= top5_col[1:5]) +
    facet_wrap(cellType~.)
  p
  ggsave(plot = p, device = "png", width=6.5, height=6, dpi=300, path = dir,
         filename = paste0(date, name))
}


# Plotting expression of markers across communities & treatment groups - density plot 
expression_density = function(data, select_ct, X, marker, xlow, xhigh, dir, name, w, h){

  p = ggplot(transform(data[which(data$top5 %in% c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2") & data$cellType %in% select_ct
                                     & data[,X] != '-Inf'),],
                       treatment = factor(treatment, levels = c("Vehicle", 'MRTX'))),
             aes(x = get(X), y = top5, fill = treatment))  + 
    geom_density_ridges(scale = 0.9, alpha = 0.8, quantile_lines = T, quantile_fun = median) + 
    scale_y_discrete(expand = c(0.05, 0)) + 
    xlim(xlow,xhigh) +  
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          strip.text = element_text(size=20, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          legend.position = "right") +
    ylab("Community") +
    scale_fill_manual(values= treat_col) 
  if(length(select_ct) == 1){
    p = p +  xlab(paste0(marker, " expression on ", select_ct)) 
  } else (
    p = p + xlab(paste0(marker, " expression")) 
  )
  p  
  ggsave(plot = p, device = "png", dpi=300, width=w, height=h, path = dir,
         filename = paste0(date, name))
}


# Plotting expression of markers to better understand tumour phenotype
tum_expr = function(data, Y, com, cols, marker, dir, name){

  p = ggplot(data[which(data$cellType == 'Tumour' & data$top5 %in% com & data$treatment == "MRTX"),],
             aes(x = factor(top5), y = get(Y), colour = factor(top5))) + 
    geom_boxplot(size = 2) + 
    theme_classic() + 
    scale_color_manual(values = cols) +
    scale_y_continuous(trans='log2', labels = label_number(accuracy = 0.0001)) + 
  xlab("") + 
    ylab(paste0("Tumour cell ", marker, " expression")) + 
    labs(color = "") + 
    theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 22),
          axis.text.x = element_blank(),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 22))
  p
  ggsave(plot = p, device = "png", width=5, height=5, dpi=300, path = dir,
         filename = paste0(date, name))
}

#########################################################################################################

#####################
#### Supp Fig 4a ####

# Supp Fig 4a 
expression_boxplot(neighb1, c("Macrophages type 1", "Macrophages type 2"), "MI_PDL1", "PD-L1", out_supp4, "_dataset1_top5_PDL1_expression_VehicleMRTX_Macrophages.png", 14)
expression_boxplot(neighb2, c("Macrophages type 1", "Macrophages type 2"), "MI_PDL1", "PD-L1", out_supp4, "_dataset2_top5_PDL1_expression_VehicleMRTX_Macrophages.png", 13.5)


#####################
#### Supp Fig 4b ####

expression_boxplot(neighb1, c("Macrophages type 1", "Macrophages type 2"), "MI_MHCcII", "MHCII", out_supp4, "_dataset1_top5_MHCII_expression_VehicleMRTX_Macrophages.png", 13)
expression_boxplot(neighb2, c("Macrophages type 1", "Macrophages type 2"), "MI_MHCcII", "MHCII", out_supp4, "_dataset2_top5_MHCII_expression_VehicleMRTX_Macrophages.png", 13)


################
#### Fig 4a ####

expression_boxplot(neighb1, c("Dendritic cells", "Dendritic cells CD103"), "MI_PDL1", "PD-L1", out_fig4, "_dataset1_top5_PDL1_expression_VehicleMRTX_DCs.png", 12.8)


#####################
#### Supp Fig 4c ####

expression_boxplot(neighb2, c("Dendritic cells", "Dendritic cells CD103"), "MI_PDL1", "PD-L1", out_supp4, "_dataset2_top5_PDL1_expression_VehicleMRTX_DCs.png", 12.8)


################
#### Fig 4b ####

expression_boxplot(neighb1, c("Dendritic cells", "Dendritic cells CD103"), "MI_CD86", "CD86", out_fig4, "_dataset1_top5_CD86_expression_VehicleMRTX_DCs.png", 13)


#####################
#### Supp Fig 4d ####

expression_boxplot(neighb1, c("Dendritic cells", "Dendritic cells CD103"),"MI_MHCcII", "MHCII", out_supp4, "_dataset1_top5_MHCII_expression_VehicleMRTX_DCs.png", 13)
expression_boxplot(neighb2, c("Dendritic cells", "Dendritic cells CD103"), "MI_MHCcII", "MHCII", out_supp4, "_dataset2_top5_MHCII_expression_VehicleMRTX_DCs.png", 13)


################
#### Fig 4c ####

# Calculated logged values of CXCL9 expression 
neighb2$log_CXCL9 = log2(neighb2$MI_CXCL9)
expression_density(neighb2, c("Dendritic cells", "Dendritic cells CD103", "Macrophages type 1", "Macrophages type 2"), "log_CXCL9", "CXCL9", -10, 5, out_fig4,
                   "_dataset2_top5_CXCL9_expression_treatment_comparison_density.png", 6, 7)

################
#### Fig 4d ####

# May need to change this to look at all 5 communities 
tum_expr(neighb1, "MI_Ki67", c("T/DC", "T/M2_1", "T/M2_2"), top5_col[3:5], "Ki67", out_fig4, "_dataset1_top5_tumour_ki67_expression_MRTXonly.png")


#####################
#### Supp Fig 4e ####

tum_expr(neighb2, "MI_Ki67", c("T/DC", "T/M2_1", "T/M2_2"), top5_col[3:5], "Ki67", out_supp4, "_dataset2_top5_tumour_ki67_expression_MRTXonly.png")

################
#### Fig 4e ####

tum_expr(neighb1, "MI_casp3", c("T/DC", "T/M2_1", "T/M2_2"), top5_col[3:5], "C-casp3", out_fig4, "_dataset1_top5_tumour_casp3_expression_MRTXonly.png")


#####################
#### Supp Fig 4f ####

tum_expr(neighb2, "MI_casp3", c("T/DC", "T/M2_1", "T/M2_2"), top5_col[3:5], "C-casp3", out_supp4, "_dataset2_top5_tumour_casp3_expression_MRTXonly.png")

################
#### Fig 4f ####

# Calculated logged values of PD-1 expression 
neighb1$log_PD1 = log2(neighb1$MI_PD1)
neighb2$log_PD1 = log2(neighb2$MI_PD1)

expression_density(neighb1, "T cells CD8", "log_PD1", "PD-1", -7, 3, out_fig4, "_dataset1_top5_CD8_PD1_expression_treatment_comparison_density.png", 5.5, 6)

#####################
#### Supp fig 4g ####

expression_density(neighb2, "T cells CD8", "log_PD1", "PD-1", -10, 5, out_supp4, "_dataset2_top5_CD8_PD1_expression_treatment_comparison_density.png", 5.5, 7)


#####################
#### Supp Fig 4h ####

expression_density(neighb1, "T cells CD4", "log_PD1", "PD-1", -7, 3, out_supp4, "_dataset1_top5_CD4_PD1_expression_treatment_comparison_density.png", 5.5, 7)
expression_density(neighb2, "T cells CD4", "log_PD1", "PD-1", -10, 5, out_supp4, "_dataset2_top5_CD4_PD1_expression_treatment_comparison_density.png", 5.5, 7)


#####################
#### Supp Fig 4i ####

expression_density(neighb1, "T reg cells", "log_PD1", "PD-1", -7, 3, out_supp4, "_dataset1_top5_Treg_PD1_expression_treatment_comparison_density.png", 5.5, 7)
expression_density(neighb2, "T reg cells", "log_PD1", "PD-1", -10, 5, out_supp4, "_dataset2_top5_Treg_PD1_expression_treatment_comparison_density.png", 5.5, 7)


#####################
#### Supp Fig 4j ####

expression_boxplot(neighb2,  c("T cells CD4", "T cells CD8", "T reg cells"), "MI_LAG3", "LAG-3", out_supp4, "_dataset2_top5_LAG3_expression_VehicleMRTX_Tcells.png", 14)


#####################
#### Supp Fig 4k ####

## Percentage of PD-1+ CD8 T cells that are LAG3+ 
p = neighb2 %>%
  dplyr::select(cellType, top5, MI_PD1, MI_LAG3, treatment) %>%
  dplyr::filter(
    cellType == "T cells CD8" &
      MI_PD1 >= 0.5 &
      top5 %in% c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2") & treatment == "MRTX"
  ) %>%
  dplyr::group_by(top5) %>%
  dplyr::mutate(pd1_cd8_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(top5) %>%
  dplyr::summarise(
    pd1_cd8_count = unique(pd1_cd8_count),
    lag3_pd1_cd8_count = sum(MI_LAG3 > 0.5),
    percent_pos = (lag3_pd1_cd8_count / pd1_cd8_count) * 100
  ) %>%
  ggplot(aes(
    x = top5,
    y = percent_pos,
    fill = top5
  )) +
  geom_col(colour = "Black", position = "dodge") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    axis.text.x = element_text(size = 22, angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  scale_colour_manual(values = treat_col) +
  xlab("") +
  ylab("PD-1+ CD8 T cells: percentage LAG3+") +
  ylim(0, 40)
p
ggsave(
  plot = p,
  device = "png",
  width = 4.5,
  height = 7,
  dpi = 300,
  path = out_supp4,
  filename = paste(
    date,
    "_dataset2_PD1CD8_percetnage_LAG3pos_community_treatment.png",
    sep = ""
  )
)



