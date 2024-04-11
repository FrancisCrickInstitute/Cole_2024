

#### Scripts for Treg paper - Figure 3 & Supplementary figure 3 ####

## Megan Cole (megan.cole@crick.ac.uk)

rm(list=ls())

# Load relevant libraries 
library(ggplot2) 
library(dplyr) 
library(tiff) 
library(reshape) 
library(tibble)

date = format(Sys.Date(), "%Y%m%d")

## Global variables ##

path = "/Users/colem/Documents/Projects/Treg_paper/Golden/Data/"

neighb1 = read.csv(paste0(path, "20231002_neighb1.csv"))
neighb2 = read.csv(paste0(path, "20231002_neighb2.csv"))

out_fig3 = "/Users/colem/Documents/Projects/Treg_paper/Figures/Figure3/"
dir.create(out_fig3)
out_supp3 = "/Users/colem/Documents/Projects/Treg_paper/Figures/Supplementary_Figure3/"
dir.create(out_supp3)

colours1 =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
              "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
              "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
              "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")

colours2 = c("B cells" = "#945931", "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Fibroblasts" = "#ABDDA4FF",
             "Leukocytes unclassified" = "#7A9F79", "Macrophages type 1" = "#336666", "Macrophages type 2" = "#4E79A7FF",
             "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF", "T cells CD8" = "#FF9D9AFF", "T cells DN" = "#FFB5CB",
             "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")

treat_col = c("Vehicle" = "#F8766D", "MRTX" = "#00BFC4")


###################
#### Functions ####
###################

# Used for Fig 3a, Supp fig 3c
stacked_bar = function(data,
                       X,
                       fill,
                       path,
                       filename,
                       position = "fill",
                       colour = NULL,
                       width = 8.5,
                       height = 8){
  
  p = ggplot(data, aes(x = as.factor(get(X)), fill = as.factor(get(fill))))+
    theme_classic()+
    coord_flip() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text=element_text(size=22)) +
    theme(legend.text = element_text(size = 16)) +
    xlab("Community")+
    labs(fill = "")
  if(position == "fill"){
    p = p + geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution")
  } else if(position == "count"){
    p = p + geom_bar(stat = "count") + ylab("Cell count")
  }
  if(!is.null(colour)){
    p = p + scale_fill_manual(values = colour)
  }
  p
  ggsave(plot = p, device = "png", width=width, height=height, dpi=300, path = path, filename = paste0(filename, ".png")
  )
}

# Cumulative flow diagram of top 5 communities across the tissue
top5_count_cumulative = function(treat, x1, x2, levels, name) {
  
  # Cumulative flow diagram - Vehicle 
  p = ggplot(neighb1[which(neighb1$scaled_Y <= 0.5 &
                             neighb1$scaled_Y >= -0.5 &
                             neighb1$domain != "n/a" &
                             neighb1$treatment == treat &
                             neighb1$top5 %in% top5),], aes(x= scaled_X)) +
    annotate("rect", xmin=-Inf, xmax=-1.1, ymin=0, ymax=Inf, alpha=0.5, fill="#F3F3F3") +
    annotate("rect", xmin=-1.1, xmax=-0.9, ymin=0, ymax=Inf, alpha=0.5, fill="#979797") +
    annotate("rect", xmin=-0.9, xmax= 0.9, ymin=0, ymax=Inf, alpha=0.5, fill="#000000") +
    annotate("rect", xmin= 0.9, xmax= 1.1, ymin=0, ymax=Inf, alpha=0.5, fill="#979797") +
    annotate("rect", xmin= 1.1, xmax= Inf, ymin=0, ymax=Inf, alpha=0.5, fill="#F3F3F3") +
    annotate("text", x = x1, y = 16000, label = "Normal", hjust = 0.5, size = 5) +
    annotate("text", x = x2, y = 16000, label = "Normal", hjust = 0.5, size = 5) +
    annotate("text", x = 0, y = 16000, label = "Tumour", hjust = 0.5, size = 5) +
    geom_density(aes(y = ..count.., fill = factor(top5, levels = levels)), colour = "grey") +
    xlim(-2,2.25) +
    ylim(0, 16000) + 
    theme_classic() +
    theme(plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18)) +
    guides(fill=guide_legend(title="Community")) +
    xlab("Cross section through tissue") +
    ylab("Counts") + 
    scale_fill_manual(values =  top5_col[1:5])
  p 
  ggsave(plot = p, device = "png", width=8, height=6, bg='white', dpi=300, path = out_fig3,
         filename = paste0(date, name))
}


# Spatial location of communities on each ROI
cell_map = function(ROI, file, files, data, clusters, path, labels){
  
  d = data[which(data[,files] == ROI),]
  
  TIFFnm = paste("/Users/colem/Documents/Fiji/", file, "/", ROI, "/mask_outline/all_cells_mask.tiff", sep = "")
  TIFFol = paste("/Users/colem/Documents/Fiji/", file, "/", ROI, "/mask_outline/Cells_outline.tiff", sep = "")
  
  # Read TIFFs of the relevant set
  TIFF = readTIFF(TIFFnm)
  TIFF2 = readTIFF(TIFFol)
  to_long_df <- function(TIFF) {
    names(TIFF)<- c(1:length(TIFF))
    TIFF <- melt(TIFF, id = c(row(TIFF), names(TIFF)))
    names(TIFF)[1] <- "y"
    names(TIFF)[2] <- "x"
    TIFF
  }  
  Outline = TIFF2
  ROI1 = TIFF
  ROI1 = to_long_df(ROI1)
  Outline = to_long_df(Outline)
  ROI1$unique_px_ID = c(1:nrow(ROI1))
  ROI2 = ROI1
  n=2
  
  ###############################################
  ## Choose the plotting of neighbour clusters ##
  for (cl in nbclusters){
    print(cl)
    cluster_xy = d[which(d$top5 == cl), c("ObjectNumber","Location_Center_X","Location_Center_Y", "top5")]
    names(cluster_xy)[names(cluster_xy) == 'top5'] <- 'cluster' 
    #If cluster cl is not found in the image/ROI, remove it from the labels list
    if (dim(cluster_xy)[1] == 0){
      print(paste("cluster ", cl, " is not found in ", ROI, sep = ""))
      next
    } else {
      cluster_xy$Location_Center_X = round(cluster_xy$Location_Center_X)
      cluster_xy$Location_Center_Y = round(cluster_xy$Location_Center_Y)
      names(cluster_xy) = c("ObjectNumber","x","y","cluster")
      colours_in_mask <- inner_join(ROI1, cluster_xy[,c("x","y","cluster")])
      min = min(unique(colours_in_mask$value))
      print(min(unique(colours_in_mask$value)))
      colours_in_mask = colours_in_mask[-(which(colours_in_mask$value == min)),]
      ROI2[which(ROI1$value %in% colours_in_mask$value), "value"] = n
      ROI2$name[ROI2$value == n] = first(cluster_xy$cluster)
      n = n+1
    }
  }
  
  ROI2$name = ifelse(ROI2$name == "1", "T/NA",
                     ifelse(ROI2$name == "2", "T/M1",
                            ifelse(ROI2$name == "3", "T/DC",
                                   ifelse(ROI2$name == "4", "T/M2_1",
                                          ifelse(ROI2$name == "5", "T/M2_2", "NA")))))
  
  background = ROI2[which(ROI2$value < 1),"unique_px_ID"]
  ROI2[which(Outline$value == 1 ),"value"] = 1
  ROI2[which(ROI2$unique_px_ID %in% background),"value"] = 0
  ROI2$name[ROI2$value == 0] = 101
  ROI2$name[ROI2$value == 1] = 102
  
  ## Subset labels based on communities present in that image
  labels = labels[which(labels$cluster %in% unique(ROI2$name)),]
  
  # Change 101 & 102 values to blank
  labels$cluster[labels$cluster == 101 | labels$cluster == 102] <- " "
  
  p = ggplot(ROI2,aes(x=x,y=-y, fill=as.factor(value))) +
    geom_raster() +
    theme_void() +
    theme(legend.title=element_blank(),
          legend.text = element_text(colour = "black")) +
    scale_fill_manual(values = alpha(labels$colour,  1), labels = labels$cluster)
  
  filename = paste(ROI, ".pdf", sep = "")
  ggsave(plot = p, device = "pdf", width=5.6, height=5, dpi=300, path = path, filename = filename)
  print("plot saved")
  return(p)
}

# T cells count of each community, as a proportion across the top5 communities as a whole
cd8_top5_perc = function(data, name) {
  
  p = data %>%
    filter(cellType == "T cells CD8",
           top5 %in% c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2")) %>% 
    select(cellType, top5, treatment) %>% 
    ggplot(aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), fill = factor(top5, levels = c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2")))) +
    geom_bar(position = "fill", colour = "Black") + 
    theme_classic() + 
    scale_fill_manual(values = top5_col[1:5]) + 
    xlab("") +
    ylab("Percentage of all CD8+ T cells") +
    labs(fill = "Community") + 
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18))
  ggsave(plot = p, device = "png", width=5, height=5, dpi=300, path = out_supp3,
         filename = paste0(date, name))
  
}

########################################################################################

#################################
#### Supplementary Figure 3a ####

# Determining which communities have the highest CD8+ T cell count - Dataset 1
sub1 = neighb1[which(neighb1$cellType == "T cells CD8"),]
sub1 <- within(sub1,
                  agglom18_average <- factor(agglom18_average,
                                             levels=names(sort(table(agglom18_average),
                                                               decreasing = TRUE))))

p = ggplot(sub1,
           aes(x = as.factor(agglom18_average), fill = as.factor(treatment)))+
  geom_bar(stat = "count") + 
  scale_fill_manual(values = treat_col) + 
  theme_classic()+
  theme(legend.text = element_text(size = 22),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        legend.title = element_blank())+
  xlab("Community") +
  ylab("Cell count")
p 
ggsave(plot = p, device = "png", width=10, height=10, dpi=300, path = out_supp3,
       filename = paste(date, "_barchart_count_CD8_Tcells_perTreatment_dataset1.png", sep = ""))



#################################
#### Supplementary Figure 3b ####

# Determining which communities have the highest CD8+ T cell count - Dataset 2 
sub2 = neighb2[which(neighb2$cellType == "T cells CD8"),]
sub2 <- within(sub2,
               agglom18_average <- factor(agglom18_average,
                                          levels=names(sort(table(agglom18_average),
                                                            decreasing = TRUE))))

p = ggplot(sub2,
           aes(x = as.factor(agglom18_average), fill = as.factor(treatment)))+
  geom_bar(stat = "count") + 
  scale_fill_manual(values = treat_col) + 
  theme_classic()+
  theme(legend.text = element_text(size = 22),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        legend.title = element_blank())+
  xlab("Community") +
  ylab("Cell count")
p 
ggsave(plot = p, device = "png", width=10, height=10, dpi=300, path = out_supp3,
       filename = paste(date, "_barchart_count_CD8_Tcells_perTreatment_dataset2.png", sep = ""))


###################
#### Figure 3a ####

# Table of comparable communities and labels - created in powerpoint 

neighb1$top5 = "Other"
neighb1[which(neighb1$agglom18_average == 18), "top5"] <- "T/NA"
neighb1[which(neighb1$agglom18_average == 10), "top5"] <- "T/M1"
neighb1[which(neighb1$agglom18_average == 11), "top5"] <- "T/DC"
neighb1[which(neighb1$agglom18_average == 16), "top5"] <- "T/M2_1"
neighb1[which(neighb1$agglom18_average == 2), "top5"] <- "T/M2_2"
unique(neighb1$top5)

neighb2$top5 = "Other"
neighb2[which(neighb2$agglom18_average == 11), "top5"] <- "T/NA"
neighb2[which(neighb2$agglom18_average == 10), "top5"] <- "T/M1"
neighb2[which(neighb2$agglom18_average == 14), "top5"] <- "T/DC"
neighb2[which(neighb2$agglom18_average == 1), "top5"] <- "T/M2_1"
neighb2[which(neighb2$agglom18_average == 2), "top5"] <- "T/M2_2"
unique(neighb2$top5)

write.csv(neighb1, paste0(path, date, "_dataset1_neighb_62communities_agglom18.csv"), row.names = F)
write.csv(neighb2, paste0(path, date, "_dataset2_neighb_254communities_agglom18.csv"), row.names = F)

top5 = c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2")
top5_col = c("T/NA" = "#F8766D", "T/M1" = "#FCBF49", "T/DC" = "#00BF7D", "T/M2_1" = "#00B0F6", "T/M2_2" = "#E76BF3", "Other" = "#E1E1E1")

neighb1$top5 = factor(neighb1$top5, levels = c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2", "Other"))
neighb2$top5 = factor(neighb2$top5, levels = c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2", "Other"))

###################

# Sub-setting for communities with highest CD8 T cell count 

stacked_bar(neighb1[which(neighb1$top5 %in% top5),],
            'top5', 'cellType', out_fig3, paste0(date, "_dataset1_top_communities_CD8TC_cellType_proportions"), colour = colours1)

stacked_bar(neighb2[which(neighb2$top5 %in% top5),],
            'top5', 'cellType', out_fig3, paste0(date, "_dataset2_top_communities_CD8TC_cellType_proportions"), colour = colours2)


#################################
#### Supplementary Figure 3c ####

## Distribution of shared cell types that contribute to top 5 communities 
shared = c("B cells", "Dendritic cells", "Dendritic cells CD103", "Fibroblasts", "Macrophages type 1",
           "Macrophages type 2", "T cells CD4", "T cell CD8", "T reg cells", "Tumour")


stacked_bar(neighb1[which(neighb1$top5 %in% top5 & neighb1$cellType %in% shared),], 'top5', 'cellType', out_supp3,
            paste0(date, "_dataset1_top_communities_CD8TC_count_shared_cellTypes"), colour = colours1)

stacked_bar(neighb2[which(neighb2$top5 %in% top5 & neighb2$cellType %in% shared),],
            'top5', 'cellType', out_supp3, paste0(date, "_dataset2_top_communities_CD8TC_count_shared_cellTypes"), colour = colours2)


#################################
#### Supplementary Figure 3d ####

#### Proportion of all CD8+ T cells in 5 selected communities - for both dataset 1 & dataset 2 ####

cd8_1 = neighb1 %>% filter(cellType == "T cells CD8") %>% select(cellType, top5) %>% mutate(dtype = "Dataset 1")
cd8_2 = neighb2 %>% filter(cellType == "T cells CD8") %>% select(cellType, top5) %>% mutate(dtype = "Dataset 2")

cd8 = rbind(cd8_1, cd8_2)

p = ggplot(cd8,  aes(x = as.factor(dtype), fill = factor(top5, levels = c("Other", "T/M2_2", "T/M2_1", "T/DC", "T/M1", "T/NA")))) +
  geom_bar(position = "fill", colour = "Black") + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() + 
  scale_fill_manual(values = top5_col) + 
  xlab("") +
  ylab("Percentage of all CD8+ T cells") + 
  labs(fill = "") + 
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))
p
ggsave(plot = p, device = "png", width=6, height=6, dpi=300, path = out_supp3,
       filename = paste(date, "_CD8_TC_proportion_top5_dataset1_dataset2.png", sep = ""))


###################
#### Figure 3b ####

#### Cumulative flow diagram of top 5 communities for dataset 1 ####

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

# Vehicle 
top5_count_cumulative("Vehicle", -1.65, 1.75, c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2"), "_dataset1_cumulative_flow_top5_Vehicle_onMRTXscale.png")
# MRTX
top5_count_cumulative("MRTX", -1.67, 1.7, c("T/M2_2", "T/M2_1", "T/DC", "T/M1", "T/NA"), "_dataset1_cumulative_flow_top5_MRTX.png")



#################################
#### Supplementary Figure 3e ####

# Treatment group split of top 5 communities 

stacked_bar(neighb1[which(neighb1$top5 %in% top5),], 'top5', 'treatment', out_supp3,
            paste0(date, "_dataset1_top_communities_CD8TC_treatment"), colour = treat_col, width = 7)

stacked_bar(neighb2[which(neighb2$top5 %in% top5),], 'top5', 'treatment',
            out_supp3, paste0(date, "_dataset2_top_communities_CD8TC_treatment"), colour = treat_col, width = 7)


###################
#### Figure 3c ####

# Spatial location of communities in dataset 1 - Vehicle & MRTX 

#############################################
nbclusters = top5
labels = as.data.frame(top5_col[1:5])
labels = rownames_to_column(labels, "cluster")
colnames(labels)[2] <- 'colour'
add = data.frame(
  cluster = c("101", "102"),
  colour = c("black", "white")
)
labels = rbind(add, labels)

## Spatial location of communities in dataset 1 - Vehicle & MRTX ##

images1 = c("20190913_BRAC3529.2d_ROI1_MRTX", "20190917_BRAC3495.3f_ROI1_Vehicle_crop1", "20190917_BRAC3495.3f_ROI1_Vehicle_crop2",
            "20190927_BRAC3529.2b_ROI1_MRTX_crop2", "20191119_BRAC3326.4e_ROI1_Vehicle_crop1", "20191121_BRAC3438.6f_ROI1_Vehicle",
            "20191121_BRAC3438.6f_ROI2_Vehicle", "20191121_BRAC3438.6f_ROI3_Vehicle", "20200130_BRAC4002.3c_ROI1_MRTX", 
            "20200130_BRAC4002.3c_ROI2_MRTX_crop1", "20200130_BRAC4002.3c_ROI2_MRTX_crop2", "20200130_BRAC4002.3c_ROI3_MRTX")


out_dir = paste0(out_fig3, "mapping_community_location_dataset1/")
dir.create(out_dir)

for(ROI in images1){
  print(ROI)
  p = cell_map(ROI, "Cont_MRTX", "ROI_name", neighb1, top5, out_dir, labels)
}


###################

## Spatial location of communities in dataset 2 - Vehicle & MRTX ##

images2 = c("BRAC3326.4e_ROI1_Vehicle", "BRAC4002.3c_ROI2_MRTX", "BRAC4002.3c_ROI3_MRTX", "BRAC3438.6f_ROI1_Vehicle",
           "BRAC3438.6f_ROI3_Vehicle","BRAC3495.3f_ROI1_Vehicle", "BRAC3529.2a_ROI1_MRTX", "BRAC3529.2b_ROI1_MRTX",
           "BRAC3529.2d_ROI3_MRTX", "BRAC3708.2d_ROI1_Vehicle", "BRAC4002.3c_ROI1_MRTX")

out_dir = paste0(out_fig3, "mapping_community_location_dataset2/")
dir.create(out_dir)

for(ROI in images2){
  print(ROI)
  p = cell_map(ROI, "TCpanel", "ROI_name", neighb2, top5, out_dir, labels)
}


#################################
#### Supplementary Figure 3f ####  

## CD8 T cells count of each community, as a proportion across the top5 communities as a whole


cd8_top5_perc(neighb1, "_dataset1_CD8TC_percentage_in_top5_perTreatment.png") 
cd8_top5_perc(neighb2, "_dataset2_CD8TC_percentage_in_top5_perTreatment.png") 


#################################
#### Supplementary Figure 3g ####

# Spatial location of community B in vehicle & MRTX groups 
nbclusters = c("T/M1")
labelsB = filter(labels, cluster %in% c("101", "102", "T/M1"))

out_dir = paste0(out_supp3, "mapping_community_location_B_1/")
dir.create(out_dir)

# dataset 1 
for(ROI in images1){
  print(ROI)
  p = cell_map(ROI, "Cont_MRTX", "ROI_name", neighb1, nbclusters, out_dir, labelsB)
}

# dataset 2
for(ROI in images2){
  print(ROI)
  p = cell_map(ROI, "TCpanel", "ROI_name", neighb2, nbclusters, out_dir, labelsB)
}









