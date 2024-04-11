

#### Figure 1 - Introduction to communities & community validation ####

## Megan Cole (megan.cole@crick.ac.uk)

rm(list=ls())

# Load relevant libraries 
library(ggplot2) 
library(dplyr)
library(Rtsne) 
library(clustree) 
library(ggdendro) 
library(tiff) 
library(reshape2) 

##########################
#### Global variables ####

date = format(Sys.Date(), "%Y%m%d")

path = "/Users/colem/Documents/Projects/Treg_paper/Files/"

out_fig1 = "/Users/colem/Documents/Projects/Treg_paper/Figures/Figure1/"
dir.create(out_dir)
out_supp1 = "/Users/colem/Documents/Projects/Treg_paper/Figures/Supplementary_Figure1/"
dir.create(out_supp)

treat_col = c("MRTX" = "#00BFC4", "Vehicle" = "#F8766D")

#### Dataset 1 setup ####
data1 = read.csv(paste0(path, "20230904_dataset1_celldata.csv"))

neighb1 = read.csv(paste0(path, "20220823_dataset1_neighb_62communities_agglom18.csv"))

agglom30 = read.csv(paste0(path, "Figure1/20220216_dataset1_average_neighbours_62clusters_agglom30.csv"))

colours1 =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
              "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
              "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
              "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")

#### Dataset 2 setup ####
neighb2 = read.csv(paste0(path, "20220823_dataset2_neighb_254communities_agglom18.csv"))

colours2 = c("B cells" = "#945931", "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Fibroblasts" = "#ABDDA4FF",
             "Leukocytes unclassified" = "#7A9F79", "Macrophages type 1" = "#336666", "Macrophages type 2" = "#4E79A7FF",
             "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF", "T cells CD8" = "#FF9D9AFF", "T cells DN" = "#FFB5CB",
             "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")

##################################################################################################################################

###################
#### Functions ####

# Spatial location of communities per ROI
cell_map = function(ROI,
                    file,
                    files,
                    data,
                    clusters,
                    path,
                    labels) {
  d = data[which(data[, files] == ROI), ]
  
  TIFFnm = paste(
    "/Users/colem/Documents/Fiji/",
    file,
    "/",
    ROI,
    "/mask_outline/all_cells_mask.tiff",
    sep = ""
  )
  TIFFol = paste(
    "/Users/colem/Documents/Fiji/",
    file,
    "/",
    ROI,
    "/mask_outline/Cells_outline.tiff",
    sep = ""
  )
  
  # Read TIFFs of the relevant set
  TIFF = readTIFF(TIFFnm)
  TIFF2 = readTIFF(TIFFol)
  to_long_df <- function(TIFF) {
    names(TIFF) <- c(1:length(TIFF))
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
  n = 2
  
  ###############################################
  ## Choose the plotting of neighbour clusters ##
  for (cl in nbclusters) {
    print(cl)
    cluster_xy = d[which(d$agglom18_average == cl), c("ObjectNumber",
                                                      "Location_Center_X",
                                                      "Location_Center_Y",
                                                      "agglom18_average")]
    names(cluster_xy)[names(cluster_xy) == 'agglom18_average'] <-
      'cluster'
    #If cluster cl is not found in the image/ROI, remove it from the labels list
    if (dim(cluster_xy)[1] == 0) {
      print(paste("cluster ", cl, " is not found in ", ROI, sep = ""))
      next
    } else {
      cluster_xy$Location_Center_X = round(cluster_xy$Location_Center_X)
      cluster_xy$Location_Center_Y = round(cluster_xy$Location_Center_Y)
      names(cluster_xy) = c("ObjectNumber", "x", "y", "cluster")
      colours_in_mask <-
        inner_join(ROI1, cluster_xy[, c("x", "y", "cluster")])
      min = min(unique(colours_in_mask$value))
      print(min(unique(colours_in_mask$value)))
      colours_in_mask = colours_in_mask[-(which(colours_in_mask$value == min)), ]
      ROI2[which(ROI1$value %in% colours_in_mask$value), "value"] = n
      ROI2$name[ROI2$value == n] = first(cluster_xy$cluster)
      n = n + 1
    }
  }
  
  background = ROI2[which(ROI2$value < 1), "unique_px_ID"]
  ROI2[which(Outline$value == 1), "value"] = 1
  ROI2[which(ROI2$unique_px_ID %in% background), "value"] = 0
  ROI2$name[ROI2$value == 0] = 101
  ROI2$name[ROI2$value == 1] = 102
  
  ## Subset labels based on communities present in that image
  labels = labels[which(labels$cluster %in% unique(ROI2$name)), ]
  
  # Change 101 & 102 values to blank
  labels$cluster[labels$cluster == 101 |
                   labels$cluster == 102] <- " "
  
  p = ggplot(ROI2, aes(
    x = x,
    y = -y,
    fill = as.factor(value)
  )) +
    geom_raster() +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(colour = "black")) +
    scale_fill_manual(values = alpha(labels$colour,  1),
                      labels = labels$cluster)
  
  filename = paste(ROI, ".pdf", sep = "")
  ggsave(
    plot = p,
    device = "pdf",
    width = 5.6,
    height = 5,
    dpi = 300,
    path = path,
    filename = filename
  )
  print("plot saved")
  return(p)
}

# Stacked bar plot  
stacked_bar = function(data,
                             path,
                             filename,
                             position,
                             colour = NULL,
                             height = 6.5,
                             width = 7.5) {
  p = ggplot(data, aes(
    x = reorder(factor(agglom18_average), cellType == "Tumour", sum),
    fill = as.factor(cellType)
  )) +
    theme_classic() +
    coord_flip() +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text = element_text(size = 16)
    ) +
    theme(legend.text = element_text(size = 16)) +
    scale_x_discrete(limits = rev) +
    xlab("Community") +
    labs(fill = "")
  if (position == "fill") {
    p = p + geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution")
  } else if (position == "count") {
    p = p + geom_bar(stat = "count") + ylab("Cell count")
  }
  if (!is.null(colour)) {
    p = p + scale_fill_manual(values = colour)
  }
  p
  ggsave(
    plot = p,
    device = "png",
    width = width,
    height = height,
    dpi = 300,
    path = path,
    filename = paste0(filename, ".png")
  )
}

##################################################################################################################################

###################
#### Figure 1a ####

nbclusters = sort(unique(neighb1$agglom18_average))

colours = c("#9c4767", "#5bb84d", "#a35cca", "#b6b236", "#586ccd", "#dc9541", "#4ca2d5", "#cf4e32", "#54be9d",
            "#c947a0", "#3c8149", "#d84169", "#99ab5b", "#9080c3", "#737129", "#da84bb", "#a16434", "#dc8074")

labels = data.frame('cluster' = nbclusters,
                       'colour' = colours)
add = data.frame(
  cluster = c("101", "102"),
  colour = c("white", "white")
)
labels = rbind(add, labels)


images1 = c("20190913_BRAC3529.2d_ROI1_MRTX", "20190917_BRAC3495.3f_ROI1_Vehicle_crop1", "20190917_BRAC3495.3f_ROI1_Vehicle_crop2",
            "20190927_BRAC3529.2b_ROI1_MRTX_crop2", "20191119_BRAC3326.4e_ROI1_Vehicle_crop1", "20191121_BRAC3438.6f_ROI1_Vehicle",
            "20191121_BRAC3438.6f_ROI2_Vehicle", "20191121_BRAC3438.6f_ROI3_Vehicle", "20200130_BRAC4002.3c_ROI1_MRTX", 
            "20200130_BRAC4002.3c_ROI2_MRTX_crop1", "20200130_BRAC4002.3c_ROI2_MRTX_crop2", "20200130_BRAC4002.3c_ROI3_MRTX")


out_dir = paste0(out_fig1, "mapping_community_location_dataset1/")
dir.create(out_dir)

for(ROI in images1){
  print(ROI)
  p = cell_map(ROI, "Cont_MRTX", "ROI_name", neighb1, nbclusters, out_dir, labels)
}


###################
#### Figure 1b ####

# Cluster tree of neighbour clustering results 
neighb1$agglom30_average = agglom30$agglomerate_30[match(neighb1$cluster, agglom30$cluster)]

# Subset data and rename columns 
subset = neighb1 %>% 
  select(B_cells:cluster, agglom18_average, agglom30_average) %>% 
  dplyr::rename(Communities62 = cluster,
                Communities30 = agglom30_average,
                Communities18 = agglom18_average)

p = clustree(subset, prefix = "Communities")
ggsave(plot = p, device = "png", width=14, height=10, dpi=300, path = out_fig1,
       filename = paste(date, "_clustree_18_30_60_clusters.png", sep = ""))



###################
#### Figure 1c ####

# tSNE of clustering results with altered k-value - illustrating the reproducibility and stability of clusters
k350 = read.csv(paste0(path, "Figure1/20220216_dataset1_k350_neighbour_clustering_40communities.csv"))

avg_neighb1 = neighb1 %>% 
  select(3:19) %>% 
  group_by(cluster) %>% 
  summarise_all(list(mean))

names(avg_neighb1) = gsub("[.]", "_", names(avg_neighb1))

k350_avg = k350 %>% 
  select(3:19) %>% 
  group_by(cluster) %>% 
  summarise_all(list(mean))

avg_neighb1$k_value = 250
k350_avg$k_value = 350

# Bind rows of averages together
comparing_k = rbind(avg_neighb1, k350_avg)

# Run tSNE & plot results 
tSNE_results = Rtsne(comparing_k[,c(2:17)], perplexity = 10, num_threads =4, verbose = TRUE, check_duplicates = FALSE)

comparing_k$tSNE1 = tSNE_results$Y[,1]
comparing_k$tSNE2 = tSNE_results$Y[,2]

p = ggplot(comparing_k, aes(x=tSNE1, y=tSNE2, color = as.factor(k_value))) +
  geom_point(size = 1.8)+ 
  theme_classic()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  labs(colour = 'k value') +
  scale_colour_manual(values=c("#00203FFF", "#ADEFD1FF"))
p
# Save plot 
ggsave(plot = p, device = "png", width=6, height=6, dpi=300, path = out_fig1,
       filename = paste(date, "_tSNE_k250_k350_comparison.png", sep = ""))


###################################
#### Supplementary Figure 1a-b ####

# a
# Table of antibodies used for dataset 1 

# b 
# Table of antibodies used for dataset 2 

###################
#### Figure 1d ####

# Stacked bar of 18 communities for dataset 1 - ordered by tumour cell count 
stacked_bar(
  neighb1,
  out_fig1,
  paste0(
    date,
    "_dataset1_stacked_bar_18communities_cellCount_orderBy_tumourCount"
  ),
  position = "count",
  colour = colours1,
  height = 6.5,
  width = 8
)



###################
#### Figure 1e ####

# Stacked bar of 18 clusters for dataset 2
stacked_bar(
  neighb2,
  out_fig1,
  paste0(
    date,
    "_dataset2_stacked_bar_18communities_cellCount_orderBy_tumourCount"
  ),
  position = "count",
  colour = colours2,
  height = 6.5,
  width = 8
)

#################################
#### Supplementary Figure 1c ####

## Plotting cell type proportions of the original 62 communities - with dendrogram attached ##

# Run clustering of 62 communities based on cell type contribution
hc <- hclust(dist(avg_neighb1[,-1]), "average") # Hierarchical clustering 
dendr <- dendro_data(hc, type="rectangle")

label = avg_neighb1$cluster

# Plot dendrogram 
p = ggdendrogram(hc, rotate = FALSE, size = 3) # Plot dendrogram 
p
ggsave(plot = p, device = "png", width=15, height=12, dpi=300, path = out_supp1,
       filename = paste(date, "_neighb1_ggdendrogram_avg_neighb_62communities.png", sep = "")) 

# Order based on dendrogram 
com_order = c(59,18,54,51,7,21,49,3,14,8,4,6,35,60,26,17,15,19,5,2,24,36,9,45,11,13,23,16,20,10,12,1,
              22,25,33,32,55,52,44,50,37,39,61,48,58,53,56,57,41,29,40,28,30,42,43,46,34,27,62,31,38,47)


p = ggplot(neighb1, aes(x = reorder(cluster, cellType == "Tumour", FUN = mean), fill = as.factor(cellType)))+
  scale_fill_manual(values = colours1) +
  geom_bar(position = "fill") +
  theme_classic()+
  scale_x_discrete(limits = factor(com_order)) + 
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20)) +
  xlab("")+
  ylab("Percentage distribution") +
  labs( fill = "")
p
ggsave(plot = p, device = "png", width=12, height=14, dpi=300, path = out_supp1,
       filename = paste(date, "_neighb1_stacked_barplot_62_communities.png", sep = ""))



#####################################
#### Supplementary Figure 1d & e ####

# d
stacked_bar(neighb1, out_supp1, paste0(date, "_dataset1_stacked_bar_18communities_proportions_orderBy_tumourCount"),  position = "fill", colour = colours1)

# e
stacked_bar(neighb2, out_supp1, paste0(date, "_dataset2_stacked_bar_18communities_proportions_orderBy_tumourCount"), position = "fill", colour = colours2)


###################
#### Figure 1f ####

# Similarity of neighbour clustering results for overlapping cell types in dataset 1 & dataset 2 - tSNE

overlap1_avg = read.csv(paste0(path, "Figure1/20220421_dataset1_overlapping_cellTypes_average_neighbours_18clusters.csv"))
overlap1_avg$dtype = "Dataset 1"

overlap2_avg = read.csv(paste0(path, "Figure1/20220421_dataset2_overlapping_cellTypes_average_neighbours_18clusters.csv"))
overlap2_avg$dtype = "Dataset 2"

avg18 = rbind(overlap1_avg, overlap2_avg)

## tSNE 
tSNE = Rtsne(avg18[,c(2:11)], perplexity = 5, num_threads = 4, verbose = TRUE, check_duplicates = FALSE)

avg18$tSNE1 = tSNE$Y[,1]
avg18$tSNE2 = tSNE$Y[,2]

p = ggplot(avg18, aes(x=tSNE1, y=tSNE2, color = as.factor(dtype))) +
  geom_point(size = 2.5)+ 
  theme_classic()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  labs(colour = "") +
  scale_colour_manual(values=c("#5F4B8BFF", "#E69A8DFF"))
p
# #Save plot 
ggsave(plot = p, device = "png", width=7, height=6, dpi=300, path = out_supp1,
       filename = paste(date, "_tSNE_dataset1_dataset2_overlapping_cellTypes_community_comparison.png", sep = ""))
# Save avg18
write.csv(avg18, paste0(path, "Figure1/", date, "_dataset1_dataset2_overlapping_cellTypes_community_agglom30_average_PCA_tSNE.csv"), row.names = F)








