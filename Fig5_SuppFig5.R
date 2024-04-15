

#### Scripts for Treg paper - Figure 5 & Supplementary figure 5 ####

## Megan Cole (megan.cole@crick.ac.uk)

rm(list=ls())

# Load libraries 
library(dplyr) 
library(ggplot2) 
library(ggsignif)
library(tiff) 
library(reshape) 
library(scales) 
library(gplots) 
library(colorspace) 
library(ggpubr) 
library(data.table) 
library(neighbouRhood) 
library(parallel)

date = format(Sys.Date(), "%Y%m%d")

## Global variables ##

path = ""

neighb1 = read.csv(paste0(path, "neighb1.csv"))
neighb2 = read.csv(paste0(path, "neighb2.csv"))

seg_neighb1 = read.csv(paste0(path, "seg_neighb1.csv"))
seg_neighb2 = read.csv(paste0(path, "seg_neighb2.csv"))

out_fig5 = paste0(path, "Figure5/")
dir.create(out_fig5)
out_supp5 = paste0(path, "Supplementary_Figure5/")
dir.create(out_supp5)

file5 = paste0(path, "Files_Figure5/")

distance_path = paste0(out_fig5, "dataset2_communityC_min_distance/")

top5 = c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2")
top5_col = c("T/NA" = "#F8766D", "T/M1" = "#FCBF49", "T/DC" = "#00BF7D", "T/M2_1" = "#00B0F6", "T/M2_2" = "#E76BF3", "Other" = "#E1E1E1")
treat_col = c("Vehicle" = "#F8766D", "MRTX" = "#00BFC4")

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

seg_neighb2$ROI_short = neighb2$ROI_short[match(seg_neighb2$cellID, neighb2$cellID)]

neighb1 = neighb1 %>%
  mutate(ROI_short = case_when(ROI_name == "20190913_BRAC3529.2d_ROI1_MRTX"~"BRAC3529.2d_MRTX",
                               ROI_name == "20200130_BRAC4002.3c_ROI2_MRTX_crop1"~"BRAC4002.3c_t2_MRTX",
                               ROI_name == "20200130_BRAC4002.3c_ROI2_MRTX_crop2"~"BRAC4002.3c_t3_MRTX",
                               ROI_name == "20200130_BRAC4002.3c_ROI3_MRTX"~"BRAC4002.3c_t4_MRTX",
                               ROI_name == "20190917_BRAC3495.3f_ROI1_Vehicle_crop1"~"BRAC3495.3f_t1_Vehicle",
                               ROI_name == "20190917_BRAC3495.3f_ROI1_Vehicle_crop2"~"BRAC3495.3f_t2_Vehicle",
                               ROI_name == "20190927_BRAC3529.2b_ROI1_MRTX_crop2"~"BRAC3529.2b_MRTX",
                               ROI_name == "20191119_BRAC3326.4e_ROI1_Vehicle_crop1"~"BRAC3326.4e_Vehicle",
                               ROI_name == "20191121_BRAC3438.6f_ROI1_Vehicle"~"BRAC3438.6f_t1_Vehicle",
                               ROI_name == "20191121_BRAC3438.6f_ROI2_Vehicle"~"BRAC3438.6f_t2_Vehicle",
                               ROI_name == "20191121_BRAC3438.6f_ROI3_Vehicle"~"BRAC3438.6f_3_Vehicle",
                               ROI_name == "20200130_BRAC4002.3c_ROI1_MRTX"~"BRAC4002.3c_t1_MRTX",
                               TRUE ~ ROI_name))
seg_neighb1$ROI_short = neighb1$ROI_short[match(seg_neighb1$cellID, neighb1$cellID)]

##################################################################################################################

###################
#### Functions ####

# Plotting distance calculation in community C 
compare_dist = function(data, ylab, dir, name){
  
  p = ggplot(transform(data, type = factor(type, levels = c("CXCL9 low", "CXCL9 high"))),
             aes(x = factor(type), y = distance)) +
    geom_boxplot() +
    # geom_signif(comparisons = list(c("CXCL9 low", "CXCL9 high")), 
    #             map_signif_level=TRUE,) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          strip.text = element_text(size=14, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          legend.position = "right") +
    ylab(paste0("Min distance to PD-1+ ", ylab)) +
    scale_y_continuous(trans='log2') +
    facet_wrap(source_cluster~.) 
  ggsave(plot = p, device = "png", width=5.5, height=5.5, dpi=300, path = dir,
         filename = paste0(date, name))
}

# Visualisation of cell types in the tissue
plot_cellmap = function(ROI, file, files, data, column, clusters, path, cellTypes){
  
  d = data[which(data[,files] == ROI),]
  
  TIFFnm = paste0("/Users/colem/Documents/Fiji/", file, "/", ROI, "/mask_outline/all_cells_mask.tiff")
  TIFFol = paste0("/Users/colem/Documents/Fiji/", file, "/", ROI, "/mask_outline/Cells_outline.tiff")
  
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
  ROI2$name = ""
  
  ###############################################
  ## Choose the plotting of neighbour clusters ##
  for (cl in nbclusters){
    cluster_xy = d[which(d[,column] == cl), c("ObjectNumber","Location_Center_X","Location_Center_Y", column)]
    #If cluster cl is not found in the image/ROI, remove it from the labels list
    if (dim(cluster_xy)[1] == 0){
      print(paste(cl, " is not found in ", ROI, sep = ""))
      next
    } else {
      cluster_xy$Location_Center_X = round(cluster_xy$Location_Center_X)
      cluster_xy$Location_Center_Y = round(cluster_xy$Location_Center_Y)
      names(cluster_xy) = c("ObjectNumber","x","y",column)
      colours_in_mask <- inner_join(ROI1, cluster_xy[,c("x","y", column)])
      min = min(unique(colours_in_mask$value))
      colours_in_mask = colours_in_mask[-(which(colours_in_mask$value == min)),]
      ROI2[which(ROI1$value %in% colours_in_mask$value), "value"] = n
      ROI2$name[ROI2$value == n] = first(cluster_xy[,column])
      n = n+1
    }
  }
  
  background = ROI2[which(ROI2$value < 1),"unique_px_ID"]
  ROI2[which(Outline$value == 1 ),"value"] = 1
  ROI2[which(ROI2$unique_px_ID %in% background),"value"] = 0
  ROI2$name[ROI2$value == 0] = 0
  ROI2$name[ROI2$value == 1] = 1
  
  ## Subset labels based on communities present in that image
  labels = labels[which(labels$clustername %in% unique(ROI2$name)),]
  
  # Only keep colours for cell types of interest
  labels = labels %>%
    mutate(colours = ifelse(clustername %in% cellTypes, colours, "#E1E1E1"))
  
  # Add background and outline colours 
  add = data.frame('clustername' = c(0,1), 'colours' = c('Black', 'White'))
  labels = rbind(add, labels)
  
  p = ggplot(ROI2,aes(x=x,y=-y, fill=as.factor(name))) +
    geom_raster() +
    theme_void() +
    theme(legend.title=element_blank(),
          legend.text = element_blank()) +
    scale_fill_manual(values = alpha(labels$colour,  1), labels = labels$clustername)
  
  filename = paste(ROI, ".pdf", sep = "")
  ggsave(plot = p, device = "pdf", width=5.6, height=5, dpi=300, path = path, filename = filename)
  print("plot saved")
}

# Plotting expression of markers across cell types and treatment groups within community C
marker_exp = function(Y, marker, name){
  
  p = ggplot(transform(neighb2[which(neighb2$top5 == "T/DC" & neighb2$cellType %in% c("T cells CD4", "T cells CD8", "T reg cells")),],
                       treatment = factor(treatment, levels = c("Vehicle", "MRTX"))),
             aes(x = factor(treatment), y = get(Y), fill = factor(treatment))) +
    geom_boxplot() +
    # geom_signif(comparisons = list(c("Vehicle", "MRTX")), 
    #             map_signif_level=TRUE,) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          strip.text = element_text(size=16, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          legend.position = "right") +
    ylab(paste0(marker, " expression")) +
    scale_y_continuous(trans='log2', labels = label_number(accuracy = 0.0001)) + 
    scale_fill_manual(values = treat_col) +
    facet_wrap(cellType~.) 
  p
  ggsave(plot = p, device = "png", width=7, height=5.5, dpi=300, path = out_fig5,
         filename = paste0(date, name))
}

# Interaction of CD8+ T cells with casp3+ tumour cells Vehicle vs MRTX - dataset 2 - community C 
close_neighbours = function(n, neighbourdata, cellType, select, name, path, w = 6) {
  
  # How often are CD8 TC and casp3+ Tumour cells neighbouring each other?
  tum = n[which(n$cellType == "Tumour" & n$MI_casp3 >= 0.5 & n$top5 %in% select),]
  cd = n[which(n$cellType == cellType & n$top5 %in% select),]
  
  # Neighbours of casp3+ tumour & CD8 T cells
  tum_neighbours = neighbourdata[which(neighbourdata$cellID %in% tum$cellID),]
  cd_neighbours = neighbourdata[which(neighbourdata$cellID %in% cd$cellID),]
  
  # Data for CD8+ T cells that have a casp3+ tumour cell in their neighbourhood
  cd_casp3Tumour_neighbours = cd_neighbours[which(cd_neighbours$cellID2 %in% tum$cellID),]
  print(unique(cd_casp3Tumour_neighbours$First.Object.top5))
  
  
  #######################################################################################
  ## Separate by ROI and treatment - only including ROIs that have casp3+ tumour cells ## 
  #######################################################################################
  paired_count = cd_casp3Tumour_neighbours %>%
    select(treatment, ROI_short, First.Object.top5, Second.Object.top5) %>%
    dplyr::group_by(treatment, First.Object.top5, Second.Object.top5) %>%
    dplyr::mutate(treatment_count=n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(treatment, ROI_short, First.Object.top5, Second.Object.top5) %>% 
    dplyr::summarize(treatment_count = unique(treatment_count),
                     ROI_count=n())
  
  paired_count = paired_count[order(paired_count$treatment, decreasing = TRUE),]
  paired_count$ID = paste(paired_count$First.Object.top5, paired_count$treatment, sep = "_")
  
  relative_count(paired_count, n, select, name, path, w)
}

relative_count = function(paired_count, n, select, name, path, w){
  
  casp3_com_treat_perc = n %>%
    dplyr::select(cellType, top5, MI_casp3, treatment) %>%
    dplyr::group_by(top5, treatment) %>%
    # Split % of casp3+ tumour cells into treatment
    dplyr:: mutate(perc_community_treat = (length(which(cellType == "Tumour" & MI_casp3 >= 0.5))/length(which(cellType == "Tumour"))*100)) %>%
    dplyr:: select(top5, treatment, perc_community_treat) %>%
    dplyr:: summarise_all(list(mean)) %>%
    # Remove rows with NaN values
    dplyr:: filter_all(any_vars(!is.na(perc_community_treat))) %>% 
    dplyr:: filter_all(any_vars(perc_community_treat > 0)) %>% 
    dplyr:: filter(top5 %in% select) %>% 
    mutate(ID = paste(top5, treatment, sep = "_"))
  
  
  # Add % of casp3+ tumour cells for each treatment group
  paired_count$perc_community_treat = casp3_com_treat_perc$perc_community_treat[match(paired_count$treatment, casp3_com_treat_perc$treatment)]
  # Calculate larger % by smaller % 
  perc_diff = max(paired_count$perc_community_treat)/min(paired_count$perc_community_treat)
  
  #########################
  ## For treatment & ROI ##
  #########################
  # Find the minimum % value from calculating % of tumour cells that are casp3+
  min_perc = min(paired_count$perc_community_treat)
  # Take ROI casp3+ count values associated with that minimum % 
  ROI_count_min_perc = paired_count$ROI_count[which(paired_count$perc_community_treat == min_perc)]
  # Times those ROI count values by perc_diff and all others by 1 
  paired_count$relative_ROI_count = ifelse(paired_count$ROI_count %in% ROI_count_min_perc,
                                           paired_count$ROI_count*perc_diff, paired_count$ROI_count*1)
  paired_count[is.na(paired_count)] <- 0
  
  # # Add 0 column per ROI and remove if there is a larger value
  df = data.frame("ROI_short" = unique(n$ROI_short), "relative_ROI_count" = rep(0, length(unique(n$ROI_short))))
  df$treatment = ifelse(grepl("MRTX", df$ROI_short), "MRTX", "Vehicle")
  df = df %>% select(treatment, ROI_short, relative_ROI_count)
  
  paired_count = paired_count %>% ungroup %>%  select(treatment, ROI_short, relative_ROI_count)
  paired_count = rbind(paired_count, df)
  paired_count = paired_count %>%
    dplyr::group_by(ROI_short) %>%
    dplyr::summarise(treatment = unique(treatment),
                     ROI_short = unique(ROI_short),
                     relative_ROI_count = max(relative_ROI_count))
  
  paired_count_mean = paired_count %>% 
    group_by(treatment) %>% 
    dplyr::summarise(treatment = unique(treatment),
                     mean_count = mean(relative_ROI_count))
  
  paired_count$ROI_short = factor(paired_count$ROI_short, levels = c( "BRAC3326.4e_Vehicle", "BRAC3438.6f_t1_Vehicle", "BRAC3438.6f_t2_Vehicle",
                                                                      "BRAC3495.3f_Vehicle", "BRAC3708.2d_Vehicle",
                                                                      "BRAC3529.2a_MRTX", "BRAC3529.2b_MRTX", "BRAC3529.2d_MRTX", 
                                                                      "BRAC4002.3c_t1_MRTX", "BRAC4002.3c_t2_MRTX", "BRAC4002.3c_t3_MRTX"
  ))
  
  p = ggplot(transform(paired_count, treatment = factor(treatment, levels = c("Vehicle", "MRTX")))) + 
    geom_col(data = transform(paired_count_mean, treatment = factor(treatment, levels = c("Vehicle", "MRTX"))),
             aes(x = factor(treatment), y = mean_count, fill = factor(treatment)),
             position = position_dodge2(width = 0, preserve = "single"), colour = "Black", fill = c("#F8766D", "#00BFC4")) + 
    geom_dotplot(binaxis='y', stackdir='center', stackgroups = TRUE, binpositions="all", dotsize = 0.8,
                 # colour = "Black", aes(x = factor(treatment), y = relative_ROI_count, fill = factor(ROI_short))) +
                 colour = "Black", fill = "Black", aes(x = factor(treatment), y = relative_ROI_count)) +
    theme_classic() + 
    xlab("") + 
    ylab("Relative count") + 
    theme(axis.title = element_text(size = 20), 
          axis.text = element_text(size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size = 18))
  p
  # Save plot
  ggsave(plot = p, device = "png", width=w, height=6, dpi=300, path = path,
         filename = paste0(date, name))
}

# Interaction of CD8+ T cells with casp3+ tumour cells in communities C, D & E - MRTX only 
close_neighbours_tumour = function(n, neighbourdata, select, name, relative, width = "7", height = "6") {
  
  # How often are CD8 TC and casp3+ Tumour cells neighbouring each other?
  tum = n[which(n$treatment == "MRTX" & n$cellType == "Tumour" & n$MI_casp3 >= 0.5 & n$top5 %in% select),]
  cd8 = n[which(n$treatment == "MRTX" & n$cellType == "T cells CD8" & n$top5 %in% select),]
  
  # Neighbours of casp3+ tumour & CD8 T cells
  tum_neighbours = neighbourdata[which(neighbourdata$cellID %in% tum$cellID),]
  cd_neighbours = neighbourdata[which(neighbourdata$cellID %in% cd8$cellID),]
  
  # Data for CD8+ T cells that have a casp3+ tumour cell in their neighbourhood
  cd8_casp3Tumour_neighbours = cd_neighbours[which(cd_neighbours$cellID2 %in% tum_neighbours$cellID),]
  unique(cd8_casp3Tumour_neighbours$First.Object.top5)
  
  paired_count = cd8_casp3Tumour_neighbours %>%
    select(First.Object.top5, ROI_short) %>%
    dplyr::group_by(First.Object.top5) %>%
    dplyr::mutate(count=n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(First.Object.top5, ROI_short) %>% 
    dplyr::summarize(count = unique(count),
                     ROI_count=n())
  
  if(relative == TRUE){
    
    relative_count_tumour(paired_count, n, select, name)
    
  } else {
    # Add 0 column per ROI and remove if there is a larger value
    df = n %>% 
      select(ROI_short, top5, treatment) %>% 
      filter(treatment == "MRTX",
             top5 %in% select) %>% 
      group_by(ROI_short, top5) %>% 
      dplyr::summarise(count = n()) %>% 
      dplyr::mutate(ROI_count = 0) %>% 
      select(top5, ROI_short, ROI_count)
    
    # Add paired_count and df together and take larger value 
    names(paired_count)[names(paired_count) == 'First.Object.top5'] <- 'top5'
    paired_count = paired_count %>% select(top5, ROI_short, ROI_count)
    paired_count = rbind(paired_count, df)
    paired_count = paired_count %>%
      dplyr::group_by(ROI_short, top5) %>%
      dplyr::summarise(top5 = unique(top5),
                       ROI_short = unique(ROI_short),
                       ROI_count = max(ROI_count))
    
    paired_count_mean = paired_count %>% 
      group_by(top5) %>% 
      dplyr::summarise(top5 = unique(top5),
                       mean_count = mean(ROI_count))
    
    paired_count$top5 = factor(paired_count$top5, levels = c("C", "D", "E"))
    
    p = ggplot(transform(paired_count, top5 = factor(top5, levels = c("C", "D", "E")))) + 
      geom_col(data = transform(paired_count_mean, top5 = factor(top5, levels = c("C", "D", "E"))),
               aes(x = factor(top5), y = mean_count, fill = factor(top5)),
               position = position_dodge2(width = 0, preserve = "single"), colour = "Black", fill = top5_col[3:5]) + 
      geom_dotplot(binaxis='y', stackdir='center', stackgroups = TRUE, binpositions="all", dotsize = 0.8,
                   #   colour = "Black", aes(x = factor(top5), y = ROI_count, fill = factor(ROI_short))) +
                   colour = "Black", fill = "Black", aes(x = factor(top5), y = ROI_count)) +
      theme_classic() + 
      xlab("") + 
      ylab("Number of interactions") + 
      theme(axis.title = element_text(size = 22), 
            axis.text = element_text(size = 22),
            legend.title = element_blank(),
            legend.text = element_text(size = 18))
    # Save plot
    ggsave(plot = p, device = "png", width=5, height=6, dpi=300, path = out_supp5,
           filename = paste0(date, name))
  }
}

# Create object relationship files ready for enrichment analysis 
object_rel_top5 = function(obj_r, obj_t){
  
  ##### MRTX only #####
  
  # Split based on top 5 communities 
  split = obj_t %>% split(f = as.factor(.$top5))
  
  
  rel_a = obj_r %>% 
    filter(cellID %in% split[["T/NA"]]$cellID) %>% 
    filter(First.Object.top5 == Second.Object.top5)
  
  rel_b = obj_r %>% 
    filter(cellID %in% split[["T/M1"]]$cellID) %>% 
    filter(First.Object.top5 == Second.Object.top5)
  
  rel_c = obj_r %>% 
    filter(cellID %in% split[["T/DC"]]$cellID) %>% 
    filter(First.Object.top5 == Second.Object.top5)
  
  rel_d = obj_r %>% 
    filter(cellID %in% split[["T/M2_1"]]$cellID) %>% 
    filter(First.Object.top5 == Second.Object.top5)
  
  rel_e = obj_r %>% 
    filter(cellID %in% split[["T/M2_2"]]$cellID) %>% 
    filter(First.Object.top5 == Second.Object.top5)
  
  print("Community relationship data separated")
  
  out = list("tna" = rel_a, "tm1" = rel_b, "tdc" = rel_c, "tm2_1" = rel_d, "tm2_2" = rel_e)
  return(out) 
}

# Run neighbouRhood enrichment 
neighbourhood_permutation = function(cdata, ndata, community, select_first, select_second, extra, out){
  
  # Convert object table to data.table format
  cdata = data.table(cdata)
  
  names(ndata) = gsub(".Object.Name", " Object Name", names(ndata))
  names(ndata) = gsub(".Image.Number", " Image Number", names(ndata))
  names(ndata) = gsub(".Object.Number", " Object Number", names(ndata))
  names(ndata) = gsub("Module.Number", "Module Number", names(ndata))
  ndata = data.table(ndata)
  
  data = prepare_tables(cdata, ndata, objname = "AllCellsMask")
  
  data_baseline = apply_labels(data[[1]], data[[2]]) %>%
    aggregate_histo()
  write.csv(data_baseline, paste0(out_stats5, date, "_dataset2_neighbouRhood_data_baseline_", community, "_MRTX", extra), row.names = F)
  
  # Calculate permutation statistics
  set.seed(12312)
  dat_perm = rbindlist(mclapply(1:1000, function(x) {
    dat_labels = shuffle_labels(data[[1]])
    apply_labels(dat_labels, data[[2]]) %>%
      aggregate_histo()
  }, mc.cores = 10), idcol = 'run')
  write.csv(dat_perm, paste0(out_stats5, date, "_dataset2_neighbouRhood_dat_perm_", community, "_MRTX", extra, ".csv"), row.names = F)
  
  # Calculate p-values
  dat_p <-
    calc_p_vals(data_baseline,
                dat_perm,
                n_perm = 1000,
                p_tresh = 0.01)
  write.csv(dat_p, paste0(out_stats5, date, "dataset2_neighbouRhood_dat_pvalues_", community, "_MRTX", extra, ".csv"), row.names = F)
  
  stats = list("data_baseline" = data_baseline, "dat_perm" = dat_perm, "dat_p" = dat_p)
  
  print("Permutation complete")
  
  # Call data prep
  data_prep(stats$data_baseline, stats$dat_perm, stats$dat_p, select_first, select_second, community, extra, out)
}

# Prep data for plotting 
data_prep = function(data_baseline, dat_perm, dat_p, select_first, select_second, community, extra, out){
  
  print("Data prep")
  dat_perm = as.data.frame(dat_perm)
  
  # Create pairing ID based on cell type pair and group 
  data_baseline$ID = paste(
    data_baseline$FirstLabel,
    data_baseline$SecondLabel,
    data_baseline$group,
    sep = "_"
  )
  dat_perm$ID = paste(
    dat_perm$FirstLabel,
    dat_perm$SecondLabel,
    dat_perm$group,
    sep = "_"
  )
  
  # Get an average of the ct values from each permutation run
  dat_perm_mean = dat_perm[,-1] %>%
    group_by(ID) %>%
    summarise_each(funs(if (is.numeric(.))
      mean(., na.rm = TRUE)
      else
        first(.)))
  
  # Subset for FirstLabel of interest
  
  data_baseline = as.data.frame(data_baseline %>% filter(FirstLabel == "T cells CD8"))
  dat_perm_mean = dat_perm_mean %>% filter(FirstLabel == "T cells CD8")
  
  # Check if dat_perm_mean and dat_baseline have the same number of objects
  data_baseline$group = as.numeric(data_baseline$group)
  dat_perm_mean$group = as.numeric(dat_perm_mean$group)
  difference = setdiff(dat_perm_mean[, c("ID", "group", "FirstLabel", "SecondLabel")], data_baseline[, c("ID", "group", "FirstLabel", "SecondLabel")])
  
  if (nrow(difference) > 0) {
    # List the values that are missing
    print(difference$ID)
    # Create a vector of zeros for the length of values that are missing
    difference$ct = 0
    # Order difference columns to match order of data_baseline, ready for binding
    difference = difference[, c("group",
                                "FirstLabel",
                                "SecondLabel",
                                "ct",
                                "ID")]
    data_baseline = rbind(data_baseline, difference)
    if (nrow(data_baseline) == nrow(dat_perm_mean)) {
      print("Dimensions of baseline and mean permutation statistics are now equal")
    } else {
      # Check if there is an interaction in baseline data not re-created in permutation data
      print("Checking if there is an interaction in baseline data that is not re-created in permutation data ")
      difference2 = as.data.frame(setdiff(data_baseline[, c("ID", "group", "FirstLabel", "SecondLabel")], dat_perm_mean[, c("ID", "group", "FirstLabel", "SecondLabel")]))
      print(difference2$ID)
      difference2$ct = 0
      difference2 = difference2[, c("ID",
                                    "group",
                                    "FirstLabel",
                                    "SecondLabel",
                                    "ct")]
      dat_perm_mean = rbind(dat_perm_mean, difference2)
      if (nrow(data_baseline) == nrow(dat_perm_mean)) {
        print("Dimensions of baseline and mean permutation statistics are now equal")
      } else {
        paste0(
          paste0(
            "Dimensions of baseline and mean permutation statistics are not equal - nrow baseline = ",
            nrow(data_baseline),
            ", nrow mean permutation = ",
            nrow(dat_perm_mean)
          )
        )
      }
    }  
  } else {
    print("Dimensions of baseline and mean permutation statistics are equal")
  }
  
  # Order baseline & mean data by ID, so it matches mean permutation data
  data_baseline = data_baseline[order(data_baseline$ID),]
  dat_perm_mean = dat_perm_mean[order(dat_perm_mean$ID),]
  
  # Check that IDs from baseline and mean permutation statistics match
  match = dat_perm_mean$ID == data_baseline$ID
  table(match)["FALSE"]
  table(match)["TRUE"]
  
  # Add baseline ct values to mean permutation data frame
  dat_perm_mean$ct_baseline = data_baseline$ct
  
  # Create pairing ID based on cell type pair and group 
  dat_p$ID = paste(
    dat_p$FirstLabel,
    dat_p$SecondLabel,
    dat_p$group,
    sep = "_"
  )
  
  # Reorder p-value data to match dat_perm_mean
  dat_p = dat_p[order(dat_p$ID),]
  
  # Check that the ordering of ID values for dat_p and dat_perm_mean are the same
  match = dat_p$ID == dat_perm_mean$ID
  table(match)["FALSE"]
  table(match)["TRUE"]
  
  # Add significance values from dat_p to dat_perm_mean
  dat_perm_mean$sig = dat_p$sig
  
  dat_perm_mean = filter(dat_perm_mean, FirstLabel %in% select_first &
                           SecondLabel %in% select_second)
  
  print("Data prep complete")
  #return(dat_perm_mean)
  
  # Call log2FC function
  log2FC(community, dat_perm_mean, select_second, extra, out)
  
}

# Calculate log2FC of enrichment and plot as dot plot 
log2FC = function(community, data, select_second, extra, out) {
  # Calculate log2 fold change between baseline and permutation ct values
  data$log2 = log2(data$ct_baseline / data$ct)
  print(length(unique(data$log2)))
  
  # Save data with calculations of difference between baseline and permutation ct values
  write.csv(
    data,
    paste0(
      out_stats5,
      date,
      "_neighbouRhood_dat_perm_mean_log2FC_",
      community,
      "_MRTX",
      extra,
      ".csv"
    ),
    row.names = F
  )
  
  for (i in 1:length(unique(data$FirstLabel))) {
    select_metacluster = data[which(data$FirstLabel == unique(data$FirstLabel)[i]),]
    print(unique(data$FirstLabel)[i])
    
    p = ggplot(select_metacluster,
               aes(
                 x = SecondLabel,
                 y = log2,
               )) +
      geom_hline(yintercept = 0) +
      geom_dotplot(
        data = select_metacluster[which(select_metacluster$sig == FALSE),],
        binaxis = 'y',
        stackdir = 'center',
        dotsize = 1.5,
        fill = "white"
      ) +
      geom_dotplot(
        data = select_metacluster[which(select_metacluster$sig == TRUE),],
        binaxis = 'y',
        stackdir = 'center',
        dotsize = 1.5
      ) +
      stat_summary(fun.data=mean_se, geom="errorbar", color="orange", width = 0.5) + 
      stat_summary(fun=mean, geom="point", color="orange", size = 0.5) + 
      theme_minimal() +
      labs(title = unique(data$FirstLabel)[i],
           x = "",
           y = "Log2FC enrichment") +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 14
        ),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14)
      ) +
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
      scale_x_discrete(limits = select_second
      )
    
    filename = paste0(
      unique(data$FirstLabel)[i],
      "_neighbours_pointPlot_log2significance_",
      community, 
      "_MRTX",
      "_",
      extra,
      ".png"
    )
    ggsave(
      plot = p,
      device = "png",
      width = 6,
      height = 4.5,
      dpi = 300,
      bg = 'White',
      path = out,
      filename = filename
    )
  }
  print("log2FC plotting complete")
  
  return(select_metacluster)
}

# T reg cells count of each community
treg_top5_count = function(data, name) {
  
  p = data %>%
    filter(cellType == "T reg cells",
           top5 %in% c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2")) %>% 
    select(cellType, top5, treatment) %>% 
    ggplot(aes(x = factor(treatment, levels = c("Vehicle", "MRTX")), fill = factor(top5, levels = c("T/NA", "T/M1", "T/DC", "T/M2_1", "T/M2_2")))) +
    geom_bar(stat = "count", colour = "Black") + 
    theme_classic() + 
    scale_fill_manual(values = top5_col[1:5]) + 
    xlab("") +
    ylab("Count of all T reg cells") +
    labs(fill = "Community") + 
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18))
  ggsave(plot = p, device = "png", width=5, height=5, dpi=300, path = out_supp5,
         filename = paste0(date, name))
}

#########################################################################################################

###################
#### Figure 5a ####

# Relationship of CXCL9 expression on DCs to distance to a CD8+ T cell 
# *NOTE: distance calculations carried out in Python* 

data_CD8 = read.csv(paste0(distance_path, "20230901_DCCXCL9_high_low_min_distance_to_CD8TC_PD1high_MRTX.csv"))
data_CD8$source_cluster = gsub("_", " ", data_CD8$source_cluster)
compare_dist(data_CD8, "T cells CD8", out_fig5, "_dataset2_comC_DC_CXCL9_high_low_minDist_CD8_PD1high.png")


################################
#### Supplementary Figure 5a ###

data_CD4 = read.csv(paste0(distance_path, "20230901_DCCXCL9_high_low_min_distance_to_CD4TC_PD1high_MRTX.csv"))
data_CD4$source_cluster = gsub("_", " ", data_CD4$source_cluster)
compare_dist(data_CD4, "T cells CD4", out_supp5, "_dataset2_comC_DC_CXCL9_high_low_minDist_CD4_PD1high.png")


################################
#### Supplementary Figure 5b ###

data_Treg = read.csv(paste0(distance_path, "20230901_DCCXCL9_high_low_min_distance_to_regTC_PD1high_MRTX.csv"))
data_Treg$source_cluster = gsub("_", " ", data_Treg$source_cluster)
compare_dist(data_Treg, "T reg cells", out_supp5, "_dataset2_comC_DC_CXCL9_high_low_minDist_Treg_PD1high.png")


###################
#### Figure 5b ####

# Visual of community C - CXCL9+ DCs interacting with PD-1+ CD8 TCs

images = c("BRAC3326.4e_ROI1_Vehicle", "BRAC4002.3c_ROI2_MRTX", "BRAC4002.3c_ROI3_MRTX", "BRAC3438.6f_ROI1_Vehicle",
           "BRAC3438.6f_ROI3_Vehicle","BRAC3495.3f_ROI1_Vehicle", "BRAC3529.2a_ROI1_MRTX", "BRAC3529.2b_ROI1_MRTX",
           "BRAC3529.2d_ROI3_MRTX", "BRAC3708.2d_ROI1_Vehicle", "BRAC4002.3c_ROI1_MRTX")

subset = neighb2 %>% 
  filter(top5 == "T/DC") %>% 
  mutate(CXCL9 = "Other") %>% 
  mutate(CXCL9 = case_when(cellType == "Dendritic cells" & MI_CXCL9 >= 0.5 ~ "Dendritic cells CXCL9 high",
                           cellType == "cDC1" & MI_CXCL9 >= 0.5 ~ "Dendritic cells CD103 CXCL9 high",
                           cellType == "T cells CD8" & MI_PD1 >= 0.5 ~ "T cells CD8 PD-1 high", 
                           TRUE ~ CXCL9))

out_dir = paste0(out_fig5, "Dataset2_CommunityC_DC_CXCL9_CD8_PD1/")
dir.create(out_dir)

nbclusters = unique(subset$CXCL9)

labels = data.frame(
  clustername = sort(unique(subset$CXCL9)),
  colours = c("#FFCC66", "#FF9933", "#E1E1E1", "#FF9D9AFF")
)

for(ROI in images){
  print(ROI)
  p = plot_cellmap(ROI,
                   "TCpanel",
                   "ROI_name",
                   subset,
                   "CXCL9",
                   nbclusters,
                   out_dir,
                   c("T cells CD8 PD-1 high", "Dendritic cells CXCL9 high", "Dendritic cells CD103 CXCL9 high"))
}


###################
#### Figure 5c #### 

# Ki67 expression 
marker_exp("MI_Ki67", "Ki67", "_dataset2_CD4TC_CD8TC_Treg_Ki67expression.png")


###################
#### Figure 5d ####

seg_neighb2$First.Object.top5 = neighb2$top5[match(seg_neighb2$cellID, neighb2$cellID)]
seg_neighb2$Second.Object.top5 = neighb2$top5[match(seg_neighb2$cellID2, neighb2$cellID)]
write.csv(seg_neighb2, paste0(path, date, "_seg_neighb2.csv"), row.names = F)

# Interaction of CD8+ T cells with casp3+ tumour cells Vehicle vs MRTX - dataset 2 - community C 
close_neighbours(neighb2, seg_neighb2, "T cells CD8", "T/DC", "_dataset2_barchart_communityC_CD8TC_casp3Tumour_interactions_mean_ROI_dots.png", out_fig5)


#################################
#### Supplementary Figure 5c ####

tumour_com = c("T/DC", "T/M2_1", "T/M2_2")

close_neighbours_tumour(neighb1, seg_neighb1, tumour_com, "_dataset1_comCDE_CD8TC_casp3Tumour_closeNeighbours_MRTX_perROI_notRelative.png", FALSE)
close_neighbours_tumour(neighb2, seg_neighb2, tumour_com, "_dataset2_comCDE_CD8TC_casp3Tumour_closeNeighbours_MRTX_perROI_notRelative.png", FALSE)


#################################
#### Supplementary Figure 5d ####

close_neighbours(neighb2, seg_neighb2, "T cells CD4", "T/DC", "_dataset2_barchart_communityC_CD4TC_casp3Tumour_interactions_mean_ROI_mean_ROI_dots.png", out_supp5, 4.5)

close_neighbours(neighb2, seg_neighb2, "T reg cells", "T/DC", "_dataset2_barchart_communityC_Treg_casp3Tumour_interactions_mean_ROI_mean_ROI_dots.png", out_supp5, 4.5)


###################################
#### Figure 5e & Supp Fig 5e-h ####

## NeighbouRhood enrichment analysis of cells in the neighbourhood of CD8 T cells ##
## MRTX only ##

# Create output folders for output statistics 
out_stats5 = paste0(file5, "out_stats/")
dir.create(out_stats5)

## Create ObjectTable 
mrtx = filter(neighb2, treatment == "MRTX")
mrtx_order = mrtx %>%
  select(ObjectNumber, ImageNumber, cellType) %>%
  dplyr::rename(label = cellType)

## Create Object relationships
rel = seg_neighb2 %>% 
  select(First.Image.Number, Second.Image.Number, First.Object.Number, Second.Object.Number,
         cellID, First.Object.top5, Second.Object.top5, cellID2, treatment) %>% 
  filter(treatment == "MRTX") %>% 
  mutate(Module.Number = 139,
         First.Object.Name = "AllCellsMask",
         Second.Object.Name = First.Object.Name,
         Module = "MeasureObjectNeighbors",
         Relationship = "Neighbors") %>% 
  select(Module, Module.Number, Relationship, First.Object.Name, First.Image.Number, First.Object.Number, Second.Object.Name,
         Second.Image.Number, Second.Object.Number, cellID, cellID2, First.Object.top5, Second.Object.top5)

# Prep object relationships 
top5_rel = object_rel_top5(rel, mrtx)

# Figure 5f - Community TDC 
select_metacluster = neighbourhood_permutation(
  mrtx_order,
  top5_rel$tdc,
  "TDC_Community",
  "T cells CD8",
  c(
    "T cells CD4",
    "T reg cells",
    "T cells CD8",
    "Dendritic cells CD103",
    "Dendritic cells",
    "B cells",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells DN",
    "Leukocytes unclassified",
    "Fibroblasts",
    "Tumour"
  ),
  "_01sig_meanse",
  out_fig5
)


# Supp Fig 5e - Community TNA
select_metacluster = neighbourhood_permutation(
  mrtx_order,
  top5_rel$tna,
  "TNA_Community",
  "T cells CD8",
  c(
    "T cells CD4",
    "T reg cells",
    "T cells CD8",
    "Dendritic cells CD103",
    "Dendritic cells",
    "B cells",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells DN",
    "Leukocytes unclassified",
    "Fibroblasts",
    "Tumour"
  ),
  "_01sig_meanse",
  out_supp5
)

# Supp Fig 5f - Community TM1
select_metacluster = neighbourhood_permutation(
  mrtx_order,
  top5_rel$tm1,
  "TM1_Community",
  "T cells CD8",
  c(
    "T cells CD4",
    "T reg cells",
    "T cells CD8",
    "Dendritic cells CD103",
    "Dendritic cells",
    "B cells",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells DN",
    "Leukocytes unclassified",
    "Fibroblasts",
    "Tumour"
  ),
  "_01sig_meanse",
  out_supp5
)

# Supp Fig 5g - Community TM2_1
select_metacluster = neighbourhood_permutation(
  mrtx_order,
  top5_rel$tm2_1,
  "TM2_1_Community",
  "T cells CD8",
  c(
    "T cells CD4",
    "T reg cells",
    "T cells CD8",
    "Dendritic cells CD103",
    "Dendritic cells",
    "B cells",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells DN",
    "Leukocytes unclassified",
    "Fibroblasts",
    "Tumour"
  ),
  "_01sig_meanse",
  out_supp5
)

# Supp Fig 5h - Community TM2_2
select_metacluster = neighbourhood_permutation(
  mrtx_order,
  top5_rel$tm2_2,
  "TM2_2_Community",
  "T cells CD8",
  c(
    "T cells CD4",
    "T reg cells",
    "T cells CD8",
    "Dendritic cells CD103",
    "Dendritic cells",
    "B cells",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells DN",
    "Leukocytes unclassified",
    "Fibroblasts",
    "Tumour"
  ),
  "_01sig_meanse",
  out_supp5
)



###################
#### Figure 5f ####

# Community 14/C and 2/E have the highest proportion of Treg cells 
heat2 = neighb2 %>% 
  filter(treatment == "MRTX") %>% 
  select(cellType, agglom18_average) %>% 
  group_by(agglom18_average) %>%
  mutate(count = n()) %>% 
  summarise(T_cells_CD8 = sum(cellType == "T cells CD8")/unique(count)*100,
            T_cells_CD4 = sum(cellType == "T cells CD4")/unique(count)*100,
            T_reg_cells = sum(cellType == "T reg cells")/unique(count)*100) 

pdf(paste0(out_fig5, date, "_dataset2_heatmap_TC_proportion_per_community_MRTX.pdf", sep = ""), width = 3, height = 6)
heatmap.2(as.matrix(heat2[,-1]), trace = "none", scale = "col", dendrogram = "none", Colv = FALSE, Rowv = FALSE, margins = c(2,6),
          density.info = "none", keysize = 1, col = diverge_hsv(10))
dev.off()


#################################
#### Supplementary Figure 5i ####

heat1 = neighb1 %>% 
  filter(treatment == "MRTX") %>% 
  select(cellType, agglom18_average) %>% 
  group_by(agglom18_average) %>%
  mutate(count = n()) %>% 
  summarise(T_cells_CD8 = sum(cellType == "T cells CD8")/unique(count)*100,
            T_cells_CD4 = sum(cellType == "T cells CD4")/unique(count)*100,
            T_reg_cells = sum(cellType == "T reg cells")/unique(count)*100) 

pdf(paste0(out_supp5, date, "_dataset1_heatmap_TC_proportion_per_community_MRTX.pdf", sep = ""), width = 3, height = 6)
heatmap.2(as.matrix(heat1[,-1]), trace = "none", scale = "col", dendrogram = "none", Colv = FALSE, Rowv = FALSE, margins = c(2,6),
          density.info = "none", keysize = 1, col = diverge_hsv(10))
dev.off()


###################
#### Figure 5g ####

# Close interactions of Tregs with DCs and CD8 TCs in community C - MRTX
subset = neighb2 %>%
  filter(top5 == "T/DC") %>%
  mutate(immune = cellType,
         immune = case_when(!(
           immune %in% c(
             "Dendritic cells",
             "Dendritic cells CD103",
             "T cells CD4",
             "T cells CD8",
             "T reg cells"
           )
         ) ~ "Other", TRUE ~ immune)) 

nbclusters = c(sort(unique(subset$immune)))

dir = paste0(out_fig5, "Dataset2_CommunityC_DCs_TCs/")
dir.create(dir)

labels = data.frame(
  clustername = sort(unique(subset$immune)),
  colours = c(
    "#FF9933",
    "#FFCC66",
    "#E1E1E1",
    "#B07AA1FF",
    "#FF9D9AFF",
    "#CC6666"
  )
)

for(ROI in images){
  print(ROI)
  p = plot_cellmap(ROI,
                   "TCpanel",
                   "ROI_name",
                   subset,
                   "immune",
                   nbclusters,
                   dir,
                   c("Dendritic cells", "Dendritic cells CD103", "T cells CD4", "T cells CD8", "T reg cells"))
}


#################################
#### Supplementary Figure 5j ####

# Count of T reg cells in the top 5 communities, split by treatment
treg_top5_count(neighb1, "_dataset1_Treg_count_in_top5_perTreatment.png")
treg_top5_count(neighb2, "_dataset2_Treg_count_in_top5_perTreatment_stat.png")





