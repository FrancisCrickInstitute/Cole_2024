
#### Scripts for Treg paper - Figure 6 & Supplementary figure 6 ####

## Megan Cole (megan.cole@crick.ac.uk)

rm(list=ls())

# Load libraries 
library(ggplot2) 
library(dplyr) 
library(ggpubr) 
library(data.table) 
library(parallel) 
library(neighbouRhood) 

date = format(Sys.Date(), "%Y%m%d")

path = "/Users/colem/Documents/Projects/Treg_paper/Golden/Data/"

neighb2 = read.csv(paste0(path, "20231002_neighb2.csv"))
seg_neighb2 = read.csv(paste0(path, "20230904_seg_neighb2.csv"))

out_fig6 = "/Users/colem/Documents/Projects/Treg_paper/Figures/Figure6/"
dir.create(out_fig6)
out_supp6 = "/Users/colem/Documents/Projects/Treg_paper/Figures/Supplementary_Figure6/"
dir.create(out_supp6)

treg_col = c("Tregs" = "#9590FF", "No Tregs" = "#00C1A3")

files6 = paste0(path, "Figure6/")
dir.create(fig6)
fig6_OT = paste0(files6, "objectTable/")
dir.create(fig6_OT)

out_stats6 = "/Users/colem/Documents/Projects/Treg_paper/Golden/Data/Figure6/out_stats/"

colours2 = c("B cells" = "#945931", "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Fibroblasts" = "#ABDDA4FF",
             "Leukocytes unclassified" = "#7A9F79", "Macrophages type 1" = "#336666", "Macrophages type 2" = "#4E79A7FF",
             "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF", "T cells CD8" = "#FF9D9AFF", "T cells DN" = "#FFB5CB",
             "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")

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

################################################################################

###################
#### Functions ####

# Create object relationships in preparation for neughbouRhood permutaiton 
create_obj_rel = function(com, seg_com){
  
  seg_com = seg_com %>% filter(treatment == "MRTX")
  seg_com$treg_pres = com$treg_pres[match(seg_com$cellID, com$cellID)]
  
  seg_com$treg_pres <- factor(seg_com$treg_pres, levels = c("Tregs", "No Tregs"))
  
  rel = seg_com %>% 
    select(First.Image.Number, Second.Image.Number, First.Object.Number, Second.Object.Number,
           cellID, First.Object.top5, Second.Object.top5, cellID2, treg_pres) %>% 
    mutate(Module.Number = 139,
           First.Object.Name = "AllCellsMask", 
           Second.Object.Name = First.Object.Name,
           Module = "MeasureObjectNeighbors",
           Relationship = "Neighbors") %>% 
    select(Module, Module.Number, Relationship, First.Object.Name, First.Image.Number, First.Object.Number, Second.Object.Name,
           Second.Image.Number, Second.Object.Number, cellID, cellID2, First.Object.top5, Second.Object.top5, treg_pres) %>% 
    group_split(treg_pres) %>% 
    setNames(levels(seg_com$treg_pres))
  
  rel$Tregs$treg_pres = as.character(rel$Tregs$treg_pres)
  rel$`No Tregs`$treg_pres = as.character(rel$`No Tregs`$treg_pres)
  
  return(rel)
}

# Run neighbourhood permutation 
neighbourhood_permutation = function(cdata, ndata, file_path, community, vs, select_first, select_second, plot_path, sig, extra, limits){
  
  # Convert object table to data.table format
  cdata = data.table(cdata)
  
  names(ndata) = gsub(".Object.Name", " Object Name", names(ndata))
  names(ndata) = gsub(".Image.Number", " Image Number", names(ndata))
  names(ndata) = gsub(".Object.Number", " Object Number", names(ndata))
  names(ndata) = gsub("Module.Number", "Module Number", names(ndata))
  ndata = data.table(ndata)
  
  ndata$`First Image Number` = as.integer(ndata$`First Image Number`)
  ndata$`Second Image Number` = as.integer(ndata$`Second Image Number`)
  
  data = prepare_tables(cdata, ndata, objname = "AllCellsMask")
  
  data_baseline = apply_labels(data[[1]], data[[2]]) %>%
    aggregate_histo()
  write.csv(data_baseline, paste0(file_path, date, "_dataset2_neighbouRhood_data_baseline_", community, "_v_", vs, "_MRTX", extra, ".csv"), row.names = F)
  
  n_perm = 1000
  
  # Calculate permutation statistics
  set.seed(12312)
  dat_perm = rbindlist(mclapply(1:n_perm, function(x) {
    dat_labels = shuffle_labels(data[[1]])
    apply_labels(dat_labels, data[[2]]) %>%
      aggregate_histo()
  }, mc.cores = 10), idcol = 'run')
  write.csv(dat_perm, paste0(file_path, date, "_dataset2_neighbouRhood_dat_perm_", community, "_v_",
                             vs, "_MRTX", extra, "2.csv"), row.names = F)
  
  # Calculate p-values
  dat_p <-
    calc_p_vals(data_baseline,
                dat_perm,
                n_perm = n_perm,
                p_tresh = sig)
  write.csv(dat_p, paste0(file_path, date, "_dataset2_neighbouRhood_dat_pvalues_", community, "_v_", 
                          vs, "_MRTX", extra, ".csv"), row.names = F)
  
  stats = list("data_baseline" = data_baseline, "dat_perm" = dat_perm, "dat_p" = dat_p)
  
  print("Permutation complete")
  #return(stats)
  
  # Call data prep
  data_prep(stats$data_baseline, stats$dat_perm, stats$dat_p, select_first, select_second, community, vs, file_path, plot_path, sig, extra, limits)
}

# Prep neighbouRhood permutation results for Log2FC plotting 
data_prep = function(data_baseline, dat_perm, dat_p, select_first, select_second, community, vs, file_path, plot_path, sig, extra, limits){
  
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
      difference2 = as.data.frame(setdiff(data_baseline[, c("ID", "group", "FirstLabel", "SecondLabel")], dat_perm_mean[, c("ID", "group", "FirstLabel", "SecondLabel")]))
      print(difference2$ID)
      difference2$ct = 0
      difference2 = difference2[, c("ID",
                                    "group",
                                    "FirstLabel",
                                    "SecondLabel",
                                    "ct")]
      data_perm_mean = rbind(dat_perm_mean, difference)
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
  
  # Order baseline data by ID, so it matches mean permutation data
  data_baseline = data_baseline[order(data_baseline$ID),]
  
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
  log2FC(community, vs, file_path, plot_path, dat_perm_mean, sig, extra, limits)
}

# Plot Log2 fold changes in enrichment from neighbouRhood analysis 
log2FC = function(community, vs, file_path, plot_path, data, sig, extra, limits) {
  # Calculate log2 fold change between baseline and permutation ct values
  data$log2 = log2(data$ct_baseline / data$ct)
  
  # Save data with calculations of difference between baseline and permutation ct values
  write.csv(
    data,
    paste0(
      file_path,
      date,
      "_neighbouRhood_dat_perm_mean_log2FC_",
      community,
      "_",
      vs,
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
      scale_x_discrete(
        limits = limits
      ) + 
      ylim(-1,1)
    
    filename = paste0(
      unique(data$FirstLabel)[i],
      "_neighbours_pointPlot_log2significance_",
      community, 
      "_v_",
      vs,
      "_MRTX",
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
      path = plot_path,
      filename = filename
    )
  }
  print("log2FC plotting complete")
  
  return(select_metacluster)
}

# Close neighbour relationships between CD8+ T cells and casp3+ tumour cells 
close_neighbours = function(n, neighbourdata, ct, name, path) {
  
  # How often are CD8 TC and casp3+ Tumour cells neighbouring each other?
  tum = n[which(n$cellType == "Tumour" & n$MI_casp3 >= 0.5),]
  cd = n[which(n$cellType == ct),]
  
  # Neighbours of casp3+ tumour & CD8 T cells
  tum_neighbours = neighbourdata[which(neighbourdata$cellID %in% tum$cellID),]
  cd_neighbours = neighbourdata[which(neighbourdata$cellID %in% cd$cellID),]
  
  # Data for CD8+ T cells that have a casp3+ tumour cell in their neighbourhood
  cd_casp3Tumour_neighbours = cd_neighbours[which(cd_neighbours$cellID2 %in% tum$cellID),]
  
  ## Separate by ROI and treatment - only including ROIs that have casp3+ tumour cells 
  paired_count = cd_casp3Tumour_neighbours %>%
    filter(treatment == "MRTX") %>% 
    select(First.Object.treg_pres, ROI_short) %>%
    dplyr::group_by(First.Object.treg_pres) %>%
    dplyr::mutate(count=n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(First.Object.treg_pres, ROI_short) %>% 
    dplyr::summarize(count = unique(count),
                     ROI_count=n())
  
  relative_count(paired_count, n, name, path)
}

relative_count = function(paired_count, n, name, path, width="6", height="6"){
  
  casp3_com_treat_perc = n %>%
    dplyr::select(cellType, MI_casp3, treg_pres) %>%
    dplyr::group_by(treg_pres) %>%
    # Split % of casp3+ tumour cells into treatment
    dplyr:: mutate(perc_community_treat = (length(which(cellType == "Tumour" & MI_casp3 >= 0.5))/length(which(cellType == "Tumour"))*100)) %>%
    dplyr:: select(MI_casp3, perc_community_treat) %>%
    dplyr:: summarise_all(list(mean)) %>%
    # Remove rows with NaN values
    dplyr:: filter_all(any_vars(!is.na(perc_community_treat))) %>% 
    dplyr:: filter_all(any_vars(perc_community_treat > 0)) 
  
  # Add % of casp3+ tumour cells for each treatment group
  paired_count$perc_community_treat = casp3_com_treat_perc$perc_community_treat[match(paired_count$First.Object.treg_pres, casp3_com_treat_perc$treg_pres)]
  # Calculate larger % by smaller % 
  perc_diff = max(paired_count$perc_community_treat)/min(paired_count$perc_community_treat)
  
  # Find the minimum % value from calculating % of tumour cells that are casp3+
  min_perc = min(paired_count$perc_community_treat)
  # Take ROI casp3+ count values associated with that minimum % 
  ROI_count_min_perc = paired_count$ROI_count[which(paired_count$perc_community_treat == min_perc)]
  # Times those ROI count values by perc_diff and all others by 1 
  paired_count$relative_ROI_count = ifelse(paired_count$ROI_count %in% ROI_count_min_perc,
                                           paired_count$ROI_count*perc_diff, paired_count$ROI_count*1)
  paired_count[is.na(paired_count)] <- 0
  names(paired_count)[names(paired_count) == 'First.Object.treg_pres'] <- 'treg_pres'
  
  # # Add 0 column per ROI and remove if there is a larger value
  df = n %>% 
    select(ROI_short, treg_pres, treatment) %>% 
    filter(treatment == "MRTX") %>% 
    group_by(ROI_short, treg_pres) %>% 
    dplyr::summarise(count = n()) %>% 
    dplyr::mutate(relative_ROI_count = 0) %>% 
    select(treg_pres, ROI_short, relative_ROI_count)
  
  # Add paired_count and df together and take larger value 
  paired_count = paired_count %>% ungroup %>% select(treg_pres, ROI_short, relative_ROI_count)
  paired_count = rbind(paired_count, df)
  paired_count = paired_count %>%
    dplyr::group_by(ROI_short, treg_pres) %>%
    dplyr::summarise(treg_pres = unique(treg_pres),
                     ROI_short = unique(ROI_short),
                     relative_ROI_count = max(relative_ROI_count))
  
  paired_count_mean = paired_count %>% 
    group_by(treg_pres) %>% 
    dplyr::summarise(treg_pres = unique(treg_pres),
                     mean_count = mean(relative_ROI_count))
  
  paired_count$ROI_short = factor(paired_count$ROI_short, levels = c("BRAC3529.2a_MRTX", "BRAC3529.2b_MRTX", "BRAC3529.2d_MRTX", 
                                                                     "BRAC4002.3c_t1_MRTX", "BRAC4002.3c_t2_MRTX", "BRAC4002.3c_t3_MRTX"
  ))
  
  paired_count$treg_pres = factor(paired_count$treg_pres, levels = c("Tregs", "No Tregs"))
  
  p = ggplot(transform(paired_count, treg_pres = factor(treg_pres, levels = c("Tregs", "No Tregs")))) + 
    geom_col(data = transform(paired_count_mean, treg_pres = factor(treg_pres, levels = c("Tregs", "No Tregs"))),
             aes(x = factor(treg_pres), y = mean_count, fill = factor(treg_pres)),
             position = position_dodge2(width = 0, preserve = "single"), colour = "Black", fill = treg_col) + 
    geom_dotplot(binaxis='y', stackdir='center', stackgroups = TRUE, binpositions="all", dotsize = 0.8,
                 #  colour = "Black", aes(x = factor(treg_pres), y = relative_ROI_count, fill = factor(ROI_short))) +
                 colour = "Black", fill = "Black", aes(x = factor(treg_pres), y = relative_ROI_count)) +
    theme_classic() + 
    xlab("") + 
    ylab("Relative count") + 
    theme(axis.title = element_text(size = 22), 
          axis.text = element_text(size = 22),
          legend.title = element_blank(),
          legend.text = element_text(size = 18))
  # Save plot
  ggsave(plot = p, device = "png", width=6, height=6, dpi=300, path = path,
         filename = paste0(date, name))
}


################################################################################

#####################
##### Figure 6a #####

## Separation of community C with and without Tregs 
comc = filter(neighb2, top5 == "T/DC")
comc$treg_pres = ifelse(comc$T_reg_cells == 0, "No Tregs", "Tregs")
comc$treg_pres <- factor(comc$treg_pres, levels = c("Tregs", "No Tregs"))

## Proportion of cell types in Treg vs no Treg neighbourhoods within community C
p = ggplot(filter(comc, treatment == "MRTX"), aes(x = factor(treg_pres, levels = c("No Tregs", "Tregs")), fill = as.factor(cellType)))+
  scale_fill_manual(values = colours2) +
  geom_bar(position = "fill") +
  theme_classic()+
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text=element_text(size=22)) +
  theme(legend.text = element_text(size = 22)) +
  xlab("")+
  ylab("Percentage distribution") +
  
  labs( fill = "")
p
ggsave(plot = p, device = "png", width=11, height=6, dpi=300, path = out_fig6,
       filename = paste(date, "_ComC_cellType_distribution_with_without_Tregs_MRTX.png", sep = ""))



#####################
#### Supp Fig 6a ####

## Percentage of CD8 T cells in Treg / No Treg neighbourhoods - community C 

comc_cd8_prop = comc %>% 
  filter(treatment == "MRTX") %>% 
  group_by(treg_pres, ROI_name) %>% 
  dplyr::mutate(count = n(),
                cd8_sum = sum(cellType == "T cells CD8")) %>% 
  ungroup() %>% 
  group_by(treg_pres, ROI_name) %>% 
  dplyr::summarise(treg_pres = unique(treg_pres),
                   ROI_name = unique(ROI_name),
                   count = unique(count),
                   cd8_sum = unique(cd8_sum),
                   cd8_prop = (cd8_sum/count)*100) 

comc_cd8_prop$treg_pres = factor(comc_cd8_prop$treg_pres, levels = c("Tregs", "No Tregs"))

p = ggbarplot(comc_cd8_prop, x = "treg_pres", y = "cd8_prop", fill = "treg_pres", levels = c("Tregs", "No Tregs"),
              add = "mean_se", width = .9) + 
  theme_classic()+
  scale_fill_manual(values = treg_col) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
        legend.position = "none") + 
  ylab("Percentage CD8 T cells") + 
  labs(fill = "") + 
  ylim(0,10)
p
ggsave(plot = p, device = "png", width=4.5, height=5, dpi=300, path = out_supp6,
       filename = paste(date, "_dataset2_treg_notreg_CD8_percentage.png", sep = ""))


#################################
#### Supplementary Figure 6b #### 

comc_treg_datp = read.csv(paste0(out_stats6, "20230821_dataset2_neighbouRhood_dat_pvalues_CommunityTDC_v_CommunityTDC_MRTX_01sig_treg_dataset2.csv"))
comc_notreg_datp = read.csv(paste0(out_stats6, "20230821_dataset2_neighbouRhood_dat_pvalues_CommunityTDC_v_CommunityTDC_MRTX_01sig_notreg_dataset2.csv"))

pmat_treg = reshape2::dcast(comc_treg_datp, 'FirstLabel ~ SecondLabel', value.var = 'sigval', fun.aggregate = sum,
                            fill=0, drop=F) 

pmat_notreg = reshape2::dcast(comc_notreg_datp, 'FirstLabel ~ SecondLabel', value.var = 'sigval', fun.aggregate = sum,
                              fill=0, drop=F)

# Make FirstLabel names = row names & remove FirstLabel column
rname_treg = pmat_treg$FirstLabel
rname_notreg = pmat_notreg$FirstLabel

pmat_treg = pmat_treg %>%
  select(-c('FirstLabel')) %>%
  as.matrix()
pmat_notreg = pmat_notreg %>%
  select(-c('FirstLabel')) %>%
  as.matrix()

row.names(pmat_treg) <- rname_treg
row.names(pmat_notreg) <- rname_notreg

# Select clusters for plotting
select = c("Dendritic cells", "Dendritic cells CD103", "Fibroblasts", "Macrophages type 1", 
           "Macrophages type 2", "T cells CD4", "T cells CD8", "Tumour")
pmat_treg = pmat_treg[select, select]
pmat_notreg = pmat_notreg[select, select]

rownames(pmat_notreg) = paste(rownames(pmat_notreg), "No Treg", sep = "-")

pmat = rbind(pmat_treg, pmat_notreg)
# Order the rows alphabetically 

pmat = pmat[ order(rownames(pmat)) , ,drop=F]
pmat = pmat[c(1,4,2,3,5:16),]

png(paste0(out_supp6, date, "_TDC_community_MRTX_Tregs_and_noTregs.png"),width=6,height=6,units="in",res=1200)
hr <- hclust(dist(pmat), method="ward.D")
heatmap.2(pmat,
          dendrogram = "none",
          Colv = FALSE,
          Rowv = FALSE,
          trace = "none",
          margins = c(15,15),
          col = diverge_hsv(11),
          density.info ='none',
          keysize = 1,
          breaks = seq(-6,5, 1) 
)
dev.off()


#################################
#### Figure 6b / Supp Fig 6c ####

## Comparing neighbourhood relationships in presence vs absence of Treg cells 

# Create object table 
comc_order = comc %>%
  select(ObjectNumber, ImageNumber, cellType, cellID) %>%
  dplyr::rename(label = cellType)

seg_comc = seg_neighb2 %>% 
  filter(cellID %in% comc$cellID)

seg_comc$First.Object.treg_pres = comc$treg_pres[match(seg_comc$cellID, comc$cellID)]
seg_comc$First.Object.treg_pres = as.character(seg_comc$First.Object.treg_pres)
seg_comc$Second.Object.treg_pres = comc$treg_pres[match(seg_comc$cellID2, comc$cellID)]
seg_comc$Second.Object.treg_pres = as.character(seg_comc$Second.Object.treg_pres)
seg_comc$Second.Object.treg_pres[is.na(seg_comc$Second.Object.treg_pres)] <- "Other"

# Create object relationships 
obj_relc = create_obj_rel(comc, seg_comc)


# Plot neighbourhood relationships 
# Community C - Treg 
select_metacluster = neighbourhood_permutation(
  comc_order[,-which(names(comc_order) == "cellID")], 
  obj_relc$Tregs,
  out_stats6, 
  "TDC_Community", 
  "TDC_Community", 
  c(
    "T cells CD4",
    "T cells CD8"
  ),
  c(
    "B cells",
    "Dendritic cells CD103",
    "Dendritic cells",
    "Fibroblasts",
    "Leukocytes unclassified",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells CD4",
    "T cells CD8",
    "T cells DN",
    "T reg cells", 
    "Tumour"
  ),
  out_fig6, 
  0.01, 
  "_01sig_treg_dataset2",
  c(
    "B cells",
    "Dendritic cells CD103",
    "Dendritic cells",
    "Fibroblasts",
    "Leukocytes unclassified",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells CD4",
    "T cells CD8",
    "T cells DN",
    "T reg cells", 
    "Tumour"
  )
)


# Community C - No Treg 
select_metacluster = neighbourhood_permutation(
  comc_order[,-which(names(comc_order) == "cellID")], 
  obj_relc$`No Tregs`,
  out_stats6, 
  "TDC_Community", 
  "TDC_Community", 
  c(
    "T cells CD4",
    "T cells CD8"
  ),
  c(
    "B cells",
    "Dendritic cells CD103",
    "Dendritic cells",
    "Fibroblasts",
    "Leukocytes unclassified",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells CD4",
    "T cells CD8",
    "T cells DN",
    "T reg cells", 
    "Tumour"
  ),
  out_fig6, 
  0.01, 
  "_01sig_notreg_dataset2", 
  c(
    "B cells",
    "Dendritic cells CD103",
    "Dendritic cells",
    "Fibroblasts",
    "Leukocytes unclassified",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells CD4",
    "T cells CD8",
    "T cells DN",
    "T reg cells", 
    "Tumour"
  )
)


#################################
#### Figure 6c & Supp fig 6d ####

# Figure 6c 
close_neighbours(comc, seg_comc, "T cells CD8", "_dataset2_barchart_TDC_community_TregNoTreg_CD8TC_casp3Tumour_interactions_mean_ROI_points.png", out_fig6)

# Supp Fig 6d 
close_neighbours(comc, seg_comc, "T cells CD4", "_dataset2_barchart_TDC_community_TregNoTreg_CD4TC_casp3Tumour_interactions_mean_ROI_points.png", out_supp6)



#########################
#### Supp Fig 6e & f ####

## Create datasets related to community T/M2_2 only  ##
comE = filter(neighb2, top5 == "T/M2_2")
comE$treg_pres = ifelse(comE$T_reg_cells == 0, "No Tregs", "Tregs")
comE$treg_pres <- factor(comE$treg_pres, levels = c("Tregs", "No Tregs"))

seg_comE = seg_neighb2 %>% 
  filter(cellID %in% comE$cellID)

seg_comE$First.Object.treg_pres = comE$treg_pres[match(seg_comE$cellID, comE$cellID)]
seg_comE$First.Object.treg_pres = as.character(seg_comE$First.Object.treg_pres)
seg_comE$Second.Object.treg_pres = comE$treg_pres[match(seg_comE$cellID2, comE$cellID)]
seg_comE$Second.Object.treg_pres = as.character(seg_comE$Second.Object.treg_pres)
seg_comE$Second.Object.treg_pres[is.na(seg_comE$Second.Object.treg_pres)] <- "Other"

# Create object table 
comE_order = comE %>%
  select(ObjectNumber, ImageNumber, cellType, cellID) %>%
  dplyr::rename(label = cellType)

# Create object relationships 
obj_relE = create_obj_rel(comE, seg_comE)

# Plot neighbourhood relationships 
# Community E - Treg 
select_metacluster = neighbourhood_permutation(
  comE_order[,-which(names(comE_order) == "cellID")], 
  obj_relE$Tregs,
  out_stats6, 
  "TM2_2_Community", 
  "TM2_2_Community", 
  c(
    "T cells CD4",
    "T cells CD8"
  ),
  c(
    "B cells",
    "Dendritic cells CD103",
    "Dendritic cells",
    "Fibroblasts",
    "Leukocytes unclassified",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells CD4",
    "T cells CD8",
    "T cells DN",
    "T reg cells", 
    "Tumour"
  ),
  out_supp6, 
  0.01, 
  "_01sig_treg_dataset2",
  c(
    "B cells",
    "Dendritic cells CD103",
    "Dendritic cells",
    "Fibroblasts",
    "Leukocytes unclassified",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells CD4",
    "T cells CD8",
    "T cells DN",
    "T reg cells", 
    "Tumour"
  )
)


# Community E - No Treg 
select_metacluster = neighbourhood_permutation(
  comE_order[,-which(names(comE_order) == "cellID")], 
  obj_relE$`No Tregs`,
  out_stats6, 
  "TM2_2_Community", 
  "TM2_2_Community", 
  c(
    "T cells CD4",
    "T cells CD8"
  ),
  c(
    "B cells",
    "Dendritic cells CD103",
    "Dendritic cells",
    "Fibroblasts",
    "Leukocytes unclassified",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells CD4",
    "T cells CD8",
    "T cells DN",
    "T reg cells", 
    "Tumour"
  ),
  out_supp6, 
  0.01, 
  "_01sig_notreg_dataset2", 
  c(
    "B cells",
    "Dendritic cells CD103",
    "Dendritic cells",
    "Fibroblasts",
    "Leukocytes unclassified",
    "Macrophages type 1",
    "Macrophages type 2",
    "T cells CD4",
    "T cells CD8",
    "T cells DN",
    "T reg cells", 
    "Tumour"
  )
)





