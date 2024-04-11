
#### Figure 5 - distance calculation ####

## Megan Cole (megan.cole@crick.ac.uk)

import pandas as pd 
import scipy
import os
from scipy.spatial import KDTree
from datetime import datetime 


data2 = pd.read_csv("/Users/colem/Documents/Projects/Treg_paper/Files/20230714_dataset2_celldata.csv", low_memory=False)
neighb2 = pd.read_csv("/Users/colem/Documents/Projects/Treg_paper/Files/20220823_dataset2_neighb_254communities_agglom18.csv")

out_dir = f"/Users/colem/Documents/Projects/Treg_paper/Files/Golden_files/Figure5/dataset2_communityC_min_distance/"
os.mkdir(out_dir)

mrtx = data2[(data2.treatment == 'MRTX')]

ROIs = set(mrtx["ROI_name"])
ROIs = list(ROIs)

##########################################################

###################
#### Functions ####

# Calculate distance between cells 
def distance_matrix(cutoff, ROI, points1, points2):
 
    tree1 = scipy.spatial.cKDTree(points1, leafsize=16)

    tree2 = scipy.spatial.cKDTree(points2,leafsize=16)
    distances = tree1.sparse_distance_matrix(tree2, cutoff, output_type='dict')
    
    # Only carry on with next steps if distances dictionary has neighbour values 
    if bool(distances) == True:
        print(f"distances found for {ROI}")
        keys = pd.DataFrame.from_dict(distances.keys())
        values = pd.DataFrame.from_dict(distances.values())
        # Give name to values dataframe 
        values.columns = ['distance']
        # Concatenate keys and values dataframes 
        neighbours = pd.concat([keys, values], axis = 1)
        # Sort data frame based on ascending order of first column 
        neighbours.sort_values([0,1], inplace = True) #Check this is sorting and not new values
        # Reset the index 
        neighbours = neighbours.reset_index(drop = True)
        # Rename column names 0 --> 'source' and 1 --> 'target'
        neighbours.rename(columns={neighbours.columns[0]: 'source', neighbours.columns[1]: 'target'}, inplace=True)
    else:
        neighbours = pd.DataFrame(distances)

    return neighbours

# Assign cell identification to distance data 
def assigning_cellID(distances, cellnames, filenames, subset1, subset2, marker1, marker2 = ""):

    # Remove any distance values below 2 pixels
    distances = distances[((distances['distance'] > 2))]

    # Link distances.source values to cellIDs in subset1
    subset1_cellID = pd.DataFrame(subset1.cellID.unique(), columns = ['source_cellID'])
    subset1_cellID = subset1_cellID[subset1_cellID.index.isin(distances.source)]

    # Link distances.target values to cellIDs in subset2
    subset2_cellID = pd.DataFrame(subset2.cellID.unique(), columns = ['cellID'])
    subset2_cellID = subset2_cellID[subset2_cellID.index.isin(distances.target)]

    # Add cellID info for source and target cells
    distances = distances.set_index(['source'])
    distj = distances.join(subset1_cellID)

    distj = distj.set_index(['target'])
    distj = distj.join(subset2_cellID)

    # Add marker info for source cells 
    distj = distj.set_index(['source_cellID'])
    subset1 = subset1.set_index(['cellID'])
    distj = distj.join(subset1)
    distj = distj.reset_index()
    # Rename new columns linking them to source cells
    distj = distj.rename(columns = {'source_cellID':'source_ID', 'index':'source_ID', cellnames:'source_cluster', 'Location_Center_X':'source_X', 'Location_Center_Y':'source_Y', marker1:f"source_{marker1}"})

    # Add marker info for target cells 
    distj = distj.set_index(['cellID', 'treatment', filenames])
    subset2 = subset2.set_index(['cellID', 'treatment', filenames])
    distj = distj.join(subset2)
    distj = distj.reset_index()

    # Rename and re-order columns 
    if marker2 == "":
        distj = distj.rename(columns = {'cellID':'target_ID', cellnames:'target_cluster', 'Location_Center_X':'target_X', 'Location_Center_Y':'target_Y'})
        distj = distj[['source_ID', 'target_ID', 'distance', 'source_X', 'source_Y', 'target_X', 'target_Y', 'source_cluster', f"source_{marker1}",
                'target_cluster', filenames, 'treatment']]
    else: 
        distj = distj.rename(columns = {'cellID':'target_ID', cellnames:'target_cluster', 'Location_Center_X':'target_X', 'Location_Center_Y':'target_Y', marker2:f"target_{marker2}"})
        distj = distj[['source_ID', 'target_ID', 'distance', 'source_X', 'source_Y', 'target_X', 'target_Y', 'source_cluster', f"source_{marker1}",
        'target_cluster', f"target_{marker2}", filenames, 'treatment']]
    
    print('Assinging cellID complete')

    # For each source_ID, take the shortest distance to each target cell type only
    distj = distj.loc[distj.groupby(['source_ID', 'target_cluster']).distance.idxmin()]

    return(distj)

# Run distance caluclation on each ROI for low and high CXCL9 expression for cell types of interest
def distance_rel(ct):

    d = pd.DataFrame()
    
    for expression in ['low', 'high']:
        print(f"{expression} CXCL9 expression")
        merged = dist_calc(expression, ct)

        if expression == 'low':
            d = merged 
        else:
            d = pd.concat([d,merged])
            d = d.reset_index(drop = True)
    # Save 
    d.to_csv(f"{out_dir}{datetime.today().strftime('%Y%m%d')}_DCCXCL9_high_low_min_distance_to_{ct}TC_PD1high_MRTX.csv") 


def dist_calc(expression, ct):

    merged = pd.DataFrame()

    for ROI in ROIs:
        print(ROI)
        subset = neighb2[(neighb2.top5 == "C")]

        data2_subset = data2[data2.cellID.isin(subset.cellID)]

        if expression == 'low':
            sub1 = data2_subset[((data2_subset['cellType'].str.contains('Dend')) & (data2_subset['MI_CXCL9'] <= 0.5) & (data2_subset['ROI_name'] == ROI))] \
                [['cellID', 'MI_PDL1', 'treatment', 'ROI_name', 'cellType', 'Location_Center_X', 'Location_Center_Y']]
        else: 
            sub1 = data2_subset[((data2_subset['cellType'].str.contains('Dend')) & (data2_subset['MI_CXCL9'] >= 0.5) & (data2_subset['ROI_name'] == ROI))] \
                [['cellID', 'MI_PDL1', 'treatment', 'ROI_name', 'cellType', 'Location_Center_X', 'Location_Center_Y']]

        sub2 = data2_subset[((data2_subset['cellType'].str.contains(ct)) & (data2_subset['MI_PD1'] >= 0.5) & (data2_subset['ROI_name'] == ROI))] \
        [['cellID', 'treatment', 'ROI_name', 'cellType', 'Location_Center_X', 'Location_Center_Y']]

    # If either sub1 or sub2 are empty, move on
        if sub1.empty == True:
            print(f"No sub1 data identified for {ROI}")
            continue

        if sub2.empty == True:
            print(f"No sub2 data identified for {ROI}")
            continue

        distances = distance_matrix(800, ROI, sub1[['Location_Center_X', 'Location_Center_Y']], sub2[['Location_Center_X', 'Location_Center_Y']])
        if distances.empty == True:
            print(f"No distances found for {ROI}")
            continue

        neighbours = assigning_cellID(distances, "cellType", "ROI_name", sub1, sub2, "MI_PDL1", marker2 = "")
        merged = pd.concat([merged, neighbours])
    
    merged['type'] = f'CXCL9 {expression}'
    return merged

##########################################################

# Fig 5a
distance_rel("CD8")

# Supp Fig 5a
distance_rel("CD4")

# Supp Fig 5b 
distance_rel("reg")


