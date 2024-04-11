
#### Agglomeration of communities ####

## Megan Cole (megan.cole@crick.ac.uk)

import pandas as pd 
import os 
from datetime import datetime 
from sklearn.cluster import AgglomerativeClustering 

# Set path 
mypath = "/Users/colem/Documents/Projects/Treg_paper/Golden/Data/"

output_dir = f"{mypath}agglomerate_communities_python/"
os.mkdir(output_dir)

#########################################################

#### Agglomerative clustering function ####

def agglomerative_clustering(k, average_neighbours, avg_noCluster, agglomerate_to):
    
    # Agglomerative clustering to chosen number of communities as determined by own biological knowledge of the data
    ac1 = AgglomerativeClustering(linkage = 'average', n_clusters = agglomerate_to)
    agglomerate1 = pd.DataFrame(ac1.fit_predict(avg_noCluster))
    # Rename column
    agglomerate1.columns = [f"agglomerate_{agglomerate_to}"]
    agglomerate1[f"agglomerate_{agglomerate_to}"] += 1
    # Add agglomerated data to average neighbours 
    average_neighbours = average_neighbours.join(agglomerate1)
   
    print(list(average_neighbours))
    
    # Save data with added columns 
    average_neighbours.to_csv(os.path.join(output_dir, f"{datetime.today().strftime('%Y%m%d')}_average_neighbours_{k[-1]+1}clusters_agglom_average.csv"), index = False)
    print('dataset saved')

    return(average_neighbours)


#########################################################

#### Dataset 1 ####
data1 = pd.read_csv(f"{mypath}20220823_dataset1_neighb_62communities_agglom18.csv")
print(list(data1))

# Create average_neighbours
data1 = data1.loc[:,"B_cells":"cluster"]
data1.columns = data1.columns.str.replace("_", " ")

# Generate average neighbours 
average_neighbours = data1.groupby('cluster', as_index=False).mean()
print(list(average_neighbours))

avg_noCluster = average_neighbours.drop('cluster',1)
print(avg_noCluster.shape)
avg_noCluster.head()

# Agglomerate to 30 communities 
average_neighbours30 = agglomerative_clustering(range(2,average_neighbours.cluster.max()), average_neighbours, avg_noCluster, 30)
# Agglomerate to 18 communites
average_neighbours30_18 = agglomerative_clustering(range(2,average_neighbours30.cluster.max()), average_neighbours30, avg_noCluster, 18)


#### Dataset 2 ####
data2 = pd.read_csv(f"{(mypath)}20220823_dataset2_neighb_254communities_agglom18.csv")
print(list(data2))

# Create average_neighbours
data2 = data2.loc[:,"B_cells":"cluster"]
data2.columns = data2.columns.str.replace("_", " ")

# Generate average neighbours 
average_neighbours2 = data2.groupby('cluster', as_index=False).mean()
print(list(average_neighbours2))

avg_noCluster2 = average_neighbours2.drop('cluster',1)
print(avg_noCluster2.shape)
avg_noCluster.head()

# Agglomerate to 18 communities 
average_neighbours2 = agglomerative_clustering(range(2,average_neighbours2.cluster.max()), average_neighbours2, avg_noCluster2, 18)
    

###########################################################################
