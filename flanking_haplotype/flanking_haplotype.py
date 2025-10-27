import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.collections import PatchCollection
#from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import seaborn as sns
from tqdm import tqdm
#from matplotlib.colors import ListedColormap
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib
##############################################################
#
# Script for Figure 2: 
# Author: Jake Marglous
#
##############################################################

###Read in combined and filtered variant table
full_data_T = pd.read_csv('gonder/gonder_flanking_region.csv', sep = ',', index_col = 0)
full_data_T = full_data_T.replace({0.5: 1, np.nan: 0})


###Convert the column labels back to int from str 
full_data_T.columns = full_data_T.columns.astype(int)
full_data_T = full_data_T.loc[:, (full_data_T != 0).any(axis=0)]

###Optionally reindex sites to set the most frequent allele within the group at each site to the ref
c580y = 1725133

def reindex_sites(variant_df_T, exclude = [c580y]):
    with open('outfile.txt', 'w') as outfile: 
        for site in variant_df_T.columns:
            if site not in exclude:
                counts = variant_df_T[site].value_counts()
                outfile.write(str(counts))
                ordered_counts = list(counts.index)
                #Swap 0 and 2 genotype labels if alt is more common at site than ref
                if 0 in ordered_counts and 2 in ordered_counts:
                    if ordered_counts.index(2) < ordered_counts.index(0):
                        variant_df_T[site] = variant_df_T[site].apply(lambda x: 2 if x == 0 else (0 if x ==2 else x))
                elif -1 in ordered_counts and 2 in ordered_counts and 0 not in ordered_counts and 1 not in ordered_counts:
                    variant_df_T[site] = variant_df_T[site].apply(lambda x: 0 if x == 2 else x)
                elif 1 in ordered_counts and 2 in ordered_counts and 0 not in ordered_counts:
                    variant_df_T[site] = variant_df_T[site].apply(lambda x: 0 if x == 2 else x)

                #replace alts with refs if site is invariant
                elif 1 not in ordered_counts and len(ordered_counts) == 1:
                    variant_df_T[site] = variant_df_T[site].apply(lambda x: 0)

                if 1 in ordered_counts and ordered_counts.index(1) == 0:
                    variant_df_T = variant_df_T.drop(site, axis = 1)
    return variant_df_T

###Use unique encoding for alt at c580y locus so it can get its own color in the plot
full_data_T[c580y] = full_data_T[c580y].apply(lambda x: 4 if x == 1 else (3 if x == 0 else x))

###This is a simple custom distance calculation I made to cluster samples by pairwise comparison over only loci where both samples have a single allele (i.e. 1 or 0 in the matrix)
###Otherwise the samples tend to cluster by missingness
def dist(df_pivot, sample_1, sample_2):

    #Get the genotypes of two samples 
    data1 = np.array(df_pivot[df_pivot.index == sample_1])[0]
    data2 = np.array(df_pivot[df_pivot.index == sample_2])[0]
    #Set the initial distance of the two samples, and the number of compared samples, to 0
    sample_dist = 0
    comp_loci = 0
    for i in range(data1.shape[0]):
        if (0<=data1[i]<=1) and (0<=data2[i]<=1):
            #If the two samples' genotypes are >=0 (i.e. either ref or alt, not missing or mixed), add their difference to the distance and add 1 to the number of compared loci
            comp_loci += 1
            sample_dist += np.abs(data1[i]-data2[i])

    if comp_loci == 0:
        return data1.shape[0]
    else:
        return sample_dist/comp_loci

def cluster_samples(df_pivot):

    ###This runs through the distance function defined above for each sample pair and gets a distance matrix, then uses the scipy distance module to 
    ###Calculate a new clustered sample order for plotting. 
    ###It also saves the clusters for plotting later. 
    dist_dict = {}
    sample_list = df_pivot.index.to_list()
    print('calculating pairwise sample distances')
    for sample_i in tqdm(range(len(sample_list))):
        sample_i_dists = []
        for sample_j in range(len(sample_list)):
            sample_i_dists.append(dist(df_pivot, sample_list[sample_i], sample_list[sample_j]))
        dist_dict[sample_i] = sample_i_dists
    dist_mat_df = pd.DataFrame.from_dict(dist_dict, orient = 'index', columns = sample_list)

    print('clustering')
    dist_matrix_reduced = sp.spatial.distance.squareform(dist_mat_df)
    clusters = linkage(dist_matrix_reduced, metric = 'euclidean')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    dn = dendrogram(clusters, ax= ax)


    ordered_samples = ax.axes.get_xticklabels()
    sample_order_indexes = []
    for i in ordered_samples:
        sample_order_indexes.append(int(i.get_text()))

    sample_order_labels = []
    for index in sample_order_indexes:
        sample_order_labels.append(sample_list[index])

    return sample_order_labels, clusters

clusters = linkage(full_data_T, method = 'complete')
dn = dendrogram(clusters)

###Reorder the samples and reorganize the dataframe according to the hierarchical clustering
new_sample_order = []
for sample in dn['ivl']:
    new_sample_order.append(list(full_data_T.index)[int(sample)])

full_data_T = full_data_T.reindex(new_sample_order)

###reverse order of rows so ayalew's sample is on top
full_data_T = full_data_T.iloc[::-1]

###Define a colormap
###(I've been plotting from a  list of colors -- this is going to use a light grey (gainsboro) 
###as the background/REF and viridis colors (purple, blue, green, yellow, etc.) for ALT alleles)
cmap = plt.colormaps['viridis']
colors = cmap(np.linspace(0, 1, 4))
vir_colors_list = colors.tolist()
cmap = ['gainsboro'] + vir_colors_list

geo_cmap = ['#E9F4FA', "#EBABBD", "#EC741E", '#023E73', '#FF7EE5', '#D3D3D3']

genotype_of_interest_color = "#0F0F10"

###Make a list of the patches (one for every (sample, locus) coordinate) to add to the plot
def make_patches(pivot_table,
                 cmap, 
                 pt_of_interest = [], 
                 pt_of_interest_color = 'magenta', 
                 index_vis_upper_lim = 4):

    sample_order = pivot_table.index
    site_list = pivot_table.columns
    
    patches = []
    for sample_i in tqdm(range(len(sample_order))):
        for site_j in range(len(site_list)):
            
            x_pad = .05
            y_pad = 0
            allele = pivot_table.loc[sample_order[sample_i], site_list[site_j]]
            if allele >= 0:
                if site_list[site_j] in pt_of_interest:
                    facecolor = pt_of_interest_color
                    
                else:
                    facecolor = cmap[int(np.min([index_vis_upper_lim, allele]))]
            else:
                facecolor = 'white'
            
            rect = Rectangle(xy = (site_j+x_pad, sample_i+y_pad), 
                            width = 1-2*x_pad, 
                            height = 1, 
                            facecolor =  facecolor,
                            edgecolor = facecolor, 
                            linewidth = .001
                            )

            if allele != 0: 
                patches.append(rect)
    
    return patches

patches_for_plot = make_patches(
    pivot_table = full_data_T, 
    cmap = geo_cmap,
    pt_of_interest = [172513], 
    pt_of_interest_color = genotype_of_interest_color
)

###Sort the samples into groups
samples = full_data_T.index
# SEA_samples = [sample for sample in samples if sample not in ['merge_GM331_580Y']]
# Eth_samples = [sample for sample in samples if sample in ['merge_GM331_580Y']]

region_rules = {r"GM": "Gonder",
                r"AC": "Gambella",
                r"MM":"Mizan",
                r"ALH":"Arbaminch",
                r"SRR":"x",
                 r"ERR11":"Amhara",
                  r"ERR10":"Oromia",
                   r"MFHC":"Metehara" }
region_color_dict = {"Gonder": "red",
                     "Gambella": "blue",
                     "Mizan": "green",
                     "Arbaminch": "purple",
                     "Amhara": "orange",
                     "Oromia": "brown",
                     "Metehara": "pink",
                     "x": "grey"}
#treatment_color_dict = {"Tf":"red","Tf0":"blue","Ts":"green"}
###Assign colors to each of the sample groups
# group_color_dict = {
#     'Pf7 SEA':'lightgrey',
#     'Eth. surveill.': 'red'
# }

###Label each of the samples by a color according to their group
metadata = pd.read_csv("./meta_data.csv")

# sample_color_dict = {}
# for sample in new_sample_order:
#     assigned = False
#     for pattern, region in region_rules.items():
#         if re.search(pattern, sample):   # check if sample matches regex
#             sample_color_dict[sample] = region_color_dict[region]
#             assigned = True
#             break
#     if not assigned:
#         sample_color_dict[sample] = "black"   # default color if no match

sample_color_dict1 = {}
for sample in new_sample_order:
    region = metadata.loc[metadata["sample"] == sample, "Treatment"].values
    if len(region) > 0:
        region_name = region[0]
        sample_color_dict1[sample] = region_color_dict.get(region_name, "black")  # default = black
    else:
        sample_color_dict1[sample] = "black" 

#regions = metadata["Treatment"].unique()

# Pick a colormap with enough distinct colors (e.g. tab20, Set3, hsv, etc.)
cmap = matplotlib.colormaps["tab20"]

# Build automatic region→color dictionary
#colors = [cmap(i / len(regions)) for i in range(len(regions))]

# Convert to hex and build dict
#region_color_dict = {region: mcolors.to_hex(colors[i]) for i, region in enumerate(regions)}
#plt.rcParams.update({'font.size': 14})

###Create figure
fig = plt.figure(figsize = [30, 10], tight_layout = True)
gs = gridspec.GridSpec(nrows = 1, ncols = 2, width_ratios=[1/8, 1], wspace = .1, figure=fig)

gs0 = gs[0].subgridspec(nrows = 1, ncols = 2, width_ratios = [1, .333], wspace=.15)

ax00 = fig.add_subplot(gs0[0])
ax0 = fig.add_subplot(gs0[1])
ax = fig.add_subplot(gs[1])

dn = dendrogram(clusters, orientation = 'left', ax = ax00, color_threshold = 0)
#ax00.yaxis.set_inverted(True)

#ax00.set_xlim(l_lim, r_lim)
ax00.set_xticks([])
ax00.axis('off')

###Draw a grey background instead of plotting default grey to reduce figure size 
background = Rectangle(xy = [0,0], width = full_data_T.shape[1], height = full_data_T.shape[0], color = geo_cmap[0])
ax.add_patch(background)

###Pull site and sample names from dataframe for easier calling later
sample_order = full_data_T.index
site_list = full_data_T.columns

###Add the variant rectangles drawn in the make_patches command 
ax.add_collection(PatchCollection(patches_for_plot, match_original = True))

###Remove labels on the x-axis for better visualization
ax.set_xlim(0, len(site_list))
ax.set_xticks([])
x = 0
height = 0
box_buffer = .15
for sample in sample_order:
    group_rect = Rectangle((x, height+box_buffer), 1, 1-box_buffer, linewidth = 1, edgecolor = 'black', facecolor = sample_color_dict1[sample])
    ax0.add_patch(group_rect)
    height+=1
    ax0.set_ylim(0,len(sample_order))
    ax0.set_xlim(0-box_buffer, 1+box_buffer)

    ax0.axis('off')

# legend_elements = []

# for item in region_color_dict.keys():
#     # legend_elements.append(Patch(facecolor=region_color_dict[item], edgecolor= 'black',
#     #                      label=item))

#     legend_elements.append(Patch(facecolor=cmap(0.0), edgecolor=cmap(0.0),
#                          label='REF'))
# # legend_elements.append(Patch(facecolor=cmap(1.0), edgecolor=cmap(1.0),
# #                          label='MIXED'))        
# legend_elements.append(Patch(facecolor=cmap(2.0), edgecolor=cmap(2.0),
#                          label='ALT')) 
# legend_elements.append(Patch(facecolor='#FF7EE5', edgecolor='#FF7EE5',
#                          label='622'))

legend_elements = [
    Patch(facecolor=geo_cmap[0], label='REF'),
    Patch(facecolor=geo_cmap[1], label='ALT'),
    Patch(facecolor=geo_cmap[3], label='R622I_REF'),
    Patch(facecolor=geo_cmap[4], label='R622I_ALT')
]

legend = ax0.legend(handles = legend_elements, loc = [-3.5, .7])

for patch in legend.get_patches()[3:]:
    patch.set_width(5)

sample_order = [item.replace('merge_', '') for item in sample_order]
sample_order = [item.replace('_580Y', '') for item in sample_order]

###Label the y-axis with sample names
ax.set_ylim(0, len(sample_order))
ax.set_yticks(np.arange(.5, len(sample_order)+.5, 1))
yticklabels = ax.set_yticklabels(sample_order)

ax.set_xticks([.5, list(full_data_T.columns).index(1725133)+.5, len(full_data_T.columns)+.5])
ax.set_xticklabels([full_data_T.columns[0], 1725133, full_data_T.columns[-1]])

###Save plot
#fig.savefig('./test_all.pdf', dpi = 500)
#fig.savefig('./test_all.svg', dpi = 500)
fig.savefig('gonder/test_622.png', dpi = 500)
   