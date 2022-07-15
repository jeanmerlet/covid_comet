import os, sys, glob, time, re
from tqdm.notebook import tqdm
import numpy as np
import pandas as pd
from scipy import signal
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Circle, Wedge
from multiprocessing import Pool

#
# parse data
#

# options
countries_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/scripts/plotting/gifs/countries_formatted.txt'
date_locs_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2021_09_30/dates_locations.tsv'
clusters_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/hip_mcl'
map_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/scripts/plotting/gifs/world_map.jpeg'
gif_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/plots/gifs'
jaccard_subset = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/data/jaccard/nz_cluster0_top_100.csv'

### Jaccard stuff ###
inf_to_get = 1.2
# get num_clusters
num_clusters = 0
with open(jaccard_subset) as in_file:
    for i, line in enumerate(in_file):
        if i == 0: continue
        num_clusters += 1
        
jaccard_clusters = pd.DataFrame(data=np.zeros((num_clusters, 3), dtype=int), columns=['inf', 'cluster', 'size'])

clusters_to_load = []
with open(jaccard_subset) as in_file:
    for i, line in enumerate(in_file):
        if i == 0: continue
        cluster_size, cluster_id, _, _ = line.strip().split(',')
        cluster_num, inf = cluster_id.split('-')
        jaccard_clusters.iloc[i-1, :] = [float(inf), int(cluster_num), int(cluster_size)]
            
jaccard_clusters = jaccard_clusters.sort_values(by='size', ascending=False)
jaccard_clusters = jaccard_clusters.loc[jaccard_clusters['inf'] == inf_to_get, :].iloc[:10, :]
print(jaccard_clusters)

# parse clusters
inf_paths = []
for path in os.listdir(clusters_dir):
    if re.match('^clusters_inflation', path):
        inf_paths.append(os.path.join(clusters_dir, path))

inf_paths.sort()

clusters = []
large_clusters_idx = []
for inf_path in inf_paths:
    _, name = os.path.split(inf_path)
    inf = re.search('\d\.\d', name).group(0)
    sizes = jaccard_clusters['size'].values[(jaccard_clusters['inf'] == float(inf)).values]
    num_sizes = len(sizes)
    if num_sizes == 0: continue
    with open(inf_path) as in_file:
        for i, line in enumerate(in_file):
            size = len(line.strip().split(' '))
            if size in sizes:
                clusters.append(line.strip().split(' '))
                large_clusters_idx.append([i, size])
                num_sizes -= 1

inf_to_get = str(inf_to_get)

# load data
countries_df = pd.read_csv(countries_path, sep=',', low_memory=False)
date_locs_df = pd.read_csv(date_locs_path, sep='\t', low_memory=False)

# build country location dictionary
countries, lats, lngs = list(countries_df['name']), list(countries_df['latitude']), list(countries_df['longitude'])
country2coord = {c: (float(lat),float(lng)) for c, lat, lng in zip(countries, lats, lngs)}
    
# options
scale = 2000
step_size = 30

# initialize book keeping
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#cluster_names = ['Cluster {0}'.format(i+1) for i in range(len(clusters)-1)] + ['Other']

# initialize map parameters
image = mpimg.imread(map_path)
nlats, nlngs = image.shape[0], image.shape[1]
lat_min, lat_max = -60, +90
lng_min, lng_max = -180, +180
lat_step = (lat_max - lat_min)/nlats
lng_step = (lng_max - lng_min)/nlngs

# cluster names from pango lineage metadata
variant_name_map = {'B.1.1.7': 'Alpha',
                    'AY.4': 'Delta (AY.4)',
                    'B.1': 'B.1',
                    'B.1.2': 'B.1.2',
                    'B.1.177': 'B.1.177'}

def plot_moving_averages(i):
    
    moving_average = moving_averages[i]

    fig = plt.figure(figsize=(15, 2*nlats/nlngs*15))
    gridspec = fig.add_gridspec(200, 100)
    
    #
    # world map
    #
    
    #ax = fig.add_subplot(3, 1, 2)
    ax = fig.add_subplot(gridspec[1:99, :])
    
    # plot world map
    plt.title('7-day moving average of COVID-19 sequences by cluster', fontsize=20)
    plt.imshow(image, aspect='auto')
    plt.xlim([0, image.shape[1]-1])
    plt.ylim([image.shape[0]-1, 0])
    plt.axis('off')
    
    # loop over countries
    for j in np.argsort(moving_average.sum(-1))[::-1]:
        values = moving_average[j].clip(1e-8)
        lat, lng = [float(c) for c in coords_u[j][1:-1].split(', ')]
        lat = nlats - int((lat - lat_min)/lat_step)
        lng = int((lng - lng_min)/lng_step)
        for k in range(len(values)):
            start = int(359*values[:k].sum()/values.sum())
            stop = int(359*values[:k+1].sum()/values.sum())
            wedge = Wedge(
                center=(lng, lat),
                r=np.sqrt(values.sum()/np.pi)*6,
                theta1=start,
                theta2=stop,
                alpha=0.9,
                color=colors[k])
            ax.add_patch(wedge)
        #circ = Circle(
        #    (lng, lat), 
        #    radius=np.sqrt(values.max()/np.pi)*8, 
        #    alpha=0.75, 
        #    color=color = colors[values.argmax()])
        #ax.add_patch(circ)
    
    #
    # histograms
    #
    
    #ax = fig.add_subplot(3, 1, 1)
    ax = fig.add_subplot(gridspec[101:149, :])
    
    # build canvas
    plt.ylim([-0.03, 10.02])
    plt.xlim([-0.01, 5.01])
    plt.axis('off')
    plt.plot([2.0, 2.0], [0, 10], 'k-', linewidth=1)
    plt.plot([3.0, 3.0], [0, 10], 'k--', linewidth=1)
    plt.plot([4.0, 4.0], [0, 10], 'k--', linewidth=1)
    plt.plot([5.0, 5.0], [0, 10], 'k--', linewidth=1)
    plt.plot([0, 0, 5, 5, 0], [0, 10, 10, 0, 0], 'k-', linewidth=1)
    plt.text(0.985, 0.15, dates_u[i], 
             transform=ax.transAxes, 
             fontsize=24, 
             bbox=dict(facecolor='white'),
             verticalalignment='center', 
             horizontalalignment='right')
    
    # loop over clusters
    for j in range(10):
        value = int(moving_average[:,j].sum())
        color = colors[j]
        plt.fill_between([2.0, value/scale+2.0], 9.9-j, 10.1-j-1, color=color)
        plt.text(2.0/5-0.01, 0.95-(j*0.1), cluster_names[j], 
                 transform=ax.transAxes, 
                 fontsize=14,
                 verticalalignment='center', 
                 horizontalalignment='right')
        plt.text(2.0/5+value/scale/5+0.01, 0.95-(j*0.1), value, 
                 transform=ax.transAxes, 
                 fontsize=14,
                 verticalalignment='center', 
                 horizontalalignment='left')
    
    #
    # total cases
    #
    
    #ax = fig.add_subplot(3, 1, 3)
    ax = fig.add_subplot(gridspec[151:199, :])
    
    for j in range(cluster_sums.shape[1]):
        plt.fill_between(np.arange(len(cluster_sums)), cluster_sums[:,-j-1], color=colors[-j-1])
    plt.plot([i, i], [y_min, y_max], 'r-', linewidth=3)
    plt.ylim([y_min, y_max])
    plt.xlim([0, len(total_cases)-1])
    plt.xticks(ticks=tick_idx, labels=tick_dates, rotation=20)
    plt.yticks([])
    
    inflation = re.search('\d\.\d', clusters_filename)[0]
    out_dir = os.path.join(gif_dir, 'inf_' + inflation + '_top3')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    plt.savefig(os.path.join(out_dir, '{0}.png'.format(str(i).zfill(len(str(len(dates_u)))))),
                bbox_inches='tight')
    plt.clf()
    plt.close("all")
        
### clusters ###

clusters_paths = []
for r, d, f in os.walk(clusters_dir):
    for clusters_path in f:
        if 'variant' not in clusters_path:
            clusters_paths.append(os.path.join(r, clusters_path))

for clusters_path in clusters_paths:
    _, clusters_filename = os.path.split(clusters_path)
    inflation = re.search('\d\.\d', clusters_filename)[0]
    #if inflation != '6.0': continue
    if inflation != inf_to_get: continue
    # parse ids, dates, and countries
    ids, dates, locs = list(date_locs_df['seq_id']), list(date_locs_df['collection_date']), list(date_locs_df['location'])
    keep = []
    for i, date in enumerate(tqdm(dates)):
        date = date.split('-')
        if len(date) < 2: continue
        if len(date) == 2: date = date + ['15']
        if int(date[0]) < 2019: continue
        if int(date[0]) == 2019 and int(date[1]) < 12: continue
        dates[i] = '-'.join(date)
        keep.append(i)

    # filter sequences
    ids = [ids[i] for i in keep]
    dates = [dates[i] for i in keep]
    locs = [locs[i].split('/')[1].strip() for i in keep]

    # build sequence id dictionary
    sequence_dict = {i: {'location':l, 'date':d} for i, d, l in zip(ids, dates, locs)}

    # make set of unique countries
    countries, coords = [c for c in country2coord], [str(country2coord[c]) for c in country2coord]
    coords_u = list(np.unique(coords))
    labels = {n:coords_u.index(c) for n, c in zip(countries, coords)}

    # make set of unique dates
    dates_u = np.unique(dates)
    dates_u.sort()
    dates_u = list(dates_u)

    #clusters = []
    with open(clusters_path, 'r') as f:
        for line in f:
            break
            line = line.split(' ')[:-1]
            clusters.append(line)

    lengths = []
    for cluster in clusters:
        break
        lengths.append(len(cluster))

    #lengths.sort()
    #lengths = lengths[::-1]
    print(f'### inflation {inflation} ###')
    for i, length in enumerate(lengths):
        break
        #if i == 15: break
        #print(f'{i}: {length}')
    #cutoff = lengths[8] - 1
    #print(f'cutoff: {cutoff}')
    # get idxs and sizes of largest clusters
    #large_clusters_idx = []
    for i, cluster in enumerate(clusters):
        break
        #if len(cluster) > cutoff:
            #print(len(cluster))
            #large_clusters_idx.append([i, len(cluster)])
    # combine small clusters into one        
    #large_clusters = [c for c in clusters if len(c)>cutoff]
    #small_clusters = sum([c for c in clusters if len(c)<cutoff], [])
    #sizes = [len(c) for c in large_clusters]
    #clusters = [large_clusters[i] for i in np.argsort(sizes)[::-1]] + [small_clusters]
    # sort idxs in same way
    #large_clusters_idx = [large_clusters_idx[i] for i in np.argsort(sizes)[::-1]]
    # build tensor
    data = np.zeros([len(dates_u), len(coords_u), len(clusters)], dtype=np.float32)
    for c_idx, cluster in enumerate(tqdm(clusters)):
        for s_idx, sequence in enumerate(tqdm(cluster)):
            if sequence in sequence_dict:
                loc, date = sequence_dict[sequence]['location'], sequence_dict[sequence]['date']
                data[dates_u.index(date), labels[loc], c_idx] += 1

    # compute moving average
    w = 7
    moving_averages = signal.convolve(data, np.ones([w,1,1])/w, mode='valid').clip(0, np.inf)
    dates_u = dates_u[w-1:]
    # initialize total cases
    total_cases = np.round(moving_averages.sum(-1).sum(-1))
    cluster_cases = np.round(moving_averages.sum(1))
    cluster_sums = np.cumsum(cluster_cases, axis=-1)
    tick_idx = np.arange(0, len(total_cases), 30)
    tick_dates = [dates_u[i][:-3] for i in tick_idx]
    y_min, y_max = 0, 1.1*total_cases.max()

    num_lineages = 3
    lineages_name = 'all_variant-ids_' + clusters_filename
    all_lineages = pd.read_csv(os.path.join(clusters_dir, lineages_name), sep='\t', index_col=None, header=None, usecols=list(range(num_lineages)), engine='python')
    cluster_names = []
    for cluster_info in large_clusters_idx:
        idx, cluster_size = cluster_info[0], cluster_info[1]
        lineages = list(all_lineages.loc[idx, :].values)
        name = []
        for lineage in lineages:
            if lineage:
                lineage_name, size = lineage.split(' ')
                pct_lineage = np.round(int(size) / cluster_size, 2)
                name.append(f'{lineage_name} ({pct_lineage})')
            else:
                name.append('none')
        name = ', '.join(name)
        cluster_names.append(name)

    cluster_names.append('Other')

    with Pool(processes=16) as pool:
        pool.map(plot_moving_averages, list(np.arange(len(moving_averages))))
