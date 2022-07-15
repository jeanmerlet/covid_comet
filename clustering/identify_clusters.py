import pandas as pd
import numpy as np
import os

meta_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2021_09_30/metadata.tsv'
clusters_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/hip_mcl'
clusters_paths = []
for r, d, f in os.walk(clusters_dir):
    for clusters_path in f:
        if 'variant' not in clusters_path:
            clusters_paths.append(os.path.join(r, clusters_path))

meta = pd.read_csv(meta_path, sep='\t', header=0, index_col=None, low_memory=False)
meta = meta.rename(columns={'Accession ID': 'seq_id', 'Pango lineage': 'pango_lineage'})
# fill sequence ids with leading 0s to match CoMet output
seq_ids = list(meta['seq_id'].values)
seq_ids = [seq_id.split('_')[-1] for seq_id in seq_ids]
seq_ids = [str(seq_id).zfill(9) + '_A' for seq_id in seq_ids]
meta['seq_id'] = seq_ids
meta = meta.loc[:, ['seq_id', 'pango_lineage']]
wuhan = pd.DataFrame(data=[['000000000_A', 'wuhan']], columns=['seq_id', 'pango_lineage'])
meta = pd.concat([wuhan, meta])
meta = meta.set_index('seq_id')

all_details = True

for clusters_path in clusters_paths:
    with open(clusters_path) as clusters_file:
        for i, cluster in enumerate(clusters_file):
            lineages = {}
            seq_ids = cluster.split(' ')
            for seq_id in seq_ids:
                if seq_id == '\n': continue
                lineage = meta.loc[seq_id, 'pango_lineage']
                if lineage in lineages.keys():
                    lineages[lineage] += 1
                else:
                    lineages[lineage] = 1
            top = max(lineages, key=lineages.get)
            _, name = os.path.split(clusters_path)
            if all_details:
                out_path = os.path.join(clusters_dir, 'all_variant-ids_' + name)
            else:
                out_path = os.path.join(clusters_dir, 'single_variant-ids_' + name)
            with open(out_path, 'a') as out_file:
                if all_details:
                    while lineages:
                        top = max(lineages, key=lineages.get)
                        out_file.write(top + ' ' + str(lineages[top]) + '\t')
                        del lineages[top]
                else:
                    out_file.write(top + '\t' + str(lineages[top]))
                out_file.write('\n')


