import pandas as pd
import numpy as np
import os, re

filtered_ids_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/combined_seq-headers.txt'
meta_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/metadata.tsv'
target_num_seqs = 250000

# make a whitelist from pre-filtered seq ids
whitelist_ids = []
with open(filtered_ids_path) as in_file:
    for i, line in enumerate(in_file):
        if i == 0: continue
        epi_id = line.strip().split('|')[-1]
        whitelist_ids.append(epi_id)

# only load cols containing epi_id, collection date, and location
usecols = [2, 3, 4]
meta = pd.read_csv(meta_path, sep='\t', header=0, index_col=None, usecols=usecols, low_memory=False)

# apply whitelist
meta = meta.set_index('Accession ID', drop=False)
meta = meta.loc[meta.index.intersection(whitelist_ids), :]
locs = meta.loc[:, 'Location'].values

# get index for sequences in USA, Canada, and Mexico
conts, countries = [], []
for loc in locs:
    cont = loc.split('/')[0].strip()
    country = loc.split('/')[1].strip()
    conts.append(cont)
    countries.append(country)

# calculate country proportions
country_choices = ['USA', 'Canada', 'Mexico']
country_counts = {}
uniq_countries, uniq_counts = np.unique(countries, return_counts=True)
for i, country in enumerate(uniq_countries):
    if country in country_choices:
        country_counts[country] = uniq_counts[i]

total_seqs = sum(country_counts.values())
num_target_seqs_per_country = {}
for country, num_seqs in country_counts.items():
    num_target_seqs_per_country[country] = int(round(num_seqs / total_seqs * target_num_seqs))

# get the index for matching countries
def get_index_by_countries(countries, country_choices):
    for i, country_choice in enumerate(country_choices):
        if i == 0:
            country_index = np.where(np.array(countries) == country_choice)[0]
            country_index = np.random.choice(country_index, size=num_target_seqs_per_country[country_choice], replace=False)
        else:
            index = np.where(np.array(countries) == country_choice)[0]
            index = np.random.choice(index, size=num_target_seqs_per_country[country_choice], replace=False)
            tmp = np.unique(np.concatenate((country_index, index)), return_counts=True)
            country_index = np.union1d(country_index, index)
    return country_index

country_index = get_index_by_countries(countries, country_choices)
meta = meta.set_index([list(range(meta.shape[0]))])
meta = meta.loc[country_index, :]
epi_ids = meta.loc[:, 'Accession ID']
out_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/usa_can_mex_epi-ids_250k.txt'
epi_ids.to_csv(out_path, header=False, index=False)
