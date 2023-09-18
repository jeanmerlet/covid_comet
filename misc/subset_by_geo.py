import pandas as pd
import numpy as np
import os, re

meta_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/metadata.tsv'

# cols containing epi_id, data, and location
usecols = [2, 3, 4]

meta = pd.read_csv(meta_path, sep='\t', header=0, index_col=None, usecols=usecols, low_memory=False)
locs = meta.loc[:, 'Location'].values

conts, countries = [], []
for loc in locs:
    cont = loc.split('/')[0].strip()
    country = loc.split('/')[1].strip()
    conts.append(cont)
    countries.append(country)

uniq_conts, counts_conts = np.unique(conts, return_counts=True)
print(f'num continents: {len(uniq_conts)}\n')

for i, cont in enumerate(uniq_conts):
    print(f'{cont}: {counts_conts[i]:,} samples')

# breakdown with a specific continent
cont_choice = 'North America'
index = np.where(np.array(conts) == cont_choice)[0]
choice_countries = np.array(countries)[index]

# further filter to specific countries
country_choices = ['USA', 'Canada', 'Mexico']

def get_uniq_country_counts(countries, country_choices=None):
    uniq_countries, counts_countries = np.unique(countries, return_counts=True)
    print(f'num countries: {len(uniq_countries)}\n')
    total = 0
    for i, country in enumerate(uniq_countries):
        if country_choices is not None:
            if country in country_choices:
                print(f'{country}: {counts_countries[i]:,} samples')
                total += counts_countries[i]
        else:
            print(f'{country}: {counts_countries[i]:,} samples')
            total += counts_countries[i]
    print(f'\ntotal samples: {total:,}')
    return uniq_countries, counts_countries

uniq_choice_countries, choice_counts = get_uniq_country_counts(choice_countries, country_choices)

# filter a country by requiring more than just country
#country = 'USA'
#country = 'Mexico'
country = 'Canada'
country_index = np.where(np.array(countries) == country)
country_locs = locs[country_index]
specific_index = []
for i, loc in enumerate(usa):
    if len(loc.split('/')) >= 3:
        specific_index.append(i)

country_locs = country_locs[specific_index]
print(f'number of filtered locs for {country}: {len(country_locs)}')
print(country_locs[:10])




