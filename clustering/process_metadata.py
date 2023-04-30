import pandas as pd
import numpy as np
import os

in_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/metadata.tsv'
out_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/dates_locations_pango.tsv'

cols_to_load = ['Accession ID', 'Collection date', 'Location', 'Pango lineage']
meta = pd.read_csv(in_path, sep='\t', header=0, index_col=None, low_memory=False, usecols=cols_to_load)
print('metadata loaded...')
# rename useful columns
meta = meta.rename(columns={'Accession ID': 'seq_id', 'Collection date': 'collection_date', 'Location': 'location', 'Pango lineage': 'pango'})

# fill sequence ids with leading 0s to match CoMet output
seq_ids = list(meta['seq_id'].values)
seq_ids = [seq_id.split('_')[-1] for seq_id in seq_ids]
seq_ids = [str(seq_id).zfill(9) + '_A' for seq_id in seq_ids]
meta['seq_id'] = seq_ids

# add 15th of month for collection dates without a day
# and remove year only dates and dates with 2019
dates = meta['collection_date'].values.astype('str')
new_dates = []
for date in dates:
    if len(date) == 7:
        new_dates.append(date + '-15')
    #elif '2019' in date:
    #    new_dates.append('2019')
    elif len(date) == 8 or len(date) == 9:
        year, month, day = date.split('-')
        if len(month) == 1:
            month = '0' + month
        if len(day) == 1:
            day = '0' + day
        new_dates.append('-'.join([year, month, day]))
    else:
        new_dates.append(date)

meta['collection_date'] = new_dates
good_date_idx = [len(date) > 4 for date in new_dates]
meta = meta.loc[good_date_idx, :]
assert np.all([len(date) == 10 for date in list(meta['collection_date'].values)])
print('dates processed...')

# drop all columns other than seq_id and collection_date and pango lineage
meta = meta.loc[:, ['seq_id', 'collection_date', 'location', 'pango']]
# add wuhan sequence
wuhan = pd.DataFrame(data=[['000000000_A', '2019-12-15', 'Asia / China / Wuhan', 'wuhan']], columns=['seq_id', 'collection_date', 'location', 'pango'])
meta = pd.concat([wuhan, meta])
# write to file
meta.to_csv(out_path, index=False, sep='\t')
print('done!')
