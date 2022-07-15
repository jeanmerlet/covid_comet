import pandas as pd
import numpy as np
import os

in_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2021_09_30/metadata.tsv'
out_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2021_09_30/dates_locations.tsv'

meta = pd.read_csv(in_path, sep='\t', header=0, index_col=None)
# rename useful columns
meta = meta.rename(columns={'Accession ID': 'seq_id', 'Collection date': 'collection_date', 'Location': 'location'})

# fill sequence ids with leading 0s to match CoMet output
seq_ids = list(meta['seq_id'].values)
seq_ids = [seq_id.split('_')[-1] for seq_id in seq_ids]
seq_ids = [str(seq_id).zfill(9) + '_A' for seq_id in seq_ids]
meta['seq_id'] = seq_ids

# add 15th of month for collection dates without a day
# and remove year only dates
dates = meta['collection_date'].values
new_dates = []
for date in dates:
    if len(date) == 7:
        new_dates.append(date + '-15')
    elif len(date) == 9:
        new_dates.append(date[:8] + '0' + date[8])
    else:
        new_dates.append(date)

meta['collection_date'] = new_dates
good_date_idx = [len(date) > 4 for date in new_dates]
meta = meta.loc[good_date_idx, :]
assert np.all([len(date) == 10 for date in list(meta['collection_date'].values)])

# drop all columns other than seq_id and collection_date
# add wuhan sequence
# write to file
meta = meta.loc[:, ['seq_id', 'collection_date', 'location']]
wuhan = pd.DataFrame(data=[['000000000_A', '2019-12-15', 'Asia / China / Wuhan']], columns=['seq_id', 'collection_date', 'location'])
meta = pd.concat([wuhan, meta])
meta.to_csv(out_path, index=False, sep='\t')
