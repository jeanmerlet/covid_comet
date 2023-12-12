import pandas as pd
import numpy as np

meta_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/metadata/metadata_2022_06_02/metadata.tsv'
meta_out_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/metadata/metadata_2022_06_02/cali-counties/cali-counties_dates_locations.tsv'
epi_out_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/metadata/metadata_2022_06_02/cali-counties/cali-counties_epi-ids.txt'

meta = pd.read_csv(meta_path, sep='\t', usecols=[2, 3, 4], index_col=0, dtype=str)

row_idx = np.full(len(meta.index.values), False)
counter = 0
for i, seq_id in enumerate(meta.index.values):
    date, loc = meta.loc[seq_id, :].values
    split_loc = loc.split('/')
    if (
           len(split_loc) > 3 and split_loc[1].strip() == 'USA' and
           split_loc[2].strip() == 'California'
       ):
        if len(meta.loc[seq_id, 'Collection date']) >= 7:
            row_idx[i] = True

meta_out = meta.loc[row_idx, :]
meta_out.index = [epi_id[8:] for epi_id in meta_out.index.values]
meta_out.columns = ['date', 'loc']
meta_out.to_csv(meta_out_path, sep='\t')

epi = meta.loc[row_idx, :].index.to_frame(index=False)
epi.to_csv(epi_out_path, header=False, index=False)
