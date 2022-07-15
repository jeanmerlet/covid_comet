import os
import pandas as pd

metadata_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/metadata.tsv'
out_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/non-human_hosts.tsv'

meta = pd.read_csv(metadata_path, sep='\t', header=0, index_col=None)
hosts = meta.loc[:, ['Accession ID', 'Host']]
hosts = hosts.loc[hosts['Host'] != 'Human', :]
hosts.to_csv(out_path, sep='\t', index=False, header=True)
