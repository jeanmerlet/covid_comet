import pandas as pd
import re

pd.set_option('display.max_rows', 500)

dates_loc_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/data/metadata/dates_location.tsv'
seqs_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/data/small_run_exploration/seqs_of_interest_10-3.fasta'

dates_loc = pd.read_csv(dates_loc_path, sep='\t', header=0, index_col=0)
dates_loc.index = [str(seq_id) for seq_id in dates_loc.index.values]

seqs = {}
line_num = 0
with open(seqs_path) as in_file:
    for line in in_file:
        line_num += 1
        if '>' in line:
            if line_num > 1:
                seq = ''.join(seq)
                seqs[epi_id] = seq
            epi_id = re.search('^>0+(\d+)00\/.*$', line).groups()[0]
            seq = []
            continue
        else:
            seq.append(line.strip())

dates_loc = dates_loc.loc[list(seqs.keys()), :]
dates_loc = dates_loc.sort_values(by=['collection_date'])
dates_loc['continent'] = [re.search('^([a-zA-Z]+)', location).groups()[0] for location in dates_loc['location']]
continuous_dates = []
for i in range(dates_loc.shape[0]):
    
datas_loc['
print(dates_loc)
