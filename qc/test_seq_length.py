import pandas as pd
import os

correct_len = 29325
in_dir = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/all-usa'
paths = [os.path.join(in_dir, path) for path in os.listdir(in_dir) if path[-4:] == '.tsv']
out_dir = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/all-usa/seq_len_check'

def check_record_seq_lens(path, correct_len):
    bad_seqs, bad_lens = [], []
    with open(path) as seq_file:
        for i, line in enumerate(seq_file):
            seq_len = len(line.split('\t'))
            if seq_len != correct_len:
                bad_seqs.append(i)
                bad_lens.append(seq_len)
    return bad_seqs, bad_lens

def write_bad_seqs(out_path, bad_seqs, bad_lens):
    df = pd.DataFrame(data=bad_lens, index=bad_seqs, columns=['seq_len'])
    df.to_csv(out_path, sep='\t')

for i, path in enumerate(paths):
    _, tail = os.path.split(path)
    out_path = os.path.join(out_dir, tail[:-4] + '_seq-len_check.tsv')
    bad_seqs, bad_lens = check_record_seq_lens(path, correct_len)
    if len(bad_seqs) > 0:
        write_bad_seqs(out_path, bad_seqs, bad_lens)
