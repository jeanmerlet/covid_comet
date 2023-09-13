import os, re

path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/unaligned/sequences_2022_06_02.fasta'
out_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/unaligned/sequences_2022_06_02'

seq_set_size = 50000
# set this to number of completed sets in case this stops running partway through, otherwise set to 0
seq_set = 196
# if the process ends early, wc -l to find the number of written lines, otherwise set to 0
lines_written = 3640356631
with open(path) as in_fasta:
    num_seqs = 0
    for i, line in enumerate(in_fasta):
        if i < lines_written: continue
        if re.search('^>.*', line):
            if num_seqs % seq_set_size == 0:
                seq_set += 1
                out_path = os.path.join(out_dir, f'seq_set_{str(seq_set).zfill(3)}.fasta')
                print(f'on sequence set {str(seq_set).zfill(3)}')
                if seq_set > 197: # set to seq_set + 1
                    out_fasta.close()
                out_fasta = open(out_path, 'w')
            num_seqs += 1
        out_fasta.write(line)

out_fasta.close()
