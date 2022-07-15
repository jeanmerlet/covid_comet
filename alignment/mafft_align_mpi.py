from mpi4py import MPI
import subprocess
import time
import os

fasta_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/unaligned/sequences_2022_06_02'
out_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02'
ref_fasta = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/unaligned/wuhan_ref.fasta'
mafft_bin = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/rna-seq_tools/alignment/mafft/andes/mafft_7.487/bin/mafft'

fasta_paths = [os.path.join(fasta_dir, fasta) for fasta in os.listdir(fasta_dir)]

def mafft_align(mafft_bin, fasta, ref_fasta, fasta_name):
    out_path = os.path.join(out_dir, 'aligned_' + fasta_name)
    with open(out_path, 'w') as out_file:
        subprocess.run([mafft_bin, '--auto', '--nomemsave', '--keeplength',
                        '--addfragments', fasta, ref_fasta], stdout=out_file) 
    
comm =  MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i, fasta in enumerate(fasta_paths):
    # distribute across ranks
    if i % size != rank: continue
    _, fasta_name = os.path.split(fasta)
    # time the first one
    if i == 0: start_time = time.time()
    mafft_align(mafft_bin, fasta, ref_fasta, fasta_name)
    # print time for the first one
    if i == 0: print(f'TOTAL TIME for first set: {time.time() - start_time}')
