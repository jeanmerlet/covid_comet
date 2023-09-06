import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

root_data_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/data/filtered_0.01N_1000D_counts/bin_size_10'
plot_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/plots/bp_dists'
primer_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/data/primers/artic_v3.tsv'

count_names = []
for r, d, f in os.walk(root_data_dir):
    for tsv in f:
        if re.search('.tsv$', tsv):
            count_names.append(os.path.join(r, tsv))

for i, count_file in enumerate(count_names):
    if i == 0:
        count_sums = pd.read_csv(count_file, sep='\t', index_col=0, header=0)
    else:
        counts = pd.read_csv(count_file, sep='\t', index_col=0, header=0)
        count_sums += counts

def totals():
    count_sums.columns = ['counts']
    count_sums = count_sums.sort_values(by=['counts'], ascending=False)
    fig, ax = plt.subplots(figsize=(12, 8))
    x = np.arange(count_sums.shape[0])
    ax.bar(x, count_sums.iloc[:, 0])
    ax.set_xticks(x)
    ax.set_xticklabels(count_sums.index.values)
    ax.set_xlabel('Nucleotide Code')
    ax.set_ylabel('(log) Count')
    ax.set_yscale('log')
    ax.set_title('Totals Across all Sequences')
    plt.savefig(os.path.join(plot_dir, 'total_counts.png'))

def get_primer_locs(path, bin_size):
    primers = pd.read_csv(path, sep='\t', header=None, index_col=None)
    forward_right = primers.iloc[::2, 2]
    reverse_left = primers.iloc[1::2, 1]
    return forward_right.values, reverse_left.values

def generate_per_nt_code_bar_plots(count_sums, bin_size=10, start=0, end=30000, primers=None):
    for i in range(count_sums.shape[1]):
        title = count_sums.columns[i].upper()
        if title != 'N': continue
        out_path = os.path.join(plot_dir, 'log_binned-' + title + f'-counts_bin-size-10_{start}-to-{end}-pos.png')
        #out_path = os.path.join(plot_dir, 'binned-' + title + f'-counts_bin-size-10_{start}-to-{end}-pos_for-slides.png')
        start = int(start / bin_size)
        end = int(end / bin_size)
        fig, ax = plt.subplots(figsize=(16, 8))
        # x axis
        x = np.arange(count_sums.shape[0])[start : end]
        bar = ax.bar(x, count_sums.iloc[:, i].values[start : end], width=1, align='edge')
        xticks = np.linspace(start, end, (end - start) // 4 + 1)
        #xticklabels = [str(int(xticks[i] * bin_size)) for i in range(len(xticks))]
        #ax.set_xlabel('Position')
        # if large numbers
        xticklabels = [str(round(xticks[i] / 100, 2)) for i in range(len(xticks))]
        #ax.set_xlabel('Position (thousands)')
        #ax.set_xticks(xticks)
        #ax.set_xticklabels(xticklabels, size=8)
        # y axis
        ax.set_ylabel('Number of sequences (log)', fontsize=28)
        ax.set_yscale('log')
        #ax.set_yticks(np.logspace(2, 7, 6))
        #ax.set_yticklabels([f'$10^{j}$' for j in range(2, 8)[2, 4, 7]], fontsize=21)
        logspace = np.logspace(1, 7, 7)
        ax.set_yticks([logspace[0], logspace[3], logspace[6]])
        ax.set_yticklabels([f'$10^{j}$' for j in [2, 4, 6]], fontsize=21)
        ax.tick_params(width=2, length=4)
        # remove ticks here
        ax.set_xticks([])
        #ax.set_yticks([])
        # other
        #ax.grid(which='both', axis='y')
        #ax.set_title(f'Binned Total for {title} Across all Sequences: {start * bin_size}bp to {end * bin_size}bp')
        if primers:
            for_locs, rev_locs = get_primer_locs(primers, 10)
            for j in range(len(for_locs)):
                if for_locs[j] > start * bin_size and for_locs[j] < end * bin_size:
                    bar[for_locs[j] // bin_size - start].set_color(u'#2ca02c')
                if rev_locs[j] > start * bin_size and rev_locs[j] < end * bin_size:
                    bar[rev_locs[j] // bin_size - start].set_color(u'#d62728')
        #out_path = os.path.join(plot_dir, 'log_binned-' + title + f'-counts_bin-size-10_{star.png')
        plt.savefig(out_path, dpi=300)
        plt.clf()

# bin_size corresponds to the size of the resolution of the bins of count_sums
# start and end correspond to actual positions in the genome (no bins)

#generate_per_nt_code_bar_plots(count_sums, bin_size=10, start=28400, end=30000, primers=primer_path)
generate_per_nt_code_bar_plots(count_sums)


