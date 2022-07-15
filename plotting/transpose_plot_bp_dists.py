import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os, re

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/transposed/binned'
plot_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/plots/bp_dists/transposed/matched_pos'
snp_pos_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2021_09_30/nt_pos.tsv'

by_pos = True
if by_pos:
    snp_pos = pd.read_csv(snp_pos_path, sep='\t', index_col=None, header=None)
    snp_pos = np.squeeze(snp_pos.values) - 342
    print(snp_pos)



count_paths = [os.path.join(data_dir, path) for path in os.listdir(data_dir)]
count_paths.sort()

for i, count_path in enumerate(count_paths):
    if i == 0:
        count_sums = pd.read_csv(count_path, sep='\t', index_col=0, header=0)
    else:
        counts = pd.read_csv(count_path, sep='\t', index_col=0, header=0)
        count_sums += counts

count_sums = count_sums.loc[:, ['-', 'n']]
if by_pos:
    count_sums = count_sums.loc[snp_pos, :]
    count_sums = count_sums.reset_index(drop=True)
    print(count_sums)

def get_matched_n_and_d_spikes(count_sums, min_count, by_pos):
    d_above = count_sums.loc[:, '-'] > min_count
    n_above = count_sums.loc[:, 'n'] > min_count
    print(np.count_nonzero(d_above))
    print(np.count_nonzero(n_above))
    print(np.count_nonzero(n_above.loc[d_above.values]))
    idx = np.logical_and(d_above, n_above)
    return count_sums.index[idx].values

def generate_per_nt_code_bar_plots(count_sums, log, grid, matched_pos, min_count, y_max, by_pos):
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    if by_pos:
        pos_name = 'by_pos'
    else:
        pos_name = ''
    if log:
        if grid:
            out_path = os.path.join(plot_dir, f'{pos_name}_min_count-{min_count}_' + 'matched_log_grid_snp_counts.png')
        else:
            out_path = os.path.join(plot_dir, f'{pos_name}_min_count-{min_count}_' + 'matched_log_snp_counts.png')
    else:
        if grid:
            out_path = os.path.join(plot_dir, f'{pos_name}_min_count-{min_count}_' + 'matched_grid_snp_counts.png')
        else:
            out_path = os.path.join(plot_dir, f'{pos_name}_min_count-{min_count}_' + 'matched_snp_counts_new.png')
    for i in range(count_sums.shape[1]):
        ax = axes[i]
        title = count_sums.columns[i].upper()
        ax.set_title(f'Summed Counts of {title}')
        print(title)
        ## x axis ##
        x = np.arange(count_sums.shape[0])
        bar = ax.bar(x, count_sums.iloc[:, i].values, width=1, align='edge')
        #xticks = np.linspace(0, count_sums.shape[0], count_sums.shape[0] // 1000 + 1)
        #xticklabels = [str(int(xticks[i] * bin_size)) for i in range(len(xticks))]
        # if large numbers
        #xticklabels = [str(round(xticks[i] / 100, 2)) for i in range(len(xticks))]
        ax.set_xlabel('Position')
        #ax.set_xticks(xticks)
        #ax.set_xticklabels(xticklabels, size=8)
        ## y axis ##
        if log:
            ax.set_ylabel('(log) Count')
            ax.set_yscale('log')
            ax.set_yticks(np.logspace(2, 7, 6))
            ax.set_yticklabels([f'$10^{j}$' for j in range(2, 8)])
        else:
            ax.set_ylabel('Count')
        y = np.arange(0, y_max + 2000, 2000)
        ax.set_yticks(y)
        ax.set_yticklabels([str(tick)[:-3] + 'k' for tick in y])
        ax.set_ylim([0, y_max])
        ## other ##
        if grid:
            ax.grid(which=grid, axis='y')
        if min_count != 'none':
            for i, pos in enumerate(matched_pos):
                if i == 0:
                    ax.axvline(pos, color='red', label='matched spike')
                else:
                    ax.axvline(pos, color='red')
            ax.legend(fontsize=14)
        plt.savefig(out_path, dpi=300)


min_count = 100000
matched_pos = get_matched_n_and_d_spikes(count_sums, min_count, by_pos)

log = False
grid = 'major'
#matched_pos = []
#min_count = 'none'
y_max = 20000
generate_per_nt_code_bar_plots(count_sums, log, grid, matched_pos, min_count, y_max, by_pos)
