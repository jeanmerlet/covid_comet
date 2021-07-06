import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

root_data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/count_filtered_0.03/counts'
plot_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/plots'

count_names = []
for r, d, f in os.walk(root_data_dir):
    for tsv in f:
        if '.tsv' in tsv:
            count_names.append(os.path.join(r, tsv))

bin_width = 0.001
bins = np.arange(bin_width, 1 + 2*bin_width, bin_width)
col_names = [str(round(x, 4)) for x in bins]
all_d_props = pd.DataFrame(data=np.zeros((len(count_names), len(col_names)), dtype=np.int), columns=col_names)
all_n_props = pd.DataFrame(data=np.zeros((len(count_names), len(col_names)), dtype=np.int), columns=col_names)

for i, count_file in enumerate(count_names):
    props = pd.read_csv(count_file, sep='\t', index_col=None, header=0)
    d_props, d_counts = np.unique(props['d_count'], return_counts=True)
    n_props, n_counts = np.unique(props['n_count'], return_counts=True)
    for n in range(len(d_props)):
        col_idx = int(d_props[n] / bin_width)
        all_d_props.iloc[i, col_idx] += d_counts[n]
    for n in range(len(n_props)):
        col_idx = int(n_props[n] / bin_width)
        all_n_props.iloc[i, col_idx] += n_counts[n]

all_d_sums = np.sum(all_d_props, axis=0)
all_n_sums = np.sum(all_n_props, axis=0)
all_prop_sums = [all_d_sums, all_n_sums]
#print(all_d_sums[100:200])
#print(all_n_sums[100:200])

start = 0
end = 11
#start = 0
#end = len(bins)

def plot_manual_bins():
    titles = ['Deletions', 'Ns']
    fig, axes = plt.subplots(2, figsize=(12, 8))
    x = np.arange(len(bins))
    xtick_labels = np.tile(np.array('', dtype='<U6'), len(x))
    for i, label in enumerate(xtick_labels):
        if i % 50 == 0:
            xtick_labels[i] = col_names[i]

    xtick_labels[-1] = '1.01'

    for i, ax in enumerate(axes):
        ax.set_title(titles[i])
        ax.bar(x[start:end], all_prop_sums[i][start:end])
        ax.set_xticks(x[start:end])
        ax.set_xticklabels(xtick_labels[start:end])
        if i == 1:
            ax.set_xlabel('Proportion Ceiling')
        ax.set_ylabel('(log) Count')
        ax.set_yscale('log')
        #ax.set_ylabel('Count')

    #plt.tight_layout()
    #plt.savefig(os.path.join(plot_dir, 'log-filter-counts_' + str(bin_width) + '-bins.png'))

# difference in deletion or N proportion of sequences
total_seq_num = np.sum(all_d_sums)
cumulative_d = []
cumulative_n = []
do_print = False
for i in range(start, end):
    num_d_skipped = np.sum(all_d_sums[:i])
    num_n_skipped = np.sum(all_n_sums[:i])
    diff_d = total_seq_num - num_d_skipped
    diff_n = total_seq_num - num_n_skipped
    cumulative_d.append(diff_d)
    cumulative_n.append(diff_n)
    if do_print:
        pct_del = str(round(bins[i] * 100, 2))
        print(f'# sequences with greater than {pct_del}% deletions: {diff_d}')
        print(f'# sequences with greater than {pct_del}% Ns: {diff_n}')
        print('')

cums = [cumulative_d, cumulative_n]
titles = ['Deletions', 'Ns']
tick_interval = 1
x_ticks = np.arange(0, len(cumulative_d), tick_interval)
plot_bins = np.arange(0, 1 + bin_width, bin_width)
x_tick_labels = [str(round(plot_bins[x] * 100, 2)) for x in range(start, end, tick_interval)]

fig, axes = plt.subplots(2, figsize=(12, 12))
for i, ax in enumerate(axes):
    ax.plot(cums[i])
    ax.set_title(titles[i])
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_ylabel('(log) Count')
    ax.set_yscale('log')
    ax.set_yticks(np.logspace(3, 6, 4))
    ax.set_yticklabels([f'$10^{j}$' for j in range(3, 7)])
    #ax.set_ylabel('Count')
    if i == 1:
        ax.set_xlabel('Percent of Sequence')
    ax.grid(which='both')

fig.suptitle('Number of Sequences with Greater than % Nucleic Acid Code')
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(os.path.join(plot_dir, 'log-diff_' + str(bin_width) + '-bins_' + str(start) + '-' + str(end) + '.png'))
#plt.savefig(os.path.join(plot_dir, 'diff_' + str(bin_width) + '-bins_' + str(start) + '-' + str(end) + '.png'))


