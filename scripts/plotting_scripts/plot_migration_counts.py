#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# input_csv = sys.argv[1]
# output_file = sys.argv[2]

input_csv = '/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2/migration_counts.csv'
output_file = '/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2/migration_counts.pdf'

data = pd.read_csv(input_csv)

p_mice = ["MMUS7351", "MMUS7317", "MMUS7358", "MMUS7323", "MMUS7354"]
pr_mice = ["MMUS7438", "MMUS7439", "MMUS7440", "MMUS7447", "MMUS7490", "MMUS7507"]

# Combine the lists to set the order
mouse_order = p_mice + pr_mice

# Compute the average number of migrations per mouse
avg_migrations_per_mouse = data.groupby('mouse')['num_migrations'].mean().reset_index()
# avg_migrations_per_mouse = data.groupby('mouse')['num_migrations'].sum().reset_index()

# Separate the data into two groups
p_mice_data = avg_migrations_per_mouse[avg_migrations_per_mouse['mouse'].isin(p_mice)]['num_migrations']
pr_mice_data = avg_migrations_per_mouse[avg_migrations_per_mouse['mouse'].isin(pr_mice)]['num_migrations']

# Compute the Mann-Whitney U test on the averages
stat, p_value = mannwhitneyu(p_mice_data, pr_mice_data)
print(f'Mann-Whitney U test p-value: {p_value}')

# Create the plot
plt.figure(figsize=(12, 6))
fs=18

# Plot individual mice data
ax1 = plt.subplot(1, 2, 1)
sns.boxplot(x='mouse', y='num_migrations', data=data, order=mouse_order, palette=['blue']*len(p_mice) + ['green']*len(pr_mice), ax=ax1, showfliers=False, boxprops=dict(alpha=0.5))
sns.stripplot(x='mouse', y='num_migrations', data=data, jitter=True, dodge=True, marker='o', alpha=1.0, color='grey', order=mouse_order, ax=ax1)
ax1.set_title('', fontsize=fs)
ax1.set_xlabel('', fontsize=fs)
ax1.set_ylabel('Number of migrations\nper CP', fontsize=fs)
ax1.tick_params(axis='both', which='major', labelsize=fs)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=-90)
# ax1.set_ylim(0, None)

# Plot overall comparison using averages
ax2 = plt.subplot(1, 2, 2)
overall_data = pd.DataFrame({
    'group': ['2P'] * len(p_mice_data) + ['2PR'] * len(pr_mice_data),
    'num_migrations': pd.concat([p_mice_data, pr_mice_data])
})
sns.boxplot(x='group', y='num_migrations', data=overall_data, palette=['blue', 'green'], ax=ax2, showfliers=False, boxprops=dict(alpha=0.5), width=0.5)
sns.stripplot(x='group', y='num_migrations', data=overall_data, jitter=True, dodge=True, marker='o', alpha=1.0, color='grey', ax=ax2)
ax2.set_title(f'Overall Comparison\nMann-Whitney U test p-value: {p_value:.6f}', fontsize=fs)
ax2.set_xlabel('', fontsize=fs)
ax2.set_ylabel('Average number of migrations\nper mouse', fontsize=fs)
ax2.tick_params(axis='both', which='major', labelsize=fs)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=-90)
# ax2.set_ylim(0, None)

# Save the plot to the output file
plt.tight_layout()
plt.savefig(output_file)
plt.close()