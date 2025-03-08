#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_csv = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(input_csv, index_col=0)

row_sums = df.sum(axis=1)
df = df.div(row_sums.replace(0, 1), axis=0)
df[row_sums == 0] = 0

plt.figure(figsize=(8, 8))
sns.heatmap(df, annot=True, cmap='Reds', fmt='.2f', cbar=False, annot_kws={"size": 18}, linewidths=0.5, linecolor='white')

plt.xlabel('Recipient Tissue', fontsize=22)
plt.ylabel('Source Tissue', fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

plt.tight_layout()
plt.savefig(outfile)
plt.close()
