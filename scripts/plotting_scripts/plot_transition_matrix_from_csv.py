#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_csv = sys.argv[1]
outfile = sys.argv[2]


df = pd.read_csv(input_csv, index_col=0)

df = df.div(df.sum(axis=1), axis=0)

plt.figure(figsize=(8, 8))
sns.heatmap(df, annot=True, cmap='Reds', fmt='.2f', cbar=False, annot_kws={"size": 18})

plt.xlabel('Recipient Tissue', fontsize=22)
plt.ylabel('Source Tissue', fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)

plt.tight_layout()
plt.savefig(outfile)
plt.close()
