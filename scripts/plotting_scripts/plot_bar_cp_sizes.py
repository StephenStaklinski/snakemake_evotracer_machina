#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


csv_file = sys.argv[1]
outfile = sys.argv[2]


df = pd.read_csv(csv_file)

# Sort dataframe by 'size' column in descending order
df = df.sort_values(by='num_asvs', ascending=False)

# Create the bubble plot
plt.figure(figsize=(10,10))

sns.barplot(x='cp', y='num_asvs', data=df, palette="tab20")

plt.xticks(rotation=-90)

# Show the plot
plt.tight_layout()
plt.savefig(outfile)
plt.close()