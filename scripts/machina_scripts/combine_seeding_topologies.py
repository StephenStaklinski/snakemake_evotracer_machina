#!/usr/bin/env python3

import sys
import os
import pandas as pd

    
input_csvs = str(sys.argv[1])
output_csv = sys.argv[2]
normalized = bool(sys.argv[3])


# Read the input CSV file
df_list = [pd.read_csv(csv) for csv in input_csvs.split(',') if csv != "" and os.path.getsize(csv) > 0]
df = pd.concat(df_list, ignore_index=True)

grouped_df = df.groupby('CP').sum().reset_index()
if normalized:
    grouped_df.iloc[:, 1:] = grouped_df.iloc[:, 1:].div(grouped_df.iloc[:, 1:].sum(axis=1), axis=0)
if len(grouped_df['CP'].unique()) == 1:
    cp = grouped_df['CP'].unique()[0]
    grouped_df = grouped_df.drop(columns=['CP'])
    grouped_df = grouped_df.mean().to_frame().T
    grouped_df['CP'] = cp
    cols = grouped_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    grouped_df = grouped_df[cols]
else:
    grouped_df = grouped_df.drop(columns=['CP'])
    grouped_df = grouped_df.mean().to_frame().T
    grouped_df['CP'] = 'All'
    cols = grouped_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    grouped_df = grouped_df[cols]

# Write the resulting DataFrame to the output CSV file
grouped_df.to_csv(output_csv, index=False)