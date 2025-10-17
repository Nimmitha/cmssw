import os
import pandas as pd
import ast
import numpy as np
import glob

import matplotlib.pyplot as plt

def load_and_combine_csvs(directory):
    """Load all CSVs and combine into a single dataframe grouped by (X, Y)."""
    csv_files = glob.glob(os.path.join(directory, "scan_*.csv"))
    all_data = []

    for file in csv_files:
        df = pd.read_csv(file)
        # Ensure list column is parsed correctly
        df['B_J1_mass_list'] = df['B_J1_mass_list'].apply(ast.literal_eval)
        all_data.append(df)

    combined_df = pd.concat(all_data, ignore_index=True)

    temp = combined_df[(combined_df['type'] == 'X1') & (combined_df['sep'] == 0)]
    temp['rate'] = temp['eventCount'] / temp['length']

    print(temp)
    print(temp['rate'].mean())
    print(temp['rate'].std())

    # Group and apply custom aggregation
    grouped = combined_df.groupby(['type', 'sep']).agg({
        'length': 'sum',                                             # simple sum
        'eventCount': 'sum',                                         # sum of event counts
        'B_J1_mass_list': lambda lists: sum(lists, []),              # concatenate lists
    }).reset_index()

    grouped['eventCount'] = grouped['B_J1_mass_list'].apply(len)  # count total events

    return grouped


# Load and combine data
combined_df = load_and_combine_csvs("output/csvs")
print(combined_df)
combined_df.to_csv("output/csvs/combined_scan_count.csv", index=False)
