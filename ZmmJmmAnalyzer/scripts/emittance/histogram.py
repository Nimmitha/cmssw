import pandas as pd
import matplotlib.pyplot as plt
import os
import glob

# --- Step 1: Read and combine all CSVs ---
folder_path = "scan_plots"
csv_files = glob.glob(os.path.join(folder_path, "scan_*.csv"))

all_data = []
for file in csv_files:
    df = pd.read_csv(file)
    # use fill and run as keys
    df['key'] = df['fill'].astype(str) + "_" + df['run'].astype(str)
    all_data.append(df)

full_df = pd.concat(all_data, ignore_index=True)

# --- Step 2: Verify sep values are consistent across types ---
sep_values = full_df.groupby('type')['sep'].unique()
if len(sep_values['X1']) != 9 or len(sep_values['Y1']) != 9:
    raise ValueError("Mismatch in number of sep values for X1 and Y1 types.")

# --- Step 3: Prepare stacked data ---
def prepare_stack_data(df, scan_type):
    df_type = df[df['type'] == scan_type]
    pivot = df_type.pivot_table(
        index='sep',
        columns='key',
        values='eventCount',
        aggfunc='sum',
        fill_value=0
    ).sort_index()
    pivot.index = pivot.index.astype(str)  # Treat 'sep' as categorical labels
    return pivot

x_stack = prepare_stack_data(full_df, 'X1')
y_stack = prepare_stack_data(full_df, 'Y1')

x_total = x_stack.sum(axis=1)
y_total = y_stack.sum(axis=1)

x_stack.to_csv("x_stack.csv")
y_stack.to_csv("y_stack.csv")

# --- Step 4: Plot stacked bar plots ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

x_stack.plot(kind='bar', stacked=False, ax=axes[0], cmap='tab10', alpha=0.9)
y_stack.plot(kind='bar', stacked=False, ax=axes[1], cmap='tab10', alpha=0.9)

axes[0].set_title('X')
axes[1].set_title('Y')
axes[0].set_xlabel('Separation (mm)')
axes[1].set_xlabel('Separation (mm)')
axes[0].set_ylabel('Total Event Count')
axes[0].grid(True, axis='y')
axes[1].grid(True, axis='y')
axes[0].legend().remove()
axes[1].legend().remove()

# axes[1].legend(title="Scan", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("bar_plot1.png", dpi=300)

fig2, ax2 = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
x_total.plot(kind='bar', stacked=False, ax=ax2[0], color='tab:blue', alpha=0.9)
y_total.plot(kind='bar', stacked=False, ax=ax2[1], color='tab:orange', alpha=0.9)
ax2[0].set_title('X Total')
ax2[1].set_title('Y Total')
ax2[0].set_xlabel('Separation (mm)')
ax2[1].set_xlabel('Separation (mm)')
ax2[0].set_ylabel('Total Event Count')
ax2[0].grid(True, axis='y')
ax2[1].grid(True, axis='y')

plt.tight_layout()
plt.savefig("bar_plot2.png", dpi=300)
