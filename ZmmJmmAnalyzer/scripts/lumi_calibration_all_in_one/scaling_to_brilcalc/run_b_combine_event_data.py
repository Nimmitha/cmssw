import pandas as pd
import glob
import os


csv_dir = 'csvs'
method = 'fit'

def combine_csv(csv_dir, files, outFileName="combined_results.csv"):
    combined_csv = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    combined_csv.sort_values(by=['run', 'lumiblock'], inplace=True)
    outFilePath = os.path.join(csv_dir, outFileName)
    combined_csv.to_csv(outFilePath, index=False)
    print(f"Combined CSV saved to: {outFilePath}")


def main():
    all_filenames = [i for i in glob.glob(os.path.join(csv_dir, f'event_data_{method}_*.csv'))]

    if not all_filenames:
        raise FileNotFoundError(f"No CSV files found in directory: {csv_dir}")

    for year in [2022, 2023, 2024, 2025]:
        files = [i for i in all_filenames if f"_{year}_" in i]
        if not files:
            print(f"No CSV files found for year {year}. Please run the fitting script first.")
            continue
        print(f"\nCombining CSV {len(files)} files for year {year}...")


        combine_csv(csv_dir, files, outFileName=f'combined_data_{method}_{year}.csv')

    print(f"\nCombining CSV {len(all_filenames)} files for all years...")
    combine_csv(csv_dir, all_filenames, outFileName=f'combined_data_{method}_Run3.csv')


if __name__ == "__main__":
    main()