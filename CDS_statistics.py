import pandas as pd
import os

# Get all files inside directory
feature_tables = os.listdir('Feature_tables')
# Prepare dataframe
columns = ["number CDS", "Total CDS length", "Mean CDS length", "Max CDS length"]
statistics_df = pd.DataFrame(columns=columns)


def count_cds(dataframe):
    statistics = [0]*4
    statistics[0] = len(dataframe.index)
    statistics[1] = dataframe['feature_interval_length'].sum()
    statistics[2] = round(dataframe['feature_interval_length'].mean(), 2)
    statistics[3] = dataframe['feature_interval_length'].max()
    return statistics


for feature_table in feature_tables:
    # Extract species
    split_file_name = feature_table.split('_')
    title = f"{split_file_name[0]}_{split_file_name[1]}"
    # Read file
    table = pd.read_csv(f'Feature_tables/{feature_table}', delimiter="\t")
    is_CDS = table['# feature'] == 'CDS'
    table = table[is_CDS]
    # Append species statistics to dataframe
    new_row = pd.DataFrame([count_cds(table)], columns=columns, index=[title])
    statistics_df = pd.concat([statistics_df, new_row])


# Save statistics file
statistics_df.to_csv("CDS_statistics.tsv", sep="\t")
