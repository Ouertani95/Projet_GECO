import pandas as pd
import matplotlib.pyplot as plt
import os

# Get all files inside directory
Eggnog_results = os.listdir('Eggnog_results')
# Prepare dataframe
columns = ["With COGs", "Without COGs", "Total"]
statistics_df = pd.DataFrame(columns=columns)


def count_cog(cog_categories_dict):
    statistics = [0]*3
    for category, count in cog_categories_dict.items():
        if category == "-":
            statistics[1] = count
        else:
            statistics[0] += count
    statistics[2] = statistics[0]+statistics[1]
    return statistics


for result in Eggnog_results:
    # Read file
    table = pd.read_csv(f'Eggnog_results/{result}', delimiter="\t")
    # Count categories
    values_count = table["COG_category"].value_counts().to_dict()
    # Extract species
    split_file_name = result.split('_')
    title = f"{split_file_name[0]}_{split_file_name[1]}"
    # Append species statistics to dataframe
    new_row = pd.DataFrame([count_cog(values_count)], columns=columns, index=[title])
    statistics_df = pd.concat([statistics_df, new_row])
    # plot pie chart and save figures
    plt.pie(values_count.values(), labels=values_count.keys())
    plt.title(title)
    plt.savefig(f"Figures/{title}_COG_distribution.png")
    plt.show()

# Save statistics file
statistics_df.to_csv("COG_statistics.tsv", sep="\t")
