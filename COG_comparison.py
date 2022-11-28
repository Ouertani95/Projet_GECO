import os

import pandas as pd
import matplotlib.pyplot as plt

# Get all files inside directory
cog_tables = os.listdir('COG_comparison')

for cog_table in cog_tables:
    # Read file
    table = pd.read_csv(f'COG_comparison/{cog_table}', delimiter="\t")
    # Count categories
    values_count = table["COG_category"].value_counts().to_dict()
    # Extract species
    split_file_name = cog_table.split('_')
    title = f"{split_file_name[0]}_{split_file_name[1]}_{split_file_name[2]}"

    plt.pie(values_count.values(), labels=values_count.keys())
    plt.title(title)
    plt.savefig(f"Figures/{title}_COG_distribution.png")
    plt.show()
