import pandas
import pandas as pd

A_hemo = 'Arcanobacterium_hemolyticum_DSM_20595_GCF_000092365.1_MM_k1jt2xw8.emapper.annotations.tsv'

# Read file
table = pd.read_csv(f'Eggnog_results/{A_hemo}', delimiter="\t")
is_general = (table["COG_category"] == "-") | (table["COG_category"] == "S") | (table["COG_category"] == "R")
table = table[~ is_general]
# Count categories
values_count = table["Description"].value_counts().to_dict()
predicted_functions = pd.DataFrame(columns=["Predicted function", "Count"])
for function, count in values_count.items():
    new_row = pandas.DataFrame([[function, count]], columns=["Predicted function", "Count"])
    predicted_functions = pd.concat([predicted_functions, new_row])
# Extract species
split_file_name = A_hemo.split('_')
title = f"{split_file_name[0]}_{split_file_name[1]}"
predicted_functions.to_csv(f"Predicted_functions_{title}.tsv", sep="\t", index=False)
