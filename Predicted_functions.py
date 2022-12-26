import os

import pandas
import pandas as pd

if  not os.path.exists('Predicted_functions'):
    os.mkdir('Predicted_functions')

eggnog_results = os.listdir('Eggnog_results/')

for result in eggnog_results:
    # Read file
    table = pd.read_csv(f'Eggnog_results/{result}', delimiter="\t")
    is_general = (table["COG_category"] == "-") | (table["COG_category"] == "S") | (table["COG_category"] == "R")
    table = table[~ is_general]
    # Count categories
    values_count = table["Description"].value_counts().to_dict()
    predicted_functions = pd.DataFrame(columns=["Predicted function", "Count"])
    for function, count in values_count.items():
        new_row = pandas.DataFrame([[function, count]], columns=["Predicted function", "Count"])
        predicted_functions = pd.concat([predicted_functions, new_row])
    # Extract species
    split_file_name = result.split('_')
    title = f"{split_file_name[0]}_{split_file_name[1]}"
    predicted_functions.to_csv(f"Predicted_functions/Predicted_functions_{title}.tsv", sep="\t", index=False)
