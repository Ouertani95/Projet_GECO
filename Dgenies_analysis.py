import os

import pandas as pd
from matplotlib import pyplot as plt

file_dict = {"1vs1":"A.haemolyticum vs A.haemolyticum",
             "1vs2":"A.haemolyticum vs A.phocisimile",
             "1vs3":"A.haemolyticum vs A.pinnipediorum",
             "1vs4":"A.haemolyticum vs A.suis",
             "1vs5":"A.haemolyticum vs Trueperella pecoris"}

if  not os.path.exists('Dgenies_analysis_results'):
    os.mkdir('Dgenies_analysis_results')

dgenies_results = []
for result in os.listdir('Dgenies_results/'):
    if result.endswith('.paf'):
        dgenies_results.append(result)

for result in dgenies_results:
    # Read file
    table = pd.read_csv(f'Dgenies_results/{result}', delimiter="\t")
    table = table.drop(table.iloc[:, 11:], axis=1)
    table = table.drop(table.iloc[:, 0:4], axis=1)
    table = table.drop(table.iloc[:, 1:5], axis=1)
    table.columns = ['side', 'match_count','total_count']
    table["side"].value_counts()\
        .to_csv(f"Dgenies_analysis_results/{result.split('.')[0]}.count",header=False)
    plt.hist(table['total_count'], bins=100, color='b')
    plt.gca().set(title=f'Distrubution des longueurs totales\n'
                        f'des blocs de synténie pour {file_dict.get(result.split("_")[1].split(".")[0])}'
                  , ylabel='Fréquence', xlabel='Longueur totale bloc de synténie')
    plt.savefig(f'Dgenies_analysis_results/{result.split(".")[0]}_total_synteny.png')
    plt.show()
    plt.hist(table['match_count'], bins=100, color='r')
    plt.gca().set(title=f'Distrubution des longueurs de match\n'
                        f'des blocs de synténie pour {file_dict.get(result.split("_")[1].split(".")[0])}'
                  , ylabel='Fréquence', xlabel='Longueur match bloc de synténie')
    plt.savefig(f'Dgenies_analysis_results/{result.split(".")[0]}_match.png')
    plt.show()