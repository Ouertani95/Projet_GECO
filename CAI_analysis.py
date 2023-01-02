import os

import pandas as pd
from matplotlib import pyplot as plt

if  not os.path.exists('CAI_analysis_results'):
    os.mkdir('CAI_analysis_results')

cai_results = os.listdir('CAI_Emboss_results/')

for result in cai_results:
    # Read file
    table = pd.read_csv(f'CAI_Emboss_results/{result}', delimiter=": ")
    table = table.drop(table.iloc[:, 0:2], axis=1)
    table.to_csv(f'CAI_analysis_results/{result}', index=False)
    name = result.split("_")[0].split(".")[1] + " " + result.split("_")[1]
    plt.hist(table, bins=100, color='g')
    plt.gca().set(title=f'Distrubution des valeurs de CAI pour {name}', ylabel='Fréquence', xlabel='CAI')
    plt.xlim(0, 1)
    plt.savefig(f'CAI_analysis_results/{result.split("_")[0].split(".")[0]}.{name}.png')
    plt.show()

