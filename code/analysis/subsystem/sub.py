import pandas as pd
import pickle
import numpy as np
from functools import reduce


def load_obj(name):
    with open('../../../obj/dicts/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


subsystem_dict = load_obj("subsystem_dict")

#
#df = pd.read_csv("../../../results/sampling_statisticalcomparison/lungmetastasis_vs_brainmetastasis.csv", delimiter=";")
df = pd.read_csv("../../../results/sampling_statisticalcomparison/lungmetastasis_vs_breastcancer.csv", delimiter=";")
#df = pd.read_csv("../../../results/sampling_statisticalcomparison/brainmetastasis_vs_breastcancer.csv", delimiter=";")
subsystems_list = subsystem_dict.keys()
used_reaction = []
print("DF", df.head(5))
print("ffff")
with open("../../../results/sampling_statisticalcomparison/subsystems_lungmetastasis_vs_breastcancer.tsv", "w") as output:
    columnnames = 'reaction;jaccardindex;KolmogorovSmirnow;KolmogorovSmirenow_pvalue;T - test_pvalue;Anova_pvalue;Wilcoxon_Mann_Whitney_pvalue;Mean;Std;Mode;Median;skew;kurt;Ttest_H;Ttest_pvalue;minFlux_lungmetastasis;minFlux_breastcancer;maxFlux_lungmetastasis;maxFlux_breastcancer;perc_brainmet/breastcancer;perc_breastcancer/brainmet\n'
    columnnames = columnnames.replace(';', '\t')
    output.write(columnnames)
    dataframe_reactions = df['reaction'].tolist()
    for subsystem in subsystems_list:
        output.write(str(subsystem) + '\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t\n')
        reaction_list = subsystem_dict.get(subsystem)
        print(reaction_list)
        for r in reaction_list:
            used_reaction.append(r)
            if r in dataframe_reactions:
                reaction_row = df.loc[df['reaction'] == r].values.flatten().tolist()
                reaction_row = [str(elem).replace(',', '.') for elem in reaction_row]
                reaction_row = reaction_row[0:21]
                output.write(str('\t'.join(reaction_row) + '\n'))
    output.write('MISC' + '\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t\n')

    used_reaction = set(used_reaction)
    missing_reactions = [x for x in dataframe_reactions if x not in used_reaction]
    for r in missing_reactions:
        if r == 'reaction':
            continue
        reaction_row = df.loc[df['reaction'] == r].values.flatten().tolist()
        reaction_row = [str(elem).replace(',', '.') for elem in reaction_row]
        reaction_row = reaction_row[0:21]
        output.write(str('\t'.join(reaction_row[0:21]) + '\n'))

output.close()