import pandas as pd
from collections import Counter
import pickle


def load_obj(name):
    with open('../obj/dicts/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


subsystem_dict = load_obj("subsystem_dict")

df = pd.read_csv("../results/fva_jaccardsimilarity/lungtumour_vs_brainmetasis.txt", delimiter=";")
df = df[df["J"] > 0.1]

subsystem_list = []
reaction_list = df["reaction"].tolist()

for reaction in reaction_list:
    reaction = reaction.strip("'")
    for key, val in subsystem_dict.items():
        if reaction in val:
            subsystem_list.append(key)

print(Counter(subsystem_list))
# {'Transport reactions': 41, 'Exchange/demand reactions': 35, 'Pentose phosphate pathway': 10,
# 'Glycolysis / Gluconeogenesis': 9, 'Purine metabolism': 8, 'Starch and sucrose metabolism': 7,
# 'Folate metabolism': 4, 'Carnitine shuttle (cytosolic)': 3, 'Bile acid biosynthesis': 3, 'Galactose metabolism': 3,
# 'Pyrimidine metabolism': 3, 'Miscellaneous': 3, 'Peptide metabolism': 2, 'ROS detoxification': 2,
# 'Nucleotide metabolism': 2, 'Fatty acid oxidation': 2, 'Fatty acid biosynthesis': 1,
# 'Fatty acid activation (cytosolic)': 1, 'Arginine and proline metabolism': 1,
# 'Glycine, serine and threonine metabolism': 1, 'Pentose and glucuronate interconversions': 1,
# 'Histidine metabolism': 1,
# 'Cysteine and methionine metabolism': 1}

df = pd.read_csv("../results/fva_jaccardsimilarity/lungtumour_vs_lungmetastasis(origin_breast).csv", delimiter=";")

# correct weird things from excel csv convert
df["J"] = df["J"].str.replace("E+", "e")
df["J"] = df["J"].str.replace("+", "-")
df["J"] = df["J"].str.replace(",", ".")
df['J'] = df['J'].astype(float)
df = df[df["J"] > 0.1]

subsystem_list = []
reaction_list = df["reaction"].tolist()

for reaction in reaction_list:
    reaction = reaction.strip("'")
    for key, val in subsystem_dict.items():
        if reaction in val:
            subsystem_list.append(key)

print(Counter(subsystem_list))
# {'Transport reactions': 12, 'Exchange/demand reactions': 7, 'Glycolysis / Gluconeogenesis': 6,
# 'Purine metabolism': 5, 'Pentose phosphate pathway': 4, 'Fatty acid oxidation': 4, 'Bile acid biosynthesis': 3,
# 'Pyrimidine metabolism': 3, 'Fatty acid biosynthesis': 2, 'Carnitine shuttle (cytosolic)': 2,
# 'Glycine, serine and threonine metabolism': 1}

df = pd.read_csv("../results/fva_jaccardsimilarity/nonlungmetastasis_vs_lungmetastasis(origin_breast).csv", delimiter=";")

# correct weird things from excel csv convert
df["J"] = df["J"].str.replace("E+", "e")
df["J"] = df["J"].str.replace("+", "-")
df["J"] = df["J"].str.replace(",", ".")
df['J'] = df['J'].astype(float)
df = df[df["J"] > 0.1]

subsystem_list = []
reaction_list = df["reaction"].tolist()

for reaction in reaction_list:
    reaction = reaction.strip("'")
    for key, val in subsystem_dict.items():
        if reaction in val:
            subsystem_list.append(key)

print(Counter(subsystem_list))
# {'Transport reactions': 27, 'Exchange/demand reactions': 24, 'Glycolysis / Gluconeogenesis': 8,
# 'Purine metabolism': 6, 'Pentose phosphate pathway': 5, 'Fatty acid oxidation': 4, 'Fatty acid biosynthesis': 2,
# 'Carnitine shuttle (cytosolic)': 2, 'Cysteine and methionine metabolism': 2, 'Folate metabolism': 2,
# 'Nucleotide metabolism': 2, 'Fatty acid activation (cytosolic)': 1, 'ROS detoxification': 1,
# 'Arginine and proline metabolism': 1, 'Pentose and glucuronate interconversions': 1, 'Histidine metabolism': 1,
# 'Oxidative phosphorylation': 1}


df = pd.read_csv("../results/fva_jaccardsimilarity/MDA_MB_231_vs_lungmetastasis(origin_breast).csv", delimiter=";")

# correct weird things from excel csv convert
df["J"] = df["J"].str.replace("E+", "e")
df["J"] = df["J"].str.replace("+", "-")
df["J"] = df["J"].str.replace(",", ".")
df['J'] = df['J'].astype(float)
df = df[df["J"] > 0.1]

subsystem_list = []
reaction_list = df["reaction"].tolist()

for reaction in reaction_list:
    reaction = reaction.strip("'")
    for key, val in subsystem_dict.items():
        if reaction in val:
            subsystem_list.append(key)

print(Counter(subsystem_list))
# {'Transport reactions': 22, 'Exchange/demand reactions': 19, 'Purine metabolism': 7,
# 'Glycolysis / Gluconeogenesis': 5, 'Pyrimidine metabolism': 3, 'Fatty acid oxidation': 3,
# 'Fatty acid activation (cytosolic)': 2, 'Carnitine shuttle (cytosolic)': 2, 'Cysteine and methionine metabolism': 2,
# 'Pentose phosphate pathway': 2, 'Nucleotide metabolism': 2, 'Fatty acid biosynthesis': 1, 'ROS detoxification': 1,
# 'Arginine and proline metabolism': 1, 'Oxidative phosphorylation': 1})

df = pd.read_csv("../results/fva_jaccardsimilarity/MDA_MB_231_vs_nonlungmetastasis(origin_breast).csv", delimiter=";")

# correct weird things from excel csv convert
df["J"] = df["J"].str.replace("E+", "e")
df["J"] = df["J"].str.replace("+", "-")
df["J"] = df["J"].str.replace(",", ".")
df['J'] = df['J'].astype(float)
df = df[df["J"] > 0.1]

subsystem_list = []
reaction_list = df["reaction"].tolist()

for reaction in reaction_list:
    reaction = reaction.strip("'")
    for key, val in subsystem_dict.items():
        if reaction in val:
            subsystem_list.append(key)

print(Counter(subsystem_list))
# {'Transport reactions': 21, 'Exchange/demand reactions': 18, 'Glycolysis / Gluconeogenesis': 5,
# 'Purine metabolism': 4, 'Fatty acid oxidation': 3, 'Carnitine shuttle (cytosolic)': 2,
# 'Cysteine and methionine metabolism': 2, 'Fatty acid biosynthesis': 1, 'Fatty acid activation (cytosolic)': 1,
# 'Sphingolipid metabolism': 1, 'ROS detoxification': 1, 'Arginine and proline metabolism': 1,
# 'Pentose phosphate pathway': 1, 'Oxidative phosphorylation': 1}