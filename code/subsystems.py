import pandas as pd
import pickle
import numpy as np
from functools import reduce


def load_obj(name):
    with open('../obj/dicts/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


subsystem_dict = load_obj("subsystem_dict")


def add_subystem(df):
    # correct weird things from excel csv convert
    df["J"] = df["J"].str.replace("E+", "e")
    df["J"] = df["J"].str.replace("+", "-")
    df["J"] = df["J"].str.replace(",", ".")
    df['J'] = df['J'].astype(float)

    reaction_list = df["reaction"].tolist()

    # add new column
    NaN = np.nan
    df["subsystem"] = NaN
    df["reaction"] = df["reaction"].str.replace("'", "")

    for reaction in reaction_list:
        reaction = reaction.strip("'")
        for key, val in subsystem_dict.items():
            if reaction in val:
                # replace NaN with metabolic subsystem
                df["subsystem"] = np.where(df["reaction"] == str(reaction), str(key), df["subsystem"])
    df_subsystem = df.drop("reaction", 1)
    df_subsystem = pd.pivot_table(df_subsystem, index=["subsystem"], values=["J"], aggfunc="sum")
    return df_subsystem



def add_subystem_2(df):
    reaction_list = df["reaction"].tolist()
    NaN = np.nan
    df["subsystem"] = NaN
    df["reaction"] = df["reaction"].str.replace("'", "")

    for reaction in reaction_list:
        reaction = reaction.strip("'")
        for key, val in subsystem_dict.items():
            if reaction in val:
                df["subsystem"] = np.where(df["reaction"] == str(reaction), str(key), df["subsystem"])
    df_subsystem = df.drop("reaction", 1)
    df_subsystem = pd.pivot_table(df_subsystem, index=["subsystem"], values=["J"], aggfunc="sum")
    return df_subsystem


df = pd.read_csv("../results/fva_jaccardsimilarity/lungtumour_vs_lungmetastasis(origin_breast).csv", delimiter=";")
lungtumour_vs_lungmetastasis_origin_breast = add_subystem(df)
lungtumour_vs_lungmetastasis_origin_breast.rename(columns={"J": "lungtumour_vs_lungmetastasis_origin_breast"},
                                                  inplace=True)

df = pd.read_csv("../results/fva_jaccardsimilarity/lungtumour_vs_brainmetasis.txt", delimiter=";")
lungtumour_vs_brainmetasis = add_subystem_2(df)
lungtumour_vs_brainmetasis.rename(columns={"J": "lungtumour_vs_brainmetasis"}, inplace=True)

df = pd.read_csv("../results/fva_jaccardsimilarity/MDA_MB_231_vs_nonlungmetastasis(origin_breast).csv", delimiter=";")
MDA_MB_231_vs_nonlungmetastasis_origin_breast = add_subystem(df)
MDA_MB_231_vs_nonlungmetastasis_origin_breast.rename(columns={"J": "MDA_MB_231_vs_nonlungmetastasis(origin_breast)"},
                                                     inplace=True)

df = pd.read_csv("../results/fva_jaccardsimilarity/MDA_MB_231_vs_lungmetastasis(origin_breast).csv", delimiter=";")
MDA_MB_231_vs_lungmetastasis_origin_breast = add_subystem(df)
MDA_MB_231_vs_lungmetastasis_origin_breast.rename(columns={"J": "MDA_MB_231_vs_lungmetastasis(origin_breast)"},
                                                  inplace=True)

df = pd.read_csv("../results/fva_jaccardsimilarity/nonlungmetastasis_vs_lungmetastasis(origin_breast).csv", ";")
nonlungmetastasis_vs_lungmetastasis_origin_breast = add_subystem(df)
nonlungmetastasis_vs_lungmetastasis_origin_breast.rename(
    columns={"J": "nonlungmetastasis_vs_lungmetastasis(origin_breast)"}, inplace=True)

data_frames = [lungtumour_vs_lungmetastasis_origin_breast, nonlungmetastasis_vs_lungmetastasis_origin_breast,
               lungtumour_vs_brainmetasis, MDA_MB_231_vs_nonlungmetastasis_origin_breast,
               MDA_MB_231_vs_lungmetastasis_origin_breast]

# merge all dataframes sum jaccard
df_merged = reduce(lambda left, right: pd.merge(left, right, on=["subsystem"], how = "outer"),
                   data_frames).fillna("void)")

df_merged.to_csv(r'/Users/s202425/Documents/GitHub/metastasis/results/fva_jaccardsimilarity/subsystems_jaccard_summed.csv')
