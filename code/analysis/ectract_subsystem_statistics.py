import csv
import pandas as pd
import pickle


def load_obj(name):
    with open('../obj/dicts/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def subsystem_df(path, subsystem, cell_type):
    subsystem_dict = load_obj("subsystem_dict")
    subsystem_list = subsystem_dict.keys()
    subsystem_list = [*subsystem_dict]
    subsystem_list.append("MISC")
    flag = False
    with open(path, 'r') as file:
        min_cell_type = 'min' + cell_type
        max_cell_type = 'max' + cell_type
        df = pd.DataFrame(columns=['reaction', min_cell_type, 'min_breast',
                                   max_cell_type, 'max_breast', 'prop1', 'prop2'], index=[1, 2])
        for row in csv.reader(file, delimiter='\t'):
            if len(row) == 0:
                continue
            if row[0] != subsystem and row[0] in subsystem_list:
                flag = False
            if flag:
                reaction = row[0]
                min1, min2 = row[15], row[16]
                max1, max2 = row[17], row[18]
                prop1_value, prop2_value = row[19], row[20]
                print(reaction, min1, min2, max1, max2, prop1_value, prop2_value)
                df2 = pd.Series([reaction, min1, min2, max1, max2, prop1_value, prop2_value], index=df.columns)
                df = df.append(df2, ignore_index=True)
            if row[0] == subsystem:
                flag = True
        file.close()
        df.dropna(subset=["reaction"], inplace=True)
        return df


df_brain = subsystem_df('../../results/sampling_statisticalcomparison/subsystems_brainmetastasis_vs_breastcancer.tsv',
                        'Exchange/demand reactions', 'brain')
df_lung = subsystem_df('../../results/sampling_statisticalcomparison/subsystems_lungmetastasis_vs_breastcancer.tsv',
                       'Exchange/demand reactions', 'lung')

outer_merged = pd.merge(df_brain, df_lung, how="outer", on=["reaction"])
print(outer_merged.head())
print(len(df_brain))
print(len(df_lung))
print(len(outer_merged))
outer_merged.to_csv('../data/subsystems_tables/exchange_demand_reactions.csv', sep = ";")
