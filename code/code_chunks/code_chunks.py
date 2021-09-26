import csv
import pandas as pd
import pickle


"""
list subsystems
"""
def list_subsystems():
    with open('/Users/s202425/Documents/GitHub/metastasis/obj/dicts/subsystem_dict.pkl', 'rb') as f:
        sub_dict =  pickle.load(f)
    return list(sub_dict.keys())

# get reactions from subsystem
def subsystem_reactions(subsystem):
    with open('/Users/s202425/Documents/GitHub/metastasis/obj/dicts/subsystem_dict.pkl', 'rb') as f:
        dict =  pickle.load(f)
    sub_list = dict.get(subsystem)
    return sub_list

"""
subsetting dataframe in subsystems for rstudio
"""
def subsystem_df(subsystem_list):
    with open("../../results/sampling_statisticalcomparison/", 'r') as file:
        df = pd.DataFrame(columns=['reaction', 'min_lung','min_brain' 'min_breast',
                                   'max_lung','max_brain', 'max_breast', 'prop_breast_lung', 'prop_breast_brain'],
                          index=[1, 2])
        for row in csv.reader(file, delimiter='\t'):
            if len(row) == 0:
                continue
            if row[0] in subsystem_list:
                reaction = row[0]
                min1, min2, min3 = row[3], row[4], row[5]
                max1, max2 ,max3 = row[6], row[7], row[8]
                prop1_value, prop2_value = row[19], row[20]
                print(reaction, min1, min2,min3, max1, max2, max3, prop2_value)
                df2 = pd.Series([reaction, min1, min2,min3, max1, max2, max3, prop1_value,
                                 prop2_value], index=df.columns)
                df = df.append(df2, ignore_index=True)
        file.close()
        df.dropna(subset=["reaction"], inplace=True)
        return df

def obsolete_subsystem_df(path, subsystem_list, cell_type):
    with open(path, 'r') as file:
        min_cell_type = 'min' + cell_type
        max_cell_type = 'max' + cell_type
        df = pd.DataFrame(columns=['reaction', min_cell_type, 'min_breast',
                                   max_cell_type, 'max_breast', 'prop1', 'prop2'], index=[1, 2])
        for row in csv.reader(file, delimiter='\t'):
            if len(row) == 0:
                continue
            if row[0] in subsystem_list:
                reaction = row[0]
                min1, min2 = row[15], row[16]
                max1, max2 = row[17], row[18]
                prop1_value, prop2_value = row[19], row[20]
                print(reaction, min1, min2, max1, max2, prop1_value, prop2_value)
                df2 = pd.Series([reaction, min1, min2, max1, max2, prop1_value,
                                 prop2_value], index=df.columns)
                df = df.append(df2, ignore_index=True)
        file.close()
        df.dropna(subset=["reaction"], inplace=True)
        return df


def dataframe_merge(df1, df2):
    outer_merged = pd.merge(df1, df2, how="outer", on=["reaction"])
    return outer_merged


def subset_merge(path1, path2, subsystem_list, celltype1, celltype2):
    df1 = subsystem_df(path1, subsystem_list, celltype1)
    df2 = subsystem_df(path2, subsystem_list, celltype2)
    df = dataframe_merge(df1, df2)
    return df


import urllib.request as urllib2
import re
import pickle
import csv

"""
load list from pickle for Rstudio
"""


# save pickle object to obj folder
def save_obj(obj, name):
    with open('../../obj/dicts/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, protocol=2)


def load_obj(name):
    with open('/Users/s202425/Documents/GitHub/metastasis/obj/dicts/brain_protein_expression_dict.pkl', 'rb') as f:
        return pickle.load(f)


# get list of low or high proteins from brain
def load_gene_list_brain(level):
    brain_gene_expression_dict = load_obj('brain_protein_expression_dict')
    if level == "low":
        list = brain_gene_expression_dict.get("low_expression")
    else:
        list = brain_gene_expression_dict.get("high_expression")
    return list

# get reactions from subsystem
def subsystem_reactions(subsystem):
    with open('/Users/s202425/Documents/GitHub/metastasis/obj/dicts/subsystem_dict.pkl', 'rb') as f:
        dict =  pickle.load(f)
    sub_list = dict.get(subsystem)
    return sub_list


# level = low/high expressed proteins translated to gene ID
def load_brain_dataframe(level):
    df = pd.read_csv('/Users/s202425/Documents/GitHub/metastasis/data/gene_expression/GSE11078.txt', delimiter="\t")
    df = df[df.ID_REF.str.contains("ENSG")]
    # Breast Cancer Met-1/8/9/10
    geneID = df.iloc[:, 0]
    # ID with brain expression data
    brain = ['GSM279958', 'GSM279969', 'GSM279971', 'GSM279972']
    # subset dataframe with all expression darta
    df_brain = df.loc[:, brain].copy()
    #mean = df_brain1.mean(axis=1)
    df_brain2 = pd.DataFrame({'geneID': geneID, 'value1': df_brain.iloc[:, 0], 'value2': df_brain.iloc[:, 1],
                              'value3': df_brain.iloc[:, 2], 'value4': df_brain.iloc[:, 3]})
    # get list of high or low expressed proteins in brain metastasis (translated to gene IDs)
    gene_list = load_gene_list_brain(level)
    # subset dataframe with list
    df_subset = df_brain2[df_brain2.geneID.isin(gene_list)]
    return df_subset

load_brain_dataframe("low")
