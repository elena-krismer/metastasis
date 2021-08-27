# convert excel table metabolites (Recon1) to human one metabolites
import pandas as pd
import csv
import cobra
import pickle

"""
translation of IDs for clustering 
h2o2[x]' -> 'MAM02041p
"""


# save pickle object to obj folder
def save_obj(obj, name):
    with open('../../obj/dicts/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, protocol=2)


def load_obj(name):
    with open('../../obj/dicts/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


biggID_new_old_dict = {}
humanone_biggID_dict = {}
oldbiggID_humanone_dict = {}
biggID_humanone_dict = {}

humanone_metabolites = []


# ----------------------------------------------------------------------------------------------------------------------
# get IDs from model and IDs from table to translate
# ----------------------------------------------------------------------------------------------------------------------

# get used metabolite IDs from table (which will be translated) mostly recon1
old_metID_table = pd.read_csv("../../data/subsystems_tables/Supp_Tables_2_ID_translation.csv", delimiter=";")
oldmet_1 = old_metID_table['IN'].tolist()
oldmet_2 = old_metID_table['OUT'].tolist()
oldID_fromtable = list(set(oldmet_1 + oldmet_2))

# get metabolite IDs from model
humanone = cobra.io.load_matlab_model("/Users/s202425/Documents/GitHub/Human-GEM/model/Human-GEM.mat")
for x in humanone.metabolites:
    humanone_metabolites.append(x.id)


# ----------------------------------------------------------------------------------------------------------------------
# create dicts
# ----------------------------------------------------------------------------------------------------------------------


# file downloaded from Bigg Database
with open("../../data/misc/bigg_models_metabolites.tsv", 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        universalBIGG_ID = row[0]
        if len(row) > 4:
            # extract column which includes various old IDs
            old_bigg_ids_list = row[5].split(";")
            for old_id in old_bigg_ids_list:
                # find oldID which was used in table for clustering
                if old_id in oldID_fromtable:
                    biggID_new_old_dict[universalBIGG_ID] = old_id
csv_file.close()


# file downloaded fom Human GEM https://github.com/SysBioChalmers/Human-GEM/blob/main/model/metabolites.tsv
with open("../../data/misc/human_gem_metabolites.tsv", 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        for metabolite in humanone_metabolites:
            if metabolite in row:
                # get compartment (last character in metabolite) MAM01965c -> c
                compartment = metabolite[-1:]
                # metabolite = metabolite[:-1]
                biggID = row[2].strip()
                # add compartment to biggID
                biggID = biggID + "[" + compartment + "]"
                humanone_biggID_dict[metabolite] = biggID
                biggID_humanone_dict[biggID] = metabolite
csv_file.close()

# manually curated metabolites
dict_manual_addition = {'glc-D[c]': 'MAM01965c', 'h2o[x]': 'MAM02040p', 'o2[x]': 'MAM02630p', 'glu-L[c]': 'MAM01974c',
                        'urate[x]': 'MAM03120p', 'h2o2[x]': 'MAM02041p', 'nadh[x]': 'MAM02553p', 'h[x]': 'MAM02039p',
                        'lac-L[c]': 'MAM02403c', 'xu5p-D[c]': 'MAM01758s', 'mi1p-D[c]': 'MAM02493c',
                        'glc-D[e]': 'MAM01965e',
                        'ala_L[c]': 'MAM01307c', 'arg_L[c]': 'MAM01365c', 'citr_L[c]': 'None', 'asp_L[c]': 'MAM01370c',
                        'asn-L[c]': 'MAM01369c', 'asp_L[m]': 'MAM01370m', 'glu-L[m]': 'MAM01974m',
                        'cys_L[c]': 'MAM01628c',
                        'gln_L[c]': 'MAM01975c', 'pro_L[m]': 'MAM02770m', 'lac-L[m]': 'MAM02403m',
                        'hcys-L[c]': 'MAM02133c',
                        'saccrp-L[m]': 'MAM02868m', 'nh4[x]': 'MAM02579p', 'nad[x]': 'MAM02552c',
                        'ser_L[c]': 'MAM02896c', 'tyr_L[c]': 'MAM03101c',
                        'mev-R[x]': 'MAM00167p', 'hdcea[e]': 'None', 'tdchola[x]': 'MAM02962p', 'ppcoa[x]': 'MAM02774p',
                        'amp[x]': 'MAM01334p',
                        'ppi[x]': 'MAM02759m', 'dgchol[x]': 'None', 'fad[x]': 'MAM01802p', 'nadp[x]': 'MAM02554p',
                        'tchola[x]': 'None',
                        'gchola[e]': 'MAM01998e', 'dolp_U[r]': 'MAM01733r', 'm4mpdol_U[c]': 'MAM01323c',
                        'doldp_U[r]': 'MAM01734r', 'fuc-L[l]': 'None'}


# combine manually curated IDs with
biggID_humanone_dict = {**biggID_humanone_dict, **dict_manual_addition}
save_obj(biggID_humanone_dict, "biggID_humanone_dict")

# ----------------------------------------------------------------------------------------------------------------------
# translation of table with IDs
# ----------------------------------------------------------------------------------------------------------------------


# translate IDs
with open("../../data/subsystems_tables/Supp_Tables_2_ID_translation.csv", 'r') as csv_file:
    with open("../../data/subsystems_tables/Supp_Tables_2_ID_translation_translated.csv", "w") as output:
        for row in csv.reader(csv_file, delimiter=';'):
            # get metabolite ID
            translation_in = biggID_humanone_dict.get(row[4])
            translation_out = biggID_humanone_dict.get(row[7])
            # replace with metabolite ID from model
            row[4] = translation_in
            row[7] = translation_out
            # write two file
            row = [str(a) for a in row]
            output.write(str(';'.join(row) + '\n'))
output.close()
csv_file.close()
