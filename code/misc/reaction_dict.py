import pickle
import csv
import pandas as pd

# save pickle object to obj folder
def save_obj(obj, name):
    with open('../../obj/dicts/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, protocol=2)

"""
bigg_reaction_details_dict = {}
# file downloaded from Bigg Database
with open("../../data/misc/bigg_models_reactions.tsv", 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        if len(row) > 0:
            universalBIGG_ID = row[0]
            reaction_details = row[1]
            bigg_reaction_details_dict[universalBIGG_ID] = reaction_details
csv_file.close()

bigg_reactions_list = bigg_reaction_details_dict.keys()

reaction_description_dict = {}
# file downloaded fom Human GEM https://github.com/SysBioChalmers/Human-GEM/blob/main/model/metabolites.tsv
with open("../../data/misc/human_gem_reactions.tsv", 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        if row[2].strip() in bigg_reactions_list:
            biggID = row[2].strip()
            humanone_reaction = row[0].strip()
            reaction_description  = bigg_reaction_details_dict.get(biggID)
            reaction_description_dict[humanone_reaction] = reaction_description
csv_file.close()

save_obj(reaction_description_dict, "reaction_description_dict")
"""
humanone_kegg_dict = {}
with open("../../data/misc/human_gem_reactions.tsv", 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        keggID = row[1].strip()
        humanone_reaction = row[0].strip()
        humanone_kegg_dict[humanone_reaction] = keggID
csv_file.close()

print(humanone_kegg_dict)

#save_obj(humanone_kegg_dict, "humanone_kegg_reaction_dict")
