import pickle
import csv
import pandas as pd

def load_obj(name):
    with open('../../obj/dicts/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
ID_dict = load_obj("humanone_old_new_reaction")

oldID_list = ID_dict.keys()
with open('../../results/sampling_statisticalcomparison/lung_brain_breast_ratio.tsv', 'r') as csv_file:
    with open('../../results/sampling_statisticalcomparison/lung_brain_breast_ratio_newIDs.tsv', 'w') as output:
        for row in csv.reader(csv_file, delimiter='\t'):
            oldID = row[0]
            if oldID in oldID_list:
                newID = ID_dict.get(oldID)
                row[0] = newID
            output.write(str('\t'.join(row) + '\n'))
csv_file.close()
output.close()