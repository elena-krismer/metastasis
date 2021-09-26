import pickle
import csv

# save pickle object to obj folder
def save_obj(obj, name):
    with open('../../obj/dicts/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f,  protocol = 2)

def load_obj(name ):
    with open('../obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

subsystem_dict = {}
# create dict with reactions belong to each subsystem
with open('../../data/misc/sub_reactions.txt', 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        if len(row) != 0:
            subsystem = row[0]
            reactions = row[1:]
            subsystem_dict[subsystem] = reactions

csv_file.close()
save_obj(subsystem_dict, "subsystem_dict")


