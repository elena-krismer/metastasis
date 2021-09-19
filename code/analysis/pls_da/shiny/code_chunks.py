import pickle


"""
list subsystems
"""
def list_subsystems():
    with open('subsystem_dict.pkl', 'rb') as f:
        sub_dict =  pickle.load(f)
    return list(sub_dict.keys())

# get reactions from subsystem
def subsystem_reactions(subsystem):
    with open('subsystem_dict.pkl', 'rb') as f:
        dict =  pickle.load(f)
    sub_list = dict.get(subsystem)
    return sub_list

# reaction description
def get_reaction_description(humanone_reaction):
    with open('reaction_description_dict.pkl', 'rb') as f:
        dict = pickle.load(f)
    description = dict.get(humanone_reaction)
    return description


#print(get_reaction_description("MAR04693"))